# R/06_viz_and_tables.R

# 06_viz_and_tables.R
library(dplyr)
library(tidyr)
library(gt)
library(scales)
library(igraph)
library(ggraph)
library(ggplot2)

#LEAGUE TABLE

post <- readRDS("output/postprocess_results.rds")
pair_results <- post$pairwise
contrast <- read.csv("output/contrast_hcc.csv", stringsAsFactors = FALSE)

df <- pair_results %>%
  mutate(
    cell_text = sprintf("%.2f (%.2f–%.2f)", est_HR, lwr_HR, upr_HR)
  )

treatments <- sort(unique(c(df$a_name, df$b_name)))
grid <- expand.grid(
  a_name = treatments,
  b_name = treatments,
  stringsAsFactors = FALSE
)

display_wide <- grid %>%
  left_join(df %>% select(a_name, b_name, cell_text),
            by = c("a_name", "b_name")) %>%
  mutate(display = if_else(a_name == b_name, "-", coalesce(cell_text, ""))) %>%
  select(a_name, b_name, display) %>%
  pivot_wider(names_from = b_name, values_from = display) %>%
  arrange(a_name)

hr_numeric_wide <- grid %>%
  left_join(df %>% select(a_name, b_name, est_HR),
            by = c("a_name", "b_name")) %>%
  mutate(est_HR = if_else(a_name == b_name, NA_real_, est_HR)) %>%
  pivot_wider(names_from = b_name, values_from = est_HR) %>%
  arrange(a_name)

gt_tbl <- display_wide %>%
  gt(rowname_col = "a_name") %>%
  tab_header(
    title = md("**League Table — Hazard Ratios**"),
    subtitle = md("Bayesian network meta-analysis (HR, 95% CrI)")
  ) %>%
  cols_align(align = "center", columns = everything()) %>%
  opt_table_outline() %>%
  tab_options(
    table.font.size = px(13),
    data_row.padding = px(6)
  )

hr_vals <- unlist(hr_numeric_wide %>% select(-a_name))
pal_fun <- col_numeric(
  palette = c("#2b83ba", "#ffffff", "#d7191c"),
  domain = range(hr_vals, na.rm = TRUE)
)

for (i in seq_len(nrow(hr_numeric_wide))) {
  for (j in seq_along(treatments)) {
    val <- hr_numeric_wide[[j + 1]][i]
    if (is.na(val)) next
    gt_tbl <- gt_tbl %>%
      tab_style(
        style = cell_fill(color = pal_fun(val)),
        locations = cells_body(columns = j + 1, rows = i)
      )
  }
}

gtsave(gt_tbl, "output/league_table.html")

#NETWORK PLOT
contrast <- contrast %>% filter(!t1 %in% c('ALL OTHERS/NOT USING THE COMPARED DRUG','CONTROL/NOT USING ANYTHING'))

edges <- contrast %>%
  select(study, t1, t2) %>%
  distinct() %>%
  group_by(t1, t2) %>%
  summarise(n_studies = n(), .groups = "drop")

edges <- edges %>%
  rowwise() %>%
  mutate(edge_id = paste(sort(c(t1, t2)), collapse = " vs ")) %>%
  ungroup() %>%
  distinct(edge_id, .keep_all = TRUE)

g <- graph_from_data_frame(d = edges, directed = FALSE)

set.seed(123)
network_plot <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = n_studies), color = "gray60", alpha = 0.7) +
  geom_node_point(size = 10, color = "maroon") +
  geom_node_text(aes(label = name), repel = TRUE, color = "black", size = 4) +
  scale_edge_width(range = c(0.5, 3)) +
  theme_void() +
  ggtitle("HCC") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

ggsave(
  filename = "output/network_hcc.pdf",
  plot = network_plot,
  width = 6,
  height = 6
)

message("Visualizations saved:")
message("- output/league_table.pdf")
message("- output/network_hcc.pdf")
