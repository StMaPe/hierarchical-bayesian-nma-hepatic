# R/05_postprocess.R

# 05_postprocess.R
library(dplyr)
library(tidyr)
library(coda)

res <- readRDS("output/jags_fit.rds")
draws_mat <- as.matrix(res$samples)

params_d <- grep("^d\\[", colnames(draws_mat), value = TRUE)

Tn <- length(params_d) + 1
n_draws <- nrow(draws_mat)

# reconstruct full d matrix (d[1] = 0 reference)
d_samples <- matrix(NA_real_, nrow = n_draws, ncol = Tn)
colnames(d_samples) <- paste0("d[", seq_len(Tn), "]")

d_samples[, 1] <- 0
for (t in 2:Tn) {
  d_samples[, t] <- draws_mat[, paste0("d[", t, "]")]
}

contrast <- read.csv("output/contrast_hcc.csv", stringsAsFactors = FALSE)
contrast <- contrast %>% filter(!t1 %in% c('ALL OTHERS/NOT USING THE COMPARED DRUG','CONTROL/NOT USING ANYTHING'))
treatments <- sort(unique(c(contrast$t1, contrast$t2)))
tmap <- tibble(
  tid = seq_along(treatments),
  treatment = treatments
)

pair_grid <- expand.grid(
  a = seq_len(Tn),
  b = seq_len(Tn),
  stringsAsFactors = FALSE
) %>%
  filter(a != b)

pair_results <- pair_grid %>%
  rowwise() %>%
  mutate(
    a_name = tmap$treatment[a],
    b_name = tmap$treatment[b],
    diffs = list(d_samples[, a] - d_samples[, b]),
    est_logHR = mean(unlist(diffs)),
    lwr_logHR = quantile(unlist(diffs), 0.025),
    upr_logHR = quantile(unlist(diffs), 0.975),
    est_HR = exp(est_logHR),
    lwr_HR = exp(lwr_logHR),
    upr_HR = exp(upr_logHR)
  ) %>%
  ungroup() %>%
  select(a_name, b_name, est_logHR, lwr_logHR, upr_logHR,
         est_HR, lwr_HR, upr_HR) %>%
  arrange(a_name, b_name)

# SUCRA
ranks_matrix <- apply(d_samples, 1, rank, ties.method = "average")
ranks_matrix <- t(ranks_matrix)

Tminus1 <- Tn - 1
SUCRA <- sapply(seq_len(Tn), function(t) {
  cum_probs <- sapply(seq_len(Tminus1), function(k) {
    mean(ranks_matrix[, t] <= k)
  })
  sum(1 - cum_probs) / Tminus1
})

SUCRA_df <- tibble(
  Treatment = tmap$treatment,
  SUCRA = SUCRA,
  SUCRA_worse = 1 - SUCRA
) %>%
  arrange(desc(SUCRA_worse))

saveRDS(
  list(
    d_samples = d_samples,
    pairwise = pair_results,
    SUCRA = SUCRA_df
  ),
  file = "output/postprocess_results.rds"
)

write.csv(pair_results, "output/nma_pairwise_results.csv", row.names = FALSE)
write.csv(SUCRA_df, "output/sucra_results.csv", row.names = FALSE)

message("Post-processing completed:")
message("- Pairwise NMA effects")
message("- SUCRA rankings")
