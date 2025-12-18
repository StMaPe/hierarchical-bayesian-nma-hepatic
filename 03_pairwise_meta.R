# R/03_pairwise_meta.R

# 03_pairwise_meta.R
library(dplyr); library(metafor)
contrast <- read.csv('output/contrast_hcc.csv', stringsAsFactors = FALSE)
contrast <- contrast %>% filter(!t1 %in% c('ALL OTHERS/NOT USING THE COMPARED DRUG','CONTROL/NOT USING ANYTHING'))
contrast2 <- contrast %>% mutate(yi = logHR, vi = var)

pair_list <- split(contrast2, contrast2$pair)
res_list <- lapply(names(pair_list), function(p) {
  dfp <- pair_list[[p]]
  dfp2 <- dfp %>% filter(!is.na(yi) & !is.na(vi) & vi>0)
  if(nrow(dfp2) >= 2) {
    m <- tryCatch(rma(yi = dfp2$yi, vi = dfp2$vi, method = 'REML'), error = function(e) NULL)
    if(!is.null(m)) return(c(pair = p, est = as.numeric(m$b), se = as.numeric(m$se)))
  } else if(nrow(dfp2)==1) {
    return(c(pair = p, est = dfp2$yi[1], se = sqrt(dfp2$vi[1])))
  }
  return(NULL)
})
res_df <- do.call(rbind, res_list)
write.csv(res_df, 'output/pairwise_reml.csv', row.names = FALSE)
message('Pairwise REML completed.')
