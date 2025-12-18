# R/02_contrast_extraction.R

# 02_contrast_extraction.R

library(dplyr); library(tidyr)
source('R/utils.R')

df <- read.csv('data/cleaned_extraction.csv', stringsAsFactors = FALSE)

# column name for the outcome text
outcome_col <- 'HAZARD.RATIO..95..CI..FOR.HCC.DEVELOPMENT'

results <- vector('list', nrow(df))
diagnostics <- tibble::tibble(row = integer(), reason = character(), outcome_text = character())

for(i in seq_len(nrow(df))) {
  row <- df[i, , drop = FALSE]
  txt <- safe_get_char(row, outcome_col)
  v <- extract_hr_ci_robust(txt)
  if(all(is.na(v))) {
    diagnostics <- bind_rows(diagnostics, tibble(row = i, reason = 'no parsable HR/CI', outcome_text = txt))
    next
  }
  t1 <- safe_get_char(row, 'COMPARISON.ASSESED..FIRST.MEDICATION..INTERVENTION.')
  t2 <- safe_get_char(row, 'COMPARISON.ASSESED..SECOND.MEDICATION..COMPARISON.')
  results[[i]] <- tibble(
    study = safe_get_char(row, '.study_id'), author = safe_get_char(row, 'AUTHOR'), year = safe_get_char(row, 'YEAR'),
    registry = safe_get_char(row, 'registry'), t1 = t1, t2 = t2,
    hr = as.numeric(v['hr']), low = as.numeric(v['low']), high = as.numeric(v['high'])
  )
}
contrast <- bind_rows(results)
# SE derived assuming symmetric CI on the log scale
contrast <- contrast %>% mutate(logHR = log(hr), se = (log(high)-log(low))/(2*1.96), var = se^2,
                                pair = paste0(t1,' v ', t2)) %>% filter(!is.na(logHR) & !is.na(se))

contrast <- contrast %>% 
  mutate( swap = t1 > t2, 
          logHR = ifelse(swap, -logHR, logHR), 
          se = se, 
          t1_new = ifelse(swap, t2, t1), 
          t2_new = ifelse(swap, t1, t2), 
          pair = if_else(t1 <= t2, 
                         paste0(t1," vs ",t2), 
                         paste0(t2," vs ",t1)) ) %>% select(-t1, -t2) %>% rename(t1 = t1_new, t2 = t2_new)

if(!dir.exists('output')) dir.create('output')
write.csv(contrast, 'output/contrast_hcc.csv', row.names = FALSE)
write.csv(diagnostics, 'output/contrast_diagnostics.csv', row.names = FALSE)
message('Contrast table and diagnostics written to output/')
