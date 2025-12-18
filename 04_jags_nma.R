# R/04_jags_nma.R

set.seed(25122004)

# 04_jags_nma.R
library(dplyr); library(rjags); library(coda); library(glue)
contrast <- read.csv('output/contrast_hcc.csv', stringsAsFactors = FALSE)
contrast <- contrast %>% filter(!t1 %in% c('ALL OTHERS/NOT USING THE COMPARED DRUG','CONTROL/NOT USING ANYTHING'))
contrast <- contrast %>% mutate(study = as.character(study), registry = as.character(registry))
contrast <- contrast %>% mutate(registry = if_else(is.na(registry) | registry=="", study, registry))

treatments <- sort(unique(c(contrast$t1, contrast$t2)))
Tn <- length(treatments)
tmap <- data.frame(treatment = treatments, tid = seq_len(Tn), stringsAsFactors = FALSE)

jags_data_df <- contrast %>% left_join(tmap, by = c('t1'='treatment')) %>% rename(t1id = tid) %>%
  left_join(tmap, by = c('t2'='treatment')) %>% rename(t2id = tid) %>%
  mutate(study_label = as.character(study), study_id = as.integer(as.factor(study_label)), db_label = as.character(registry), db_id = as.integer(as.factor(db_label)))

N <- nrow(jags_data_df); S <- length(unique(jags_data_df$study_id)); K <- length(unique(jags_data_df$db_id))

jags_data <- list(N = N, yi = jags_data_df$logHR, vi = jags_data_df$var, t1 = jags_data_df$t1id, t2 = jags_data_df$t2id, study = jags_data_df$study_id, db = jags_data_df$db_id, S = S, K = K, T = Tn)

# read model file
model_file <- 'models/nma_hierarchical_db_study.jags'
model_text <- paste(readLines(model_file), collapse = '\n')

inits_fun <- function(){ list(tau_db = runif(1,0,1), tau_study = runif(1,0,1), u = rnorm(K,0,0.1), s = rnorm(S,0,0.1)) }
params <- c(paste0('d[', 2:Tn, ']'), 'tau_db', 'tau_study')

run_jags_until_converge <- function(n.chains = 4, iter = 4000, burnin = 1000, adapt = 1000, max_extra = 5000, step_increase = 500){
  total_iter <- iter
  repeat {
    message(glue('Running JAGS: chains={n.chains}, iter={total_iter} (warmup={burnin}), adapt={adapt}'))
    jm <- jags.model(textConnection(model_text), data = jags_data, inits = replicate(n.chains, inits_fun(), simplify = FALSE), n.chains = n.chains, n.adapt = adapt)
    update(jm, burnin)
    samp <- coda.samples(jm, variable.names = params, n.iter = total_iter - burnin)
    gel <- tryCatch(gelman.diag(samp, multivariate = FALSE), error = function(e) NULL)
    max_psrf <- if(!is.null(gel)) max(gel$psrf[, 'Point est.'], na.rm = TRUE) else Inf
    message('Max PSRF: ', round(max_psrf,3))
    if(max_psrf <= 1.1) return(list(mod = jm, samples = samp, gelman = gel))
    if(total_iter + step_increase > iter + max_extra) { warning('No convergence under limits'); return(list(mod = jm, samples = samp, gelman = gel)) }
    total_iter <- total_iter + step_increase
  }
}

res <- run_jags_until_converge()
saveRDS(res, 'output/jags_fit.rds')
message('JAGS finished and saved to output/jags_fit.rds')
