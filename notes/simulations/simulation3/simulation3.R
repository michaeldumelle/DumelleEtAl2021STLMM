# simulation 3

write <- TRUE # set equal to TRUE if you want to write

library(parallel)

path <- "C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/STLMM/DumelleEtAl2021STLMM/notes/simulations/"
source(paste0(path, "conduct_simulations.R"))

# this is fixed across simulations
n_sim <- 2000
n <- list(n_s = 35, n_t = 35, n_m = 1)
cors <- list(s_cor = "exponential", t_cor = "tent")
beta <- list(beta0 = 0, beta1 = 0, beta2 = 0, beta3 = 0)
error <- "component_squared"

# this changes simulation to simulation
set.seed(3) # for simulation 3
seeds <- as.list(sample.int(.Machine$integer.max, n_sim))
covparams <- list(s_de = 4, s_ie = 4, t_de = 4, t_ie = 4,
                  st_de = 10, st_ie = 4, s_range = 0.5 * sqrt(2), t_range = 0.5 * 1)



cl <- makeCluster(detectCores()) # clusterExport(cl, varlist = c()) exports objects to each cluster
clusterEvalQ(cl, library(purrr))
output <-  parLapply(cl, seeds, conduct_simulations, n = n, covparams = covparams, cors = cors, beta = beta, error = error)
stopCluster(cl)


library(dplyr)
fixed_full <- lapply(output, function(x) x$fixed) %>%
  do.call(rbind, .) %>%
  mutate(trial = rep(1:n_sim, each = nrow(.) / n_sim))

fixed <- fixed_full %>%
  filter(beta != "beta0") %>%
  group_by(stcov, estmethod, beta) %>%
  summarize(typeone = mean(abs(z) > 1.96),
            rmse = sqrt(mean((est)^2)),
            mbias = mean(est)) %>%
  arrange(stcov)

print(fixed, n = Inf)

predictions_full <- lapply(output, function(x) x$prediction) %>%
  do.call(rbind, .) %>%
  mutate(trial = rep(1:n_sim, each = nrow(.) / n_sim))

predictions <- predictions_full %>%
  group_by(stcov, estmethod) %>%
  summarize(coverage = mean(abs(z) <= 1.96),
            rmspe = sqrt(mean((response - est)^2)),
            mbias = mean(response - est))

print(predictions, n = Inf)


covparams_full <- lapply(output, function(x) x$covparams) %>%
  do.call(rbind, .) %>%
  mutate(trial = rep(1:n_sim, each = nrow(.) / n_sim))

covparams <- covparams_full %>%
  group_by(stcov, estmethod, covparam) %>%
  summarize(lowq = quantile(value, 0.1),
            q1 = quantile(value, 0.25),
            median_result = median(value),
            mean_result = mean(value),
            q3 = quantile(value, 0.75),
            highq = quantile(value, 0.9)) %>%
  arrange(covparam, estmethod)

print(covparams, n = Inf)

objectives_full <- lapply(output, function(x) x$objectives) %>%
  do.call(rbind, .) %>%
  mutate(trial = rep(1:n_sim, each = nrow(.) / n_sim))

objectives <- objectives_full %>%
  group_by(stcov, estmethod, quantity) %>%
  summarize(mean_result = mean(result)) %>%
  arrange(quantity, estmethod)

print(objectives, n = Inf)

seeds_full <- lapply(output, function(x) x$simseed) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(trial = rep(1:n_sim, each = nrow(.) / n_sim))

colnames(seeds_full) <- c("seed", "trial")

seeds <- seeds_full


if (write) {
  library(readr)
  write_csv(fixed_full, paste0(path, "simulation3/simulation3_fixed_full.csv"))
  write_csv(fixed, paste0(path, "simulation3/simulation3_fixed.csv"))
  write_csv(predictions_full, paste0(path, "simulation3/simulation3_predictions_full.csv"))
  write_csv(predictions, paste0(path, "simulation3/simulation3_predictions.csv"))
  write_csv(objectives_full, paste0(path, "simulation3/simulation3_objectives_full.csv"))
  write_csv(objectives, paste0(path, "simulation3/simulation3_objectives.csv"))
  write_csv(covparams_full, paste0(path, "simulation3/simulation3_covparams_full.csv"))
  write_csv(covparams, paste0(path, "simulation3/simulation3_covparams.csv"))
  write_csv(seeds, paste0(path, "simulation3/simulation3_seeds.csv"))
}
