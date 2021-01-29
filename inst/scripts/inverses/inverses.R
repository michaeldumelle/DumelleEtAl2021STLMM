# inverse times

write <- TRUE # set equal to TRUE if you want to write

set.seed(1)

library(microbenchmark)
library(tidyverse)

path <- "C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/STLMM/DumelleEtAl2021STLMM/notes/inverses/"
source(paste0(path, "conduct_inverses.R"))

n_st_sizes <- seq(1000, 15000, by = 1000)
sizes <- round(sqrt(n_st_sizes))
n_m <- 1
n_rep <- 100
n <- list(n_s = sizes, n_t = sizes, n_m = n_m, n_rep = n_rep)

#output <- pmap_dfr(.l = n, .f = conduct_inverses)

conduct_small <- function(index, n) {
  output <- conduct_inverses(n_s = n$n_s[index], n_t = n$n_t[index], n_m = n$n_m, n_rep = n$n_rep)
  write_csv(output, paste0(path, "conduct_small/", "output", index, ".csv"))
  return(NULL)
}
index <- 1:length(sizes)
map(index, conduct_small, n)


stempsv_full <- output %>%
  filter(algorithm == "stempsv_fast")

stempsv_times <- stempsv_full %>%
  group_by(n_st) %>%
  summarize(mean_time = mean(time))

ggplot(stempsv_times, aes(x = n_st, y = mean_time)) +
  geom_point() +
  geom_line()

inverse_full <- output %>%
  filter(algorithm != "stempsv_fast")

inverse_times <- inverse_full %>%
  group_by(n_st, algorithm) %>%
  summarize(mean_time = mean(time))

ggplot(inverse_times, aes(x = n_st, y = mean_time, col = algorithm)) +
  geom_line() +
  geom_point()

inverse_ratios <- inverse_full %>%
  pivot_wider(names_from = algorithm, values_from = time) %>%
  group_by(n_st) %>%
  summarize(psratio = mean(cholesky) / mean(productsum),
            sweratio = mean(cholesky) / mean(sum_with_error),
            pratio = mean(cholesky) / mean(product)
            ) %>%
  pivot_longer(c(psratio, sweratio, pratio), names_to = "algorithm", values_to = "mean_ratio")
# may be able to do grouped mutate here without pivoting

ggplot(inverse_ratios, aes(x = n_st, y = mean_ratio, col = algorithm)) +
  geom_point() +
  geom_line()

if (write) {
  library(readr)
  write_csv(stempsv_full, paste0(path, "stempsv_full.csv"))
  write_csv(stempsv_times, paste0(path, "stempsv_times.csv"))
  write_csv(inverse_full, paste0(path, "inverse_full.csv"))
  write_csv(inverse_times, paste0(path, "inverse_times.csv"))
  write_csv(inverse_ratios, paste0(path, "inverse_ratios.csv"))
}
