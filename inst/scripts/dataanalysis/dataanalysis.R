# conduct the data analysis

write <- TRUE #set to TRUE if you don't want to write csv's

path <- "C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/STLMM/DumelleEtAl2021STLMM/notes/dataanalysis/"
source(paste0(path, "conduct_dataanalysis.R"))

for (f in list.files("C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/STLMM/DumelleEtAl2021STLMM/R", pattern="*.R")) {
  source(paste("C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/STLMM/DumelleEtAl2021STLMM/R", f, sep = "/"))
}



or_train <- read.csv("C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/STLMM/DumelleEtAl2021STLMM/notes/dataanalysis/or_train.csv", stringsAsFactors = FALSE)
or_test <- read.csv("C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/STLMM/DumelleEtAl2021STLMM/notes/dataanalysis/or_test.csv", stringsAsFactors = FALSE)


set.seed(1)
dataanalysis <-  lapply(as.list(1), conduct_dataanalysis)


# return the fixed effect table
fixed <- dataanalysis[[1]]$fixed
print(fixed)

# return a summary of predictions
library(dplyr)
predictions_full <- dataanalysis[[1]]$prediction
predictions <- predictions_full %>%
  group_by(stcov, estmethod) %>%
  summarize(
    coverage = mean(abs(z) <= 1.96),
    mspe = sqrt(mean((response - est)^2)),
    mean_bias = mean(response - est)
  )
print(predictions)

# return a summary of covariance parameters
covparams <- dataanalysis[[1]]$covparams
print(covparams)


# return a summary of objectives
objectives <- dataanalysis[[1]]$objectives
print(objectives)


if (write) {
  library(readr)
  write_csv(fixed, paste0(path, "dataanalysis_fixed.csv"))
  write_csv(predictions_full, paste0(path, "dataanalysis_predictions_full.csv"))
  write_csv(predictions, paste0(path, "dataanalysis_predictions.csv"))
  write_csv(covparams, paste0(path, "dataanalysis_covparams.csv"))
  write_csv(objectives, paste0(path, "dataanalysis_objectives.csv"))
}
