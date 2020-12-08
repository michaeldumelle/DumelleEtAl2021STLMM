# conduct the data analysis

write <- FALSE #set to TRUE if you don't want to write csv's

path <- "C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/STLMM/DumelleEtAl2021STLMM/notes/dataanalysis/"
source(paste0(path, "conduct_dataanalysis.R"))

cl <- makeCluster(1) # clusterExport(cl, varlist = c()) exports objects to each cluster
clusterEvalQ(cl, library(purrr))
dataanalysis <-  parLapply(cl, "null_argument_for_parLapply", conduct_dataanalysis) #need an argument to run parLapply
stopCluster(cl)

# return the fixed effect table
fixed <- dataanalysis[[1]]$fixed
print(fixed)

# return a summary of predictions
library(dplyr)
predictions <- dataanalysis[[1]]$prediction %>%
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
  write.csv(fixed, paste0(path, "dataanalysis_fixed.csv"))
  write.csv(predictions, paste0(path, "dataanalysis_predictions.csv"))
  write.csv(covparams, paste0(path, "dataanalysis_covparams.csv"))
  write.csv(objectives, paste0(path, "dataanalysis_objectives.csv"))
}
