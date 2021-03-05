# Preliminaries ---------------------------------------------------------------

# set to TRUE if you don't want to write csv's
write <- TRUE

# load the package
library(DumelleEtAl2021STLMM)

# load dplyr for summarizing later
# summarizing can be done with aggregate(formula, data, ...) in base R
library(dplyr)

# load the data
data("or_data")

# subset into training and test
or_train <- subset(or_data, TYPE == "TRAIN")
or_test <- subset(or_data, TYPE == "TEST")

# Save the conduct data analysis function -------------------------------------
# conduct data anlaysis
conduct_dataanalysis <- function(trial) {


  create_models <- function(data, siminitial) {


    ps_reml_mod <- stlmm(
      formula = TMAX ~ ELEVATION + TIMES + PRCP,
      data = data,
      xcoord = "LONGITUDE",
      ycoord = "LATITUDE",
      tcoord = "TIMES",
      stcov = "productsum",
      estmethod = "reml",
      s_cor = s_cor,
      t_cor = t_cor,
      initial = siminitial$ps_reml,
      condition = 0
    )

    ps_svwls_mod <- stlmm(
      formula = TMAX ~ ELEVATION + TIMES + PRCP,
      data = data,
      xcoord = "LONGITUDE",
      ycoord = "LATITUDE",
      tcoord = "TIMES",
      stcov = "productsum",
      estmethod = "svwls",
      s_cor = s_cor,
      t_cor = t_cor,
      initial = siminitial$ps_svwls,
      condition = 0
    )

    swe_reml_mod <- stlmm(
      formula = TMAX ~ ELEVATION + TIMES + PRCP,
      data = data,
      xcoord = "LONGITUDE",
      ycoord = "LATITUDE",
      tcoord = "TIMES",
      stcov = "sum_with_error",
      estmethod = "reml",
      s_cor = s_cor,
      t_cor = t_cor,
      initial = siminitial$swe_reml,
      condition = 0
    )

    swe_svwls_mod <- stlmm(
      formula = TMAX ~ ELEVATION + TIMES + PRCP,
      data = data,
      xcoord = "LONGITUDE",
      ycoord = "LATITUDE",
      tcoord = "TIMES",
      stcov = "sum_with_error",
      estmethod = "svwls",
      s_cor = s_cor,
      t_cor = t_cor,
      initial = siminitial$swe_svwls,
      condition = 0
    )

    p_reml_mod <- stlmm(
      formula = TMAX ~ ELEVATION + TIMES + PRCP,
      data = data,
      xcoord = "LONGITUDE",
      ycoord = "LATITUDE",
      tcoord = "TIMES",
      stcov = "product",
      estmethod = "reml",
      s_cor = s_cor,
      t_cor = t_cor,
      initial = siminitial$p_reml,
      condition = 0
    )

    p_svwls_mod <- stlmm(
      formula = TMAX ~ ELEVATION + TIMES + PRCP,
      data = data,
      xcoord = "LONGITUDE",
      ycoord = "LATITUDE",
      tcoord = "TIMES",
      stcov = "product",
      estmethod = "svwls",
      s_cor = s_cor,
      t_cor = t_cor,
      initial = siminitial$p_svwls,
      condition = 0
    )

    ols_mod <- lm(TMAX ~ ELEVATION + TIMES + PRCP, data = data)



    models <- list(
      ps_reml_mod = ps_reml_mod,
      ps_svwls_mod = ps_svwls_mod,
      swe_reml_mod = swe_reml_mod,
      swe_svwls_mod = swe_svwls_mod,
      p_reml_mod = p_reml_mod,
      p_svwls_mod = p_svwls_mod,
      ols_mod = ols_mod
    )

    return(models)
  }
  s_cor <- "exponential"
  t_cor <- "exponential"

  s_de_initial <- 40
  s_ie_initial <- 2
  t_de_initial <- 2
  t_ie_initial <- 2
  st_de_initial <- 7
  st_ie_initial <- 7
  total_var_initial <- sum(s_de_initial, s_ie_initial, t_de_initial,
                           t_ie_initial, st_de_initial, st_ie_initial)
  s_range_initial <- 1000
  t_range_initial <- 4
  swe_scale <- total_var_initial / (total_var_initial - st_de_initial)

  siminitial <- list(ps_reml = make_covparam_object(s_de = s_de_initial, s_ie = s_ie_initial,
                                                    t_ie = t_ie_initial, t_de = t_de_initial,
                                                    st_de = st_de_initial, st_ie = st_ie_initial,
                                                    s_range = s_range_initial, t_range = t_range_initial,
                                                    stcov = "productsum", estmethod = "reml"
  ),
  swe_reml = make_covparam_object(s_de = swe_scale * s_de_initial, s_ie = swe_scale * s_ie_initial,
                                  t_ie = swe_scale * t_de_initial, t_de = swe_scale * t_ie_initial,
                                  st_ie = swe_scale * st_ie_initial, s_range = s_range_initial, t_range = t_range_initial,
                                  stcov = "sum_with_error", estmethod = "reml"
  ),
  p_reml = make_covparam_object(st_de = total_var_initial,
                                v_s = s_ie_initial / (s_ie_initial + s_de_initial),
                                v_t = t_ie_initial / (t_ie_initial + t_de_initial),
                                s_range = s_range_initial, t_range = t_range_initial,
                                stcov = "product", estmethod = "reml"
  ),
  ps_svwls = make_covparam_object(s_de = s_de_initial, s_ie = s_ie_initial,
                                  t_ie = t_ie_initial, t_de = t_de_initial,
                                  st_de = st_de_initial, st_ie = st_ie_initial,
                                  s_range = s_range_initial, t_range = t_range_initial,
                                  stcov = "productsum", estmethod = "svwls"
  ),
  swe_reml = make_covparam_object(s_de = swe_scale * s_de_initial, s_ie = swe_scale * s_ie_initial,
                                  t_ie = swe_scale * t_de_initial, t_de = swe_scale * t_ie_initial,
                                  st_ie = swe_scale * st_ie_initial, s_range = s_range_initial, t_range = t_range_initial,
                                  stcov = "sum_with_error", estmethod = "svwls"
  ),
  p_reml = make_covparam_object(st_de = total_var_initial,
                                v_s = s_ie_initial / (s_ie_initial + s_de_initial),
                                v_t = t_ie_initial / (t_ie_initial + t_de_initial),
                                s_range = s_range_initial, t_range = t_range_initial,
                                stcov = "product", estmethod = "svwls"
  )
  )

  models <- create_models(or_train, siminitial)

  ## create predictions function
  create_predictions <- function(data_m, models) {
    ps_reml_predictions <- predict(models$ps_reml_mod, data_m, "prediction", se.fit = TRUE, predcov = FALSE)

    ps_svwls_predictions <- predict(models$ps_svwls_mod, data_m, "prediction", se.fit = TRUE, predcov = FALSE)

    swe_reml_predictions <- predict(models$swe_reml_mod, data_m, "prediction", se.fit = TRUE, predcov = FALSE)

    swe_svwls_predictions <- predict(models$swe_svwls_mod, data_m, "prediction", se.fit = TRUE, predcov = FALSE)

    p_reml_predictions <- predict(models$p_reml_mod, data_m, "prediction", se.fit = TRUE, predcov = FALSE)

    p_svwls_predictions <- predict(models$p_svwls_mod, data_m, "prediction", se.fit = TRUE, predcov = FALSE)

    ols_predictions <- predict(models$ols_mod, data_m, interval = "prediction", se.fit = TRUE, df = Inf)

    predictions <- list(
      ps_reml_predictions = ps_reml_predictions,
      ps_svwls_predictions = ps_svwls_predictions,
      swe_reml_predictions = swe_reml_predictions,
      swe_svwls_predictions = swe_svwls_predictions,
      p_reml_predictions = p_reml_predictions,
      p_svwls_predictions = p_svwls_predictions,
      ols_predictions = ols_predictions
    )
  }

  ## create predictions
  predictions <- create_predictions(data_m = or_test, models = models)

  get_fixed <- function(models) {

    fixed_ps_reml <- data.frame(
      stcov = "productsum",
      estmethod = "reml",
      beta = c("beta0", "beta1", "beta2", "beta3"),
      est = unname(models$ps_reml_mod$Coefficients),
      se = sqrt(diag(models$ps_reml_mod$CovCoefficients))
    )

    fixed_ps_svwls <- data.frame(
      stcov = "productsum",
      estmethod = "svwls",
      beta = c("beta0", "beta1", "beta2", "beta3"),
      est = unname(models$ps_svwls_mod$Coefficients),
      se = sqrt(diag(models$ps_svwls_mod$CovCoefficients))
    )

    fixed_swe_reml <- data.frame(
      stcov = "sum_with_error",
      estmethod = "reml",
      beta = c("beta0", "beta1", "beta2", "beta3"),
      est = unname(models$swe_reml_mod$Coefficients),
      se = sqrt(diag(models$swe_reml_mod$CovCoefficients))
    )

    fixed_swe_svwls <- data.frame(
      stcov = "sum_with_error",
      estmethod = "svwls",
      beta = c("beta0", "beta1", "beta2", "beta3"),
      est = unname(models$swe_svwls_mod$Coefficients),
      se = sqrt(diag(models$swe_svwls_mod$CovCoefficients))
    )

    fixed_p_reml <- data.frame(
      stcov = "product",
      estmethod = "reml",
      beta = c("beta0", "beta1", "beta2", "beta3"),
      est = unname(models$p_reml_mod$Coefficients),
      se = sqrt(diag(models$p_reml_mod$CovCoefficients))
    )

    fixed_p_svwls <- data.frame(
      stcov = "product",
      estmethod = "svwls",
      beta = c("beta0", "beta1", "beta2", "beta3"),
      est = unname(models$p_svwls_mod$Coefficients),
      se = sqrt(diag(models$p_svwls_mod$CovCoefficients))
    )

    fixed_ols <- data.frame(
      stcov = "ind",
      estmethod = "ols",
      beta = c("beta0", "beta1", "beta2", "beta3"),
      est = unname(coef(models$ols_mod)),
      se = unname(sqrt(diag(vcov(models$ols_mod))))
    )

    fixed <- list(
      fixed_ps_reml,
      fixed_ps_svwls,
      fixed_swe_reml,
      fixed_swe_svwls,
      fixed_p_reml,
      fixed_p_svwls,
      fixed_ols
    )

    fixed <- do.call(rbind, fixed)

    fixed$z <- fixed$est / fixed$se

    return(fixed)
  }

  fixed <- get_fixed(models = models)

  get_prediction <- function(data_m, predictions) {

    prediction_ps_reml <- data.frame(
      stcov = "productsum",
      estmethod = "reml",
      response = data_m$TMAX,
      est = unname(predictions$ps_reml_predictions$fit),
      se = predictions$ps_reml_predictions$se.fit
    )

    prediction_ps_svwls <- data.frame(
      stcov = "productsum",
      estmethod = "svwls",
      response = data_m$TMAX,
      est = unname(predictions$ps_svwls_predictions$fit),
      se = predictions$ps_svwls_predictions$se.fit
    )

    prediction_swe_reml <- data.frame(
      stcov = "sum_with_error",
      estmethod = "reml",
      response = data_m$TMAX,
      est = unname(predictions$swe_reml_predictions$fit),
      se = predictions$swe_reml_predictions$se.fit
    )

    prediction_swe_svwls <- data.frame(
      stcov = "sum_with_error",
      estmethod = "svwls",
      response = data_m$TMAX,
      est = unname(predictions$swe_svwls_predictions$fit),
      se = predictions$swe_svwls_predictions$se.fit
    )

    prediction_p_reml <- data.frame(
      stcov = "product",
      estmethod = "reml",
      response = data_m$TMAX,
      est = unname(predictions$p_reml_predictions$fit),
      se = predictions$p_reml_predictions$se.fit
    )

    prediction_p_svwls <- data.frame(
      stcov = "product",
      estmethod = "svwls",
      response = data_m$TMAX,
      est = unname(predictions$p_svwls_predictions$fit),
      se = predictions$p_svwls_predictions$se.fit
    )


    prediction_ols <- data.frame(
      stcov = "ind",
      estmethod = "ols",
      response = data_m$TMAX,
      est = unname(predictions$ols_predictions$fit[, "fit"]),
      se = unname(sqrt(predictions$ols_predictions$se.fit^2 + predictions$ols_predictions$residual.scale^2))
    )

    prediction <- list(
      prediction_ps_reml,
      prediction_ps_svwls,
      prediction_swe_reml,
      prediction_swe_svwls,
      prediction_p_reml,
      prediction_p_svwls,
      prediction_ols
    )

    prediction <- do.call(rbind, prediction)

    prediction$z <- (prediction$est - prediction$response) / prediction$se

    return(prediction)

  }

  get_covparams <- function(models) {

    ps_reml_covparams <- data.frame(
      stcov = "productsum",
      estmethod = "reml",
      covparam = names(unclass(models$ps_reml_mod$CovarianceParameters)),
      value = unclass(models$ps_reml_mod$CovarianceParameters)
    )

    ps_svwls_covparams <- data.frame(
      stcov = "productsum",
      estmethod = "svwls",
      covparam = names(unclass(models$ps_svwls_mod$CovarianceParameters)),
      value = unclass(models$ps_svwls_mod$CovarianceParameters)
    )

    swe_reml_covparams <- data.frame(
      stcov = "sum_with_error",
      estmethod = "reml",
      covparam = names(unclass(models$swe_reml_mod$CovarianceParameters)),
      value = unclass(models$swe_reml_mod$CovarianceParameters)
    )

    swe_svwls_covparams <- data.frame(
      stcov = "sum_with_error",
      estmethod = "svwls",
      covparam = names(unclass(models$swe_svwls_mod$CovarianceParameters)),
      value = unclass(models$swe_svwls_mod$CovarianceParameters)
    )

    p_reml_covparams <- data.frame(
      stcov = "product",
      estmethod = "reml",
      covparam = names(unclass(models$p_reml_mod$CovarianceParameters)),
      value = unclass(models$p_reml_mod$CovarianceParameters)
    )

    p_svwls_covparams <- data.frame(
      stcov = "product",
      estmethod = "svwls",
      covparam = names(unclass(models$p_svwls_mod$CovarianceParameters)),
      value = unclass(models$p_svwls_mod$CovarianceParameters)
    )

    covparams <- list(
      ps_reml_covparams = ps_reml_covparams,
      ps_svwls_covparams = ps_svwls_covparams,
      swe_reml_covparams = swe_reml_covparams,
      swe_svwls_covparams = swe_svwls_covparams,
      p_reml_covparams = p_reml_covparams,
      p_svwls_covparams = p_svwls_covparams
    )

    covparams <- do.call(rbind, covparams)

    return(covparams)

  }

  covparams <- get_covparams(models = models)

  prediction <- get_prediction(data_m = or_test, predictions = predictions)

  get_objectives <- function(models) {

    ps_reml_objective <- data.frame(
      stcov = "productsum",
      estmethod = "reml",
      quantity = c("value", "counts", "convergence", "stemp_seconds", "optim_seconds"),
      result = models$ps_reml_mod$Objective[-3] # this is the counts.gradient output
    )


    ps_svwls_objective <- data.frame(
      stcov = "productsum",
      estmethod = "svwls",
      quantity = c("value", "counts", "convergence", "stemp_seconds", "optim_seconds"),
      result = models$ps_svwls_mod$Objective[-3] # this is the counts.gradient output
    )

    swe_reml_objective <- data.frame(
      stcov = "sum_with_error",
      estmethod = "reml",
      quantity = c("value", "counts", "convergence", "stemp_seconds", "optim_seconds"),
      result = models$swe_reml_mod$Objective[-3] # this is the counts.gradient output
    )


    swe_svwls_objective <- data.frame(
      stcov = "sum_with_error",
      estmethod = "svwls",
      quantity = c("value", "counts", "convergence", "stemp_seconds", "optim_seconds"),
      result = models$swe_svwls_mod$Objective[-3] # this is the counts.gradient output
    )

    p_reml_objective <- data.frame(
      stcov = "product",
      estmethod = "reml",
      quantity = c("value", "counts", "convergence", "stemp_seconds", "optim_seconds"),
      result = models$p_reml_mod$Objective[-3] # this is the counts.gradient output
    )


    p_svwls_objective <- data.frame(
      stcov = "product",
      estmethod = "svwls",
      quantity = c("value", "counts", "convergence", "stemp_seconds", "optim_seconds"),
      result = models$p_svwls_mod$Objective[-3] # this is the counts.gradient output
    )

    #storing the average stempsv time as the exact quantity is computed 3 times
    stempsv_average <- mean(c(ps_svwls_objective$result[ps_svwls_objective$quantity == "stemp_seconds"],
                              swe_svwls_objective$result[swe_svwls_objective$quantity == "stemp_seconds"],
                              p_svwls_objective$result[p_svwls_objective$quantity == "stemp_seconds"]))
    ps_svwls_objective$result[ps_svwls_objective$quantity == "stemp_seconds"] = stempsv_average
    swe_svwls_objective$result[swe_svwls_objective$quantity == "stemp_seconds"] = stempsv_average
    p_svwls_objective$result[p_svwls_objective$quantity == "stemp_seconds"] = stempsv_average

    objectives <- list(
      ps_reml_objective = ps_reml_objective,
      ps_svwls_objective = ps_svwls_objective,
      swe_reml_objective = swe_reml_objective,
      swe_svwls_objective = swe_svwls_objective,
      p_reml_objective = p_reml_objective,
      p_svwls_objective = p_svwls_objective
    )

    objectives <- do.call(rbind, objectives)

    return(objectives)
  }

  objectives <- get_objectives(models = models)



  output <- list(fixed = fixed, prediction = prediction, covparams = covparams, objectives = objectives, stempsv = models$ps_svwls_mod$stempsv)
}




# Run the data analysis function ----------------------------------------------
dataanalysis <-  lapply(as.list(1), conduct_dataanalysis)

# Summarize the output --------------------------------------------------------

# return the fixed effect table
fixed <- dataanalysis[[1]]$fixed
fixed$pvalue <- 2 * pnorm(abs(fixed$z), lower.tail = FALSE)


# return a summary of predictions
predictions_full <- dataanalysis[[1]]$prediction
predictions <- predictions_full %>%
  group_by(stcov, estmethod) %>%
  summarize(
    coverage = mean(abs(z) <= 1.96),
    mspe = sqrt(mean((response - est)^2)),
    mean_bias = mean(response - est)
  )


# return a summary of covariance parameters
covparams <- dataanalysis[[1]]$covparams



# return a summary of objectives
objectives <- dataanalysis[[1]]$objectives


# Write the Output ------------------------------------------------------------
if (write) {
  # load readr if needed
  library(readr)
  write_csv(fixed, "inst/output/dataanalysis/fixed.csv")
  write_csv(predictions, "inst/output/dataanalysis/predictions.csv")
  write_csv(covparams, "inst/output/dataanalysis/covparams.csv")
  write_csv(objectives, "inst/output/dataanalysis/objectives.csv")
}
