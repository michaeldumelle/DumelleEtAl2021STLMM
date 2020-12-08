# conduct data anlaysis
conduct_dataanalysis <- function(null_argument_for_parLapply) {

  for (f in list.files("C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/STLMM/DumelleEtAl2021STLMM/R", pattern="*.R")) {
    source(paste("C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/STLMM/DumelleEtAl2021STLMM/R", f, sep = "/"))
  }
  or_train <- read.csv("C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/STLMM/DumelleEtAl2021STLMM/notes/dataanalysis/or_train.csv", stringsAsFactors = FALSE)
  or_test <- read.csv("C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/STLMM/DumelleEtAl2021STLMM/notes/dataanalysis/or_test.csv", stringsAsFactors = FALSE)




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
      initial = siminitial$ps_reml
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
      initial = siminitial$ps_svwls
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
      initial = siminitial$swe_reml
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
      initial = siminitial$swe_svwls
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
      initial = siminitial$p_reml
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
      initial = siminitial$p_svwls
    )

    ols_mod <- lm(TMAX ~ ELEVATION + TIMES + PRCP, data = data)

    ## cholesky models
    chol_ps_reml_mod <- stlmm(
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
      chol = TRUE
    )

    chol_swe_reml_mod <- stlmm(
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
      chol = TRUE
    )

    chol_p_reml_mod <- stlmm(
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
      chol = TRUE
    )




    models <- list(
      ps_reml_mod = ps_reml_mod,
      ps_svwls_mod = ps_svwls_mod,
      swe_reml_mod = swe_reml_mod,
      swe_svwls_mod = swe_svwls_mod,
      p_reml_mod = p_reml_mod,
      p_svwls_mod = p_svwls_mod,
      ols_mod = ols_mod,
      chol_ps_reml_mod = chol_ps_reml_mod,
      chol_swe_reml_mod = chol_swe_reml_mod,
      chol_p_reml_mod = chol_p_reml_mod
    )

    return(models)
  }
  s_cor <- "exponential"
  t_cor <- "exponential"
  siminitial <- list(ps_reml = NULL, ps_svwls = NULL, swe_reml = NULL, swe_svwls = NULL, p_reml = NULL, p_svwls = NULL)

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

    chol_ps_reml_objective <- data.frame(
      stcov = "productsum",
      estmethod = "reml_chol",
      quantity = c("value", "counts", "convergence", "stemp_seconds", "optim_seconds"),
      result = models$chol_ps_reml_mod$Objective[-3] # this is the counts.gradient output
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

    chol_swe_reml_objective <- data.frame(
      stcov = "sum_with_error",
      estmethod = "reml_chol",
      quantity = c("value", "counts", "convergence", "stemp_seconds", "optim_seconds"),
      result = models$chol_swe_reml_mod$Objective[-3] # this is the counts.gradient output
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

    chol_p_reml_objective <- data.frame(
      stcov = "product",
      estmethod = "reml_chol",
      quantity = c("value", "counts", "convergence", "stemp_seconds", "optim_seconds"),
      result = models$chol_p_reml_mod$Objective[-3] # this is the counts.gradient output
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
      chol_ps_reml_objective = chol_ps_reml_objective,
      ps_svwls_objective = ps_svwls_objective,
      swe_reml_objective = swe_reml_objective,
      chol_swe_reml_objective = chol_swe_reml_objective,
      swe_svwls_objective = swe_svwls_objective,
      p_reml_objective = p_reml_objective,
      chol_p_reml_objective = chol_p_reml_objective,
      p_svwls_objective = p_svwls_objective
    )

    objectives <- do.call(rbind, objectives)

    return(objectives)
  }

  objectives <- get_objectives(models = models)



  output <- list(fixed = fixed, prediction = prediction, covparams = covparams, objectives = objectives)
}


