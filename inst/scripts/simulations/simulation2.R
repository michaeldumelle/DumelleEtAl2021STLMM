# Simulation 2 Preliminaries --------------------------------------------------

# set to TRUE if you don't want to write csv's
write <- FALSE

# load the package
library(DumelleEtAl2021STLMM)

# load dplyr for summarizing later
# summarizing can be done with aggregate(formula, data, ...) in base R
library(dplyr)

# fix some stuff across simulations
## set number of trials
n_sim <- 2000
## sample sizes for spatial, temporal, missing
n <- list(n_s = 35, n_t = 35, n_m = 1)
## spatial and temporal correlation forms
cors <- list(s_cor = "exponential", t_cor = "tent")
## beta parameter values
beta <- list(beta0 = 0, beta1 = 0, beta2 = 0, beta3 = 0)


# Save the conduct simulations function ---------------------------------------
conduct_simulations <- function(simseed = NULL, n, covparams, cors, beta, error) {


  if (is.null(simseed)) {
    simseed <- sample.int(.Machine$integer.max, 1)
  }
  #running the simulations
  set.seed(simseed)

  ## sample sizes
  n_s <- n$n_s
  n_t <- n$n_t
  n_st <- n_s * n_t

  ## predictors
  pred_s <- rep(rnorm(n_s), times = n_t)
  pred_t <- rep(rnorm(n_t), each = n_s)
  pred_st <- rep(rnorm(n_st))

  ## mean
  beta0 <- beta$beta0
  beta1 <- beta$beta1
  beta2 <- beta$beta2
  beta3 <- beta$beta3
  mu <- beta0 + beta1 * pred_s + beta2 * pred_t + beta3 * pred_st

  ## coordinates
  xcoord <- rep(runif(n_s), times = n_t)
  ycoord <- rep(runif(n_s), times = n_t)
  #tcoord <- rep(runif(n_t), each = n_s)
  tcoord <- rep(seq(0, 1, length.out = n_t), each = n_s)

  ## data
  data <- data.frame(
    xcoord = xcoord,
    ycoord = ycoord,
    tcoord = tcoord,
    pred_s = pred_s,
    pred_t = pred_t,
    pred_st = pred_st,
    mu = mu
  )

  data <- data[with(data, order(tcoord, ycoord)), ]

  ## create large distance matrices
  h_s_large <- make_h(data$xcoord, data$ycoord, distmetric = "euclidean")
  h_t_large <- make_h(data$tcoord, distmetric = "euclidean")

  ## creating correlation forms
  s_cor <- cors$s_cor
  t_cor <- cors$t_cor

  ## setting covariance parameter values
  s_de <- covparams$s_de
  s_ie <- covparams$s_ie
  t_de <- covparams$t_de
  t_ie <- covparams$t_ie
  st_de <- covparams$st_de
  st_ie <- covparams$st_ie
  s_range <- covparams$s_range
  t_range <- covparams$t_range
  total_var <- sum(s_de + s_ie + t_de + t_ie + st_de + st_ie)

  ## setting the true covariance parameter vector
  covparams <- make_covparam_object(
    s_de = s_de,
    s_ie = s_ie,
    t_de = t_de,
    t_ie = t_ie,
    st_de = st_de,
    st_ie = st_ie,
    s_range = s_range,
    t_range = t_range,
    stcov = "productsum"
  )

  ## simulating the response
  stcovariance <- make_stcovariance(
    covparam_object = covparams,
    h_s_large = h_s_large,
    h_t_large = h_t_large,
    s_cor = s_cor,
    t_cor = t_cor
  )

  #if (error == "normal") {
  # data$response <- as.vector(strnorm(stcovariance, data$mu, size = 1))
  #}
  data$response <- as.vector(strnorm(
    object = covparams,
    mu = data$mu,
    size = 1,
    xcoord = "xcoord",
    ycoord = "ycoord",
    tcoord = "tcoord",
    s_cor = s_cor,
    t_cor = t_cor,
    error = error,
    data = data)
  )






  ## simulating missing values
  n_m <- n$n_m
  data$observed <- sample(c(rep(TRUE, n_st - n_m), rep(FALSE, n_m)), size = n_st, replace = F)


  ## make intial value function
  create_siminitial <- function(s_de, s_ie, t_de, t_ie, st_de, st_ie, s_range, t_range, sd_scaling) {

    sd_scaling <- 10
    s_de_initial <- pmax(0.1, s_de + rnorm(1, mean = 0, sd = s_de / sd_scaling))
    s_ie_initial <- pmax(0.1, s_ie + rnorm(1, mean = 0, sd = s_ie / sd_scaling))
    t_de_initial <- pmax(0.1, t_de + rnorm(1, mean = 0, sd = t_de / sd_scaling))
    t_ie_initial <- pmax(0.1, t_ie + rnorm(1, mean = 0, sd = t_ie / sd_scaling))
    st_de_initial <- pmax(0.1, st_de + rnorm(1, mean = 0, sd = st_de / sd_scaling))
    st_ie_initial <- pmax(0.1, st_ie + rnorm(1, mean = 0, sd = st_ie / sd_scaling))
    s_range_initial <- pmax(0.1, s_range + rnorm(1, mean = 0, sd = s_range / sd_scaling))
    t_range_initial <- pmax(0.1, t_range + rnorm(1, mean = 0, sd = s_range / sd_scaling))

    total_var_initial <- sum(s_de_initial + s_ie_initial + t_de_initial + t_ie_initial + st_de_initial + st_ie_initial)

    ps_reml <- initial(
      s_de = s_de_initial,
      s_ie = s_ie_initial,
      t_de = t_de_initial,
      t_ie = t_ie_initial,
      st_de = st_de_initial,
      st_ie = st_ie_initial,
      s_range = s_range_initial,
      t_range = t_range_initial,
      stcov = "productsum",
      estmethod = "reml"
    )

    ps_svwls <- initial(
      s_de = s_de_initial,
      s_ie = s_ie_initial,
      t_de = t_de_initial,
      t_ie = t_ie_initial,
      st_de = st_de_initial,
      st_ie = st_ie_initial,
      s_range = s_range_initial,
      t_range = t_range_initial,
      stcov = "productsum",
      estmethod = "svwls"
    )

    swe_scale <- total_var_initial / (total_var_initial - st_de)
    if (swe_scale == Inf) swe_scale <- 0.1

    swe_reml <- initial(
      s_de = swe_scale * s_de_initial,
      s_ie = swe_scale * s_ie_initial,
      t_de = swe_scale * t_de_initial,
      t_ie = swe_scale * t_ie_initial,
      st_ie = swe_scale * st_ie_initial,
      s_range = s_range_initial,
      t_range = t_range_initial,
      stcov = "sum_with_error",
      estmethod = "reml"
    )

    swe_svwls <- initial(
      s_de = swe_scale * s_de_initial,
      s_ie = swe_scale * s_ie_initial,
      t_de = swe_scale * t_de_initial,
      t_ie = swe_scale * t_ie_initial,
      st_ie = swe_scale * st_ie_initial,
      s_range = s_range_initial,
      t_range = t_range_initial,
      stcov = "sum_with_error",
      estmethod = "svwls"
    )

    p_reml <- initial(
      st_de = total_var_initial,
      v_s = s_ie_initial / (s_ie_initial + s_de_initial),
      v_t = t_ie_initial / (t_ie_initial + t_de_initial),
      s_range = s_range_initial,
      t_range = t_range_initial,
      stcov = "product",
      estmethod = "reml"
    )

    p_svwls <- initial(
      st_de = total_var_initial,
      v_s = s_ie_initial / (s_ie_initial + s_de_initial),
      v_t = t_ie_initial / (t_ie_initial + t_de_initial),
      s_range = s_range_initial,
      t_range = t_range_initial,
      stcov = "product",
      estmethod = "svwls"
    )

    siminitial <- list(
      ps_reml = ps_reml,
      ps_svwls = ps_svwls,
      swe_reml = swe_reml,
      swe_svwls = swe_svwls,
      p_reml = p_reml,
      p_svwls = p_svwls
    )
    return(siminitial)
  }

  ## construct initial values
  siminitial <- create_siminitial(s_de = s_de,
                                  s_ie = s_ie,
                                  t_de = t_de,
                                  t_ie = t_ie,
                                  st_de = st_de,
                                  st_ie = st_ie,
                                  s_range = s_range,
                                  t_range = t_range,
                                  sd_scaling = 10)

  ## subset the observed data
  data_o <-  subset(data, observed)

  ## subset the "missing" data
  data_m <- subset(data, !observed)

  ## create models function
  create_models <- function(data_o, siminitial) {


    ps_reml_mod <- stlmm(
      formula = response ~ pred_s + pred_t + pred_st,
      data = data_o,
      xcoord = "xcoord",
      ycoord = "ycoord",
      tcoord = "tcoord",
      stcov = "productsum",
      estmethod = "reml",
      s_cor = s_cor,
      t_cor = t_cor,
      initial = siminitial$ps_reml
    )

    ps_svwls_mod <- stlmm(
      formula = response ~ pred_s + pred_t + pred_st,
      data = data_o,
      xcoord = "xcoord",
      ycoord = "ycoord",
      tcoord = "tcoord",
      stcov = "productsum",
      estmethod = "svwls",
      s_cor = s_cor,
      t_cor = t_cor,
      initial = siminitial$ps_svwls
    )

    swe_reml_mod <- stlmm(
      formula = response ~ pred_s + pred_t + pred_st,
      data = data_o,
      xcoord = "xcoord",
      ycoord = "ycoord",
      tcoord = "tcoord",
      stcov = "sum_with_error",
      estmethod = "reml",
      s_cor = s_cor,
      t_cor = t_cor,
      initial = siminitial$swe_reml
    )

    swe_svwls_mod <- stlmm(
      formula = response ~ pred_s + pred_t + pred_st,
      data = data_o,
      xcoord = "xcoord",
      ycoord = "ycoord",
      tcoord = "tcoord",
      stcov = "sum_with_error",
      estmethod = "svwls",
      s_cor = s_cor,
      t_cor = t_cor,
      initial = siminitial$swe_svwls
    )

    p_reml_mod <- stlmm(
      formula = response ~ pred_s + pred_t + pred_st,
      data = data_o,
      xcoord = "xcoord",
      ycoord = "ycoord",
      tcoord = "tcoord",
      stcov = "product",
      estmethod = "reml",
      s_cor = s_cor,
      t_cor = t_cor,
      initial = siminitial$p_reml
    )

    p_svwls_mod <- stlmm(
      formula = response ~ pred_s + pred_t + pred_st,
      data = data_o,
      xcoord = "xcoord",
      ycoord = "ycoord",
      tcoord = "tcoord",
      stcov = "product",
      estmethod = "svwls",
      s_cor = s_cor,
      t_cor = t_cor,
      initial = siminitial$p_svwls
    )

    ols_mod <- lm(response ~ pred_s + pred_t + pred_st, data = data_o)

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

  ## fit models
  models <- create_models(data_o = data_o, siminitial = siminitial)

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
  predictions <- create_predictions(data_m = data_m, models = models)

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
      response = data_m$response,
      est = unname(predictions$ps_reml_predictions$fit),
      se = predictions$ps_reml_predictions$se.fit
    )

    prediction_ps_svwls <- data.frame(
      stcov = "productsum",
      estmethod = "svwls",
      response = data_m$response,
      est = unname(predictions$ps_svwls_predictions$fit),
      se = predictions$ps_svwls_predictions$se.fit
    )

    prediction_swe_reml <- data.frame(
      stcov = "sum_with_error",
      estmethod = "reml",
      response = data_m$response,
      est = unname(predictions$swe_reml_predictions$fit),
      se = predictions$swe_reml_predictions$se.fit
    )

    prediction_swe_svwls <- data.frame(
      stcov = "sum_with_error",
      estmethod = "svwls",
      response = data_m$response,
      est = unname(predictions$swe_svwls_predictions$fit),
      se = predictions$swe_svwls_predictions$se.fit
    )

    prediction_p_reml <- data.frame(
      stcov = "product",
      estmethod = "reml",
      response = data_m$response,
      est = unname(predictions$p_reml_predictions$fit),
      se = predictions$p_reml_predictions$se.fit
    )

    prediction_p_svwls <- data.frame(
      stcov = "product",
      estmethod = "svwls",
      response = data_m$response,
      est = unname(predictions$p_svwls_predictions$fit),
      se = predictions$p_svwls_predictions$se.fit
    )


    prediction_ols <- data.frame(
      stcov = "ind",
      estmethod = "ols",
      response = data_m$response,
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

    prediction$z <- (prediction$est - prediction$respons) / prediction$se

    return(prediction)

  }

  prediction <- get_prediction(data_m = data_m, predictions = predictions)

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

  output <- list(
    fixed = fixed,
    prediction = prediction,
    covparams = covparams,
    objectives = objectives,
    simseed = simseed
  )

  return(output)

}

# Run the conduct simulations function for simulation 2 -----------------------
set.seed(2)
# generate seeds for reproducibility
seeds <- as.list(sample.int(.Machine$integer.max, n_sim))
#set covariance parameter vector
covparams <- list(s_de = 18, s_ie = 0, t_de = 10, t_ie = 0,
                  st_de = 0, st_ie = 2, s_range = 0.5 * sqrt(2), t_range = 0.5 * 1)
error <- "normal"

# load parallel is using parallelization
library(parallel)
# number of cores (change if have many)
n_clust <- 1
# make the cluster object
cl <- makeCluster(n_clust)
clusterEvalQ(cl, library(DumelleEtAl2021STLMM))
# run the simulation
output <-  parLapply(cl, seeds, conduct_simulations, n = n, covparams = covparams, cors = cors, beta = beta, error = error)
# stop the cluster
stopCluster(cl)

# if not using parallel uncomment
# output <- lapply(seeds, conduct_simulations, n = n, covparams = covparams, cors = cors, beta = beta, error = error)


# Summarize the output --------------------------------------------------------
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


# Write the Output ------------------------------------------------------------
if (write) {
  # load readr if needed
  library(readr)
  write_csv(fixed, "inst/output/simulations/simulation2/fixed.csv")
  write_csv(predictions, "inst/output/simulations/simulation2/predictions.csv")
  write_csv(objectives,  "inst/output/simulations/simulation2/objectives.csv")
  write_csv(covparams, "inst/output/simulations/simulation2/covparams.csv")
  write_csv(seeds, "inst/output/simulations/simulation2/seeds.csv")
}
