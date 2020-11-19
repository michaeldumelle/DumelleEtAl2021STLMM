library(parallel)
library(tidyverse)

fsims <- function(nosimseed) {
  # running the simulations
  #set.seed(simseed)
  library(purrr)
  for (f in list.files("C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/STLMM/DumelleEtAl2021STLMM/R", pattern="*.R")) {
    source(paste("C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/STLMM/DumelleEtAl2021STLMM/R", f, sep = "/"))
  }

  ## sample sizes
  n_s <- 30
  n_t <- 30
  n_st <- n_s * n_t

  ## predictors
  pred_s <- rep(rnorm(n_s), times = n_t)
  pred_t <- rep(rnorm(n_t), each = n_s)
  pred_st <- rep(rnorm(n_st))

  ## mean
  beta0 <- 0
  beta1 <- 0
  beta2 <- 0
  beta3 <- 0
  mu <- beta0 + beta1 * pred_s + beta2 * pred_t + beta3 * pred_st

  ## coordinates
  xcoord <- rep(runif(n_s), times = n_t)
  ycoord <- rep(runif(n_s), times = n_t)
  tcoord <- rep(runif(n_t), each = n_s)

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

  ## create large distance matrices
  h_s_large <- make_h(data$xcoord, data$ycoord, distmetric = "euclidean")
  h_t_large <- make_h(data$tcoord, distmetric = "euclidean")

  ## creating correlation forms
  s_cor <- "exponential"
  t_cor <- "tent"

  ## setting covariance parameter values
  s_de <- 5
  s_ie <- 5
  t_de <- 5
  t_ie <- 5
  st_de <- 5
  st_ie <- 5
  s_range <- sqrt(2)
  t_range <- 1
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

  data$response <- as.vector(strnorm(stcovariance, data$mu, size = 1))


  ## simulating missing values
  n_m <- 25
  data$observed <- sample(c(rep(TRUE, n_st - n_m), rep(FALSE, n_m)), size = n_st, replace = F)

  datafail <<- data
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
  data_o <- data %>%
    dplyr::filter(observed)

  ## subset the "missing" data
  data_m <- data %>%
    dplyr::filter(!observed)

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

  output <- list(fixed = fixed, prediction = prediction)

  return(output)

}
n_sim <- 1000
cl <- makeCluster(45)
clusterExport(cl, "fsim")
clusterEvalQ(cl, fsim)
test = parLapply(cl, 1:n_sim, fsim)
stopCluster(cl)

testfixed <- lapply(test, function(x) x$fixed)
testfixed <- do.call(rbind, testfixed)
summary(testfixed$se)
fixed <- testfixed %>%
  filter(beta != "beta0") %>%
  group_by(stcov, estmethod, beta) %>%
  summarize(typeone = mean(abs(z) > 1.96),
            mse = sqrt(mean((est)^2)),
            bias = mean(est)) %>%
  arrange(stcov)
print(fixed, n = Inf)
testpreds <- lapply(test, function(x) x$prediction)
testpreds <- do.call(rbind, testpreds)
summary(testpreds$se)

pred <- testpreds %>%
  group_by(stcov, estmethod) %>%
  summarize(coverage = mean(abs(z) <= 1.96),
            mspe = sqrt(mean((response - est)^2)),
            bias = mean(response - est))
pred
