r2plo <- function(covparam_object, ...) {
  # calling the appropriate estimation method / stcov type generic (class/sub)
  UseMethod("r2plo", object = covparam_object)
}

# the regular to profiled log odds for semivariogram-weighted least squares
r2plo.svwls <- function(covparam_object, ...) {
  UseMethod("r2plo.svwls", object = covparam_object)
}

# the svwls productsum r2plo
r2plo.svwls.productsum <- function(covparam_object, max_options) {

  # giving covparam_object a new name
  params <- covparam_object

  # storing the spatial variance components
  s_vc <- c(params[["s_de"]], params[["s_ie"]])

  # computing the overall spatial variance
  svar <- sum(s_vc)

  # storing the temporal variance components
  t_vc <- c(params[["t_de"]], params[["t_ie"]])

  # computing the overall temporal variance
  tvar <- sum(t_vc)

  # storing the spatio-temporal variance components
  st_vc <- c(params[["st_de"]], params[["st_ie"]])

  # computing the overall spatio-temporal variance
  stvar <- sum(st_vc)

  # storing the overall variances
  vc <- c(s_vc, t_vc, st_vc)

  # computing the proportion of main effect variance
  lambda <- (svar + tvar) / (svar + tvar + stvar)

  # computing the proportion of main effect variance that is spatial
  alpha <- svar / (svar + tvar)

  # computing the proportion of spatial variance that is independent error
  n_s <- params[["s_ie"]] / svar

  # computing the proportion of temporal variance that is independent error
  n_t <- params[["t_ie"]] / tvar

  # computing the proportion of spatial-temporal that is independent error
  n_st <- params[["st_ie"]] / stvar

  # storing the profiled parameters
  pparm <- c(lambda = lambda, alpha = alpha, n_s = n_s, n_t = n_t, n_st = n_st)

  # storing the proportion of variance/spatial range/temporal range relative
  # to their maximums
  pparm <- c(
    pparm,
    var_prop = pmin(1, (svar + tvar + stvar) / max_options$max_v),
    srange_prop = params[["s_range"]] / max_options$max_s_range,
    trange_prop = params[["t_range"]] / max_options$max_t_range
  )

  # storing the logit of the profiled parameters
  pparm <- log(pparm / (1 - pparm))

  # returning the logit of the profiled parameters
  return(pparm)
}

# the svwls sum with error r2plo
r2plo.svwls.sum_with_error <- function(covparam_object, max_options) {

  # giving covparam_object a new name
  params <- covparam_object

  # storing the spatial variance components
  s_vc <- c(params[["s_de"]], params[["s_ie"]])

  # computing the overall spatial variance
  svar <- sum(s_vc)

  # storing the temporal variance components
  t_vc <- c(params[["t_de"]], params[["t_ie"]])

  # computing the overall temporal variance
  tvar <- sum(t_vc)

  # computing the spatio-temporal independent error
  st_vc <- params[["st_ie"]]

  # restoring it to make it clear this is the only parameter
  stvar <- sum(st_vc)

  # storing the overall variances
  vc <- c(s_vc, t_vc, st_vc)

  # computing the proportion of main effect variance
  lambda <- (svar + tvar) / (svar + tvar + stvar)

  # computing the proportion of main effect variance that is spatial
  alpha <- svar / (svar + tvar)

  # computing the proportion of spatial variance that is independent error
  n_s <- params[["s_ie"]] / svar

  # computing the proportion of temporal variance that is independent error
  n_t <- params[["t_ie"]] / tvar

  # storing the profiled parameters
  pparm <- c(lambda = lambda, alpha = alpha, n_s = n_s, n_t = n_t)

  # storing the proportion of variance/spatial range/temporal range relative
  # to their maximums
  pparm <- c(
    pparm,
    var_prop = pmin(1, (svar + tvar + stvar) / max_options$max_v),
    srange_prop = params[["s_range"]] / max_options$max_s_range,
    trange_prop = params[["t_range"]] / max_options$max_t_range
  )

  # storing the logit of the profiled parameters
  pparm <- log(pparm / (1 - pparm))

  # returning the logit of the profiled parameters
  return(pparm)
}

# the svwls product r2plo
r2plo.svwls.product <- function(covparam_object, max_options) {

  # giving covparam_object a new name
  params <- covparam_object

  # storing the profiled parameters
  pparm <- c(v_s = params[["v_s"]], v_t = params[["v_t"]])

  # storing the proportion of variance/spatial range/temporal range relative
  # to their maximums
  pparm <- c(
    pparm,
    var_prop = pmin(1, params[["st_de"]] / max_options$max_v),
    srange_prop = params[["s_range"]] / max_options$max_s_range,
    trange_prop = params[["t_range"]] / max_options$max_t_range
  )

  # storing the logit of the profiled parameters
  pparm <- log(pparm / (1 - pparm))

  # returning the logit of the profiled parameters
  return(pparm)
}

# the regular to profiled log odds for semivariogram-weighted least squares
r2plo.reml <- function(covparam_object, ...) {
  UseMethod("r2plo.reml", object = covparam_object)
}

# the reml productsum r2plo
r2plo.reml.productsum <- function(covparam_object, max_options) {
  # giving covparam_object a new name
  params <- covparam_object

  # storing the spatial variance components
  s_vc <- c(params[["s_de"]], params[["s_ie"]])

  # computing the overall spatial variance
  svar <- sum(s_vc)

  # storing the temporal variance components
  t_vc <- c(params[["t_de"]], params[["t_ie"]])

  # computing the overall temporal variance
  tvar <- sum(t_vc)

  # storing the spatio-temporal variance components
  st_vc <- c(params[["st_de"]], params[["st_ie"]])

  # computing the overall spatio-temporal variance
  stvar <- sum(st_vc)

  # storing the overall variances
  vc <- c(s_vc, t_vc, st_vc)

  # computing the proportion of main effect variance
  lambda <- (svar + tvar) / (svar + tvar + stvar)

  # computing the proportion of main effect variance that is spatial
  alpha <- svar / (svar + tvar)

  # computing the proportion of spatial variance that is independent error
  n_s <- params[["s_ie"]] / svar

  # computing the proportion of temporal variance that is independent error
  n_t <- params[["t_ie"]] / tvar

  # computing the proportion of spatial-temporal that is independent error
  n_st <- params[["st_ie"]] / stvar

  # storing the profiled parameters
  pparm <- c(lambda = lambda, alpha = alpha, n_s = n_s, n_t = n_t, n_st = n_st)

  # storing the proportion of spatial range/temporal range relative
  # to their maximums
  pparm <- c(
    pparm,
    srange_prop = params[["s_range"]] / max_options$max_s_range,
    trange_prop = params[["t_range"]] / max_options$max_t_range
  )

  # storing the logit of the profiled parameters
  pparm <- log(pparm / (1 - pparm))

  # returning the logit of the profiled parameters
  return(pparm)
}

# the reml sum with error r2plo
r2plo.reml.sum_with_error <- function(covparam_object, max_options) {

  # giving covparam_object a new name
  params <- covparam_object

  # storing the spatial variance components
  s_vc <- c(params[["s_de"]], params[["s_ie"]])

  # computing the overall spatial variance
  svar <- sum(s_vc)

  # storing the temporal variance components
  t_vc <- c(params[["t_de"]], params[["t_ie"]])

  # computing the overall temporal variance
  tvar <- sum(t_vc)

  # computing the spatio-temporal independent error
  st_vc <- params[["st_ie"]]

  # restoring it to make it clear this is the only parameter
  stvar <- sum(st_vc)

  # storing the overall variances
  vc <- c(s_vc, t_vc, st_vc)

  # computing the proportion of main effect variance
  lambda <- (svar + tvar) / (svar + tvar + stvar)

  # computing the proportion of main effect variance that is spatial
  alpha <- svar / (svar + tvar)

  # computing the proportion of spatial variance that is independent error
  n_s <- params[["s_ie"]] / svar

  # computing the proportion of temporal variance that is independent error
  n_t <- params[["t_ie"]] / tvar

  # storing the profiled parameters
  pparm <- c(lambda = lambda, alpha = alpha, n_s = n_s, n_t = n_t)

  # storing the proportion of variance/spatial range/temporal range relative
  # to their maximums
  pparm <- c(
    pparm,
    srange_prop = params[["s_range"]] / max_options$max_s_range,
    trange_prop = params[["t_range"]] / max_options$max_t_range
  )

  # storing the logit of the profiled parameters
  pparm <- log(pparm / (1 - pparm))

  # returning the logit of the profiled parameters
  return(pparm)
}

# the reml product r2plo
r2plo.reml.product <- function(covparam_object, max_options) {

  # giving covparam_object a new name
  params <- covparam_object

  # storing the profiled parameters
  pparm <- c(v_s = params[["v_s"]], v_t = params[["v_t"]])

  # storing the proportion of variance/spatial range/temporal range relative
  # to their maximums
  pparm <- c(
    pparm,
    srange_prop = params[["s_range"]] / max_options$max_s_range,
    trange_prop = params[["t_range"]] / max_options$max_t_range
  )

  # storing the logit of the profiled parameters
  pparm <- log(pparm / (1 - pparm))

  # returning the logit of the profiled parameters
  return(pparm)
}
