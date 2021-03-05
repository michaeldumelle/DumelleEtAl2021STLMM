covest <- function(par, covest_object, ...){
  UseMethod("covest", object = covest_object)
}

# semivariogram weighted least squares optimization
covest.svwls <- function(par, covest_object, data_object){

  # transform profiled variance parameters to regular
  plo2r <- plo2r.svwls(par, covest_object)

  # copy the class of covest_object to us the appropriate generic
  class(plo2r) <- class(covest_object)

  # make the spatio-temopral semivariogram
  theo_sv <- make_stsemivariogram(
    covparam_object = plo2r,
    h_s_large = covest_object$stempsv$h_s_avg,
    h_t_large = covest_object$stempsv$h_t_avg,
    s_cor = covest_object$s_cor,
    t_cor = covest_object$t_cor
  )
  # create the weights used in the weighted least squares optimization
  wts <- switch(
    covest_object$weights,
      "cressie" = weights_cressie(sv = covest_object$stempsv, theo_sv = theo_sv),
      stop("choose valid weights")
    )

  # create the objective function
  sumsq <- (covest_object$stempsv$gammahat - theo_sv)^2

  # return the objective function
  return(sum(wts * sumsq))
}

# create cressie (1985) weights
weights_cressie <- function(sv, theo_sv) {
  wts <- sv$n / theo_sv^2
  return(wts)
}

# reml optimization
covest.reml <- function(par, covest_object, invert_object){

  # transform profiled variance parameters to regular
  ## the overall variance is 1 because it has been profiled
  plo2r <- plo2r.reml(par, covest_object, ov_var = 1)

  # change the covariance parameter vecotr in invert_object
  invert_object$covparams <- plo2r

  # compute the inverse
  invert_output <- invert(invert_object)

  # compute minus twice the negative log likelihood
  m2ll <- minus2loglik.reml(invert_object = invert_object, invert_output = invert_output)

  # return minus twice the negative log likelihood
  return(m2ll)
}


