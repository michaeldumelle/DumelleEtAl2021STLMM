covest_wrapper <- function(covest_object, data_object){
  UseMethod("covest_wrapper", object = covest_object)
}

# wrapper for the weighted least squarse optimization
covest_wrapper.svwls <- function(covest_object, data_object){

  # performing the optimization

  # timing
  optim_start <- Sys.time()
  covest_output <- optim(
    par = covest_object$initial_plo,
    fn = covest.svwls,
    covest_object = covest_object,
    data_object = data_object,
    method = covest_object$optim_options$method,
    control = covest_object$optim_options$control
  )
  optim_end <- Sys.time()
  covest_output$optim_seconds <- as.numeric(optim_end - optim_start, units = "secs")

  # set reasonable bounds on what the exponentiated plo parameters can be - this was a frustrating bug to find
  # covest_output$par[covest_output$par > 7] <- 7
  # covest_output$par[covest_output$par < -7] <- -7

  # transforming the optimized parameter values to regular values
  covest_output$par_r <- plo2r.svwls(par = covest_output$par, covest_object = covest_object)

  # a warning if convergence did not occur
  if (covest_output$convergence != 0) {
    warning("covariance parameter convergence may not have been achieved - consider
            setting new initial values, lowering the relative tolerance, or increasing
            the maximum iterations")
  }

  # returning the covariance parameter output
  return(covest_output)


}

covest_wrapper.reml <- function(covest_object, data_object){

  # newly storing the initial values (profiled log odds)
  initial_plo_noclass <- covest_object$initial_plo

  # giving this the same class as the covest_object
  class(covest_object$initial_plo) <- class(covest_object)

  # making the inverse object - many are from the data object
  invert_object <- make_invert_object(
    covparam_object = covest_object$initial_plo,
    chol = covest_object$chol,
    co = NULL,
    condition = covest_object$condition,
    h_s_large = data_object$h_s_large,
    h_t_large = data_object$h_t_large,
    h_s_small = data_object$h_s_small,
    h_t_small = data_object$h_t_small,
    logdet = covest_object$logdet,
    m_index = data_object$m_index,
    o_index = data_object$o_index,
    s_cor = covest_object$s_cor,
    t_cor = covest_object$t_cor,
    xo = data_object$ordered_xo,
    yo = data_object$ordered_yo
  )

  # performing the optimization

  # timing
  optim_start <- Sys.time()
  covest_output <- optim(
    par = initial_plo_noclass,
    fn = covest.reml,
    covest_object = covest_object,
    invert_object = invert_object,
    method = covest_object$optim_options$method,
    control = covest_object$optim_options$control
  )
  optim_end <- Sys.time()
  covest_output$optim_seconds <- as.numeric(optim_end - optim_start, units = "secs")


  # saving the profiled covariance parameter output
  invert_object$covparams <- plo2r.reml(covest_output$par, covest_object = covest_object, ov_var = 1)

  # computing the inverse
  invert_output <- invert(invert_object)

  # computing the overall variance
  ov_var <- varest.reml(invert_object = invert_object, invert_output = invert_output)

  # storing the regular covariance parameter output
  covest_output$par_r <- plo2r.reml(covest_output$par, covest_object = covest_object, ov_var = ov_var)

  # a warning if convergence did not occur
  if (covest_output$convergence != 0) {
    warning("covariance parameter convergence may not have been achieved - consider
            setting new initial values, lowering the relative tolerance, or increasing
            the maximum iterations")
  }

  # returning the covariance parameter output
  return(covest_output)


}


