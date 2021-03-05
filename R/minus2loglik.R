minus2loglik <- function(invert_object, invert_output) {

  # appropriate generic for likelihood based estimation
  UseMethod("minus2loglik", object = invert_object)
}

# compute the restricted log likelihood
minus2loglik.reml <- function(invert_object, invert_output) {

  # storing the sample size
  n <- nrow(invert_output$sigmainv_o)

  # computing the current beta estimate
  betaest_output <- betaest(
    xo = invert_object$xo,
    sigmainv_xyo = invert_output$sigmainv_o,
    condition = invert_object$condition,
    return_estlist = TRUE
  )

  # computing sigma inverse times the residual vector
  siginv_r <- betaest_output$estlist$sigmainv_yo -
    betaest_output$estlist$sigmainv_xo %*% betaest_output$betahat

  # computing the quadratic residual function
  r_siginv_r <- t(invert_object$yo - invert_object$xo %*% betaest_output$betahat) %*%
    siginv_r

  # storing n - p
  nminusp <- n - betaest_output$estlist$p

  # computing the restricted likelihood
  m2ll <- invert_output$logdet +
    (nminusp) * (log(r_siginv_r) + 1 + log(2 * pi / nminusp)) +
    betaest_output$estlist$ldet_cicb

  # returning the restricted likelihood
  return(m2ll)
}
