#' Title
#'
#' @param xo
#' @param sigmainv_xyo
#' @param condition
#' @param return_estlist
#'
#' @import stats
#' @return
#' @export
#'
#' @examples
betaest <- function(xo, sigmainv_xyo, condition, return_estlist = FALSE){

  # number of columns of X (of the observed data)
  ncol_xo <- ncol(xo)

  # a sequence from 1 to the number of columns in X
  xo_dims <- seq.int(1, ncol_xo)

  # Taking sigma_inverse %*% X
  sigmainv_xo <- sigmainv_xyo[, xo_dims, drop = FALSE]

  # Taking sigma_inverse %*% Y
  sigmainv_yo <- sigmainv_xyo[, ncol_xo + 1, drop = FALSE]

  # Inverse of cov beta hat
  invcov_betahat <- t(xo) %*% sigmainv_xo

  # Adding diagonal stability
  diag(invcov_betahat) <- diag(invcov_betahat) + condition

  # The cholesky of this matrix
  chol_invcov_betahat <- chol(invcov_betahat)

  # Finding cov beta hat
  cov_betahat <- chol2inv(chol_invcov_betahat)

  # Computing beta hat
  betahat <- cov_betahat %*% t(xo) %*% sigmainv_yo

  # Returning betahat and the covariance in a list
  betaest_output <- list(betahat = betahat, cov_betahat = cov_betahat)

  # returning relevant output necessary for estimation
  if (return_estlist) {
    betaest_output$estlist <- list(
      # the log determinant of the inverse of cov beta hat
      ldet_cicb = 2 * sum(log(diag(chol_invcov_betahat))),
      sigmainv_xo = sigmainv_xo,
      sigmainv_yo = sigmainv_yo,
      p = ncol_xo
    )
  }
  return(betaest_output)
}
