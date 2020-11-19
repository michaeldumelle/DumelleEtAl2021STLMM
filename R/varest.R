#' Title
#'
#' @param invert_object
#' @param invert_output
#' @import stats
#' @return
#' @export
#'
#' @examples
varest <- function(invert_object, invert_output) {

  # calling the appropriate generic
  UseMethod("varest", object = invert_object)
}

# compute the overall variance for reml estimation
varest.reml <- function(invert_object, invert_output) {

  # store the sample size
  n <- nrow(invert_output$sigmainv_o)

  # compute the beta estimate and return estmation information (via return_estlist)
  betaest_output <- betaest(
    xo = invert_object$xo,
    sigmainv_xyo = invert_output$sigmainv_o,
    condition = invert_object$condition,
    return_estlist = TRUE
  )

  # compute sigma inverse times the residual vector
  siginv_r <- betaest_output$estlist$sigmainv_yo - betaest_output$estlist$sigmainv_xo %*% betaest_output$betahat

  # compute the quadratic residual function
  r_siginv_r <- t(invert_object$yo - invert_object$xo %*% betaest_output$betahat) %*% siginv_r

  # compute the overall variance
  ws_l2 <- as.vector(r_siginv_r / (n - betaest_output$estlist$p))

  # return the overall variance
  return(ws_l2) # thanks Russ and John!
}
