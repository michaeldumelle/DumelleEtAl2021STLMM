#' Title
#'
#' @param stlmm_object
#' @param type
#' @import stats
#' @return
#' @export
#'
#' @examples
residuals.stlmm <- function(object, type = "raw", ...){

  # calling the function for the raw residuals
  if (type == "raw") {

  # computing the raw residuals
  residuals <- object$model$Response - object$model$FixedDesignMatrix %*% object$Coefficients
  }

  # returning the type of the residuals as an attribute
  attr(residuals, "type") <- type

  # returning the residuals
  return(residuals)
}
