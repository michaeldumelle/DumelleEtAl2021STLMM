residuals.stlmm <- function(stlmm_object, type = "raw"){

  # calling the function for the raw residuals
  if (type == "raw") {

  # computing the raw residuals
  residuals <- stlmm_object$model$Response - stlmm_object$model$FixedDesignMatrix %*% stlmm_object$Coefficients
  }

  # returning the type of the residuals as an attribute
  attr(residuals, "type") <- type

  # returning the residuals
  return(residuals)
}
