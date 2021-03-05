#' Summarize a Spatio-Temporal Linear Mixed Model
#'
#' @param object A \code{stlmm} object.
#' @param ... Additional arguments.
#'
#' @name summary
#'
#' @method summary stlmm
#'
#' @return A list containing several objects
#'   \describe{
#'     \item{\code{Call}}{The original function call.}
#'     \item{\code{FixedEffects}}{Fixed effects estimates and standard errors.}
#'     \item{\code{CovarianceParameters}}{Covariance parameter estimates.}
#'     \item{\code{CovarianceForms}}{The spatial, temopral, and spatio-temporal correlation forms.}
#'     \item{\code{Residuals}}{Raw residuals.}
#'     \item{\code{ObjectiveFn}}{Objective function value.}
#'   }
#'
#' @export
summary.stlmm <- function(object, ...) {

  # store the formula
  call <- object$formula

  # store the coefficients
  regcoefs <- as.vector(object$Coefficients)

  # store the covariance of beta hat
  regvar <- object$CovCoefficients

  # store the number of fixed effect predictors
  p <- ncol(object$model$FixedDesignMatrix)

  # store the number of observations
  n <- nrow(object$model$FixedDesignMatrix)

  # store the standard errors
  sereg <- sqrt(diag(as.matrix(regvar)))

  # store a vector of t statistics
  tvec <- regcoefs / sereg

  # store a vector of p-values appropriately rounded
  pvec <- round(100000 * (1 - pt(abs(regcoefs / sereg), df = n - p)) * 2) / 100000

  # save the fixed effect data frame
  fixed.effect.estimates <- data.frame(
    Estimate = regcoefs,
    Std.Error = sereg,
    t.value = tvec,
    prob.t = pvec
  )

  # assign row names
  rownames(fixed.effect.estimates) <- object$NamesCoefficients

  # save the covariance parameters as a list and then data frame
  covmodels <- as.list(object$CovarianceParameters)
  covmodelout <- data.frame(covmodels, stringsAsFactors = FALSE)

  # save the covariance information as a list and then data frame
  covinfo <- as.list(object$CovarianceForms)
  covinfoout <- data.frame(covinfo, stringsAsFactors = FALSE)

  # save the residuals
  resid_vec <- object$Residuals

  # save the objective functions
  objective <- object$Objective

  # store the final summary output
  output <- structure(
    list(
      Call = call,
      FixedEffects = fixed.effect.estimates,
      CovarianceParameters = covmodelout,
      CovarianceForms = covinfoout,
      Residuals = resid_vec,
      ObjectiveFn = objective
    ),
    class = "summary.stlmm"
  )

  # return the summary output
  return(output)
}
