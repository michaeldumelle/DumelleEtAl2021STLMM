#' Print a Spatio-Temporal Linear Mixed Model Summary
#'
#' @param x A \code{stlmm} summary object.
#'
#' @param digits The number of printing digits.
#'
#' @param signif.stars The significance stars on p-values of the fixed effects.
#'
#' @param ... Additional arguments.
#'
#' @name print
#'
#' @method print summary.stlmm
#'
#' @return A printed summary of the spatio-temporal linear mixed model fit.
#'
#' @export
print.summary.stlmm <- function(x,
                                digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"),
                                ...) {

  # pasting the formula call
  cat("\nCall:\n", paste(deparse(x$Call), sep = "\n", collapse = "\n"), "\n", sep = "")

  # pasting the objective function
  cat("\nObjective Function:\n")
  print(x$Objective)

  # pasting the residual summary
  cat("\nResiduals:\n")
  resQ <- c(
    min(x$Residuals),
    quantile(x$Residuals,
      p = c(0.25, 0.5, 0.75),
      na.rm = TRUE
    ),
    max(x$Residuals)
  )
  names(resQ) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(resQ, digits = digits)

  # pasting the coefficient summary
  cat("\nCoefficients:\n")
  coefs <- x$FixedEffects
  colnames(coefs) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  # pasting the covariance parameter summary
  cat("\nCovariance Parameters:\n")
  print(x$CovarianceParameters)

  # pasting the covariance form summary
  cat("\nCovariance Forms:\n")
  print(x$CovarianceForms)
}






print.stlmm <- function(x, ...) {
  print(summary(x, ...))
}
