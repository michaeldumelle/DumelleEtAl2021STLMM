#' Fit a Spatio-Temporal Linear Mixed Model
#'
#' @param data A data object containing all necessary variables.
#'
#' @param formula A formula of the form \code{y ~ x}, where \code{y} is the response variable
#'   and \code{x} are the predictor variables.
#'
#' @param xcoord A character vector specifying the column name of the x-coordinate
#'   variable in \code{data}.
#'
#' @param ycoord A character vector specifying the column name of the y-coordinate
#'   variable in \code{data}.
#'
#' @param tcoord A character vector specifying the column name of the t-coordinate (time)
#'   variable in \code{data}.
#'
#' @param stcov The spatio-temporal covariance type
#'  \describe{
#'    \item{\code{product}}{The product LMM}
#'    \item{\code{sum_with_error}}{The sum-with-error LMM}
#'    \item{\code{productsum}}{The product sum LMM}
#'  }
#'
#' @param estmethod The estimation method
#'  \describe{
#'    \item{\code{reml}}{Restricted Maximum Likelihood}
#'    \item{\code{svwls}}{Semivariogram Weighted Least Squares}
#'  }
#'
#' @param s_cor The spatial correlation
#'   \describe{
#'     \item{\code{exponential}}{The exponential correlation (the default).}
#'     \item{\code{spherical}}{The spherical correlation.}
#'     \item{\code{gaussian}}{The Gaussian correlation.}
#'   }
#'
#' @param t_cor The temporal correlation
#'   \describe{
#'     \item{\code{exponential}}{The exponential correlation (the default).}
#'     \item{\code{spherical}}{The spherical correlation.}
#'     \item{\code{gaussian}}{The Gaussian correlation.}
#'     \item{\code{tent}}{The tent (linear with sill) correlation.}
#'   }
#'
#' @param chol Should the Cholesky decomposition be used? If \code{FALSE},
#'   efficient inversion algorithms are implemented. Defaults to \code{FALSE}.
#'
#' @param condition A small number added to the diagonals of matrices before
#'   inverting them to prevent ill-conditioning (defaults to \code{1e-4}).
#'
#' @param logdet Should the log determinant be returned? (defaults to \code{FALSE}).
#'
#' @param weights Weights when \code{estmethod = "svwls"} (defaults to
#'   \code{"cressie"}) for the Cressie's weighted least squares weights).
#'
#' @param initial Initial values for the parameters. Must be made with
#'   \code{make_covparam_object()} (defaults to even spread of across variance
#'   parameters, a total variance matching the sample variance of OLS residuals,
#'   and ranges equaling half the maximum distance in the domain).
#'
#' @param optim_options A list containing additional options to pass to
#'   \link[stats]{optim}.
#'
#' @param h_options A list containing options to compute distances if
#'   \code{response}, \code{xcoord}, \code{ycoord}, and \code{tcoord} are
#'   provided. Named arguments are
#'   \describe{
#'     \item{\code{h_t_distmetric}}{The temporal distance matrix (defaults to
#'     \code{"euclidean"}).}
#'     \item{\code{h_s_distmetric}}{The spatial distance matrix (defaults to
#'     \code{"euclidean"}).}
#'  }
#'
#' @param max_options A list containing additonal options for placing upper
#' bounds on the total variance, spatial range, and temporal range. This can
#' be helpful for numerical stability in optimization.
#'   \describe{
#'     \item{\code{max_v}}{The maximum total variance (defaults to four times the
#'      variance of OLS residuals)}
#'     \item{\code{max_s_range}}{The maximum spatial range variance (defaults to
#'     four times the maximum observed spatial distance)}
#'     \item{\code{max_t_range}}{The maximum temporal range variance (defaults to
#'     four times the maximum observed temporal distance)}
#'   }
#'
#' @param stempsv_options A list containing additional options for the
#'   empirical spatio-temporal semivariogram. Named arguments are
#'   \describe{
#'     \item{\code{n_s_lag}}{The number of spatial distance classes (defaults to 16).}
#'     \item{\code{n_t_lag}}{The number of temporal distance classes (defaults to 16).}
#'     \item{\code{h_s_max}}{The maximum spatial distance. Deafaults to half the
#'       maximum distance in the spatial domain.}
#'     \item{\code{h_t_max}}{The maximum temporal distance. Deafaults to half the
#'       maximum distance in the temporal domain.}
#'   }
#'
#' @param ... Additonal arguments.
#'
#' @return A list containing several objects
#'   \describe{
#'     \item{\code{CovarianceParameters}}{Estimated covariance parameters.}
#'     \item{\code{Coefficients}}{Fixed effect estimates.}
#'     \item{\code{NamesCoefficients}}{Names of the fixed effect estimates.}
#'     \item{\code{CovCoefficients}}{The covariance matrix of the fixed effect estimates.}
#'     \item{\code{Objective}}{A list containing optimization information.}
#'     \item{\code{CovarianceForms}}{The spatial, temopral, and spatio-temporal correlation forms.}
#'     \item{\code{formula}}{The model formula.}
#'     \item{\code{model}}{A list containing the fixed effect design matrix and response vector.}
#'     \item{\code{data_object}}{An ordered data object.}
#'     \item{\code{invert_object}}{An inverse object.}
#'     \item{\code{coord_names}}{The names of the coordinate vectors.}
#'     \item{\code{coords}}{The coordinate vectors.}
#'     \item{\code{h_options}}{Returning the \code{h_options} argument.}
#'     \item{\code{stempsv_options}}{Returning the \code{stempsv_options} argument.}
#'     \item{\code{stempsv}}{The empirical spatio-temporal semivariogram (if \code{estmethod = "svwls"})}
#'     \item{\code{optim_options}}{Returning the \code{stempsv_options} argument.}
#'     \item{\code{max_options}}{Returning the \code{max_options} argument.}
#'     \item{\code{chol}}{Returning the \code{chol} argument.}
#'     \item{\code{condition}}{Returning the \code{condition} argument.}
#'     \item{\code{residuals}}{Raw residuals.}
#'   }
#'
#' @export
stlmm <- function(data, formula, ...) {
  UseMethod("stlmm", object = data)
}

#' @name stlmm
#' @method stlmm data.frame
#' @export
stlmm.data.frame <- function(data, formula, xcoord, ycoord = NULL, tcoord, stcov,
                             estmethod = "reml", s_cor = "exponential", t_cor = "exponential", chol = FALSE, condition = 1e-4,
                             logdet = FALSE, weights = "cressie", initial = NULL,
                             optim_options = NULL, h_options = NULL,
                             max_options = NULL, stempsv_options = NULL, ...) {

  # create the data object
  data_object <- make_data_object(
    formula = formula,
    xcoord = xcoord,
    ycoord = ycoord,
    tcoord = tcoord,
    data = data,
    h_options = h_options
  )

  # create the covest object
  covest_object <- make_covest_object(
    initial = initial,
    estmethod = estmethod,
    stcov = stcov,
    data_object = data_object,
    condition = condition,
    chol = chol,
    s_cor = s_cor,
    t_cor = t_cor,
    weights = weights,
    max_options = max_options,
    optim_options = optim_options,
    stempsv_options = stempsv_options
  )




  # estimate the profiled covariance parameters
  covest_output <- covest_wrapper(covest_object = covest_object, data_object = data_object)

  # give the estimate parameters the same class as the original covest_object
  class(covest_output$par_r) <- class(covest_object)

  # create the invert object
  invert_object <- make_invert_object(
    covparam_object = covest_output$par_r,
    chol = chol,
    condition = condition,
    h_s_large = data_object$h_s_large,
    h_t_large = data_object$h_t_large,
    h_s_small = data_object$h_s_small,
    h_t_small = data_object$h_t_small,
    logdet = logdet,
    m_index = data_object$m_index,
    o_index = data_object$o_index,
    s_cor = s_cor,
    t_cor = t_cor,
    xo = data_object$ordered_xo,
    yo = data_object$ordered_yo
  )

  # compute the inverse covariance matrix (times the design matrix and response)
  invert_output <- invert(invert_object)


  # estimate the fixed effects
  betaest_output <- betaest(
    xo = data_object$ordered_xo,
    sigmainv_xyo = invert_output$sigmainv_o,
    condition = condition,
    return_estlist = FALSE # don't need the extra output required for reml optimization
  )

  # return the relevant summary output
  stlmm_object <- structure(
    list(
      CovarianceParameters = covest_output$par_r,
      Coefficients = betaest_output$betahat,
      NamesCoefficients = colnames(data_object$original_xo),
      CovCoefficients = betaest_output$cov_betahat,
      Objective = c(
        value = covest_output$value,
        counts = covest_output$count,
        convergence = covest_output$convergence,
        stempsv_seconds = covest_object$stempsv_seconds,
        optim_seconds = covest_output$optim_seconds
      ),
      CovarianceForms = c(
        stcov = stcov,
        s_cor = s_cor,
        t_cor = t_cor
      ),
      formula = formula,
      model = list(
        FixedDesignMatrix = data_object$original_xo,
        Response = data_object$original_yo
      ),
      data_object = data_object,
      invert_output = invert_output,
      coordnames = list(
        xcoord = xcoord,
        ycoord = ycoord,
        tcoord = tcoord
      ),
      coords = list(
        xcoord = data_object$ordered_data_o[[xcoord]],
        ycoord = data_object$ordered_data_o[[ycoord]],
        tcoord = data_object$ordered_data_o[[tcoord]]
      ),
      h_options = data_object$h_options,
      stempsv_options = covest_object$stempsv_options,
      stempsv = covest_object$stempsv,
      optim_options = covest_object$optim_options,
      max_options = covest_object$max_options,
      chol = chol,
      condition = condition
    ),
    class = "stlmm" # give the object class "stlmm"
  )

  # compute the residuals
  stlmm_object$Residuals <- residuals(object = stlmm_object)

  # return the stlmm_object
  return(stlmm_object)
}
