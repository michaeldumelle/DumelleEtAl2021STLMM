#' Make the Inversion Object.
#'
#' @param covparam_object A covariance parameter object from \code{make_covparam_object()}.
#'
#' @param chol Should the Cholesky decomposition be used? If \code{FALSE},
#'   efficient inversion algorithms are implemented. Defaults to \code{FALSE}.
#'
#' @param co The covariance at prediction locations (if specified)
#'
#' @param condition A small number added to the diagonals of matrices before
#'   inverting them to prevent ill-conditioning (defaults to \code{1e-4}).
#'
#' @param h_s_large A spatial distance matrix of all spatio-temporal observations (if specified)
#'
#' @param h_t_large A temporal distance matrix of all spatio-temopral observations (if specified)
#'
#' @param h_s_small A spatial distance matrix of all spatial locations (if specified)
#'
#' @param h_t_small A temporal distance matrix of all temporal locations (if specified)
#'
#' @param logdet Should the log determinant be returned? (defaults to \code{FALSE}).
#'
#' @param m_index An index of missing values (from the space time cube).
#'
#' @param o_index An index of observed values (from the space time cube).
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
#' @param xo The fixed effects design matrix.
#'
#' @param yo A response vector.
#'
#' @return A list of relevant inversion information.
#'
#' @export
make_invert_object <- function(covparam_object,
                               chol,
                               co = NULL,
                               condition,
                               h_s_large = NULL,
                               h_t_large = NULL,
                               h_s_small = NULL,
                               h_t_small = NULL,
                               logdet,
                               m_index = NULL,
                               o_index = NULL,
                               s_cor,
                               t_cor,
                               xo,
                               yo
                               ) {

  # making the right multiplication matrix
  xyc_o <- cbind(xo, yo, co)

  # error message if small distance matrices are not provided and
  # cholesky is not being used
  if (!chol & (is.null(h_s_small) | is.null(h_t_small))){
    stop("If not using Cholesky decomposition, h_s_small and h_t_small must be provided")
  }

  # error message if large distance matrices are not provided and
  # cholesky is being used
  if (chol & (is.null(h_s_large) | is.null(h_t_large))){
    stop("If using Cholesky decomposition, h_s_large and h_t_large must be provided")
  }
    # storing the spatial dimension
    n_s <- nrow(h_s_small)

    # storing the temporal dimension
    n_t <- nrow(h_t_small)

    # creating the invert_object
    invert_object <- structure(list(
        covparams = covparam_object,
        chol = chol,
        condition = condition,
        logdet = logdet,
        h_s_small = h_s_small,
        h_t_small = h_t_small,
        h_s_large = h_s_large,
        h_t_large = h_t_large,
        o_index = o_index,
        m_index = m_index,
        n_s = n_s,
        n_t = n_t,
        s_cor = s_cor,
        t_cor = t_cor,
        xo = xo,
        yo = yo,
        co = co,
        xyc_o = cbind(xo, yo, co)),
      class = class(covparam_object)
    )
  return(invert_object)
}
