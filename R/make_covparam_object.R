#' Make a Covariance Parameter Object
#'
#' @param s_de The spatial dependent variance (spatial partial sill).
#'
#' @param s_ie The spatial independent variance (spatial nugget).
#'
#' @param t_de The temporal dependent variance (temporal partial sill).
#'
#' @param t_ie The temporal independent variance (temporal nugget).
#'
#' @param st_de The spatio-temporal dependent variance (spatio-temporal partial sill).
#'
#' @param st_ie The spatio-temporal independent variance (spatio-temporal nugget).
#'
#' @param v_s The proportion of spatial dependent variance
#'   (if \code{estmethod = "product"}).
#'
#' @param v_t The proportion of temporal dependent variance
#'   (if \code{estmethod = "product"}).
#'
#' @param s_range The spatial effective range (the spatial distance at which
#'   the correlation equals 0.05 (for non-compact) or 0 (for compact correlations)
#'
#' @param t_range The spatial effective range (the spatial distance at which
#'   the correlation equals 0.05 (for non-compact) or 0 (for compact correlations)
#'
#' @param estmethod The estimation method
#'  \describe{
#'    \item{\code{reml}}{Restricted Maximum Likelihood}
#'    \item{\code{svwls}}{Semivariogram Weighted Least Squares}
#'  }
#'
#' @param stcov The spatio-temporal covariance type
#'  \describe{
#'    \item{\code{product}}{The product LMM}
#'    \item{\code{sum_with_error}}{The sum-with-error LMM}
#'    \item{\code{productsum}}{The product sum LMM}
#'  }
#'
#' @return A named vector with covariance parameters having class equal to
#' the \code{estmethod} argument (if provided) and the \code{stcov} argument.
#'
#' @seealso [initial()]
#'
#' @export
make_covparam_object <- function(s_de,
                                 s_ie,
                                 t_de,
                                 t_ie,
                                 st_de,
                                 st_ie,
                                 v_s,
                                 v_t,
                                 s_range,
                                 t_range,
                                 estmethod = NULL,
                                 stcov) {

    # conditional call for type of spatio-temporal covariance
    covparam_object <- switch(
      stcov,
      "productsum" = covparam_object_productsum(
        s_de = s_de,
        s_ie = s_ie,
        t_de = t_de,
        t_ie = t_ie,
        st_de = st_de,
        st_ie = st_ie,
        s_range = s_range,
        t_range = t_range
      ),
      "sum_with_error" = covparam_object_sum_with_error(
        s_de = s_de,
        s_ie = s_ie,
        t_de = t_de,
        t_ie = t_ie,
        st_ie = st_ie,
        s_range = s_range,
        t_range = t_range
      ),
      "product" = covparam_object_product(
        st_de = st_de,
        v_s = v_s,
        v_t = v_t,
        s_range = s_range,
        t_range = t_range
      ),
      stop("Use a valid error structure")
    )

    # giving the object the appropriate estimation method and stcov class
    covparam_object <- structure(covparam_object, class = c(estmethod, stcov))

    # returning the object
    return(covparam_object)
}

# make product sum covariance parameter object
covparam_object_productsum <- function(s_de,
                                       s_ie,
                                       t_de,
                                       t_ie,
                                       st_de,
                                       st_ie,
                                       s_range,
                                       t_range
                                       ) {
  # create the covariance parameter vector
  cov_vec <- c(s_de = s_de,
               s_ie = s_ie,
               t_de = t_de,
               t_ie = t_ie,
               st_de = st_de,
               st_ie = st_ie,
               s_range = s_range,
               t_range = t_range
               )

  # return the covariance parameter vector
  return(cov_vec)
}

# make sum with error covariance parameter object
covparam_object_sum_with_error <- function(s_de,
                                           s_ie,
                                           t_de,
                                           t_ie,
                                           st_ie,
                                           s_range,
                                           t_range
                                           ) {
  # create the covariance parameter vector
  cov_vec <- c(s_de = s_de,
               s_ie = s_ie,
               t_de = t_de,
               t_ie = t_ie,
               st_ie = st_ie,
               s_range = s_range,
               t_range = t_range
               )

  # return the covariance parameter vector
  return(cov_vec)
}

# make product covariance parameter object
covparam_object_product <- function(st_de,
                                    v_s,
                                    v_t,
                                    s_range,
                                    t_range
                                    ) {
  # create the covariance parameter vector
  cov_vec <- c(st_de = st_de,
               v_s = v_s,
               v_t = v_t,
               s_range = s_range,
               t_range = t_range
               )

  # return the covariance parameter vector
  return(cov_vec)
}



