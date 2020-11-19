#' Title
#'
#' @param s_de
#' @param s_ie
#' @param t_de
#' @param t_ie
#' @param st_de
#' @param st_ie
#' @param v_s
#' @param v_t
#' @param s_range
#' @param t_range
#' @param estmethod
#' @param stcov
#' @import stats
#' @return
#' @export
#'
#' @examples
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



