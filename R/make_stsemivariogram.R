#' Make a Semivariogram Matrix
#'
#' @inheritParams make_stcovariance
#'
#' @return A semivariogram matrix
#'
#' @export
make_stsemivariogram <- function(covparam_object,
                                 h_s_large,
                                 h_t_large,
                                 s_cor,
                                 t_cor
                                 ) {

  # call the appropriate generic
  UseMethod("make_stsemivariogram", object = covparam_object)
}

#' @name make_stsemivariogram
#'
#' @method make_stsemivariogram productsum
#'
#' @export make_stsemivariogram.productsum
#' @export
make_stsemivariogram.productsum <- function(covparam_object,
                                            h_s_large,
                                            h_t_large,
                                            s_cor,
                                            t_cor
                                            ) {

  # make the productsum semivariogram
  # taking the variance parameters from the covparam_object
  variances <- c(covparam_object[c("s_de", "s_ie", "t_de", "t_ie", "st_de", "st_ie")])

  # computing the semivariogram
  gamma <- sum(variances) -
    make_stcovariance.productsum(
      covparam_object = covparam_object,
      h_s_large = h_s_large,
      h_t_large = h_t_large,
      s_cor = s_cor,
      t_cor = t_cor
    )

  # returning the semivariogram
  return(gamma)
}

#' @name make_stsemivariogram
#'
#' @method make_stsemivariogram sum_with_error
#'
#' @export make_stsemivariogram.sum_with_error
#' @export
make_stsemivariogram.sum_with_error <- function(covparam_object,
                                                h_s_large,
                                                h_t_large,
                                                s_cor,
                                                t_cor) {

  # taking the variance parameters from the covparam_object
  variances <- c(covparam_object[c("s_de", "s_ie", "t_de", "t_ie", "st_ie")])

  # computing the semivariogram
  gamma <- sum(variances) -
    make_stcovariance.sum_with_error(
      covparam_object = covparam_object,
      h_s_large = h_s_large,
      h_t_large = h_t_large,
      s_cor = s_cor,
      t_cor = t_cor
    )

  # returning the semivariogram
  return(gamma)
}

#' @name make_stsemivariogram
#'
#' @method make_stsemivariogram product
#'
#' @export make_stsemivariogram.product
#' @export
make_stsemivariogram.product <- function(covparam_object, h_s_large, h_t_large,
                                             s_cor, t_cor){

  # taking the variance parameters from the covparam_object
  variances <- c(covparam_object[c("st_de")])

  # computing the semivariogram
  gamma <- sum(variances) -
    make_stcovariance.product(
      covparam_object = covparam_object,
      h_s_large = h_s_large, h_t_large = h_t_large,
      s_cor = s_cor, t_cor = t_cor
    )

  # returning the semivariogram
  return(gamma)
}
