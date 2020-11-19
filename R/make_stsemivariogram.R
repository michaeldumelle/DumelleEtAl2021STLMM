#' Title
#'
#' @param covparam_object
#' @param h_s_large
#' @param h_t_large
#' @param s_cor
#' @param t_cor
#' @import stats
#' @return
#' @export
#'
#' @examples
make_stsemivariogram <- function(covparam_object,
                                 h_s_large,
                                 h_t_large,
                                 s_cor,
                                 t_cor
                                 ) {

  # call the appropriate generic
  UseMethod("make_stsemivariogram", object = covparam_object)
}

# make the productsum semivariogram
make_stsemivariogram.productsum <- function(covparam_object,
                                            h_s_large,
                                            h_t_large,
                                            s_cor,
                                            t_cor
                                            ) {

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

# make the sum with error semivariogram
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

# making the product semivariogram
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
