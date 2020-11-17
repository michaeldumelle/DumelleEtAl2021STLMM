make_stcovariance <- function(covparam_object,
                              h_s_large,
                              h_t_large,
                              s_cor,
                              t_cor
                              ) {

  # call the appropriate generic
  UseMethod("make_stcovariance", object = covparam_object)
}

make_stcovariance.productsum <- function(covparam_object,
                                         h_s_large,
                                         h_t_large,
                                         s_cor,
                                         t_cor
                                         ) {
  # creating the spatial correlation matrix
  r_s_large <- make_r(
    h = h_s_large,
    range = covparam_object[["s_range"]],
    structure = s_cor
  )

  # creating the temporal correlation matrix
  r_t_large <- make_r(
    h = h_t_large,
    range = covparam_object[["t_range"]],
    structure = t_cor
  )

  # creating the spatio-temporal correlation matrix
  r_st <- r_s_large * r_t_large

  # make the spatial sigma matrix
  sigma_s <- make_sigma(
    de = covparam_object[["s_de"]],
    r_mx = r_s_large,
    ie = covparam_object[["s_ie"]]
  )

  # make the temporal sigma matrix
  sigma_t <- make_sigma(
    de = covparam_object[["t_de"]],
    r_mx = r_t_large,
    ie = covparam_object[["t_ie"]]
  )

  # make the spatio-temporal sigma matrix
  sigma_st <- make_sigma(
    de = covparam_object[["st_de"]],
    r_mx = r_st,
    ie = covparam_object[["st_ie"]]
  )

  # making the overall covariance matrix
  sigma <- sigma_s + sigma_t + sigma_st

  # returning the overall covariance matrix
  return(sigma)
}

make_stcovariance.sum_with_error <- function(covparam_object,
                                             h_s_large,
                                             h_t_large,
                                             s_cor,
                                             t_cor
                                             ) {
  # creating the spatial correlation matrix
  r_s_large <- make_r(
    h = h_s_large,
    range = covparam_object[["s_range"]],
    structure = s_cor
  )

  # creating the temporal correlation matrix
  r_t_large <- make_r(
    h = h_t_large,
    range = covparam_object[["t_range"]],
    structure = t_cor
  )

  # creating the spatio-temporal correlation matrix
  r_st <- r_s_large * r_t_large

  # make the spatial sigma matrix
  sigma_s <- make_sigma(
    de = covparam_object[["s_de"]],
    r_mx = r_s_large,
    ie = covparam_object[["s_ie"]]
  )

  # make the temporal sigma matrix
  sigma_t <- make_sigma(
    de = covparam_object[["t_de"]],
    r_mx = r_t_large,
    ie = covparam_object[["t_ie"]]
  )

  # making the spatio-temporal matrix
  sigma_st <- make_sigma(
    de = 0,
    r_mx = r_st,
    ie = covparam_object[["st_ie"]]
  )

  # make the overall covariance matrix
  sigma <- sigma_s + sigma_t + sigma_st

  # return the overall covariance matrix
  return(sigma)
}



make_stcovariance.product <- function(covparam_object,
                                      h_s_large,
                                      h_t_large,
                                      s_cor,
                                      t_cor) {
  # make the spatial scaled correlation matrix
  # (e = 1 because the variance must be 1 as these are
  # correlation matrices)
  r_s <- make_sigma(
    r_mx = make_r(
      h = h_s_large,
      range = covparam_object[["s_range"]],
      structure = s_cor
    ),
    v_ie = covparam_object[["v_s"]],
    e = 1,
    scale = TRUE
  )

  # make the temporal scaled correlation matrix
  # (e = 1 because the variance must be 1 as these are
  # correlation matrices)
  r_t <- make_sigma(
    r_mx = make_r(
      h = h_t_large,
      range = covparam_object[["t_range"]],
      structure = t_cor
    ),
    v_ie = covparam_object[["v_t"]],
    e = 1,
    scale = TRUE
  )

  # make the spatio-temporal scaled correlation matrix
  r_st <- r_s * r_t

  # make overall covariance matrix (the independent
  # error must be zero for the product covariance)
  sigma_st <- make_sigma(
    de = covparam_object[["st_de"]],
    r_mx = r_st,
    ie = 0
  )

  # making it explicit that sigma is the same as sigma_st
  sigma <- sigma_st

  # returning the overall covariance
  return(sigma)
}
