initial <- function(s_de,
                    s_ie,
                    t_de,
                    t_ie,
                    st_de,
                    st_ie,
                    v_s,
                    v_t,
                    s_range,
                    t_range,
                    estmethod,
                    stcov) {

  # make the initial object covariance parameter object
  # it is the same as make_covparam_object() but forces
  # the user to specify an estimation method
  initialcov <- make_covparam_object(
    s_de = s_de,
    s_ie = s_ie,
    t_de = t_de,
    t_ie = t_ie,
    st_de = st_de,
    st_ie = st_ie,
    v_s = v_s,
    v_t = v_t,
    s_range = s_range,
    t_range = t_range,
    estmethod = estmethod,
    stcov = stcov
  )

  # returning the initial covaiance parameter vector
  return(initialcov)
}


# old way to get initial values that did not seem to work

# get_varinitial <- function(stempsv, stcov){
#   inf_inf <- max(stempsv$gammahat)
#   zeroplus_inf <- stempsv$gammahat[stempsv$h_s_avg == min(stempsv$h_s_avg[stempsv$h_s_avg > 0]) &
#                                      stempsv$h_t_avg == max(stempsv$h_t_avg)]
#   zero_inf <- stempsv$gammahat[stempsv$h_s_avg == 0 & stempsv$h_t_avg == max(stempsv$h_t_avg)]
#   inf_zeroplus <- stempsv$gammahat[stempsv$h_t_avg == min(stempsv$h_t_avg[stempsv$h_t_avg > 0]) &
#                                      stempsv$h_s_avg == max(stempsv$h_s_avg)]
#   inf_zero <- stempsv$gammahat[stempsv$h_t_avg == 0 & stempsv$h_s_avg == max(stempsv$h_s_avg)]
#   zeroplus_zeroplus <- stempsv$gammahat[stempsv$h_s_avg == min(stempsv$h_s_avg[stempsv$h_s_avg > 0]) &
#                                           stempsv$h_t_avg == min(stempsv$h_t_avg[stempsv$h_t_avg > 0])]
#   s_de <- inf_inf - zeroplus_inf
#   s_ie <- zeroplus_inf - zero_inf
#   t_de <- inf_inf - inf_zeroplus
#   t_ie <- inf_zeroplus - inf_zero
#   st_de <- inf_inf - zeroplus_zeroplus - s_de - t_de
#   st_ie <- inf_inf - (s_de + s_ie + t_de + t_ie + st_de)
#   varparams <- pmax(1e-10, c(s_de, s_ie, t_de, t_ie, st_de, st_ie))
#   numzero <- sum(varparams == 1e-10)
#   varparams <- (inf_inf - numzero * 1e-2) * varparams / (sum(varparams))
#   varparams <- pmax(1e-2, varparams)
#   names(varparams) <- c("s_de", "s_ie", "t_de", "t_ie", "st_de", "st_ie")
#   varparams <- rescale_varinitial(varparams = varparams, stcov = stcov)
#   return(varparams)
# }

# rescale_varinitial <- function(varparams, stcov){
#   varparams <- switch(stcov,
#                       "productsum" = rescale_varinitial_productsum(varparams = varparams),
#                       "sum_with_error" = rescale_varinitial_sum_with_error(varparams = varparams),
#                       "product" = rescale_varinitial_product(varparams = varparams),
#                       stop("Not a valid stcov type for varinitial"))
#   return(varparams)
# }
#
# rescale_varinitial_productsum <- function(varparams){
#   varparams <- c(varparams, v_s = NA, v_t = NA)
#   return(varparams)
# }
#
# rescale_varinitial_sum_with_error <- function(varparams){
#   varsum <- sum(varparams)
#   varparams["st_de"] <- 0
#   varparams <- varsum * varparams / sum(varparams)
#   varparams <- c(varparams, v_s = NA, v_t = NA)
#   return(varparams)
# }
#
# rescale_varinitial_product <- function(varparams){
#   v_s <- varparams[["s_ie"]] / (varparams[["s_ie"]] + varparams[["s_de"]])
#   v_t <- varparams[["t_ie"]] / (varparams[["t_ie"]] + varparams[["t_de"]])
#   varsum <- sum(varparams)
#   varparams[which(names(varparams) != "st_de")] <- 0
#   varparams <- varsum * varparams / sum(varparams)
#   varparams <- c(varparams, v_s = v_s, v_t = v_t)
#   return(varparams)
# }
