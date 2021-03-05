plo2r <- function(par, covest_object, ...){
  # calling the appropriate estimation method / stcov type generic (class/sub)
  UseMethod("plo2r", object = covest_object)
}

# the profiled log odds to regular for semivariogram-weighted least squares
plo2r.svwls <- function(par, covest_object){

  # calling the appropriate estimation method generic
  UseMethod("plo2r.svwls", object = covest_object)
}

# the svwls productsum plo2r
plo2r.svwls.productsum <- function(par, covest_object){

  # set reasonable bounds on what the exponentiated plo parameters can be - this was a frustrating bug to find
  par[par > 7] = 7
  par[par < -7] = -7

  # inverse logit the parameter vector
  invlogit <- exp(par) / (1 + exp(par))

  # store lambda, alpha, n_s, n_t, n_st
  lambda <- invlogit[["lambda"]] # proportion of main effect variance
  alpha <- invlogit[["alpha"]] # proportion of spatial main effect variance
  n_s <- invlogit[["n_s"]] # proportion of spatial main efffect ind error
  n_t <- invlogit[["n_t"]] # proportion of temporal main effect ind error
  n_st <- invlogit[["n_st"]] # proportion of interaction ind error

  # transform the plo parameters to regular parameters
  s_de <- lambda * alpha * (1 - n_s)
  s_ie <- lambda * alpha * n_s
  t_de <- lambda * (1 - alpha) * (1 - n_t)
  t_ie <- lambda * (1 - alpha) * n_t
  st_de <- (1 - lambda) * (1 - n_st)
  st_ie <- (1 - lambda) * n_st
  # overall variance capped by max_v
  ov_var <- covest_object$max_options$max_v * invlogit[["var_prop"]]

  # storing the parameter vector
  rparm <- c(ov_var * c(
    s_de = s_de,
    s_ie = s_ie,
    t_de = t_de,
    t_ie = t_ie,
    st_de = st_de,
    st_ie = st_ie
    ),
    # range capped by max_s_range
    s_range = covest_object$max_options$max_s_range * invlogit[["srange_prop"]],
    # range capped by max_t_range
    t_range = covest_object$max_options$max_t_range * invlogit[["trange_prop"]]
  )

  # returning the parameters
  return(rparm)
}


plo2r.svwls.sum_with_error <- function(par, covest_object){

  # set reasonable bounds on what the exponentiated plo parameters can be - this was a frustrating bug to find
  par[par > 7] = 7
  par[par < -7] = -7

  # inverse logit the parameter vector
  invlogit <- exp(par) / (1 + exp(par))

  # store lambda, alpha, n_s, n_t
  lambda <- invlogit[["lambda"]]
  alpha <- invlogit[["alpha"]]
  n_s <- invlogit[["n_s"]]
  n_t <- invlogit[["n_t"]]
  # n_st is automatically 1

  # transform the plo parameters to regular parameters
  s_de <- lambda * alpha * (1 - n_s)
  s_ie <- lambda * alpha * n_s
  t_de <- lambda * (1 - alpha) * (1 - n_t)
  t_ie <- lambda * (1 - alpha) * n_t
  st_ie <- (1 - lambda)

  # overall variance capped by max_v
  ov_var <- covest_object$max_options$max_v * invlogit[["var_prop"]]

  # storing the parameter vector
  rparm <- c(ov_var * c(
    s_de = s_de,
    s_ie = s_ie,
    t_de = t_de,
    t_ie = t_ie,
    st_ie = st_ie
    ),
    # range parameter capped by max_s_range
    s_range = covest_object$max_options$max_s_range * invlogit[["srange_prop"]],
    # range parameter capepd by max_t_range
    t_range = covest_object$max_options$max_t_range * invlogit[["trange_prop"]]
  )
  # returning the parameters
  return(rparm)
}

plo2r.svwls.product <- function(par, covest_object){

  # set reasonable bounds on what the exponentiated plo parameters can be - this was a frustrating bug to find
  par[par > 7] = 7
  par[par < -7] = -7

  # inverse logit the parameter vector
  invlogit <- exp(par) / (1 + exp(par))

  # store v_s and v_t
  v_s <- invlogit[["v_s"]]
  v_t <- invlogit[["v_t"]]

  # overall variance capped by max_v
  ov_var <- covest_object$max_options$max_v * invlogit[["var_prop"]]

  # storing the parameter vector
  rparm <- c(
    v_s = v_s,
    v_t = v_t,
    st_de = ov_var,
    # range parameter capped by max_s_range
    s_range = covest_object$max_options$max_s_range * invlogit[["srange_prop"]],
    # range parameter capped by max_t_range
    t_range = covest_object$max_options$max_t_range * invlogit[["trange_prop"]]
  )

  # return the variance parameters
  return(rparm)
}


# profiled log odds generic for reml
plo2r.reml <- function(par, covest_object, ...){
  UseMethod("plo2r.reml", object = covest_object)
}

plo2r.reml.productsum <- function(par, covest_object, ov_var){

  # set reasonable bounds on what the exponentiated plo parameters can be - this was a frustrating bug to find
  par[par > 7] = 7
  par[par < -7] = -7

  # inverse logit the parameter vector
  invlogit <- exp(par) / (1 + exp(par))

  # store lambda, alpha, n_s, n_t, n_st
  lambda <- invlogit[["lambda"]] # proportion of main effect variance
  alpha <- invlogit[["alpha"]] # proportion of spatial main effect variance
  n_s <- invlogit[["n_s"]] # proportion of spatial main efffect ind error
  n_t <- invlogit[["n_t"]] # proportion of temporal main effect ind error
  n_st <- invlogit[["n_st"]] # proportion of interaction ind error

  # transform the plo parameters to regular parameters
  s_de <- lambda * alpha * (1 - n_s)
  s_ie <- lambda * alpha * n_s
  t_de <- lambda * (1 - alpha) * (1 - n_t)
  t_ie <- lambda * (1 - alpha) * n_t
  st_de <- (1 - lambda) * (1 - n_st)
  st_ie <- (1 - lambda) * n_st

  # storing the parameter vector
  rparm <- c(ov_var * c(
    s_de = s_de,
    s_ie = s_ie,
    t_de = t_de,
    t_ie = t_ie,
    st_de = st_de,
    st_ie = st_ie
  ),
  # range capped by max_s_range
  s_range = covest_object$max_options$max_s_range * invlogit[["srange_prop"]],
  # range capped by max_t_range
  t_range = covest_object$max_options$max_t_range * invlogit[["trange_prop"]]
  )

  # returning the parameters
  return(rparm)
}


plo2r.reml.sum_with_error <- function(par, covest_object, ov_var){

  # set reasonable bounds on what the exponentiated plo parameters can be - this was a frustrating bug to find
  par[par > 7] = 7
  par[par < -7] = -7

  # inverse logit the parameter vector
  invlogit <- exp(par) / (1 + exp(par))

  # store lambda, alpha, n_s, n_t
  lambda <- invlogit[["lambda"]]
  alpha <- invlogit[["alpha"]]
  n_s <- invlogit[["n_s"]]
  n_t <- invlogit[["n_t"]]
  # n_st is automatically 1

  # transform the plo parameters to regular parameters
  s_de <- lambda * alpha * (1 - n_s)
  s_ie <- lambda * alpha * n_s
  t_de <- lambda * (1 - alpha) * (1 - n_t)
  t_ie <- lambda * (1 - alpha) * n_t
  st_ie <- (1 - lambda)

  # storing the parameter vector
  rparm <- c(ov_var * c(
    s_de = s_de,
    s_ie = s_ie,
    t_de = t_de,
    t_ie = t_ie,
    st_ie = st_ie
  ),
  # range parameter capped by max_s_range
  s_range = covest_object$max_options$max_s_range * invlogit[["srange_prop"]],
  # range parameter capepd by max_t_range
  t_range = covest_object$max_options$max_t_range * invlogit[["trange_prop"]]
  )
  # returning the parameters
  return(rparm)
}

plo2r.reml.product <- function(par, covest_object, ov_var){

  # set reasonable bounds on what the exponentiated plo parameters can be - this was a frustrating bug to find
  par[par > 7] = 7
  par[par < -7] = -7

  # inverse logit the parameter vector
  invlogit <- exp(par) / (1 + exp(par))

  # store v_s and v_t
  v_s <- invlogit[["v_s"]]
  v_t <- invlogit[["v_t"]]

  # storing the parameter vector
  rparm <- c(
    v_s = v_s,
    v_t = v_t,
    st_de = ov_var,
    # range parameter capped by max_s_range
    s_range = covest_object$max_options$max_s_range * invlogit[["srange_prop"]],
    # range parameter capped by max_t_range
    t_range = covest_object$max_options$max_t_range * invlogit[["trange_prop"]]
  )

  # return the variance parameters
  return(rparm)
}
