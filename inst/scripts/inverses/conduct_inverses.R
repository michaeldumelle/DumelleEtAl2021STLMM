conduct_inverses <- function(n_s, n_t, n_m, n_rep) {



  # loading required functions
  for (f in list.files("C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/STLMM/DumelleEtAl2021STLMM/R", pattern="*.R")) {
    source(paste("C:/Users/mdumelle/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/STLMM/DumelleEtAl2021STLMM/R", f, sep = "/"))
  }

  ## sample sizes

  n_st <- n_s * n_t

  print(n_st)

  ## mean
  mu <- 0

  ## coordinates
  xcoord <- rep(runif(n_s), times = n_t)
  ycoord <- rep(runif(n_s), times = n_t)
  #tcoord <- rep(runif(n_t), each = n_s)
  tcoord <- rep(seq(0, 1, length.out = n_t), each = n_s) * (n_t - 1) + 1


  ## data
  data <- data.frame(
    xcoord = xcoord,
    ycoord = ycoord,
    tcoord = tcoord,
    mu = mu
  )

  data <- data[with(data, order(tcoord, ycoord)), ]

  ## creating correlation forms
  s_cor <- "exponential"
  t_cor <- "tent"

  ## setting covariance parameter values
  s_de <- 1
  s_ie <- 1
  t_de <- 1
  t_ie <- 1
  st_de <- 1
  st_ie <- 1
  s_range <- 1
  t_range <- 1
  total_var <- sum(s_de + s_ie + t_de + t_ie + st_de + st_ie)

  ## setting the true covariance parameter vector
  covparams <- make_covparam_object(
    s_de = s_de,
    s_ie = s_ie,
    t_de = t_de,
    t_ie = t_ie,
    st_de = st_de,
    st_ie = st_ie,
    s_range = s_range,
    t_range = t_range,
    stcov = "productsum"
  )

  data$response <- as.vector(strnorm(
    object = covparams,
    mu = data$mu,
    size = 1,
    xcoord = "xcoord",
    ycoord = "ycoord",
    tcoord = "tcoord",
    s_cor = s_cor,
    t_cor = t_cor,
    error = "normal",
    data = data)
  )


  ## make intial value function
  create_siminitial <- function(s_de, s_ie, t_de, t_ie, st_de, st_ie, s_range, t_range) {


    s_de_initial <- s_de
    s_ie_initial <- s_ie
    t_de_initial <- t_de
    t_ie_initial <- t_ie
    st_de_initial <- st_de
    st_ie_initial <- st_ie
    s_range_initial <- s_range
    t_range_initial <- t_range

    total_var_initial <- sum(s_de_initial + s_ie_initial + t_de_initial + t_ie_initial + st_de_initial + st_ie_initial)

    ps_initial <- make_covparam_object(
      s_de = s_de_initial,
      s_ie = s_ie_initial,
      t_de = t_de_initial,
      t_ie = t_ie_initial,
      st_de = st_de_initial,
      st_ie = st_ie_initial,
      s_range = s_range_initial,
      t_range = t_range_initial,
      stcov = "productsum"
    )


    swe_scale <- total_var_initial / (total_var_initial - st_de)

    if (swe_scale == Inf) {

      swe_initial <- make_covparam_object(
        s_de = total_var_initial / 5,
        s_ie = total_var_initial / 5,
        t_de = total_var_initial / 5,
        t_ie = total_var_initial / 5,
        st_ie = total_var_initial / 5,
        s_range = s_range_initial,
        t_range = t_range_initial,
        stcov = "sum_with_error"
      )

    } else {

      swe_initial <- make_covparam_object(
        s_de = swe_scale * s_de_initial,
        s_ie = swe_scale * s_ie_initial,
        t_de = swe_scale * t_de_initial,
        t_ie = swe_scale * t_ie_initial,
        st_ie = swe_scale * st_ie_initial,
        s_range = s_range_initial,
        t_range = t_range_initial,
        stcov = "sum_with_error",
        estmethod = "reml"
      )
    }

    p_initial <- make_covparam_object(
      st_de = total_var_initial,
      v_s = s_ie_initial / (s_ie_initial + s_de_initial),
      v_t = t_ie_initial / (t_ie_initial + t_de_initial),
      s_range = s_range_initial,
      t_range = t_range_initial,
      stcov = "product"
    )

    siminitial <- list(
      ps_initial = ps_initial,
      swe_initial = swe_initial,
      p_initial = p_initial
    )
    return(siminitial)
  }

  ## construct initial values
  siminitial <- create_siminitial(s_de = s_de,
                                  s_ie = s_ie,
                                  t_de = t_de,
                                  t_ie = t_ie,
                                  st_de = st_de,
                                  st_ie = st_ie,
                                  s_range = s_range,
                                  t_range = t_range
  )

  ## simulating missing values
  data$observed <- sample(c(rep(TRUE, n_st - n_m), rep(FALSE, n_m)), size = n_st, replace = F)

  data_object <- make_data_object(response ~ 1, "xcoord", "ycoord", "tcoord", subset(data, observed), h_options = NULL)

  create_invert_objects <- function(data_object, siminitial) {

    ps_invert_object <- make_invert_object(covparam_object = siminitial$ps_initial,
                                      chol = FALSE,
                                      condition = 1e-4,
                                      h_s_small = data_object$h_s_small,
                                      h_t_small = data_object$h_t_small,
                                      logdet = FALSE,
                                      o_index = data_object$o_index,
                                      m_index = data_object$m_index,
                                      s_cor = "exponential",
                                      t_cor = "tent",
                                      xo = data_object$ordered_xo,
                                      yo = data_object$ordered_yo
                                      )

    swe_invert_object <- make_invert_object(covparam_object = siminitial$swe_initial,
                                           chol = FALSE,
                                           condition = 1e-4,
                                           h_s_small = data_object$h_s_small,
                                           h_t_small = data_object$h_t_small,
                                           logdet = FALSE,
                                           o_index = data_object$o_index,
                                           m_index = data_object$m_index,
                                           s_cor = "exponential",
                                           t_cor = "tent",
                                           xo = data_object$ordered_xo,
                                           yo = data_object$ordered_yo
    )

    p_invert_object <- make_invert_object(covparam_object = siminitial$p_initial,
                                           chol = FALSE,
                                           condition = 1e-4,
                                           h_s_small = data_object$h_s_small,
                                           h_t_small = data_object$h_t_small,
                                           logdet = FALSE,
                                           o_index = data_object$o_index,
                                           m_index = data_object$m_index,
                                           s_cor = "exponential",
                                           t_cor = "tent",
                                           xo = data_object$ordered_xo,
                                           yo = data_object$ordered_yo
    )

    chol_invert_object <- make_invert_object(covparam_object = siminitial$ps_initial,
                                           chol = TRUE,
                                           condition = 1e-4,
                                           h_s_large = data_object$h_s_large,
                                           h_t_large = data_object$h_t_large,
                                           logdet = FALSE,
                                           o_index = data_object$o_index,
                                           m_index = data_object$m_index,
                                           s_cor = "exponential",
                                           t_cor = "tent",
                                           xo = data_object$ordered_xo,
                                           yo = data_object$ordered_yo
    )

    return(list(ps_invert_object = ps_invert_object,
                swe_invert_object = swe_invert_object,
                p_invert_object = p_invert_object,
                chol_invert_object = chol_invert_object
                )
           )
  }

  invert_objects <- create_invert_objects(data_object = data_object, siminitial = siminitial)




  ps_times <- microbenchmark::microbenchmark(productsum = invert.productsum(invert_objects$ps_invert_object),
                                             times = n_rep)

  ps_times$rep <- 1:n_rep

  swe_times <- microbenchmark::microbenchmark(sum_with_error = invert.sum_with_error(invert_objects$swe_invert_object),
                                              times = n_rep)
  swe_times$rep <- 1:n_rep

  p_times <- microbenchmark::microbenchmark(product = invert.product(invert_objects$p_invert_object),
                                            times = n_rep)
  p_times$rep <- 1:n_rep


  chol_times <- microbenchmark::microbenchmark(cholesky = invert.productsum(invert_objects$chol_invert_object),
                                              times = n_rep)
  chol_times$rep <- 1:n_rep

  h_response <- make_h(coord1 = data_object$ordered_yo, distmetric = "euclidean")^2

  # stempsv_times <- microbenchmark::microbenchmark(stempsv = stempsv(h_response = h_response,
  #                                                                   h_s_large = data_object$h_s_large,
  #                                                                   h_t_large = data_object$h_t_large),
  #                                                                   times = n_rep)
  # stempsv_times$rep <- 1:n_rep

  stempsv_fast_times <- microbenchmark::microbenchmark(stempsv_fast = stempsv(h_response = h_response,
                                                                    h_s_large = data_object$h_s_large,
                                                                    h_t_large = data_object$h_t_large),
                                                  times = n_rep)
  stempsv_fast_times$rep <- 1:n_rep

  invert_times <- as.data.frame(rbind(ps_times, swe_times, p_times, chol_times, stempsv_fast_times)) #, stempsv_times))

  invert_times$expr <- as.character(invert_times$expr)
  invert_times$time <- invert_times$time * 1e-9
  invert_times$n_st <- n_s * n_t - 1

  colnames(invert_times) <- c("algorithm", "time", "rep", "n_st")

  return(invert_times)

}
