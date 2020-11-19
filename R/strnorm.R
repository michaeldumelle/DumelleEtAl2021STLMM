#' Title
#'
#' @param object
#' @param mu
#' @param size
#' @param condition
#' @param ...
#' @import stats
#' @return
#' @export
#'
#' @examples
strnorm <- function(object, mu, size, condition, ...) {

  # dispatching the appropriate generic
  UseMethod("strnorm", object = object)
}

# simulating a random variable if a covariance matrix is provided
strnorm.matrix <- function(object, mu, size, condition = 1e-4) {

  # record the sample size
  n_st <- nrow(object)

  # return an error if the length of mu does not equal the sample size
  if (length(mu) != 1 & length(mu) != n_st) {
    stop("choose mu as a vector having size nst or a single scalar")
  }

  # taking the lower triangular cholesky decomposition
  # (chol() returns the upper triangular)
  chol_siginv <- t(chol(object))

  # returning the simulated random vector as many times as "size" indicates
  return(
    vapply(
      1:size,
      FUN = function(x) mu + chol_siginv %*% rnorm(n_st),
      double(n_st)
    )
  )
}

# the deafult simulation method
# the user does not have to provide a covariance matrix
strnorm.default <- function(object,
                            mu,
                            size,
                            xcoord,
                            ycoord = NULL,
                            tcoord,
                            data,
                            s_cor,
                            t_cor,
                            chol = FALSE,
                            condition = 1e-4,
                            h_options = NULL
                            ) {
  # setting default h options if none are provided
  if (is.null(h_options)){
    h_options = list(
      h_large = TRUE,
      h_t_distmetric = "euclidean",
      h_s_distmetric = "euclidean"
    )
  }

  # simulation the random vector by cholesky decompositing
  # the full covariance matrix
  if (chol){

    # compute the large spatial distance matrices in 1d
    if (is.null(ycoord)){

      # compute the large spatial distance matrix
      h_s_large <- make_h(
        coord1 = data[[xcoord]],
        distmetric = h_options$h_s_distmetric
      )

      # compute the large temporal distance matrix
      h_t_large <- make_h(
        coord1 = data[[tcoord]],
        distmetric = h_options$h_t_distmetric
      )

    } else {  # compute the large spatial distance matrices in 2d

      # compute the large spatial distance matrix
      h_s_large <- make_h(
        coord1 = data[[xcoord]],
        coord2 = data[[ycoord]],
        distmetric = h_options$h_s_distmetric
      )

      # compute the large temporal distance matrix
      h_t_large <- make_h(
        coord1 = data[[tcoord]],
        distmetric = h_options$h_t_distmetric
      )
    }

    # make the spatio-temporal covariance
    st_covariance <- make_stcovariance(
      covparam_object = object,
      h_s_large = h_s_large,
      h_t_large = h_t_large,
      s_cor = s_cor,
      t_cor = t_cor
    )

    # simulate the random vector by calling strnorm.matrix
    strnorm_sim <- strnorm.matrix(
      object = st_covariance,
      mu = mu,
      size = size,
      condition = condition
    )
  } else {
    # rename the object as covparam_object (should provide a check for this later)
    covparam_object <- object

    # storing an original key of provided spatio-temporal locations
    data$original_key <- seq.int(1, nrow(data))

    # order the data by space within time
    spint <- storder(
      data = data,
      xcoord = xcoord,
      ycoord = ycoord,
      tcoord = tcoord,
      h_options = h_options
    )

    # compute the small spatial correlation matrix
    r_s_small <- make_r(
      h = spint$h_s_small,
      range = covparam_object[["s_range"]],
      structure = s_cor
    )

    # adding a small diagonal constant for inversion stability
    diag(r_s_small) <- diag(r_s_small) + condition

    # compute the small temporal correlation matrix
    r_t_small <- make_r(
      h = spint$h_t_small,
      range = covparam_object[["t_range"]],
      structure = t_cor
    )

    # adding a small diagonal constant for inversion stability
    diag(r_t_small) <- diag(r_t_small) + condition


    # simulating the random vector
    strnorm_sim <- strnorm_small(
      covparam_object = covparam_object,
      mu = mu,
      size = size,
      r_s_small = r_s_small,
      r_t_small = r_t_small
    )

    # removing the spatio-temporal observations not provided
    strnorm_sim <- strnorm_sim[spint$ordered_data_dense$observed, , drop = FALSE ]

    # ordering by the original data
    strnorm_sim <- strnorm_sim[order(spint$ordered_data_o$original_key), , drop = FALSE ]

    # adding the mean
    strnorm_sim <- mu + strnorm_sim

    # return the random vector
    return(strnorm_sim)
  }
}

# the function that actually simulates the random vector in strnorm.deafult
# and object is a covariance parameter object
strnorm_small <- function(covparam_object, mu, size, r_s_small, r_t_small) {

  # calling the appropriate generic
  UseMethod("strnorm_small", object = covparam_object)
}


strnorm_small.productsum <- function(covparam_object, mu, size, r_s_small, r_t_small) {

  # computing the lower spatial cholesky
  chol_r_s_small <- t(chol(r_s_small))

  # computing the lower temporal cholesky
  chol_r_t_small <- t(chol(r_t_small))

  # storing the number of unique spatial locations
  n_s <- nrow(chol_r_s_small)

  # storing the number of unique temporal locations
  n_t <- nrow(chol_r_t_small)

  # storing the number of possible spatio-temporal locations
  n_st <- n_s * n_t

  # multiply the spatial cholesky by zs on the left
  zs_chol_r_s_small <- multiply_z(
    mx = chol_r_s_small,
    z_type = "spatial",
    n_s = n_s,
    n_t = n_t,
    side = "left"
  )

  # multiply the temporal cholesky by zt on the left
  zt_chol_r_t_small <- multiply_z(
    mx = chol_r_t_small,
    z_type = "temporal",
    n_s = n_s,
    n_t = n_t,
    side = "left"
  )

  # computing the lower spatio-temporal cholesky
  chol_r_st <- kronecker(chol_r_t_small, chol_r_s_small)

  # simulate the random error vector
  strnorm_small_sim <- vapply(
    1:size,
    function(x) {
      # spatial dependent error simulation
      s_de_sim <- zs_chol_r_s_small %*% rnorm(n_s, sd = sqrt(covparam_object[["s_de"]]))

      # spatial independent error simulation
      s_ie_sim <- multiply_z(
        mx = rnorm(n_s, sd = sqrt(covparam_object[["s_ie"]])),
        z_type = "spatial",
        n_s = n_s,
        n_t = n_t,
        side = "left"
      )

      # temporal dependent random error simulation
      t_de_sim <- zt_chol_r_t_small %*% rnorm(n_t, sd = sqrt(covparam_object[["t_de"]]))

      # temporal independent error simulation
      t_ie_sim <- multiply_z(
        mx = rnorm(n_t, sd = sqrt(covparam_object[["t_ie"]])),
        z_type = "temporal",
        n_s = n_s,
        n_t = n_t,
        side = "left"
      )

      # spatio-temporal dependent error simulation
      st_de_sim <- chol_r_st %*% rnorm(n_st, sd = sqrt(covparam_object[["st_de"]]))

      # spatio-temporal independent error simulation
      st_ie_sim <- rnorm(n_st, sd = sqrt(covparam_object[["st_ie"]]))

      # simulating the random error vector
      randomerror_sim <- s_de_sim + s_ie_sim + t_de_sim + t_ie_sim + st_de_sim + st_ie_sim
      return(randomerror_sim)
    },
    double(n_st)
  )

  # return the random error vector
  return(strnorm_small_sim)
}

strnorm_small.sum_with_error <- function(covparam_object, mu, size, r_s_small, r_t_small){

  # storing the lower spatial cholesky
  chol_r_s_small <- t(chol(r_s_small))

  # storing the lower temporal cholesky
  chol_r_t_small <- t(chol(r_t_small))

  # storing the number of unique spatial locations
  n_s <- nrow(chol_r_s_small)

  # storing the number of unique temporal locations
  n_t <- nrow(chol_r_t_small)

  # storing the number of possible spatio-temporal locations
  n_st <- n_s * n_t

  # multiply the spatial cholesky by zs on the left
  zs_chol_r_s_small <- multiply_z(
    mx = chol_r_s_small,
    z_type = "spatial",
    n_s = n_s,
    n_t = n_t,
    side = "left"
  )

  # multiply the spatial cholesky by zs on the left
  zt_chol_r_t_small <- multiply_z(
    mx = chol_r_t_small,
    z_type = "temporal",
    n_s = n_s,
    n_t = n_t,
    side = "left"
  )

  # simulate the random error vector
  strnorm_small_sim <- vapply(
    1:size,
    function(x) {

      # spatial dependent error simulation
      s_de_sim <- zs_chol_r_s_small %*% rnorm(n_s, sd = sqrt(covparam_object[["s_de"]]))

      # spatial independent error simulation
      s_ie_sim <- multiply_z(
        mx = rnorm(n_s, sd = sqrt(covparam_object[["s_ie"]])),
        z_type = "spatial",
        n_s = n_s,
        n_t = n_t,
        side = "left"
      )

      # temporal dependent error simulation
      t_de_sim <- zt_chol_r_t_small %*% rnorm(n_t, sd = sqrt(covparam_object[["t_de"]]))

      # temporal independent error simulation
      t_ie_sim <- multiply_z(
        mx = rnorm(n_t, sd = sqrt(covparam_object[["t_ie"]])),
        z_type = "temporal",
        n_s = n_s,
        n_t = n_t,
        side = "left"
      )

      # spatio-temporal independent error simulation
      st_ie_sim <- rnorm(n_st, sd = sqrt(covparam_object[["st_ie"]]))

      # simulating the random error vector
      randomerror_sim <- s_de_sim + s_ie_sim + t_de_sim + t_ie_sim + st_ie_sim
      return(randomerror_sim)
    },
    double(n_st)
  )

  # return the random error vector
  return(strnorm_small_sim)
}


strnorm_small.product <- function(covparam_object, mu, size, r_s_small, r_t_small){

  # store the lower spatial cholesky
  chol_r_s_small <- t(chol(r_s_small))

  # store the lower temporal cholesky
  chol_r_t_small <- t(chol(r_t_small))

  # storing the number of unique spatial locations
  n_s <- nrow(chol_r_s_small)

  # storing the number of unique temporal locations
  n_t <- nrow(chol_r_t_small)

  # storing the number of possible spatio-temporal locations
  n_st <- n_s * n_t

  # store the lower spatio-temporal cholesky
  chol_r_st <- kronecker(chol_r_t_small, chol_r_s_small)

  # simulate the random error vector
  strnorm_small_sim <- vapply(
    1:size,
    function(x) {

      # spatio-temporal dependent error simulation
      st_de_sim <- chol_r_st %*% rnorm(n_st, sd = sqrt(covparam_object[["st_de"]]))

      # simulating the random error vector
      randomerror_sim <- st_de_sim
      return(randomerror_sim)
    },
    double(n_st)
  )

  # return the random error vector
  return(strnorm_small_sim)
}
