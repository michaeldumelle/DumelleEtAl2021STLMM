#' Title
#'
#' @param response
#' @param xcoord
#' @param ycoord
#' @param tcoord
#' @param h_options
#' @param h_response
#' @param h_s_large
#' @param h_t_large
#' @param stempsv_options
#' @import stats
#' @importFrom purrr pmap_dfr
#' @return
#' @export
#'
#' @examples
stempsv <- function(response,
                         xcoord,
                         ycoord = NULL,
                         tcoord,
                         h_options = NULL,
                         h_response = NULL,
                         h_s_large = NULL,
                         h_t_large = NULL,
                         stempsv_options = NULL
) {

  # setting default options if none are given
  if (is.null(stempsv_options)) {
    stempsv_options <- list(n_s_lag = 16, n_t_lag = 16, h_s_max = NULL, h_t_max = NULL)
  }

  # setting default h options if none are given
  if (is.null(h_options)) {
    h_options = list(h_t_distmetric = "euclidean", h_s_distmetric = "euclidean")
  }

  # creating a large spatial distance matrix if not provided
  if (is.null(h_s_large)) {
    h_s_large <- make_h(coord1 = xcoord, coord2 = ycoord, distmetric = h_options$h_s_distmetric)
  }

  # creating a large temporal distance matrix if not provided
  if (is.null(h_t_large)){
    h_t_large <- make_h(coord1 = tcoord, distmetric = h_options$h_t_distmetric)
  }

  # creating a large squared response matrix if not provided
  if (is.null(h_response)) {
    h_response <- make_h(coord1 = response, distmetric = "euclidean")^2
  }

  # only storing the upper triangular portions of the spatial, temporal,
  # and squared distance matrices
  h_s_large <- h_s_large[upper.tri(h_s_large, diag = F)]
  h_t_large <- h_t_large[upper.tri(h_t_large, diag = F)]
  h_response <- h_response[upper.tri(h_response, diag = F)]

  # setting a max for the spatial disance if none is provided
  # (as the max distance over 2)
  if (is.null(stempsv_options$h_s_max)) {
    stempsv_options$h_s_max <- max(h_s_large) / 2
  }

  # setting a max for the temporal disance if none is provided
  # (as the max distance over 2)
  if (is.null(stempsv_options$h_t_max)) {
    stempsv_options$h_t_max <- max(h_t_large) / 2
  }

  hmax_index <- (h_s_large <= stempsv_options$h_s_max) & (h_t_large <= stempsv_options$h_t_max)
  h_s_large <- h_s_large[hmax_index]
  h_t_large <- h_t_large[hmax_index]
  h_response <- h_response[hmax_index]

  s_lags <- c(-0.1, seq(0, stempsv_options$h_s_max, length.out = stempsv_options$n_s_lag))
  h_s_index <- cut(h_s_large, s_lags, include.lowest = T)
  t_lags <- c(-0.1, seq(0, stempsv_options$h_t_max, length.out = stempsv_options$n_t_lag))
  h_t_index <- cut(h_t_large, t_lags, include.lowest = T)
  st_index <- interaction(h_s_index, h_t_index)


  h_s_avg <- tapply(h_s_large, st_index, mean)
  h_t_avg <- tapply(h_t_large, st_index, mean)
  gammahat <- tapply(h_response, st_index, FUN = function(x) mean(x) / 2)
  n <- tapply(h_response, st_index, length)

  return_output <- na.omit(data.frame(gammahat, n, h_s_avg, h_t_avg))
  attr(return_output, "na.action") <- NULL
  row.names(return_output) <- NULL
  return(return_output)
}

# stempsv_slow <- function(response,
#                     xcoord,
#                     ycoord = NULL,
#                     tcoord,
#                     h_options = NULL,
#                     h_response = NULL,
#                     h_s_large = NULL,
#                     h_t_large = NULL,
#                     stempsv_options = NULL
#                     ) {
#
#   # setting default options if none are given
#   if (is.null(stempsv_options)) {
#     stempsv_options <- list(n_s_lag = 16, n_t_lag = 16, h_s_max = NULL, h_t_max = NULL)
#   }
#
#   # setting default h options if none are given
#   if (is.null(h_options)) {
#     h_options = list(h_t_distmetric = "euclidean", h_s_distmetric = "euclidean")
#   }
#
#   # creating a large spatial distance matrix if not provided
#   if (is.null(h_s_large)) {
#     h_s_large <- make_h(coord1 = xcoord, coord2 = ycoord, distmetric = h_options$h_s_distmetric)
#   }
#
#   # creating a large temporal distance matrix if not provided
#   if (is.null(h_t_large)){
#     h_t_large <- make_h(coord1 = tcoord, distmetric = h_options$h_t_distmetric)
#   }
#
#   # creating a large squared response matrix if not provided
#   if (is.null(h_response)) {
#     h_response <- make_h(coord1 = response, distmetric = "euclidean")^2
#   }
#
#   # only storing the upper triangular portions of the spatial, temporal,
#   # and squared distance matrices
#   h_s_large <- h_s_large[upper.tri(h_s_large, diag = F)]
#   h_t_large <- h_t_large[upper.tri(h_t_large, diag = F)]
#   h_response <- h_response[upper.tri(h_response, diag = F)]
#
#   # setting a max for the spatial disance if none is provided
#   # (as the max distance over 2)
#   if (is.null(stempsv_options$h_s_max)) {
#     stempsv_options$h_s_max <- max(h_s_large) / 2
#   }
#
#   # setting a max for the temporal disance if none is provided
#   # (as the max distance over 2)
#   if (is.null(stempsv_options$h_t_max)) {
#     stempsv_options$h_t_max <- max(h_t_large) / 2
#   }
#
#   # creating the base upper limits of a vector of spatial lags
#   s_lags_upper <- seq(0, stempsv_options$h_s_max, length.out = stempsv_options$n_s_lag)
#
#   # creating the base lower limits of a vector of spatial lags
#   # the -0.1 is provided so that the (0, 0] endpoint becomes
#   # [0, 0] because the code will be lower < distance <= upper
#   s_lags_lower <- c(-0.1, s_lags_upper[-stempsv_options$n_s_lag])
#
#   # repeating these lags to for each time lag
#   s_lags_upper <- rep(s_lags_upper, times = stempsv_options$n_t_lag)
#   s_lags_lower <- rep(s_lags_lower, times = stempsv_options$n_t_lag)
#
#   # finding the mid point of each lag set
#   h_s_mid <- pmax(0, rowMeans(cbind(s_lags_lower, s_lags_upper)))
#
#   # creating the base upper limits of a vector of temporal lags
#   t_lags_upper <- rep(seq(0, stempsv_options$h_t_max, length.out = stempsv_options$n_t_lag))
#
#   # creating the base lower limits of a vector of temporal lags
#   # the -0.1 is provided so that the (0, 0] endpoint becomes
#   # [0, 0] because the code will be lower < distance <= upper
#   t_lags_lower <- rep(c(-0.1, t_lags_upper[-stempsv_options$n_t_lag]))
#
#   # repeating these lags to for each spatial lag
#   t_lags_upper <- rep(t_lags_upper, each = stempsv_options$n_s_lag)
#   t_lags_lower <- rep(t_lags_lower, each = stempsv_options$n_s_lag)
#
#   # finding the mid point of each lag set
#   h_t_mid <- pmax(0, rowMeans(cbind(t_lags_lower, t_lags_upper)))
#
#   # computing the semivariogram
#   output <- pmap_dfr(
#     list(
#     s_lag_lower = s_lags_lower,
#     s_lag_upper = s_lags_upper,
#     t_lag_lower = t_lags_lower,
#     t_lag_upper = t_lags_upper,
#     h_s_mid = h_s_mid, h_t_mid = h_t_mid
#     ),
#     .f = compute_stempsv, h_s_large, h_t_large, h_response
#   )
#
#   # only returning the lag sets having at least one observation
#   return_output <- output[output[["n"]] > 0, , drop = FALSE]
#
#   # semivariogram returned if there is at least one lag set having an observation
#   if (nrow(return_output) > 0) {
#     return(return_output)
#   } else { # an error if there are no distance classes provided having observations
#     stop("No semivariogram bins meet distance requirements: Choose larger values for h_s_max and h_t_max")
#   }
# }
#
# compute_stempsv <- function(s_lag_lower,
#                             s_lag_upper,
#                             t_lag_lower,
#                             t_lag_upper,
#                             h_s_mid,
#                             h_t_mid,
#                             h_s_large,
#                             h_t_large, h_response
#                             ) {
#
#   # subsetting by the spatial distance in the appropriate lag set
#   h_s <- (h_s_large > s_lag_lower) & (h_s_large <= s_lag_upper)
#
#   # taking the average of the spatial distances in the lag set
#   h_s_avg <- mean(h_s_large[h_s])
#
#   # subsetting by the temporal distance in the appropriate lag set
#   h_t <- (h_t_large > t_lag_lower) & (h_t_large <= t_lag_upper)
#
#   # taking the average of the temporal distances in the lag set
#   h_t_avg <- mean(h_t_large[h_t])
#
#   # taking the appropriate squared differences in the response
#   sqrdifs <- h_response[h_s & h_t]
#
#   # number of unique pairs of squared differences (could have double this
#   # and it would no affect optimization)
#   n_gammahat <- length(sqrdifs)
#
#   # computing the semivariance
#   gammahat <- mean(sqrdifs)/2
#
#   # returning a data frame with relevant information
#   return(data.frame(
#     n = n_gammahat,
#     gammahat = gammahat,
#     h_s_mid = h_s_mid,
#     h_s_avg = h_s_avg,
#     h_t_mid = h_t_mid,
#     h_t_avg = h_t_avg
#     )
#   )
# }
# these two are slightly different in the way that they calculate distances --
# the new method calculates the averages after bin assignment (which is better)
# and the former does it BEFORE -- so it is possible in the new way to have
# the same spatial bins at different temporal bins having different average
# distances when not every spatial location is observed at ever temporal location
