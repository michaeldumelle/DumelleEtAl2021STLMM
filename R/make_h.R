#' Create a Distance Matrix
#'
#' @param coord1 Coordinate 1.
#' @param coord2 Coordinate 2.
#' @param distmetric Distance metric. Defaults to \code{"euclidean"}.
#'
#' @return A distance matrix
#'
#' @export
make_h <- function(coord1, coord2 = NULL, distmetric = "euclidean") {

  # show the available distance metrics
  distmetric <- match.arg(distmetric)

  # calling the appropriate distance calculation
  switch(
    distmetric,
    euclidean = eucdist(coord1, coord2),
    stop("invalid distance metric")
  )
}



# compute the euclidean distance
eucdist <- function(coord1, coord2 = NULL) {
  if (is.null(coord2)) {
    # euclidean distance if 1d
    eucdist_1d <- sqrt(outer(coord1, coord1, sqr_dif))

    # return the distance
    return(eucdist_1d)
  } else {

    # euclidean distance if 2d
    eucdist_2d <- sqrt(outer(coord1, coord1, sqr_dif) + outer(coord2, coord2, sqr_dif))

    # return the distance
    return(eucdist_2d)
  }
}

# a squared difference function for outer to call
sqr_dif <- function(a, b) {
  (a - b)^2
}
