make_newdata_h <- function(newdata_coord1, data_coord1, newdata_coord2 = NULL, data_coord2 = NULL, distmetric = "euclidean") {

  # show the available distance metrics
  distmetric <- match.arg(distmetric)

  # calling the appropriate distance calculation
  switch(
    distmetric,
    euclidean = eucdist_newdata(newdata_coord1, data_coord1, newdata_coord2, data_coord2),
    stop("invalid distance metric")
  )
}

# compute the euclidean distance
eucdist_newdata <- function(newdata_coord1, data_coord1, newdata_coord2, data_coord2) {
  if (is.null(newdata_coord2) & is.null(data_coord2)) {

    # euclidean distance if 1d
    eucdist_1d <- sqrt(outer(newdata_coord1, data_coord1, sqr_dif))

    # return the distance
    return(eucdist_1d)
  } else {

    # euclidean distance if 2d
    eucdist_2d <- sqrt(outer(newdata_coord1, data_coord1, sqr_dif) +
      outer(newdata_coord2, data_coord2, sqr_dif))

    # return the distance
    return(eucdist_2d)
  }
}



# this will be rewritten to merge with make_h
