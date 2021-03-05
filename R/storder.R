storder <- function(data, xcoord, ycoord = NULL, tcoord, h_options){

  # find unique temporal coordinates
  key_t <- unique(data[, tcoord, drop = FALSE])

  # order the unique temporal coodrinates
  key_t <- key_t[order(key_t[[tcoord]]), , drop = FALSE]

  # compute the ordered small temporal distance matrix
  h_t_small <- make_h(
    coord1 = key_t[[tcoord]],
    distmetric = h_options$h_t_distmetric
  )

  # record the number of unique temporal observations
  n_t <- nrow(key_t)

  # create an index for the unique temporal locations
  key_t$tindex <- seq.int(1, n_t)

  # find unique spatial coordinates
  key_s <- unique(data[, c(xcoord, ycoord), drop = FALSE])

  # compute the small distance matrix in 1d
  if (is.null(ycoord)) {

    # order the unique spatial coodrinates
    key_s <- key_s[order(key_s[[xcoord]]), , drop = FALSE]

    # compute the ordered small spatial distance matrix
    h_s_small <- make_h(
      coord1 = key_s[[xcoord]],
      distmetric = h_options$h_s_distmetric
    )

  } else {   # compute the small distance matrix in 1d

    # order the unique spatial coodrinates
    key_s <- key_s[order(key_s[[ycoord]], key_s[[xcoord]]), , drop = FALSE]

    # compute the ordered small spatial distance matrix
    h_s_small <- make_h(
      coord1 = key_s[[xcoord]],
      coord2 = key_s[[ycoord]],
      distmetric = h_options$h_s_distmetric
    )
  }

  # record the number of unique spatial observations
  n_s <- nrow(key_s)

  # create an index for the unique temporal locations
  key_s$sindex <- seq.int(1, n_s)

  # merge the data and the spatial data key
  # and then merge that data and the temporal data key
  data <- merge(merge(data, key_s), key_t)

  # create a grid containing every combination of spatial and temporal indices
  full_grid <- expand.grid(sindex = key_s$sindex, tindex = key_t$tindex)

  # create an overall index
  full_grid$index <- seq.int(1, n_t * n_s)

  # merge the grid and data
  data <- merge(full_grid, data, all = TRUE)

  # order the data by the ordered index (by space within time)
  data <- data[order(data$index), , drop = FALSE]

  # create a new variable, observed, which is logical indicating whether
  # the index value was observed in the original data
  data$observed <- !(is.na(data[[tcoord]]) & is.na(data[[xcoord]]))

  # find indicies of the observed data and "missing" data - which
  # are index values having no observation
  o_index <- data$index[data$observed]
  m_index <- data$index[!data$observed]

  # create a subsetted data frame of only the observed values
  ordered_data_o <- data[o_index, , drop = FALSE]

  # compute large distance matrices if needed
  if (h_options$h_large) {

    # compute the large distance matrices in 1d
    if (is.null(ycoord)) {

      # save the times
      hdist_start <- Sys.time()
      # compute the large spatial distance matrix
      h_s_large <- make_h(
        coord1 = ordered_data_o[[xcoord]],
        distmetric = h_options$h_s_distmetric
      )

      # compute the large temporal distance matrix
      h_t_large <- make_h(
        coord1 = ordered_data_o[[tcoord]],
        distmetric = h_options$h_t_distmetric
      )
      hdist_end <- Sys.time()
      hdist_seconds <- as.numeric(hdist_end - hdist_start, units = "secs")

    } else { # compute the large distance matrices in 2d

      # compute the large spatial distance matrix
      # save the times
      hdist_start <- Sys.time()
      h_s_large <- make_h(
        coord1 = ordered_data_o[[xcoord]],
        coord2 = ordered_data_o[[ycoord]],
        distmetric = h_options$h_s_distmetric
      )

      # compute the large temporal distance matrix
      h_t_large <- make_h(
        coord1 = ordered_data_o[[tcoord]],
        distmetric = h_options$h_t_distmetric
      )
      hdist_end <- Sys.time()
      hdist_seconds <- as.numeric(hdist_end - hdist_start, units = "secs")

    }
  } else { # set the large distance matrices equal to NULL if not requested
    h_s_large <- NULL
    h_t_large <- NULL
  }

  #return the relevant output
  return(list(
    ordered_data_dense = data,
    ordered_data_o = ordered_data_o,
    h_s_small = h_s_small,
    h_t_small = h_t_small,
    n_s = n_s,
    n_t = n_t,
    o_index = o_index,
    m_index = m_index,
    hdist_seconds = hdist_seconds,
    h_s_large = h_s_large,
    h_t_large = h_t_large,
    key_s = key_s,
    key_t = key_t)
  )
}
