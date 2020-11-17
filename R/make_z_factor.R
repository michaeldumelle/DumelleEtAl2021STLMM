make_z_factor <- function(factor_index) {

  # if factor_index is not a character vector, store it as one
  if (!is.character(factor_index)) {
    factor_index <- as.character(factor_index)
  }

  # make the fixed effect spatial design matrix
  z <- stats::model.matrix(~ factor_index - 1)

  # return the matrix
  return(z)
}

# remove attributes - useful to not carry over
# various attributes from matrix sums
# attributes(z) <- list(dim = dim(z))
