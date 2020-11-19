#' Title
#'
#' @param de
#' @param r_mx
#' @param ie
#' @param v_ie
#' @param e
#' @param scale
#' @import stats
#' @return
#' @export
#'
#' @examples
make_sigma <- function(de, r_mx, ie, v_ie, e = 1, scale = FALSE) {

  # a warning message if the diagonal bug (being slightly greater than 1) appears
  if (any(r_mx > 1)) {
    stop("Bug - Diagonal of correlation matrix different from one")
  }

  # If the scale = FALSE, the standard variance parameters are used
  if (!scale) {

    # create the sigma matrix (equal to the weighted average of
    # dependent correlation and independent correlation)
    sigma <- de * r_mx + ie * (r_mx == 1)

    # return the sigma matrix
    return(sigma)
  }
  # If the scale = TRUE, the scaled variance parameters are used
  else if (scale) {

    # create the proportion of independent error as 1 - the proportion
    # of dependent error
    v_de <- 1 - v_ie

    # create the sigma matrix (equal to the weighted average of
    # dependent correlation and independent correlation)
    sigma <- e * v_de * r_mx + e * v_ie * (r_mx == 1)

    # return the sigma matrix
    return(sigma)
  } else {

    # a warning if scale is not appropriately specified
    stop("scale must be TRUE or FALSE")
  }
}

# could make this a generic based on class of r_mx
