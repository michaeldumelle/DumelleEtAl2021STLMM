make_r <- function(h, range, structure = c("exponential", "spherical", "gaussian", "tent")) {

  # call the appropriate correlation function
  switch(
    structure,
    exponential = r_exp(h, range),
    spherical = r_sph(h, range),
    gaussian = r_gau(h, range),
    tent = r_tent(h, range),
    stop("Choose a valid covariance structure")
  )
}

# the correlation functions are parameterized in terms of their effective ranges

# exponential correlation
r_exp <- function(h, range) {

  # compute the exponential correlation
  r <- exp(-(3 * (h / range)))

  # return the exponential correlation
  return(r)
}

# spherical correlation
r_sph <- function(h, range) {

  # compute the spherical correlation
  r <- (1 - (3 / 2) * (h / range) + (1 / 2) * (h / range)^3) * (h <= range)

  # return the spherical correlation
  return(r)
}

# gaussian correlation
r_gau <- function(h, range) {

  # compute the gaussian correlation
  r <- exp(-(3 * (h / range)^2))

  # return the gaussian correlation
  return(r)
}

# tent correlation (only valid in 1d)
r_tent <- function(h, range) {

  # compute the tent correlation
  r <- (1 - h / range) * (h <= range)

  # return the tent correlation
  return(r)
}
