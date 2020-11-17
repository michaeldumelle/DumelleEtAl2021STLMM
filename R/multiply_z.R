multiply_z <- function(mx, z_type, n_s, n_t, side = c("right", "left", "p_right", "p_left", "pz_z", "z_pz")) {

  # show the appropriate multiplication side arguments
  side <- match.arg(side)

  # performing the appropriate multiplication
  switch(side,
         right = multiply_z_r(mx = mx, z_type = z_type, n_s = n_s, n_t = n_t),
         left = multiply_z_l(mx = mx, z_type = z_type, n_s = n_s, n_t = n_t),
         p_right = multiply_zp_r(mx = mx, z_type = z_type, n_s = n_s, n_t = n_t),
         p_left = multiply_zp_l(mx = mx, z_type = z_type, n_s = n_s, n_t = n_t),
         pz_z = multiply_zp_z(z_type = z_type, n_s = n_s, n_t = n_t),
         z_pz = multiply_z_zp(z_type = z_type, n_s = n_s, n_t = n_t),
         stop("Please choose an approporiate multiplication function for Z")
         )
}

# multiplying by z on the right (AZ)
multiply_z_r <- function(mx, z_type, n_s, n_t) {

  # is mx is a 1d vector, turn it into a matrix
  if (is.vector(mx)){
    mx <- matrix(mx, ncol = 1)
  }

  # compute the sample size of a dense data rectangle
  n_st <- n_s * n_t

  # find the row dimension of the matrix A
  n_row <- nrow(mx)

  # multiplying by Z_s
  if (z_type == "spatial") {

    # iterating through the appropriate row sums
    a_zs <- vapply(
      1:n_s,
      function(a) rowSums(mx[, seq(a, n_st, by = n_s)]),
      double(n_row)
    )

    # returning the appropriate product
    return(a_zs)

  } else if (z_type == "temporal") { # multiplying by Z_t

    # iterating through the appropriate row sums
    a_zt <- vapply(
      1:n_t,
      function(a) rowSums(mx[, seq(n_s * (a - 1) + 1, n_s * a, by = 1), drop = FALSE]),
      double(n_row)
    )

    # returning the appropriate product
    return(a_zt)

  } else {
    # returning an error message if no proper z product is chosen
    stop("inappropriate z type")
  }
}

# multiplying by z transpose on the left (Z'A)
# if A is symmetric this equals AZ
multiply_zp_l <- function(mx, z_type, n_s, n_t) {

  # is mx is a 1d vector, turn it into a matrix
  if (is.vector(mx)){
    mx <- matrix(mx, ncol = 1)
  }

  # compute the sample size of a dense data rectangle
  n_st <- n_s * n_t

  # find the column dimension of the matrix A
  n_col <- ncol(mx)

  # multiplying by Z_s
  if (z_type == "spatial") {

    # iterating through the appropriate column sums
    zps_a <- t(vapply(
      1:n_s,
      function(a) colSums(mx[seq(a, n_st, by = n_s), , drop = FALSE]),
      double(n_col))
    )

    # returning the appropriate product
    return(zps_a)

  } else if (z_type == "temporal") { #   # multiplying by Z_s

    # iterating through the appropriate column sums
    zpt_a <- t(vapply(
      1:n_t,
      function(a) colSums(mx[seq(n_s * (a - 1) + 1, n_s * a, by = 1), , drop = FALSE]),
      double(n_col))
    )

    # returning the appropriate product
    return(zpt_a)

  } else {
    # returning an error message if no proper z product is chosen
    stop("innapropriate Z type")
  }
}

# multiplying by z on the left (ZA)
multiply_z_l <- function(mx, z_type, n_s, n_t) {

  # is mx is a 1d vector, turn it into a matrix
  if (is.vector(mx)){
    mx <- matrix(mx, ncol = 1)
  }

  # compute the sample size of a dense data rectangle
  n_st <- n_s * n_t

  # multiplying by Z_s
  if (z_type == "spatial") {

    # computing the appropriate sums
    zs_a <- mx[rep(seq(1, n_s), times = n_t), , drop = FALSE]

    # returning the appropriate product
    return(zs_a)

  } else if (z_type == "temporal") {

    # computing the appropriate sums
    zt_a <- mx[rep(seq(1, n_t), each = n_s), , drop = FALSE]

    # returning the appropriate product
    return(zt_a)

  } else {
    stop("innapropriate Z type")
  }
}

# multiplying by z transpose on the right (AZ')
# if A is symmetric this equals ZA
multiply_zp_r <- function(mx, z_type, n_s, n_t){

  # is mx is a 1d vector, turn it into a matrix
  if (is.vector(mx)){
    mx <- matrix(mx, ncol = 1)
  }

  # compute the sample size of a dense data rectangle
  n_st <- n_s * n_t

  # multiplying by Z_s
  if (z_type == "spatial"){

    # computing the appropriate sums
    a_zps <- mx[ , rep(seq(1, n_s), times = n_t), drop = FALSE]

    # returning the appropriate product
    return(a_zps)


  } else if (z_type == "temporal"){ # multiplying by Z_t

    # computing the appropriate sums
    a_zpt <- mx[ , rep(seq(1, n_t), each = n_s), drop = FALSE]

    # returning the appropriate product
    return(a_zpt)

  } else {
    stop("innapropriate Z type")
  }
}

# multiplying z transpose and z (of the same type)
multiply_zp_z <- function(z_type, n_s, n_t){

  # compute the sample size of a dense data rectangle
  n_st <- n_s * n_t

  # multiplying Z_s
  if (z_type == "spatial"){

    # computing the appropriate sums
    zps_zs <- diag(n_t, nrow = n_s, ncol = n_s)

    # returning the appropriate product
    return(zps_zs)

  } else if (z_type == "temporal") {   # multiplying Z_t

    # computing the appropriate sums
    zpt_zt <- diag(n_s, nrow = n_t, ncol = n_t)

    # returning the appropriate product
    return(zpt_zt)

  } else {
    stop("inappropriate z type")
  }
}

# multiplying z and z transpose (of the same type)
multiply_z_zp <- function(z_type, n_s, n_t){

  # compute the sample size of a dense data rectangle
  n_st <- n_s * n_t

  # multiplying by Z_s
  if (z_type == "spatial"){

    # computing the appropriate sums
    zs_zps <- kronecker(matrix(1, nrow = n_t, ncol = n_t), diag(1, nrow = n_s, ncol = n_s))

    # returning the appropriate product
    return(zs_zps)

  } else if (z_type == "temporal") {

    # computing the appropriate sums
    zt_zpt <- kronecker(diag(1, nrow = n_t, ncol = n_t), matrix(1, nrow = n_s, ncol = n_s))

    # returning the appropriate product
    return(zt_zpt)

  } else {
    stop("inappropriate z type")
  }
}
