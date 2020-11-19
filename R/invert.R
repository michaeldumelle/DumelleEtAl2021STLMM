#' Title
#'
#' @param invert_object
#' @import stats
#'
#' @return
#' @export
#'
#' @examples
invert <- function(invert_object) {
  UseMethod("invert", object = invert_object)
}

# invert a product sum covariance matrix
invert.productsum <- function(invert_object) {

  # invert using a the standard cholesky decomposition approach
  if (invert_object$chol) {

    # make the covariance matrix
    sigma <- make_stcovariance.productsum(
      covparam_object = invert_object$covparams,
      h_s_large = invert_object$h_s_large,
      h_t_large = invert_object$h_t_large,
      s_cor = invert_object$s_cor,
      t_cor = invert_object$t_cor
    )

    # adding condition number stability
    diag(sigma) <- diag(sigma) + invert_object$condition

    # finding the Cholesky decomposition
    chol_sigma <- chol(sigma)

    # computing the inverse
    siginv <- chol2inv(chol_sigma)

    # multiplying this inverse on the right
    siginv_o <- siginv %*% invert_object$xyc_o

    # computing the log determinant if requested
    if (invert_object$logdet){
      logdet <- 2 * sum(log(diag(chol_sigma)))
    } else {
      logdet <- NULL
    }
  } else {

    # saving the required objects to invert the matrix - this needs to be cleaned up
    # at some point
    r_s_small <-  make_r(h = invert_object$h_s_small,
                         range = invert_object$covparams[["s_range"]],
                         structure = invert_object$s_cor)
    r_t_small <- make_r(h = invert_object$h_t_small,
                        range = invert_object$covparams[["t_range"]],
                        structure = invert_object$t_cor)
    s_de <- invert_object$covparams[["s_de"]]
    s_ie <- invert_object$covparams[["s_ie"]]
    t_de <- invert_object$covparams[["t_de"]]
    t_ie <- invert_object$covparams[["t_ie"]]
    st_de <- invert_object$covparams[["st_de"]]
    st_ie <- invert_object$covparams[["st_ie"]]
    xyc_o <- invert_object$xyc_o
    condition <- invert_object$condition
    logdet <- invert_object$logdet
    o_index <- invert_object$o_index
    m_index <- invert_object$m_index
    n_s <- invert_object$n_s
    n_t <- invert_object$n_t
    n_st <- n_s * n_t

    dense <- length(o_index) == n_st




    ## adding diagonal tolerances for invertibility stability - for the correlation matrices, they are
    ## also rescaled so the diagonal is 1
    r_s_small <- r_s_small / (1 + condition)
    diag(r_s_small) <- 1
    r_t_small <- r_t_small / (1 + condition)
    diag(r_t_small) <- 1
    s_ie <- s_ie + condition
    t_ie <- t_ie + condition
    st_ie <- st_ie + condition

    ## finding eigendecompositions
    r_s_small_eigen <- eigen(r_s_small)
    r_t_small_eigen <- eigen(r_t_small)

    ## creating w matrix
    w <- kronecker(r_t_small_eigen$vectors, r_s_small_eigen$vectors)
    # creating v matrix
    v <- st_de * kronecker(r_t_small_eigen$values, r_s_small_eigen$values) + st_ie
    # creating inverse square root of v matrix
    vinvroot <- 1 / sqrt(v)
    # creating inverse of w * sqrt(v)
    t_w_vinvroot <- as.vector(vinvroot) * t(w)
    #creating the transpose
    w_vinvroot <- t(t_w_vinvroot)


    # storing the output we will need for the iterative smw
    c_t <- chol(make_sigma(de = t_de, r_mx = r_t_small, ie = t_ie))
    c_s <- chol(make_sigma(de = s_de, r_mx = r_s_small, ie = s_ie))

    ist_zt <- w_vinvroot %*% multiply_z(t_w_vinvroot, "temporal", n_s, n_t, "right")
              # st x st %*% (st x st * st x t) = st x t
    tr_ist_zt <- t(ist_zt)
    c_mt <- chol(chol2inv(c_t) + multiply_z(ist_zt, "temporal", n_s, n_t, "p_left"))
                   #t(multiply_z(mx = t(ist_zt), z_type = "temporal", n_s = n_s)))
              # (t x t + tr(tr(st x t) * st x t) = t x t
    ic_mt <- chol2inv(c_mt)
              # t x t

    istpt_zs <- w_vinvroot %*% multiply_z(mx = t_w_vinvroot, z_type = "spatial", n_s = n_s, n_t = n_t, side = "right") -
      ist_zt %*% (ic_mt %*% multiply_z(mx = tr_ist_zt, z_type = "spatial", n_s = n_s, n_t = n_t, side = "right"))
      # st x st * (st x st * st x s) - st x t (t x t * (tr(st x t) * st x s)) =
      # st x s - st x t * t x s = st x s
    tr_istpt_zs <- t(istpt_zs)
    c_ms <- chol(chol2inv(c_s) + multiply_z(mx = istpt_zs, z_type = "spatial", n_s = n_s, n_t = n_t, side = "p_left"))
                   #t(multiply_z(mx = tr_istpt_zs, z_type = "spatial", n_s = n_s)))
      # s x s + tr(tr(st x s) * st x s) = s x s
    ic_ms <- chol2inv(c_ms)
      # s x s

    # now to implement algorithm

    if (dense) {
      siginv_o <- w_vinvroot %*% (t_w_vinvroot %*% xyc_o) -
      # st x st * (st x st * st x p) = st x p
        ist_zt %*% (ic_mt %*% (tr_ist_zt %*% xyc_o)) -
        # st x t * (t x t * (t x st * st x p)) = st x p
        istpt_zs %*% (ic_ms %*% (tr_istpt_zs %*% xyc_o))
        # st x s * (s x s * (s x st * st x p)) = st x p
    } else {
      d_oo <- w_vinvroot[o_index, o_index, drop = FALSE] %*% (t_w_vinvroot[o_index, o_index, drop = FALSE] %*% xyc_o) +
        w_vinvroot[o_index, m_index, drop = FALSE] %*% (t_w_vinvroot[m_index, o_index, drop = FALSE] %*% xyc_o) -
        ist_zt[o_index, , drop = FALSE] %*%
        (ic_mt %*% (tr_ist_zt[ , o_index, drop = FALSE] %*% xyc_o)) -
        istpt_zs[o_index, , drop = FALSE] %*%
        (ic_ms %*% (tr_istpt_zs[ , o_index, drop = FALSE] %*% xyc_o))

      d_om <- w_vinvroot[o_index, o_index, drop = FALSE] %*% t_w_vinvroot[o_index, m_index, drop = FALSE]  +
        w_vinvroot[o_index, m_index, drop = FALSE] %*% t_w_vinvroot[m_index, m_index, drop = FALSE] -
        ist_zt[o_index, , drop = FALSE] %*%
        (ic_mt %*% tr_ist_zt[ , m_index, drop = FALSE]) -
        istpt_zs[o_index, , drop = FALSE] %*%
        (ic_ms %*% tr_istpt_zs[ , m_index, drop = FALSE])

      d_mm <- w_vinvroot[m_index, o_index, drop = FALSE] %*% t_w_vinvroot[o_index, m_index, drop = FALSE]  +
        w_vinvroot[m_index, m_index, drop = FALSE] %*% t_w_vinvroot[m_index, m_index, drop = FALSE] -
        ist_zt[m_index, , drop = FALSE] %*%
        (ic_mt %*% tr_ist_zt[ , m_index, drop = FALSE]) -
        istpt_zs[m_index, , drop = FALSE] %*%
        (ic_ms %*% tr_istpt_zs[ , m_index, drop = FALSE])

      # return the correct object
      c_mm <- chol(d_mm)
      siginv_o <- d_oo - d_om %*% (chol2inv(c_mm) %*% (t(d_om) %*% xyc_o))
    }

    if (logdet) {
      logdet <- sum(log(v)) +
        2 * sum(log(diag(c_t))) + 2 * sum(log(diag(c_mt))) +
        2 * sum(log(diag(c_s))) + 2 * sum(log(diag(c_ms)))

      if (!dense)
        logdet <- logdet + 2 * sum(log(diag(c_mm)))
    } else {
      logdet <- NULL
    }
  }

  output <- list(sigmainv_o = siginv_o, logdet = logdet)
  output_non_null <- output[!unlist(lapply(output, is.null))]
  return(output_non_null)
}

invert.sum_with_error <- function(invert_object) {

  if (invert_object$chol) {

    # make the covariance matrix
    sigma <- make_stcovariance.sum_with_error(
      covparam_object = invert_object$covparams,
      h_s_large = invert_object$h_s_large,
      h_t_large = invert_object$h_t_large,
      s_cor = invert_object$s_cor,
      t_cor = invert_object$t_cor
    )

    # adding condition number stability
    diag(sigma) <- diag(sigma) + invert_object$condition

    # finding the Cholesky decomposition
    chol_sigma <- chol(sigma)

    # computing the inverse
    siginv <- chol2inv(chol_sigma)

    # multiplying this inverse on the right
    siginv_o <- siginv %*% invert_object$xyc_o

    # computing the log determinant if requested
    if (invert_object$logdet){
      logdet <- 2 * sum(log(diag(chol_sigma)))
    } else {
      logdet <- NULL
    }
  } else {

    # saving the required objects to invert the matrix - this needs to be cleaned up
    # at some point
    r_s_small <-  make_r(h = invert_object$h_s_small,
                         range = invert_object$covparams[["s_range"]],
                         structure = invert_object$s_cor)
    r_t_small <- make_r(h = invert_object$h_t_small,
                        range = invert_object$covparams[["t_range"]],
                        structure = invert_object$t_cor)
    s_de <- invert_object$covparams[["s_de"]]
    s_ie <- invert_object$covparams[["s_ie"]]
    t_de <- invert_object$covparams[["t_de"]]
    t_ie <- invert_object$covparams[["t_ie"]]
    st_ie <- invert_object$covparams[["st_ie"]]
    xyc_o <- invert_object$xyc_o
    condition <- invert_object$condition
    logdet <- invert_object$logdet
    o_index <- invert_object$o_index
    m_index <- invert_object$m_index
    n_s <- invert_object$n_s
    n_t <- invert_object$n_t
    n_st <- n_s * n_t


    dense <- length(o_index) == n_st




    ## adding diagonal tolerances for invertibility stability - for the correlation matrices, they are
    ## also rescaled so the diagonal is 1
    r_s_small <- r_s_small / (1 + condition)
    diag(r_s_small) <- 1
    r_t_small <- r_t_small / (1 + condition)
    diag(r_t_small) <- 1
    s_ie <- s_ie + condition
    t_ie <- t_ie + condition
    st_ie <- st_ie + condition

    # storing the output we will need for the iterative smw
    c_t <- chol(make_sigma(de = t_de, r_mx = r_t_small, ie = t_ie))
    c_s <- chol(make_sigma(de = s_de, r_mx = r_s_small, ie = s_ie))

    c_mt <- chol(chol2inv(c_t) + multiply_z(z_type = "temporal", n_s = n_s, n_t = n_t, side = "pz_z")/st_ie)
    ic_mt <- chol2inv(c_mt)
    istpt <- - multiply_z(multiply_z(ic_mt, z_type = "temporal", n_s = n_s, n_t = n_t, side = "p_right"),
                          z_type = "temporal", n_s = n_s, n_t = n_t, side = "left") / (st_ie^2)
    diag(istpt) <- diag(istpt) + 1/st_ie


    istpt_zs <- multiply_z(mx = istpt, z_type = "spatial", n_s = n_s, n_t = n_t, side = "right")
    c_ms <- chol(chol2inv(c_s) + multiply_z(mx = istpt_zs, z_type = "spatial", n_s = n_s, n_t = n_t, side = "p_left"))
    ic_ms <- chol2inv(c_ms)


    siginv <- istpt - istpt_zs %*% (ic_ms %*% t(istpt_zs))

    if (dense) {
      siginv_o <- siginv %*% xyc_o
    } else {
      c_mm <- chol(siginv[m_index, m_index])
      siginv_o <- (siginv[o_index, o_index] - siginv[o_index, m_index] %*% (chol2inv(c_mm) %*% siginv[m_index, o_index])) %*% xyc_o
    }


    if (logdet) {
      logdet <- n_st * (log(st_ie)) +
        2 * sum(log(diag(c_t))) + 2 * sum(log(diag(c_mt))) +
        2 * sum(log(diag(c_s))) + 2 * sum(log(diag(c_ms)))

      if (!dense)
        logdet <- logdet + 2 * sum(log(diag(c_mm)))
    } else {
      logdet <- NULL
    }
  }

  output <- list(sigmainv_o = siginv_o, logdet = logdet)
  output_non_null <- output[!unlist(lapply(output, is.null))]
  return(output_non_null)
}







invert.product <- function(invert_object) {

  if (invert_object$chol) {

    # make the covariance matrix
    sigma <- make_stcovariance.product(
      covparam_object = invert_object$covparams,
      h_s_large = invert_object$h_s_large,
      h_t_large = invert_object$h_t_large,
      s_cor = invert_object$s_cor,
      t_cor = invert_object$t_cor
    )

    # adding condition number stability
    diag(sigma) <- diag(sigma) + invert_object$condition

    # finding the Cholesky decomposition
    chol_sigma <- chol(sigma)

    # computing the inverse
    siginv <- chol2inv(chol_sigma)

    # multiplying this inverse on the right
    siginv_o <- siginv %*% invert_object$xyc_o

    # computing the log determinant if requested
    if (invert_object$logdet){
      logdet <- 2 * sum(log(diag(chol_sigma)))
    } else {
      logdet <- NULL
    }
  } else {

    # saving the required objects to invert the matrix - this needs to be cleaned up
    # at some point
    r_s_small <-  make_r(h = invert_object$h_s_small,
                         range = invert_object$covparams[["s_range"]],
                         structure = invert_object$s_cor)
    r_t_small <- make_r(h = invert_object$h_t_small,
                        range = invert_object$covparams[["t_range"]],
                        structure = invert_object$t_cor)
    st_de <- invert_object$covparams[["st_de"]]
    v_s <- invert_object$covparams[["v_s"]]
    v_t <- invert_object$covparams[["v_t"]]
    xyc_o <- invert_object$xyc_o
    condition <- invert_object$condition
    logdet <- invert_object$logdet
    o_index <- invert_object$o_index
    m_index <- invert_object$m_index
    n_s <- invert_object$n_s
    n_t <- invert_object$n_t
    n_st <- n_s * n_t


    dense <- length(o_index) == n_st

    r_s_small <- r_s_small / (1 + condition)
    diag(r_s_small) <- 1
    r_t_small <- r_t_small / (1 + condition)
    diag(r_t_small) <- 1

    scale_r_s_small <- make_sigma(r_mx = r_s_small, v_ie = v_s, e = 1, scale = TRUE)
    c_scale_r_s_small  <- chol(scale_r_s_small)
    scale_r_t_small <- make_sigma(r_mx = r_t_small, v_ie = v_t, e = 1, scale = TRUE)
    c_scale_r_t_small  <- chol(scale_r_t_small)

    siginv <- kronecker(chol2inv(c_scale_r_t_small), chol2inv(c_scale_r_s_small)) / st_de

    if (dense) {
      siginv_o <- siginv %*% xyc_o
    } else {
      c_mm <- chol(siginv[m_index, m_index])
      siginv_o <- (siginv[o_index, o_index] - siginv[o_index, m_index] %*% (chol2inv(c_mm) %*% siginv[m_index, o_index])) %*% xyc_o
    }


    if (logdet){
      logdet <- n_st * log(st_de) +
        n_s * 2 * sum(log(diag(c_scale_r_t_small))) +
        n_t * 2 * sum(log(diag(c_scale_r_s_small)))
      if (!dense){
        logdet <- logdet + 2 * sum(log(diag(c_mm)))
      }
    } else {
      logdet <- NULL
    }
  }
  output <- list(sigmainv_o = siginv_o, logdet = logdet)
  output_non_null <- output[!unlist(lapply(output, is.null))]
  return(output_non_null)
}



