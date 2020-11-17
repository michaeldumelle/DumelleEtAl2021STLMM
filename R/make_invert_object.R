make_invert_object <- function(covparam_object,
                               chol,
                               co = NULL,
                               condition,
                               h_s_large = NULL,
                               h_t_large = NULL,
                               h_s_small = NULL,
                               h_t_small = NULL,
                               logdet,
                               m_index = NULL,
                               o_index = NULL,
                               s_cor,
                               t_cor,
                               xo,
                               yo
                               ) {

  # making the right multiplication matrix
  xyc_o <- cbind(xo, yo, co)

  # error message if small distance matrices are not provided and
  # cholesky is not being used
  if (!chol & (is.null(h_s_small) | is.null(h_t_small))){
    stop("If not using Cholesky decomposition, h_s_small and h_t_small must be provided")
  }

  # error message if large distance matrices are not provided and
  # cholesky is being used
  if (chol & (is.null(h_s_large) | is.null(h_t_large))){
    stop("If using Cholesky decomposition, h_s_large and h_t_large must be provided")
  }
    # storing the spatial dimension
    n_s <- nrow(h_s_small)

    # storing the temporal dimension
    n_t <- nrow(h_t_small)

    # creating the invert_object
    invert_object <- structure(list(
        covparams = covparam_object,
        chol = chol,
        condition = condition,
        logdet = logdet,
        h_s_small = h_s_small,
        h_t_small = h_t_small,
        h_s_large = h_s_large,
        h_t_large = h_t_large,
        o_index = o_index,
        m_index = m_index,
        n_s = n_s,
        n_t = n_t,
        s_cor = s_cor,
        t_cor = t_cor,
        xo = xo,
        yo = yo,
        co = co,
        xyc_o = cbind(xo, yo, co)),
      class = class(covparam_object)
    )
  return(invert_object)
}
