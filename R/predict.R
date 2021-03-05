#' Predict
#'
#' Compute best linear unbiased predictions (Kriging).
#'
#' @param object A model object of class \code{stlmm}.
#'
#' @param newdata A data frame containing columns whose names match the names
#'   of the x-coordinate, y-coordinate, t-coordinate, and predictor variables
#'   in \code{object}.
#'
#' @param interval The interval type. \code{"none"} implies point estimates,
#'   \code{"confidence"} implies point estimates whose standard errors are related
#'   to the mean. \code{"predction"} implies point estimates whose standard errors
#'   are related to a prediction.
#'
#' @param se.fit Should the standard error of the point estimate be returned?
#'   Defaults to \code{TRUE}.
#'
#' @param predcov Should the appropriate full covariance matrix of predictions
#' be returned? Deafults to \code{FALSE}.
#'
#' @param ... Additional arguments
#'
#' @return A list containing the point estimates, standard errors (if requested), and
#'   prediction covariance matrix (if requested).
#'
#' @export
predict.stlmm <- function(object,
                          newdata,
                          interval = c("none", "confidence", "prediction"),
                          se.fit = TRUE,
                          predcov = FALSE,
                          ...) {

  # show the interval types in documentation
  interval <- match.arg(interval)

  # calling the appropriate interval
  pred_output <- switch(
    interval,
    "none" = predict.stlmm_none(
      object = object,
      newdata = newdata,
      ...
    ),
    "confidence" = predict.stlmm_confidence(
      object = object,
      newdata = newdata,
      se.fit = se.fit,
      predcov = predcov,
      ...
    ),
    "prediction" = predict.stlmm_prediction(
      object = object,
      newdata = newdata,
      se.fit = se.fit,
      predcov = predcov,
      ...
    ),
    stop("Must choose confidence or prediction interval")
  )

  # returning the interval argument in the prediction output
  pred_output$interval <- interval

  # returning the prediction output
  return(pred_output)
}


predict.stlmm_none <- function(object,
                               newdata,
                               ...) {

  # creating the model frame
  newdata_stmodel_frame <- model.frame(object$formula[-2], newdata,
    na.action = stats::na.omit
  )
  # creating the fixed design matrix
  # -2 is to remove response from formula
  newdata_xo <- model.matrix(object$formula[-2], newdata_stmodel_frame)

  # computing the estimate
  fit <- newdata_xo %*% object$Coefficients

  # storing the output as a list
  pred_output <- list(fit = fit)

  # returning the output
  return(pred_output)
}

predict.stlmm_confidence <- function(object,
                                     newdata,
                                     se.fit,
                                     predcov,
                                     ...) {
  # creating the model frame
  newdata_stmodel_frame <- model.frame(object$formula[-2], newdata,
    na.action = stats::na.omit
  )

  # creating the fixed design matrix
  # -2 is to remove response from formula
  newdata_xo <- model.matrix(object$formula[-2], newdata_stmodel_frame)

  # computing the estimate
  fit <- newdata_xo %*% object$Coefficients


  if (predcov) {
    # computing the full covariance matrix of the estimates
    predcov <- newdata_xo %*% object$CovCoefficients %*% t(newdata_xo)
    if (se.fit) {
      # computing the standard errors
      se.fit <- sqrt(diag(predcov))
    } else {
      # setting the standard errors to NULL if not requested
      se.fit <- NULL
    }
  } else {
    # setting the full covariance matrix NULL if requested
    predcov <- NULL
    if (se.fit) {
      # the number of predictions
      npred <- nrow(newdata_xo)

      # computing the standard errors
      var.fit <- vapply(
        1:npred,
        function(x) newdata_xo[x, , drop = FALSE] %*% object$CovCoefficients %*% t(newdata_xo[x, , drop = FALSE]),
        double(1)
      )
      se.fit <- sqrt(var.fit)
    } else {
      # setting the standard errors to NULL if not requested
      se.fit <- NULL
    }
  }

  # storing output in a list
  pred_output <- list(fit = fit, se.fit = se.fit, predcov = predcov)

  # removing the non-NULL elements
  pred_output <- pred_output[!vapply(pred_output, function(x) is.null(x), logical(1))]

  # returnign the output
  return(pred_output)
}

predict.stlmm_prediction <- function(object,
                                     newdata,
                                     se.fit,
                                     predcov,
                                     ...) {

  # creating the model frame
  newdata_stmodel_frame <- model.frame(object$formula[-2], newdata,
    na.action = stats::na.omit
  )
  # creating the fixed design matrix
  # -2 is to remove response from formula
  newdata_xo <- model.matrix(object$formula[-2], newdata_stmodel_frame)

  # make spatial distance matrices in 1d
  if (is.null(object$coordnames$ycoord)) {

    # distances between new observations and data
    newdata_data_h_s_large <- make_newdata_h(
      newdata_coord1 = newdata[[object$coordnames$xcoord]],
      data_coord1 = object$coords$xcoord,
      distmetric = object$h_options$h_s_distmetric
    )

    # distances between new obsevations and new observations
    newdata_h_s_large <- make_h(
      coord1 = newdata[[object$coordnames$xcoord]],
      distmetric = object$h_options$h_s_distmetric
    )
  } else { # make spatial distance matrices in 2d

    # distances between new observations and data
    newdata_data_h_s_large <- make_newdata_h(
      newdata_coord1 = newdata[[object$coordnames$xcoord]],
      data_coord1 = object$coords$xcoord,
      newdata_coord2 = newdata[[object$coordnames$ycoord]],
      data_coord2 = object$coords$ycoord,
      distmetric = object$h_options$h_s_distmetric
    )
    # distances between new observations and new observations
    newdata_h_s_large <- make_h(
      coord1 = newdata[[object$coordnames$xcoord]],
      coord2 = newdata[[object$coordnames$ycoord]],
      distmetric = object$h_options$h_s_distmetric
    )
  }

  # make temporal distance matrices
  # distances between new observations and data
  newdata_data_h_t_large <- make_newdata_h(
    newdata_coord1 = newdata[[object$coordnames$tcoord]],
    data_coord1 = object$coords$tcoord,
    distmetric = object$h_options$h_t_distmetric
  )

  # distances between new observations and new observations
  newdata_h_t_large <- make_h(
    coord1 = newdata[[object$coordnames$tcoord]],
    distmetric = object$h_options$h_t_distmetric
  )


  # make covariance matrix between new observations and data
  newdata_data_stcovariance <- make_stcovariance(
    covparam_object = object$CovarianceParameters,
    h_s_large = newdata_data_h_s_large,
    h_t_large = newdata_data_h_t_large,
    s_cor = object$CovarianceForms[["s_cor"]],
    t_cor = object$CovarianceForms[["t_cor"]]
  )


  # make the object to invert and multiply by covariance on the right
  invert_object <- make_invert_object(
    covparam_object = object$CovarianceParameters,
    chol = object$chol,
    condition = object$condition,
    co = t(newdata_data_stcovariance),
    h_s_small = object$data_object$h_s_small,
    h_t_small = object$data_object$h_t_small,
    h_s_large = object$data_object$h_s_large,
    h_t_large = object$data_object$h_t_large,
    logdet = FALSE,
    m_index = object$data_object$m_index,
    o_index = object$data_object$o_index,
    s_cor = object$CovarianceForms[["s_cor"]],
    t_cor = object$CovarianceForms[["t_cor"]],
    xo = NULL,
    yo = NULL
  )

  # compute the inverse
  invert_output <- invert(invert_object)

  # store sigmainverse %*% X
  invxo <- object$invert_output$sigmainv_o[, 1:(ncol(object$invert_output$sigmainv_o) - 1), drop = FALSE]

  # store cov(new, data) %*% sigmainverse %*% X
  newdata_invxo <- newdata_data_stcovariance %*% invxo

  # computing the BLUP estimates
  fit <- newdata_xo %*% object$Coefficients + # beta hat estimates
    newdata_data_stcovariance %*% object$invert_output$sigmainv_o[, ncol(object$invert_output$sigmainv_o), drop = FALSE] -
    newdata_invxo %*% object$Coefficients # blup residuals



  if (predcov) {
    # computing the covariance of new and new
    newdata_stcovariance <- make_stcovariance(
      covparam_object = object$CovarianceParameters,
      h_s_large = newdata_h_s_large,
      h_t_large = newdata_h_t_large,
      s_cor = object$CovarianceForms[["s_cor"]],
      t_cor = object$CovarianceForms[["t_cor"]]
    )

    # storing Xo - cov(new, data) %*% sigmainverse %*% X
    H <- newdata_xo - newdata_invxo

    # computing the prediction covariance matrix
    predcov <- newdata_stcovariance - newdata_data_stcovariance %*% invert_output$sigmainv_o +
      H %*% object$CovCoefficients %*% t(H)

    if (se.fit) {
      # computing the standard errors
      se.fit <- sqrt(diag(predcov))
    } else {
      # setting the standard errors to NULL if not requested
      se.fit <- NULL
    }
  } else {
    # setting the predcov matrix to NULL if not requested
    predcov <- NULL
    if (se.fit) {
      # storing the parameter names
      vparm_names <- c(
        "s_de", "s_ie", "t_de", "t_ie",
        "st_de", "st_ie"
      )

      # computing the overall variance by taking only the variance parameters in the
      # covariance parameter object
      varsum <- sum(object$CovarianceParameters[names(object$CovarianceParameters) %in% vparm_names])

      # computing the standard error
      se.fit <- vapply(
        1:nrow(newdata_xo),
        function(x) {
          H <- newdata_xo[x, , drop = FALSE] - newdata_invxo[x, , drop = FALSE]
          predvar <- varsum - newdata_data_stcovariance[x, , drop = FALSE] %*% invert_output$sigmainv_o[, x, drop = FALSE] +
            H %*% object$CovCoefficients %*% t(H)
          se.fit <- sqrt(predvar)
        },
        double(1)
      )
    } else {
      # setting the standard errors to NULL if not requested
      se.fit <- NULL
    }
  }

  # storing output in a list
  pred_output <- list(fit = fit, se.fit = se.fit, predcov = predcov)

  # returnign the non-NULL elements
  pred_output <- pred_output[!vapply(pred_output, function(x) is.null(x), logical(1))]

  # returning the output
  return(pred_output)
}
