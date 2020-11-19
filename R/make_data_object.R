#' Title
#'
#' @param formula
#' @param xcoord
#' @param ycoord
#' @param tcoord
#' @param data
#' @param h_options
#' @import stats
#' @return
#' @export
#'
#' @examples
make_data_object <- function(formula, xcoord, ycoord, tcoord, data, h_options) {

  # setting a default for h_options
  if (is.null(h_options)) {
    h_options = list(h_large = TRUE,
                     h_t_distmetric = "euclidean",
                     h_s_distmetric = "euclidean")
  }

  # restoring the original data
  original_data <- data

  # making the original model frame
  original_stmodel_frame <- model.frame(formula, original_data, na.action = stats::na.omit)

  # creating the fixed design matrix
  original_xo <- model.matrix(formula, original_stmodel_frame)

  # creating the response column
  original_yo <- model.response(original_stmodel_frame)

  # order the data by space within time
  spint <- storder(
    data = original_data,
    xcoord = xcoord,
    ycoord = ycoord,
    tcoord = tcoord,
    h_options = h_options
  )

  # create the model frame using the provided formula
  ordered_stmodel_frame <- model.frame(formula, spint$ordered_data_o,
                                       na.action = stats::na.omit)

  # creating the fixed design matrix
  ordered_xo <- model.matrix(formula, ordered_stmodel_frame)

  # creating the response column
  ordered_yo <- model.response(ordered_stmodel_frame)

  # creating the data object
  data_object <- list(
    formula = formula,
    original_data = original_data,
    original_xo = original_xo,
    original_yo = original_yo,
    ordered_data_dense = spint$ordered_data_dense,
    ordered_data_o = spint$ordered_data_o,
    h_s_small = spint$h_s_small,
    h_t_small = spint$h_t_small,
    n_s = spint$n_s,
    n_t = spint$n_t,
    o_index = spint$o_index,
    m_index = spint$m_index,
    h_s_large = spint$h_s_large,
    h_t_large = spint$h_t_large,
    key_s = spint$key_s,
    key_t = spint$key_t,
    ordered_xo = ordered_xo,
    ordered_yo = ordered_yo,
    h_options = h_options
    )

  # returning the data object
  return(data_object)
}
