#' Checks if model type is valid
#' Will \link{stop} if not.
#' @inheritParams default_params_doc
#' @return nothing. Will \link{stop} if not
#' @author Giovanni Laudanno, Richèl J.C. Bilderbeek
#' @examples
#' check_model_type(get_model_type_names()[1])
#' check_model_type(get_model_type_names()[2])
#' @export
check_model_type <- function(
  model_type
) {
  model_types <- pirouette::get_model_type_names()
  out <- rep(NA, length(model_types))
  for (l in seq_along(model_types)) {
    out[l] <- paste("'", model_types[l], "'", sep = "")
  }
  if (!(model_type %in% model_types)) {
    stop(
      paste0(
        "'model_type' must be among the following: ",
        paste(out, collapse = ", "), "."
      )
    )
  }
  invisible(model_type)
}
