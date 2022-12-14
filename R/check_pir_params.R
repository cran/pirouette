#' Checks if the argument is a valid \link{pirouette} parameter set.
#'
#' Will \link{stop} if not.
#' A valid \link{pirouette} parameter set
#' can be created by \link{create_pir_params}.
#' @inheritParams default_params_doc
#' @return nothing. Will \link{stop} if not
#' @author Giovanni Laudanno, Richèl J.C. Bilderbeek
#' @examples
#' if (beautier::is_on_ci()) {
#'   check_pir_params(create_test_pir_params())
#' }
#' @export
check_pir_params <- function(
  pir_params
) {
  pirouette::check_pir_params_names(pir_params)
  pirouette::check_pir_params_data_types(pir_params)

  # Cannot call 'will_measure_evidence', as that function calls
  # 'check_pir_params', resulting in infinite recursion
  will_measure_evidence <- FALSE
  for (experiment in pir_params$experiments) {
    if (isTRUE(experiment$inference_conditions$do_measure_evidence)) {
      will_measure_evidence <- TRUE
    }
  }
  evidence_filename <- pir_params$evidence_filename
  if (will_measure_evidence) {
    file_extenstion <- substr(
      basename(evidence_filename),
      nchar(basename(evidence_filename)) - 3,
      nchar(basename(evidence_filename))
    )
    if (file_extenstion != ".csv") {
      stop("'evidence_filename' must be a csv filename")
    }
  }
}

#' Checks if the \code{pir_params} has all the named elements needed
#'
#' Will \link{stop} if not.
#' A valid \link{pirouette} parameter set
#' can be created by \link{create_pir_params}.
#' @inheritParams default_params_doc
#' @return nothing. Will \link{stop} if not
#' @author Richèl J.C. Bilderbeek
#' @export
check_pir_params_names <- function(pir_params) {
  argument_names <- c(
    "alignment_params",
    "twinning_params",
    "experiments",
    "error_measure_params",
    "evidence_filename",
    "verbose"
  )
  for (arg_name in argument_names) {
    if (!arg_name %in% names(pir_params)) {
      stop(
        "'", arg_name, "' must be an element of an 'pir_params'.\n",
        "Tip: use 'create_pir_params'"
      )
    }
  }
}

#' Checks if the \code{pir_params} elements are all of the right data type.
#'
#' Will \link{stop} if not.
#' A valid \link{pirouette} parameter set
#' can be created by \link{create_pir_params}.
#' @inheritParams default_params_doc
#' @return nothing. Will \link{stop} if not
#' @author Richèl J.C. Bilderbeek
#' @export
check_pir_params_data_types <- function(pir_params) {
  tryCatch(
    pirouette::check_alignment_params(pir_params$alignment_params),
    error = function(e) {
      msg <- paste0(
        "'alignment_params' must be a set of alignment parameters.\n",
        "Tip: use 'create_alignment_params'\n",
        "Error message: ", e$message, "\n",
        "Actual value: ", pir_params$alignment_params
      )
      stop(msg)
    }
  )
  # Cannot use the 'has_twinning' function, as this will check
  # the pir_params, leading to infinite recursion
  has_twinning <- !beautier::is_one_na(pir_params$twinning_params)

  if (has_twinning) {
    tryCatch(
      pirouette::check_twinning_params(pir_params$twinning_params),
      error = function(e) {
        msg <- paste0(
          "'twinning_params' must be NA or a set of twinning parameters.\n",
          "Tip: use 'create_twinning_params'\n",
          "Error message: ", e$message, "\n",
          "Actual value: ", pir_params$twinning_params
        )
        stop(msg)
      }
    )
  }
  tryCatch(
    pirouette::check_error_measure_params(pir_params$error_measure_params),
    error = function(e) {
      msg <- paste0(
        "'error_measure_params' must be a set of error measurement ",
        "parameters.\n",
        "Tip: use 'create_error_measure_params'\n",
        "Error message: ", e$message, "\n",
        "Actual value: ", pir_params$error_measure_params
      )
      stop(msg)
    }
  )
  tryCatch(
    pirouette::check_experiments(pir_params$experiments),
    error = function(e) {
      msg <- paste0(
        "'experiments' must be one experiment or a list of one or more ",
        "experiments. \n",
        "Tip: use a list of 'create_experiment'. \n",
        "Error message: ", e$message, " \n",
        "Actual value: ", pir_params$experiments
      )
      stop(msg)
    }
  )
  # Cannot call 'will_measure_evidence', as that function calls
  # 'check_pir_params', resulting in infinite recursion
  will_measure_evidence <- FALSE
  for (experiment in pir_params$experiments) {
    if (isTRUE(experiment$inference_conditions$do_measure_evidence)) {
      will_measure_evidence <- TRUE
    }
  }
  if (will_measure_evidence) {
    if (!is.character(pir_params$evidence_filename)) {
      stop(
        "'evidence_filename' must be a string ",
        "if there is an evidence estimation"
      )
    }
    if (has_twinning) {
      if (!is.character(pir_params$twinning_params$twin_evidence_filename)) {
        stop(
          "'twinning_params$twin_evidence_filename' must be a string ",
          "if there is an evidence estimation"
        )
      }

    }
  } else {
    if (!beautier::is_one_na(pir_params$evidence_filename)) {
      stop(
        "'evidence_filename' must be NA ",
        "if there is no evidence estimation"
      )
    }
    if (has_twinning) {
      if (
        !beautier::is_one_na(pir_params$twinning_params$twin_evidence_filename)
      ) {
        stop(
          "'twinning_params$twin_evidence_filename' must be NA ",
          "if there is no evidence estimation"
        )
      }
    }
  }
  if (!beautier::is_one_bool(pir_params$verbose)) {
    stop("'verbose' must be one boolean")
  }
  invisible(pir_params)
}
