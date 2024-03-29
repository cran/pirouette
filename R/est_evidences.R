#' Estimate the evidences
#' @inheritParams default_params_doc
#' @return a data frame with evidences. Returns NULL if there
#'   are no experiments that have their evidence measured.
#' @author Richèl J.C. Bilderbeek
#' @examples
#' if (rappdirs::app_dir()$os != "win" &&
#'   beautier::is_on_ci() &&
#'   beastier::is_beast2_installed() &&
#'   is_beast2_ns_pkg_installed()
#' ) {
#'   fasta_filename <- system.file(
#'     "extdata", "alignment.fas", package = "pirouette"
#'   )
#'
#'   # Create a single one candidate experiment
#'   experiments <- list(create_test_cand_experiment())
#'
#'   # Collect the evidences
#'   est_evidences(
#'     fasta_filename = fasta_filename,
#'     experiments = experiments
#'   )
#' }
#' @export
est_evidences <- function(
  fasta_filename,
  experiments,
  evidence_filename,
  verbose = FALSE
) {
  pirouette::check_beast2_installed()
  beautier::check_file_exists(
    fasta_filename,
    filename_description = "fasta_filename"
  )

  # Must use a different FASTA file name, else the evidence estimation
  # will overwrite normal inference files
  #
  # Do not forget: this is only a marginal likelihood estimation, not
  # the actual BEAST2 run, so the files in the experiments (that determine
  # where actual run results are stored) should not be created.
  # Instead, some temporary files that BEAST2 uses in a nested sampling
  # run should be created and deleted afterwards.
  #
  evidence_fasta_filename <- pirouette::to_evidence_filename(fasta_filename)
  # Create the subsubsubfolder for target, do not warn if it already exists
  dir.create(
    dirname(evidence_fasta_filename),
    recursive = TRUE, showWarnings = FALSE
  )
  file.copy(from = fasta_filename, to = evidence_fasta_filename)
  beautier::check_file_exists(
    evidence_fasta_filename,
    "evidence_fasta_filename"
  )
  if (isTRUE(verbose)) {
    message(
      "Copied FASTA file from '", fasta_filename,
      "' to '", evidence_fasta_filename, "'"
    )
  }

  pirouette::check_experiments(experiments)

  # Collect inference models and BEAST2 optionses
  inference_models <- list()
  beast2_optionses <- list()
  i <- 1
  for (experiment in experiments) {
    if (experiment$inference_conditions$do_measure_evidence) {
      testthat::expect_true(
        beautier::is_nested_sampling_mcmc(experiment$est_evidence_mcmc)
      )
      inference_models[[i]] <- experiment$inference_model
      # Overwrite Nested Sampling MCMC
      inference_models[[i]]$mcmc <- experiment$est_evidence_mcmc
      beast2_optionses[[i]] <- experiment$beast2_options
      # Change the filenames, else good BEAST2 runs may be deleted
      beast2_optionses[[i]]$input_filename <- to_evidence_filename(
        beast2_optionses[[i]]$input_filename
      )
      beast2_optionses[[i]]$output_state_filename <- to_evidence_filename(
        beast2_optionses[[i]]$output_state_filename
      )
      # Overwrite BEAST2 bin path
      beast2_optionses[[i]]$beast2_path <- experiment$beast2_bin_path
      i <- i + 1
    }
  }
  testthat::expect_true(length(inference_models) == length(beast2_optionses))
  if (length(inference_models) == 0) {
    return(NULL)
  }

  testthat::expect_true(length(inference_models) > 0)
  beautier::check_inference_models(inference_models)
  beastier::check_beast2_optionses(beast2_optionses)
  pirouette::check_is_ns_beast2_pkg_installed()

  if (verbose) {
    for (i in seq_along(beast2_optionses)) {
      message(
        paste(
          i,
          beast2_optionses[[i]]$input_filename,
          beast2_optionses[[i]]$output_state_filename
        )
      )
    }
  }

  marg_liks <- mcbette::est_marg_liks(
    fasta_filename = evidence_fasta_filename,
    inference_models = inference_models,
    beast2_optionses =  beast2_optionses,
    verbose = verbose
  )
  if (verbose == TRUE) {
    message(marg_liks)
  }

  # Create a sub-sub-folder for the evidence file to be put in,
  # do not warn if it already exists
  dir.create(
    dirname(evidence_filename),
    recursive = TRUE, showWarnings = FALSE
  )
  utils::write.csv(
    x = marg_liks,
    file = evidence_filename
  )
  beautier::check_file_exists(evidence_filename)

  file.remove(evidence_fasta_filename)

  # Delete files
  delete_beast2_state_files(
    beast2_optionses = beast2_optionses,
    verbose = verbose
  )
  sum_marg_liks <- sum(marg_liks$weight)
  tolerance <- 0.1
  testthat::expect_lt(abs(1.0 - sum_marg_liks), tolerance)
  marg_liks
}
