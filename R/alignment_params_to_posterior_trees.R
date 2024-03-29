#' Creates a posterior of phylogenies from a known alignment.
#'
#' Ignores the jointly-estimated posterior estimates
#' @inheritParams default_params_doc
#' @return a list of phylogenies in the posterior,
#'   as a \link[ape]{multiphylo}
#' @author Richèl J.C. Bilderbeek
#' @examples
#' if (beautier::is_on_ci() && beastier::is_beast2_installed()) {
#'
#'   # Check cleanup by other functions
#'   beastier::check_empty_beaustier_folders()
#'
#'   alignment_params <- create_test_alignment_params()
#'   create_tral_file(
#'     phylogeny = ape::read.tree(text = "((A:1, B:1):1, C:2);"),
#'     alignment_params = alignment_params
#'   )
#'
#'   alignment_params_to_posterior_trees(
#'     alignment_params = alignment_params,
#'     experiment = create_test_experiment()
#'   )
#'
#'   # Cleanup
#'   beastier::remove_beaustier_folders()
#'   beastier::check_empty_beaustier_folders()
#' }
#' @noRd
alignment_params_to_posterior_trees <- function(# nolint indeed a long name
  alignment_params,
  experiment,
  verbose = FALSE
) {
  pirouette::check_alignment_params(alignment_params)
  pirouette::check_experiment(experiment)
  beautier::check_file_exists(alignment_params$fasta_filename)
  testthat::expect_true(
    !beautier::is_nested_sampling_mcmc(experiment$inference_model$mcmc)
  )

  bbt_out <- babette::bbt_run_from_model(
    fasta_filename = alignment_params$fasta_filename,
    inference_model = experiment$inference_model,
    beast2_options = experiment$beast2_options
  )
  if (verbose) {
    message(
      paste0(
        "Saved BEAST2 input file to '",
        experiment$beast2_options$input_filename, "'"
      )
    )
    message(
      paste0(
        "Saved BEAST2 output log file to '",
        experiment$inference_model$mcmc$tracelog$filename, "'"
      )
    )
    message(
      paste0(
        "Saved BEAST2 output trees file to '",
        experiment$inference_model$mcmc$treelog$filename, "'"
      )
    )
    message(
      paste0(
        "Saved BEAST2 output state file to '",
        experiment$beast2_options$output_state_filename, "'"
      )
    )
  }
  beautier::check_file_exists(experiment$beast2_options$input_filename)

  # tracelog must be uninitialized, won't do it here
  testthat::expect_true(
    !beautier::is_one_na(experiment$inference_model$mcmc$tracelog$filename)
  )
  # Well, if you'd really want to, this is how:
  # if (beautier::is_one_na(
  #  experiment$inference_model$mcmc$tracelog$filename
  # )) {
  #   experiment$inference_model$mcmc$tracelog$filename <- paste0(
  #       beautier::get_alignment_id(
  #       alignment_params$fasta_filename
  #     ), ".log"
  #   )
  # }
  #
  # treelog must be initialized, won't do it here
  testthat::expect_true(
    stringr::str_count(
      string = experiment$inference_model$mcmc$treelog$filename,
      pattern = "\\$"
    ) == 0
  )
  # Well, if you'd really want to, this is how:
  # replace '$(tree)' in the filename
  # by 'beautier::get_alignment_id(alignment_params$fasta_filename)

  trees <- c(bbt_out[[grep(x = names(bbt_out), pattern = "trees")]])

  # Check the number of trees
  mcmc <- experiment$inference_model$mcmc
  testthat::expect_true(!beautier::is_nested_sampling_mcmc(mcmc))
  expected_n_trees <- 1 + (mcmc$chain_length / mcmc$treelog$log_every)
  testthat::expect_true(length(trees) == expected_n_trees)

  trees
}
