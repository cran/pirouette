## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(pirouette)

## -----------------------------------------------------------------------------
true_phylogeny <- ape::read.tree(text = "(((A:1, B:1):1, C:2):1, D:3);")
ape::plot.phylo(true_phylogeny, main = "The 'true' phylogeny")

## -----------------------------------------------------------------------------
alignment_params <- create_alignment_params(
  root_sequence = create_blocked_dna(length = 20)
)

## -----------------------------------------------------------------------------
experiment <- create_test_gen_experiment()
experiments <- list(experiment)

## -----------------------------------------------------------------------------
pir_params <- create_pir_params(
  alignment_params = alignment_params,
  experiments = experiments
)

## -----------------------------------------------------------------------------
errors <- NULL

if (beastier::is_beast2_installed()) {
  errors <- pir_run(
    phylogeny = true_phylogeny,
    pir_params = pir_params
  )
} else {
  errors <- create_test_pir_run_output(
    add_twin = FALSE,
    add_best = FALSE
  )
}

## -----------------------------------------------------------------------------
knitr::kable(utils::head(errors))

## ----fig.width=7, fig.height=7------------------------------------------------
pir_plot(errors)

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win") {
  experiment_yule <- create_test_cand_experiment(
    inference_model = create_test_inference_model(
      tree_prior = create_yule_tree_prior()
    )
  )
  experiment_bd <- create_test_cand_experiment(
    inference_model = create_test_inference_model(
      tree_prior = create_bd_tree_prior()
    )
  )
  # Use the same files to work on, as only one will actually run an experiment
  experiment_bd$beast2_options <- experiment_yule$beast2_options
  experiment_bd$inference_model$mcmc <- experiment_yule$inference_model$mcmc
  experiment_bd$errors_filename <- experiment_yule$errors_filename

  experiments <- list(experiment_yule, experiment_bd)
  check_experiments(experiments)

  pir_params <- create_pir_params(
    alignment_params = alignment_params,
    experiments = experiments,
    evidence_filename = get_temp_evidence_filename()
  )
}

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" &&
    is_beast2_installed() &&
    is_beast2_ns_pkg_installed()
) {
  errors <- pir_run(
    phylogeny = true_phylogeny,
    pir_params = pir_params
  )
} else {
  errors <- create_test_pir_run_output(add_best = TRUE)
}

## -----------------------------------------------------------------------------
knitr::kable(utils::head(errors))

## ----fig.width=7, fig.height=7------------------------------------------------
pir_plot(errors)

## -----------------------------------------------------------------------------
experiment_yule <- create_test_gen_experiment(
  inference_model = create_test_inference_model(
    tree_prior = create_yule_tree_prior()
  )
)
if (rappdirs::app_dir()$os != "win") {
  experiment_bd <- create_test_cand_experiment(
    inference_model = create_test_inference_model(
      tree_prior = create_bd_tree_prior()
    )
  )
  experiments <- list(experiment_yule, experiment_bd)

  pir_params <- create_pir_params(
    alignment_params = create_test_alignment_params(),
    experiments = experiments,
    evidence_filename = get_temp_evidence_filename()
  )

} else {
  experiments <- list(experiment_yule)

  pir_params <- create_pir_params(
    alignment_params = create_test_alignment_params(),
    experiments = experiments
  )
}


## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" &&
    is_beast2_installed() &&
    is_beast2_ns_pkg_installed()
) {
  errors <- pir_run(
    phylogeny = true_phylogeny,
    pir_params = pir_params
  )
} else {
  errors <- create_test_pir_run_output(add_best = TRUE)
}

## -----------------------------------------------------------------------------
knitr::kable(utils::head(errors))

## ----fig.width=7, fig.height=7------------------------------------------------
pir_plot(errors)

