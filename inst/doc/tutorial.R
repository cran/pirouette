## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(pirouette)

## ----fig.width=7, fig.height=7------------------------------------------------
set.seed(42)
phylogeny <- create_exemplary_dd_tree(n_taxa = 6)
ape::plot.phylo(main = "The True Phylogeny", phylogeny)

## -----------------------------------------------------------------------------
alignment_params <- pirouette::create_alignment_params(
  root_sequence = create_blocked_dna(length = 20),
  rng_seed = 314
)
print(alignment_params$root_sequence)

## ----fig.width=7, fig.height=7------------------------------------------------
alignment <- create_true_alignment(
  true_phylogeny = phylogeny,
  alignment_params = alignment_params
)
ape::image.DNAbin(
  alignment,
  main = "A Possibly-True Alignment",
  legend = FALSE
)

## -----------------------------------------------------------------------------
experiment <- create_experiment(
  inference_conditions = create_inference_conditions(
    model_type = "generative",
    run_if = "always",
    do_measure_evidence = FALSE
  ),
  inference_model = beautier::create_inference_model(
    tree_prior = beautier::create_bd_tree_prior(),
    mcmc = beautier::create_test_mcmc()
  ),
  beast2_options = beastier::create_beast2_options(rng_seed = 314)
)
experiments <- list(experiment)

## -----------------------------------------------------------------------------
df <- NULL
if (rappdirs::app_dir()$os != "win" && beastier::is_beast2_installed()) {

  df <- pirouette::pir_run(
    phylogeny = phylogeny,
    pir_params = create_pir_params(
      alignment_params = alignment_params,
      experiments = experiments
    )
  )
} else {
  df <- create_test_pir_run_output()
}
knitr::kable(df)

## ----fig.width=7, fig.height=7------------------------------------------------
pirouette::pir_plot(df)

## ----fig.width=7, fig.height=7------------------------------------------------
ape::plot.phylo(main = "The True Phylogeny", phylogeny)

## ----fig.width=7, fig.height=7------------------------------------------------
ape::image.DNAbin(
  alignment,
  main = "A Possibly-True Alignment",
  legend = FALSE
)

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" && beastier::is_beast2_installed()) {
  mcmc <- beautier::create_test_mcmc()
  experiment_yule <- pirouette::create_test_cand_experiment()
  experiment_bd <- pirouette::create_test_cand_experiment(
    inference_model = beautier::create_inference_model(
      tree_prior = beautier::create_bd_tree_prior(),
    )
  )
  # Candidat experiments must produce the same files
  experiment_yule$beast2_options <- experiment_bd$beast2_options
  experiment_yule$inference_model$mcmc <- experiment_bd$inference_model$mcmc
  experiment_yule$errors_filename <- experiment_bd$errors_filename

  experiments <- list(experiment_yule, experiment_bd)
  check_experiments(experiments)

  pir_params <- create_pir_params(
    alignment_params = alignment_params,
    experiments = experiments,
    evidence_filename = get_temp_evidence_filename()
  )
}

## ----fig.width=7, fig.height=7------------------------------------------------
if (rappdirs::app_dir()$os != "win" &&
    beastier::is_beast2_installed() &&
    mauricer::is_beast2_ns_pkg_installed()) {
  df <- pirouette::pir_run(
    phylogeny = phylogeny,
    pir_params = pir_params
  )
  knitr::kable(df)
  pirouette::pir_plot(df)
}

