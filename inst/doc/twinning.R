## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(pirouette)
library(ggplot2)

## -----------------------------------------------------------------------------
phylogeny <- ape::read.tree(text = "((A:1, B:1):1, (C:1, D:1):1);")
ape::plot.phylo(phylogeny, main = "The True Phylogeny")

## -----------------------------------------------------------------------------
twinning_params <- pirouette::create_twinning_params()

## -----------------------------------------------------------------------------
alignment_params <- pirouette::create_alignment_params(
  root_sequence = pirouette::create_blocked_dna(length = 100),
  rng_seed = 314
)

## -----------------------------------------------------------------------------
experiment_1 <- create_test_gen_experiment(
  inference_model = beautier::create_test_inference_model(
    tree_prior = beautier::create_bd_tree_prior()
  )
)
if (rappdirs::app_dir()$os != "win") {
  experiment_2 <- create_test_cand_experiment(
    inference_model = beautier::create_test_inference_model(
      tree_prior = beautier::create_yule_tree_prior()
    )
  )
  experiment_1$inference_conditions$do_measure_evidence <- TRUE

  experiments <- list(experiment_1, experiment_2)

  twinning_params$twin_evidence_filename <- get_temp_evidence_filename()

  pir_params <- create_pir_params(
    alignment_params = alignment_params,
    experiments = experiments,
    twinning_params = twinning_params,
    evidence_filename = get_temp_evidence_filename()
  )

} else {
  experiments <- list(experiment_1)

  pir_params <- create_pir_params(
    alignment_params = alignment_params,
    experiments = experiments,
    twinning_params = twinning_params
  )
}

## -----------------------------------------------------------------------------
df <- NULL
if (rappdirs::app_dir()$os != "win" &&
    beastier::is_beast2_installed() &&
    mauricer::is_beast2_ns_pkg_installed()
) {
  df <- pir_run(
    phylogeny = phylogeny,
    pir_params = pir_params
  )
} else {
  df <- create_test_pir_run_output()
}

## -----------------------------------------------------------------------------
knitr::kable(df)

## ----fig.width=7, fig.height=7------------------------------------------------
pirouette::pir_plot(df)

## -----------------------------------------------------------------------------
library(ggplot2)

## -----------------------------------------------------------------------------
ape::plot.phylo(phylogeny, main = "The True Phylogeny")

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" &&
    beastier::is_beast2_installed() &&
    mauricer::is_beast2_ns_pkg_installed()
) {
  check_file_exists(pir_params$alignment_params$fasta_filename)
  ape::image.DNAbin(
    ape::read.FASTA(file = pir_params$alignment_params$fasta_filename)
  )
}

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" &&
    beastier::is_beast2_installed() &&
    mauricer::is_beast2_ns_pkg_installed()
) {
  # BEAUti offers the '$(tree)' shorthand notation.
  # Here, do what BEAUti does...
  treelog_filename <-
    pir_params$experiments[[1]]$inference_model$mcmc$treelog$filename

  treelog_filename <- gsub(
    x = treelog_filename,
    pattern = "\\$\\(tree\\)",
    replacement = beautier::get_alignment_id(
      pir_params$alignment_params$fasta_filename
    )
  )
  check_file_exists(treelog_filename)
  babette::plot_densitree(tracerer::parse_beast_trees(treelog_filename))
}

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" &&
    beastier::is_beast2_installed() &&
    mauricer::is_beast2_ns_pkg_installed()
) {
  knitr::kable(readr::read_csv(pir_params$evidence_filename))
}

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" &&
    beastier::is_beast2_installed() &&
    mauricer::is_beast2_ns_pkg_installed()
) {
  experiment <- pir_params$experiments[[2]]
  trees_filename <- experiment$inference_model$mcmc$treelog$filename
  check_file_exists(trees_filename)
  babette::plot_densitree(tracerer::parse_beast_trees(trees_filename))
}

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" &&
    beastier::is_beast2_installed() &&
    mauricer::is_beast2_ns_pkg_installed()
) {
    # A tracelog's filename is set to NA by default.
    # Here, do what BEAUti does...
    tracelog_filename <-
      pir_params$experiments[[1]]$inference_model$mcmc$tracelog$filename
    if (is.na(tracelog_filename)) {
      pir_params$experiments[[1]]$inference_model$mcmc$tracelog$filename <-
      paste0(
        beautier::get_alignment_id(pir_params$alignment_params$fasta_filename),
        ".log"
      )
    }

  check_file_exists(
    pir_params$experiments[[1]]$inference_model$mcmc$tracelog$filename
  )
  df <- tracerer::parse_beast_tracelog_file(
    pir_params$experiments[[1]]$inference_model$mcmc$tracelog$filename
  )
  ggplot(data = df, aes(x = Sample, y = likelihood)) + geom_line()
}

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" &&
    beastier::is_beast2_installed() &&
    mauricer::is_beast2_ns_pkg_installed()
) {
  experiment <- pir_params$experiments[[2]]
  log_filename <- experiment$inference_model$mcmc$tracelog$filename
  check_file_exists(log_filename)
  df <- tracerer::parse_beast_tracelog_file(log_filename)
  ggplot(data = df, aes(x = Sample, y = likelihood)) + geom_line()
}

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" &&
    beastier::is_beast2_installed() &&
    mauricer::is_beast2_ns_pkg_installed()
) {
  ape::plot.phylo(
    ape::read.tree(pir_params$twinning_params$twin_tree_filename),
    main = "The Twin Tree"
  )
}

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" &&
    beastier::is_beast2_installed() &&
    mauricer::is_beast2_ns_pkg_installed()
) {
  check_file_exists(pir_params$twinning_params$twin_alignment_filename)
  ape::image.DNAbin(
    ape::read.FASTA(file = pir_params$twinning_params$twin_alignment_filename)
  )
}

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" &&
    beastier::is_beast2_installed() &&
    mauricer::is_beast2_ns_pkg_installed()
) {
  trees_filename <- to_twin_filename(
    pir_params$experiments[[1]]$inference_model$mcmc$treelog$filename
  )
  check_file_exists(trees_filename)
  babette::plot_densitree(
    tracerer::parse_beast_trees(trees_filename)
  )
}

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" &&
    beastier::is_beast2_installed() &&
    mauricer::is_beast2_ns_pkg_installed()
) {
  twin_evidence_filename <- pir_params$twinning_params$twin_evidence_filename
  check_file_exists(twin_evidence_filename)
  knitr::kable(readr::read_csv(twin_evidence_filename))
}

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" &&
    beastier::is_beast2_installed() &&
    mauricer::is_beast2_ns_pkg_installed()
) {
  trees_filename <- to_twin_filename(
    pir_params$experiments[[2]]$inference_model$mcmc$treelog$filename
  )
  check_file_exists(trees_filename)
  babette::plot_densitree(tracerer::parse_beast_trees(trees_filename))
}

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" &&
    beastier::is_beast2_installed() &&
    mauricer::is_beast2_ns_pkg_installed()
) {
  log_filename <- to_twin_filename(
    pir_params$experiments[[2]]$inference_model$mcmc$tracelog$filename
    )
  check_file_exists(log_filename)
  df <- tracerer::parse_beast_tracelog_file(log_filename)
  ggplot(data = df, aes(x = Sample, y = likelihood)) + geom_line()
}

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" &&
    beastier::is_beast2_installed() &&
    mauricer::is_beast2_ns_pkg_installed()
) {
  log_filename <- to_twin_filename(
    pir_params$experiments[[2]]$inference_model$mcmc$tracelog$filename
  )
  check_file_exists(log_filename)
  df <- tracerer::parse_beast_tracelog_file(log_filename)
  ggplot(data = df, aes(x = Sample, y = likelihood)) + geom_line()
}

