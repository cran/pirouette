## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(pirouette)

## -----------------------------------------------------------------------------
ext_rates <- seq(0.0, 0.4, length.out = 3)

## -----------------------------------------------------------------------------
spec_rate <- 1.0

## ----fig.width=7, fig.height=7------------------------------------------------
set.seed(42)
crown_age <- 4.0
n_taxa <- 6
phylogenies <- list()
for (i in seq_along(ext_rates)) {
  ext_rate <- ext_rates[i]
  phylogeny <- create_exemplary_dd_tree(
    n_taxa = n_taxa,
    crown_age = crown_age,
    extinction_rate = ext_rate
  )
  phylogenies[[i]] <- phylogeny
}
ape::plot.phylo(phylogenies[[1]])
ape::plot.phylo(phylogenies[[2]])
ape::plot.phylo(phylogenies[[3]])

## -----------------------------------------------------------------------------
alignment_params <- create_alignment_params(
  root_sequence = create_blocked_dna(length = 20)
)
experiment <- create_experiment(
  inference_conditions = create_inference_conditions(
    model_type = "generative",
    run_if = "always",
    do_measure_evidence = FALSE
  ),
  inference_model = beautier::create_inference_model(
    tree_prior = beautier::create_bd_tree_prior(),
    mcmc = beautier::create_test_mcmc()
  )
)
experiments <- list(experiment)
pir_params <- create_pir_params(
  alignment_params = alignment_params,
  experiments = experiments
)
errors <- list()
if (rappdirs::app_dir()$os != "win" && beastier::is_beast2_installed()) {
  for (i in seq_along(ext_rates)) {
    phylogeny <- phylogenies[[i]]
    df <- pirouette::pir_run(
      phylogeny = phylogeny,
      pir_params = pir_params
    )
    errors[[i]] <- df
  }
}

## -----------------------------------------------------------------------------
if (rappdirs::app_dir()$os != "win" && beastier::is_beast2_installed()) {
  first_error_col <- which(colnames(errors[[1]]) == "error_1")
  last_error_col <- ncol(errors[[1]])
  n_errors <- 1 + last_error_col - first_error_col
  df <- data.frame(
    ext_rate = as.factor(rep(ext_rates, each = n_errors)),
    idx = seq(1, n_errors),
    error = NA
  )

  for (i in seq_along(ext_rates)) {
    this_df <- errors[[i]]
    nltts <- this_df[1, first_error_col:last_error_col]
    to_row_index <- 1 + (i * n_errors) - n_errors
    df$error[to_row_index:(to_row_index + n_errors - 1)] <- t(as.numeric(nltts))
  }
}

## ----fig.width=7, fig.height=7------------------------------------------------
if (rappdirs::app_dir()$os != "win" && beastier::is_beast2_installed()) {
  ggplot2::ggplot(
    df, ggplot2::aes(x = ext_rate, y = error)
  ) + ggplot2::geom_boxplot()
  ggplot2::ggplot(
    df, ggplot2::aes(x = ext_rate, y = error)
  ) + ggplot2::geom_violin()
}

