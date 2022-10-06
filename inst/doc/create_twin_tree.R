## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE------------------------------------------------------------
library(babette)
library(pirouette)
library(phangorn)
library(pryr)

## -----------------------------------------------------------------------------
phylogeny <- ape::read.tree(text = "(((A:1, B:1):1, C:2):1, D:3);")
ape::plot.phylo(phylogeny)

## -----------------------------------------------------------------------------
bd_twin <- create_twin_tree(
  phylogeny = phylogeny,
  twinning_params = create_twinning_params(
    sim_twin_tree_fun = get_sim_bd_twin_tree_fun()
  )
)

## -----------------------------------------------------------------------------
yule_twin <- create_twin_tree(
  phylogeny = phylogeny,
  twinning_params = create_twinning_params(
    sim_twin_tree_fun = create_sim_yule_twin_tree_fun()
  )
)

## -----------------------------------------------------------------------------
copy_twin <- create_twin_tree(
  phylogeny = phylogeny,
  twinning_params = create_twinning_params(
    sim_twin_tree_fun = create_copy_twtr_from_true_fun()
  )
)

## ----fig.width=7, fig.height=7------------------------------------------------
plot_densitree(
  phylos = c(bd_twin, yule_twin, copy_twin),
  col = c("red", "blue", "gray"),
  width = 4,
  consensus = phylogeny,
  alpha = 0.5
)

## ----fig.width=7, fig.height=7------------------------------------------------
plot_densitree(
  phylos = c(phylogeny, yule_twin),
  col = c("red", "green"),
  width = 4,
  consensus = phylogeny,
  alpha = 0.5
)

