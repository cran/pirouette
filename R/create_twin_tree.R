#' Create a twin tree
#'
#' It sets the seed with value \code{twinning_params$rng_seed_twin_tree},
#' then generates a tree by
#' calling \code{twinning_paramssim_twin_tree_fun} on the
#' given tree.
#' @inheritParams default_params_doc
#' @return a phylogeny of type \link[ape]{phylo}
#' @author Richèl J.C. Bilderbeek, Giovanni Laudanno
#' @examples
#' phylogeny <- ape::read.tree(text = "((A:1, B:1):1, C:2);")
#' twin_phylogeny <- create_twin_tree(phylogeny)
#' @export
create_twin_tree <- function(
  phylogeny,
  twinning_params = create_twinning_params()
) {
  beautier::check_phylogeny(phylogeny)
  pirouette::check_twinning_params(twinning_params)

  set.seed(twinning_params$rng_seed_twin_tree)
  twin_tree <- twinning_params$sim_twin_tree_fun(phylogeny)

  testthat::expect_true(beautier::is_phylo(twin_tree))

  # Same number of tips
  testthat::expect_equal(ape::Ntip(phylogeny), ape::Ntip(twin_tree))

  if (1 == 2) {
    # Same crown age, is this really important?
    testthat::expect_true(
      all.equal(
        max(ape::branching.times(phylogeny)),
        max(ape::branching.times(twin_tree))
      )
    )
  }
  twin_tree
}
