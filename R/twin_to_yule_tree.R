#' Create a twin tree from a phylogeny using a Yule process
#' @inheritParams default_params_doc
#' @return a twin Yule tree of class \link[ape]{phylo},
#' @author Richèl J.C. Bilderbeek, Giovanni Laudanno
#' @examples
#' phylogeny <- ape::read.tree(text = "((A:2, B:2):1, C:3);")
#' yule_tree <- sim_yule_twin_tree(phylogeny)
#' @seealso
#' Use \link{sim_bd_twin_tree} to simulate a Birth-Death twin tree.
#' Use \link{create_sim_yule_twin_tree_fun} to get a partially
#' evaluated function to use in the \code{twinning_params} (as
#' created by \link{create_twinning_params})
#' @export
sim_yule_twin_tree <- function(
  true_phylogeny,
  method = "random_tree",
  n_replicates = 1e4
) {
  phylogeny <- true_phylogeny
  age <- beautier::get_crown_age(phylogeny)
  phylo_brts <- sort(
    pirouette::convert_tree2brts(phylogeny),
    decreasing = TRUE
  )
  n_tips <- ape::Ntip(phylogeny)
  soc <- 1 + n_tips - length(phylo_brts)
  testthat::expect_true(all(soc == 1 | soc == 2))
  difference <- (log(n_tips) - log(soc)) / age
  mu <- 0
  lambda <- mu + difference

  sink(tempfile())
  yule_pars <- DDD::bd_ML(
    brts = sort(phylo_brts, decreasing = TRUE),
    cond = 1, #conditioning on stem or crown age # nolint
    initparsopt = c(lambda),
    idparsopt = c(1),
    idparsfix = c(2, 3, 4),
    parsfix = c(mu, 0, 0),
    missnumspec = 0,
    tdmodel = 0,
    btorph = 1,
    soc = soc
  )
  sink()

  lambda_yule <- as.numeric(unname(yule_pars[1]))
  mu_yule <- mu
  testthat::expect_true(beautier::is_one_double(lambda_yule))

  # generate bd branching times from the inferred parameters
  yule_brts0 <- pirouette::create_twin_branching_times(
    phylogeny = phylogeny,
    lambda = lambda_yule,
    mu = mu_yule,
    n_replicates = n_replicates,
    method = method
  )

  pirouette::combine_brts_and_topology(
    brts = yule_brts0,
    tree = phylogeny
  )
}
