#' Adapter function to simulate the twin alignment
#' using the \code{lns} site model
#' @inheritParams default_params_doc
#' @return an alignment of type \code{DNAbin}
#' @seealso
#' Use \link{sim_twal_with_uns_nsm}
#' to simulate using an unlinked node substitution model.
#' Use \link{sim_tral_with_lns_nsm}
#' to simulate a true alignment.
#' @author Richèl J.C. Bilderbeek
#' @export
sim_twal_with_lns_nsm <- function(
  twin_phylogeny,
  true_alignment = "irrelevant",
  root_sequence,
  subst_matrix = rep(1, 6),
  branch_mutation_rate = 1.0,
  node_mutation_rate = 1.0,
  base_frequencies = rep(0.25, 4),
  node_time = 0.001
) {
  beautier::check_phylogeny(twin_phylogeny)
  pirouette::check_reconstructed_phylogeny(twin_phylogeny)
  pirouette::check_root_sequence(root_sequence)
  nodesub_result <- nodeSub::sim_linked(
    phy = twin_phylogeny,
    Q = subst_matrix,
    rate = branch_mutation_rate,
    node_mut_rate_double = node_mutation_rate,
    bf = base_frequencies,
    rootseq = strsplit(root_sequence, split = "")[[1]],
    l = nchar(root_sequence),
    node_time = node_time
  )
  alignment_phydat <- nodesub_result$alignment
  alignment <- ape::as.DNAbin(alignment_phydat)
  alignment
}
