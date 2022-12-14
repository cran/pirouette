#' Get the number of taxa of an alignment
#' @inheritParams default_params_doc
#' @return the number of taxa
#' @author Richèl J.C. Bilderbeek
#' @examples
#' get_alignment_n_taxa(
#'   alignment = ape::as.DNAbin(
#'     x = list(species_1 = strsplit("aaaa", split = "")[[1]])
#'   )
#' )
#' @export
get_alignment_n_taxa <- function(alignment, verbose = FALSE) {
  pirouette::check_alignment(alignment)
  nrow(as.matrix(alignment))
}
