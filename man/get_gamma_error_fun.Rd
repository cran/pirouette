% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_gamma_error_fun.R
\name{get_gamma_error_fun}
\alias{get_gamma_error_fun}
\title{Get an error function that uses the difference in gamma statistic.}
\usage{
get_gamma_error_fun()
}
\value{
an error function
}
\description{
Get an error function that uses the difference in gamma statistic.
}
\note{
the gamma statistic can be found in Pybus and Harvey, 2000,
'Testing macro-evolutionary models using incomplete molecular phylogenies.',
Proc. R. Soc. Lond. B, 267, 2267-2272.
}
\examples{
error_fun <- get_gamma_error_fun()

phylogeny <- ape::read.tree(text = "((A:1.5, B:1.5):1.5, C:3.0);")

tree_1 <- ape::read.tree(text = "((A:1.0, B:1.0):2.0, C:3.0);")
tree_2 <- ape::read.tree(text = "((A:2.0, B:2.0):1.0, C:3.0);")
trees <- c(tree_1, tree_2)

error_fun(phylogeny, c(phylogeny))
error_fun(phylogeny, c(tree_1))
error_fun(phylogeny, c(tree_2))
}
\author{
Richèl J.C. Bilderbeek
}
