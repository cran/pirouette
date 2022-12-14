% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pir_save.R
\name{pir_save}
\alias{pir_save}
\title{Save all output from \link{pir_run}}
\usage{
pir_save(phylogeny, pir_params, pir_out, folder_name)
}
\arguments{
\item{phylogeny}{a phylogeny of class \link[ape]{phylo}}

\item{pir_params}{the parameters of \link[pirouette]{pirouette}.
They are created by \link{create_pir_params}.}

\item{pir_out}{the output of \link{pir_run}}

\item{folder_name}{name of the main folder}
}
\value{
nothing
}
\description{
Save all output from \link{pir_run}
}
