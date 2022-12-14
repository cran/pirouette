% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_std_pir_params.R
\name{create_std_pir_params}
\alias{create_std_pir_params}
\title{Create a standard \code{pir_params}}
\usage{
create_std_pir_params(
  folder_name = rappdirs::user_cache_dir(),
  rng_seed = 314,
  crown_age = 10,
  sequence_length = 1000,
  mutation_rate = 1/crown_age,
  os = rappdirs::app_dir()$os
)
}
\arguments{
\item{folder_name}{name of the main folder}

\item{rng_seed}{a random number generator seed}

\item{crown_age}{the fixed crown age of the posterior. Set to NA
to let it be estimated}

\item{sequence_length}{the length of each DNA sequence in an alignment}

\item{mutation_rate}{the mutation rate per base pair per time unit.
Use \link{check_mutation_rate} to check if a mutation rate is valid.}

\item{os}{name of the operating system, can be \code{mac}, \code{unix}
or \code{win}. Use \link[beastier]{check_os} if the operating system
is valid.}
}
\value{
a `pir_params`, as can be checked by \link{check_pir_params}
}
\description{
Create a standard \code{pir_params}, as used in Bilderbeek, Laudanno
and Etienne.
}
\details{
Create a standard \code{pir_params}, as used in Bilderbeek, Laudanno
and Etienne,
by calling \link{create_pir_params} with these settings:

\itemize{
  \item alignment_params default alignment parameters,
    in which the alignment is simulated using the Jukes-Cantor
    nucleotide substitution model and a strict clock
    (as created by \link{create_alignment_params})
  \item twinning_params default twinning parameters,
    in which the twin tree is simulated using a Birth-Death process
    (using \link{get_sim_bd_twin_tree_fun}),
    the twin alignment is simulated using the Jukes-Cantor
    nucleotide substitution model and a strict clock
    and has an equal amount of nucleotide substitutions as the
    true alignment (using \link{get_sim_twal_same_n_muts_fun})
  \item experiments
    a list of a generative and multiple candidate models.
    The generative model is the default generative model,
    which uses JC69, strict, and Yule,
    as created by \link{create_gen_experiment}.
    The candidate models are all other
    (that is, excluding the generative model)
    birth-death (including Yule) models, which are
    all nucleotide substitution models (see \link{create_site_models}),
    all clock models (see \link{create_clock_models}),
    and the Yule and BD model
    (see \link{create_yule_tree_prior} and \link{create_bd_tree_prior})
  \item error_measure_params the default error measurement parameters,
    (as created by \link{create_error_measure_params})
    which uses the nLTT statistic to determine the difference between
    two (or more) trees
}
}
