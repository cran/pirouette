% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_inference_conditions.R
\name{create_inference_conditions}
\alias{create_inference_conditions}
\title{Create the inference conditions}
\usage{
create_inference_conditions(
  model_type = "generative",
  run_if = "always",
  do_measure_evidence = FALSE,
  os = rappdirs::app_dir()$os
)
}
\arguments{
\item{model_type}{type of inference model supplied for an experiment.
Possible values:
\itemize{
  \item \code{generative}: the inference model is (or is assumed to be)
    the inference model underlying the phylogeny
  \item \code{candidate}: the inference model is a candidate model,
    that competes with other models for having the most
    evidence (aka highest marginal likelihood)
}}

\item{run_if}{the condition for an experiment's inference model to be run.
Possible values:
\itemize{
  \item \code{always}: always
  \item \code{best_candidate}: if the inference model is the
    candidate model with the most evidence (aka highest marginal
    likelihood)
}}

\item{do_measure_evidence}{boolean to indicate if the
evidence (aka marginal likelihood) of an experiment must be
measured}

\item{os}{name of the operating system, can be \code{mac}, \code{unix}
or \code{win}. Use \link[beastier]{check_os} if the operating system
is valid.}
}
\value{
the inference conditions
}
\description{
Create the parameters to determine how to choose
a model for the inference
}
\examples{
if (beautier::is_on_ci()) {

  # Create the inference conditions parameter set
  if (rappdirs::app_dir()$os != "win") {
    # it does not work on Windows

    # Model type can be 'generative' or 'candidate'
    model_type <- "candidate"
    # Run condition can be 'always' or 'best_candidate'
    run_if <- "best_candidate"
    # Evidence (aka marginal likelihood) can be measured yes or no
    do_measure_evidence <- TRUE
  } else {
    # Model type can be 'generative' or 'candidate'
    model_type <- "generative"
    # Run condition can be 'always' or 'best_candidate'
    run_if <- "always"
    # Evidence (aka marginal likelihood) can be measured yes or no
    do_measure_evidence <- FALSE
  }

  inference_conditions <- create_inference_conditions(
    model_type = model_type,
    run_if = run_if,
    do_measure_evidence = do_measure_evidence
  )

  # Using the inference conditions, create a testing candidate experiment
  experiment <- create_test_cand_experiment(
    inference_conditions = inference_conditions
  )

  evidence_filename <- NA
  if (do_measure_evidence) evidence_filename <- get_temp_evidence_filename()

  # Use the experiment to create the full pirouette parameter set
  pir_params <- create_pir_params(
    alignment_params = create_test_alignment_params(),
    experiments = list(experiment),
    evidence_filename = evidence_filename
  )

  # Run that experiment on a continuous integration service,
  # only when BEAST2 is installed
  if (beautier::is_on_ci() &&
    beastier::is_beast2_installed() &&
    is_beast2_ns_pkg_installed()
  ) {
    pir_out <- pir_run(
      phylogeny = ape::read.tree(text = "((A:2, B:2):1, C:3);"),
      pir_params = pir_params
    )
    pir_plot(pir_out)
  }
}
}
\author{
Giovanni Laudanno, Richèl J.C. Bilderbeek
}
