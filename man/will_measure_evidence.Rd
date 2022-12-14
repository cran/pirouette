% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/will_measure_evidence.R
\name{will_measure_evidence}
\alias{will_measure_evidence}
\title{Determine if there is an experiment in which the evidence will
be measured.}
\usage{
will_measure_evidence(pir_params)
}
\arguments{
\item{pir_params}{the parameters of \link[pirouette]{pirouette}.
They are created by \link{create_pir_params}.}
}
\value{
TRUE if yes, FALSE otherwise
}
\description{
Determine if there is an experiment in which the evidence will
be measured.
}
\examples{
if (beautier::is_on_ci()) {

  if (rappdirs::app_dir()$os != "win") {
    pir_params <- create_test_pir_params_setup(has_candidate = TRUE)

    # Yes, the evidence will be measured on Linux with candidate models
    will_measure_evidence(pir_params)
  }

  pir_params <- create_test_pir_params_setup(has_candidate = FALSE)
  # No, the evidence will be measured without candidate models
  will_measure_evidence(pir_params)

}
}
