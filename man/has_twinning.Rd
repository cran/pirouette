% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/has_twinning.R
\name{has_twinning}
\alias{has_twinning}
\title{Determine if these \code{pir_params} use twinning}
\usage{
has_twinning(pir_params)
}
\arguments{
\item{pir_params}{the parameters of \link[pirouette]{pirouette}.
They are created by \link{create_pir_params}.}
}
\value{
TRUE if the pir_params uses twinnning
}
\description{
Determine if these \code{pir_params} use twinning
}
\examples{
if (beautier::is_on_ci()) {

  pir_params <- create_test_pir_params()
  # Returns FALSE
  has_twinning(pir_params)

  pir_params <- create_test_pir_params(
    twinning_params = create_twinning_params()
  )
  # Returns TRUE
  has_twinning(pir_params)

}
}
