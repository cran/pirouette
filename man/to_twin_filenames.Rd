% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/to_twin_filenames.R
\name{to_twin_filenames}
\alias{to_twin_filenames}
\title{Convert multiple filenames to their twin equivalent}
\usage{
to_twin_filenames(filenames)
}
\arguments{
\item{filenames}{the paths to files, may be relative or absolute paths}
}
\value{
a character vector with filenames
}
\description{
Convert multiple filenames to their twin equivalent
}
\examples{
filenames <- c("a.csv", "b.xml")
# c("a_twin.csv", "b_twin.xml")
to_twin_filenames(filenames)
}
\author{
Richèl J.C. Bilderbeek
}
