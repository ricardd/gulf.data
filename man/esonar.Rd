% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/esonar.R
\name{esonar}
\alias{esonar}
\alias{esonar.default}
\title{\strong{eSonar} Data}
\usage{
esonar(x, ...)

\method{esonar}{default}(x, ...)
}
\arguments{
\item{x}{An \code{esonar} object, data file, survey year or keyword search term.}

\item{header}{\strong{eSonar} file header information to be assigned as metadata.}
}
\description{
Create \strong{eSonar} trawl acoustic monitoring probe data objects.
}
\section{Methods (by class)}{
\itemize{
\item \code{esonar(default)}: Create an \code{esonar} object.

}}
\examples{
# eSonar files for the 2020 snow crab survey:
locate.esonar(year = 2020)

# Locate files with a specific tow ID from snow crab survey 2018-2020:
locate.esonar("GP001", year = 2018:2020)

# Working example:
x <- read.esonar("GP001", year = 2020)
describe(x)  # Description of file contents.
header(x)    # File header information.

}
\seealso{
\link{scanmar}, \link{netmind}, \link{notus}
}
