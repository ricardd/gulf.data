% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scanmar.R
\name{scanmar}
\alias{scanmar}
\alias{scanmar.default}
\title{\strong{Scanmar} Data}
\usage{
scanmar(x, ...)

\method{scanmar}{default}(x, ...)
}
\arguments{
\item{x}{An \code{scanmar} object, data file, survey year or keyword search term.}

\item{header}{\strong{Scanmar} file header information to be assigned as metadata.}
}
\description{
Create \strong{Scanmar} trawl acoustic monitoring probe data objects.
}
\section{Methods (by class)}{
\itemize{
\item \code{scanmar(default)}: Create an \code{scanmar} object.

}}
\examples{
# Scanmar files for the 2020 snow crab survey:
locate.scanmar(year = 1990)

# Locate files with a specific tow ID from snow crab survey 2018-2020:
locate.scanmar(1990, tow.id = 223)

# Working example:
locate.scanmar(1990, tow.id = 223)

}
\seealso{
\link{locate.probe}, \link{read.scanmar}
}
