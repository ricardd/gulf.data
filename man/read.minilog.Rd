% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.minilog.R, R/read.star.oddi.R
\name{read.minilog}
\alias{read.minilog}
\alias{read.minilog.default}
\alias{read.star.oddi.default}
\title{Read VEMCO Minilog Data}
\usage{
read.minilog(x, ...)

\method{read.minilog}{default}(x, file, offset = 0, echo = FALSE, ...)

\method{read.star.oddi}{default}(x, file, offset = 0, repeats = FALSE, verbose = FALSE, ...)
}
\arguments{
\item{x}{Survey year or file name.}

\item{...}{Other parameters passed onto \code{locate} functions or used to subset data.}

\item{file}{File name(s).}

\item{offset}{Numeric value specifying the offset time (in minutes) to include as a corrective in the data time stamps.}

\item{echo}{Logical value specifying whether to report files being read to the R console.}

\item{year}{Survey year(s).}

\item{survey}{Survey type, as determined by the \link{survey.scsset} function.}

\item{tow.id}{Numeric value or character string specifying the ID corresponding to a particular tow sampling station.}
}
\description{
VEMCO Minilog probes are used to measure water temperature and depth. 
             These functions reads the ASCII formatted Minilog data or header information.
}
\section{Methods (by class)}{
\itemize{
\item \code{read.minilog(default)}: Read a Minilog data file.

}}
\section{Functions}{
\itemize{
\item \code{read.star.oddi(default)}: Read a Star Oddi data file.

}}
\examples{
# Read snow crab survey Minilog data from 1997:
x <- read.minilog(1997)
x <- read.minilog(tow.id = 300)

# Read minilog headers:
files <- locate.minilog(year = 1999, project = "scs")
read.minilog.header(files)
}
