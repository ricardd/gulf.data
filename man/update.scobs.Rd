% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update.scobs.R
\name{update.scobs}
\alias{update.scobs}
\title{Update Snow Crab Observer Data}
\usage{
\method{update}{scobs}(
  year,
  zone,
  path = options("gulf.path")[[1]]$snow.crab$observer,
  type = c("sea", "port"),
  source = "oracle",
  Rfile = TRUE,
  csv = TRUE,
  ...
)
}
\arguments{
\item{year}{Year.}

\item{zone}{Fishing zone identifier.}

\item{path}{Data path location.}

\item{type}{Sampling type ('sea' or 'port').}

\item{Rfile}{Logical value specifying whether output to Rdata format.}

\item{csv}{Logical value specifying whether to write to CSV data format.}

\item{...}{Not used.}

\item{password}{Oracle password.}
}
\description{
Update snow crab observer data from Oracle database download.
}
