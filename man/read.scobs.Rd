% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.scobs.R
\name{read.scobs}
\alias{read.scobs}
\title{Read Observer Data}
\usage{
read.scobs(source = "R", ...)
}
\arguments{
\item{source}{Data source ('r', 'ascii', 'csv', or 'oracle').}

\item{year}{Year.}

\item{file}{File name(s).}

\item{path}{Data file path.}

\item{cfvn}{Canadian Fishing Vessel Number(s).}

\item{type}{Sampling type ('sea' or 'port').}

\item{trip.number}{Trip number identifier.}

\item{database}{Oracle database name.}

\item{username}{Oracle user name.}

\item{password}{Oracle password.}
}
\description{
Function to access and read at-sea and port observer data.
}
