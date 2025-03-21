% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mishmash.R
\name{flash_lfsr_LF_E}
\alias{flash_lfsr_LF_E}
\title{Get Local False Sign Rate (lfsr) for Flash Model}
\usage{
flash_lfsr_LF_E(fl, fit_E, verbose, nsamp)
}
\arguments{
\item{fl}{A \code{flash} object.}

\item{fit_E}{Boolean indicating if E was fit in the flash model. This is
equivalent to \code{var.type = NULL} in \code{flashier::flash}.}

\item{verbose}{boolean indicating if messages should be printed about lfsr.}

\item{nsamp}{Number of samples to use to calculate lfsr.}
}
\value{
A matrix with the lfsr for each element of the data matrix.
}
\description{
Get Local False Sign Rate (lfsr) for Flash Model
}
