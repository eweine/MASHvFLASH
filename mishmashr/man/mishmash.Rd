% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mishmash.R
\name{mishmash}
\alias{mishmash}
\title{Fit Mash and Flash models}
\usage{
mishmash(mash_data, mash_params = list(), flash_params = list())
}
\arguments{
\item{mash_data}{Output from \code{mashr::mash_set_data}.}

\item{mash_params}{List of parameters to be passed to
\code{mashr::mash}}

\item{flash_params}{List of parameters to be passed to
\code{flashier::flash}.}
}
\value{
A list with components
\itemize{
  \item flash - a flash object (what is returned by \code{flashier::flash})
  \item mash - a list of mash outputs (what is returned by \code{mashr::mash})
}
}
\description{
Fit Mash and Flash models
}
\examples{

# Generate 10 x 100 matrix with rank 3.
n <- 10
p <- 100
k <- 3
A <- matrix(
  data = rnorm(n = n * k, sd = 2), nrow = n, ncol = k
)

cov <- A \%*\% t(A) / (n - 1)

X <- t(MASS::mvrnorm(n = p, mu = rep(0, n), Sigma = cov))
Y <- X + matrix(rnorm(n * p, 0, 1), n, p)

#set mash data
mash_data <- mashr::mash_set_data(Y)

Ulist <- mashr::cov_canonical(mash_data)

#call mishmash with specified parameters for mash and flash
mm <- mishmash(
  mash_data,
  mash_params = list(Ulist = Ulist, algorithm.version = 'R', posterior_samples = 10),
  flash_params = list(greedy.Kmax = 5, ebnm.fn = ebnm::ebnm_unimodal_symmetric)
)

}
