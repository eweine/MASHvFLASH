#' Sample Posterior Effect Sizes from Flash Model
#'
#' @param flash flash object
#' @param nsamp number of samples
#'
#' @return List of samples
#'
flash_sample_fx <- function(flash, nsamp) {

  LF_samp <- flash$sampler(nsamp)
  fx_samp <- lapply(LF_samp, function(x) {x[[1]] %*% t(x[[2]])})
  return(fx_samp)

}

#' Get posterior mean of E from flash model
#'
#' @param fl A \code{flash} object.
#' @param fit_E Boolean indicating if E was fit in the flash model. This is
#' equivalent to \code{var.type = NULL} in \code{flashier::flash}.
#'
#' @return
#'
flash_get_E_pm <- function(fl, fit_E) {

  if (!fit_E || flashier:::get.n.factors(fl$flash.fit) == 0) {

    return(
      matrix(data = 0, nrow = nrow(fl$flash.fit$Y), ncol = ncol(fl$flash.fit$Y))
    )

  } else {

    resid <- fl$flash.fit$Y - fl$L.pm %*% t(fl$F.pm)
    E_pm <- resid * (1 / fl$flash.fit$given.S2) / (1 / fl$flash.fit$given.S2 + fl$flash.fit$tau)
    return(E_pm)

  }

}

#' Get Local False Sign Rate (lfsr) for Flash Model
#'
#' @param fl A \code{flash} object.
#' @param fit_E Boolean indicating if E was fit in the flash model. This is
#' equivalent to \code{var.type = NULL} in \code{flashier::flash}.
#' @param verbose boolean indicating if messages should be printed about lfsr.
#' @param nsamp Number of samples to use to calculate lfsr.
#'
#' @return A matrix with the lfsr for each element of the data matrix.
#'
flash_lfsr_LF_E <- function(fl, fit_E, verbose, nsamp) {

  if (nsamp == 0) {

    return(NULL)

  }

  n <- nrow(fl$flash.fit$Y)
  p <- ncol(fl$flash.fit$Y)

  fx_samp <- flash_sample_fx(fl, nsamp)

  if (fit_E) {

    if (verbose) {

      cat("Generating Posterior Samples for lfsr...\n")
      pb <- txtProgressBar(min = 0, max = nsamp, initial = 0, style = 3)

    }

    # sample from posterior of E | LF'
    for (i in 1:nsamp) {

      if (verbose) {

        setTxtProgressBar(pb, i)

      }

      resid <- fl$flash.fit$Y - fx_samp[[i]]
      E_mean_mat <- resid * (1 / fl$flash.fit$given.S2) / (1 / fl$flash.fit$given.S2 + fl$flash.fit$tau)
      E_sd_mat <- sqrt(1 / (1 / fl$flash.fit$given.S2 + fl$flash.fit$tau))
      E_samp <- matrix(
        data = rnorm(
          n = n * p, mean = as.vector(E_mean_mat), sd = as.vector(E_sd_mat)
        ),
        nrow = n,
        ncol = p
      )

      fx_samp[[i]] <- fx_samp[[i]] + E_samp

    }
    if(verbose) {

      close(pb)

    }

  }

  pos_samp <- lapply(fx_samp, function(x) {x > 0})
  neg_samp <- lapply(fx_samp, function(x) {x < 0})

  num_pos <- Reduce('+', pos_samp)
  num_neg <- Reduce('+', neg_samp)
  num_zero <- matrix(data = nsamp, nrow = n, ncol = p) - num_pos - num_neg

  lfsr_neg <- (num_pos + num_zero) / nsamp
  lfsr_pos <- (num_neg + num_zero) / nsamp

  lfsr <- pmin(lfsr_neg, lfsr_pos)
  return(lfsr)

}

#' Run Flash with mash data object
#'
#' @param mash_data Output from \code{mashr::mash_set_data}.
#' @param ... Additional arguments to be passed to \code{flashier::flash}.
#' @param lfsr_nsamp number of samples for computation of lfsr. Set to 0
#' to skip computation of lfsr.
#'
#' @return A \code{flash} object. For more details see the documentation for
#' \code{flashier::flash}.
#' @export
#'
#' @examples
#'
#' # Generate 10 x 500 matrix with rank 3.
#' n <- 10
#' p <- 500
#' k <- 3
#' A <- matrix(
#'    data = rnorm(n = n * k, sd = 2), nrow = n, ncol = k
#' )
#'
#' cov <- A %*% t(A) / (n - 1)
#'
#' X <- t(MASS::mvrnorm(n = p, mu = rep(0, n), Sigma = cov))
#' Y <- X + matrix(rnorm(n * p, 0, 1), n, p)
#'
#' # set mash data
#' mash_data <- mashr::mash_set_data(Y)
#'
#' # call flash
#' fl <- mashy_flash(mash_data, backfit = TRUE)
#'
mashy_flash <- function(mash_data, lfsr_nsamp = 1000, ...) {

  fl <- flashier::flash(
    data = mash_data$Bhat,
    S = mash_data$Shat,
    ...
  )

  optional_args <- list(...)
  fit_E <- TRUE
  if ('var.type' %in% names(optional_args)) {

    if(is.null(optional_args$var.type)) {

      fit_E <- FALSE

    }

  }

  verbose <- TRUE
  if ('verbose' %in% names(optional_args)) {

    if(optional_args$verbose == 0) {

      verbose <- FALSE

    }

  }

  lfsr <- flash_lfsr_LF_E(fl, fit_E, verbose, lfsr_nsamp)

  if (!is.null(lfsr)) {

    fl[['LF_E.lfsr']] <-lfsr

  }

  fl[['E.pm']] <- flash_get_E_pm(fl, fit_E)
  return(fl)

}

#' Get Posterior Mean of Flash Model
#'
#' Computes posterior mean of LF' + E. Note that this will in general
#' give different results than \code{flashier::flash}, as that function doesn't
#' add the posterior mean of E.
#'
#' @param fl flash object
#'
#' @return matrix with posterior mean of each element of data matrix
#' @export
#'
#' @examples
#'
#' # Generate 10 x 500 matrix with rank 3.
#' n <- 10
#' p <- 500
#' k <- 3
#' A <- matrix(
#'    data = rnorm(n = n * k, sd = 2), nrow = n, ncol = k
#' )
#'
#' cov <- A %*% t(A) / (n - 1)
#'
#' X <- t(MASS::mvrnorm(n = p, mu = rep(0, n), Sigma = cov))
#' Y <- X + matrix(rnorm(n * p, 0, 1), n, p)
#'
#' # set mash data
#' mash_data <- mashr::mash_set_data(Y)
#'
#' # call flash
#' fl <- mashy_flash(mash_data, backfit = TRUE)
#'
#' # get posterior mean
#' pm <- flash_get_pm(fl)
#'
flash_get_pm <- function(fl) {

  if (flashier:::get.n.factors(fl$flash.fit) == 0) {
    warning("Flash object does not have any factors.")
    return(matrix(data = 0, nrow = nrow(fl$flash.fit$Y), ncol = ncol(fl$flash.fit$Y)))
  }

  if(!("E.pm" %in% names(fl))) {

    warning("Posterior mean of E not found in fl. returning fitted(fl)
            Please fit fl with mishmashr::mashy_flash to get E[LF' + E | Y]")

    fl$L.pm %*% t(fl$F.pm)

  }

  fl$L.pm %*% t(fl$F.pm) + fl$E.pm

}

#' Fit Mash and Flash models
#'
#' @param mash_data Output from \code{mashr::mash_set_data}.
#' @param mash_params List of parameters to be passed to
#' \code{mashr::mash}
#' @param flash_params List of parameters to be passed to
#' \code{flashier::flash}.
#'
#' @return A list with components
#' \itemize{
#'   \item flash - a flash object (what is returned by \code{flashier::flash})
#'   \item mash - a list of mash outputs (what is returned by \code{mashr::mash})
#' }
#'
#' @examples
#'
#' # Generate 10 x 100 matrix with rank 3.
#' n <- 10
#' p <- 100
#' k <- 3
#' A <- matrix(
#'   data = rnorm(n = n * k, sd = 2), nrow = n, ncol = k
#' )
#'
#' cov <- A %*% t(A) / (n - 1)
#'
#' X <- t(MASS::mvrnorm(n = p, mu = rep(0, n), Sigma = cov))
#' Y <- X + matrix(rnorm(n * p, 0, 1), n, p)
#'
#' #set mash data
#' mash_data <- mashr::mash_set_data(Y)
#'
#' Ulist <- mashr::cov_canonical(mash_data)
#'
#' #call mishmash with specified parameters for mash and flash
#' mm <- mishmash(
#'   mash_data,
#'   mash_params = list(Ulist = Ulist, algorithm.version = 'R', posterior_samples = 10),
#'   flash_params = list(greedy.Kmax = 5, ebnm.fn = ebnm::ebnm_unimodal_symmetric)
#' )
#'
mishmash <- function(mash_data, mash_params = list(), flash_params = list()) {

  fit_list <- list()
  mash_params[['data']] <- mash_data
  flash_params[['mash_data']] <- mash_data
  fit_list[['flash']] <- do.call(mashy_flash, flash_params)
  fit_list[['mash']] <- do.call(mashr::mash, mash_params)
  return(fit_list)

}

