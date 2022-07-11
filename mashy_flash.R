#' Sample Posterior Effect Sizes from Flash Model
#'
#' @param flash flash object
#' @param nsamp number of samples
#'
#' @return List of samples
#'
#' @examples
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
#' @export
#'
#' @examples
flash_get_E_pm <- function(fl, fit_E) {

  if (!fit_E) {

    return(
      matrix(data = 0, nrow = nrow(fl$flash.fit$Y), ncol = ncol(fl$flash.fit$Y))
    )

  } else {

    resid <- fl$flash.fit$Y - fl$L.pm %*% t(fl$F.pm)
    E_pm <- resid * (1 / fl$flash.fit$given.S2)
    return(E_pm)

  }

}

#' Get Local False Sign Rate (lfsr) for Flash Model
#'
#' @param fl A \code{flash} object.
#' @param fit_E Boolean indicating if E was fit in the flash model. This is
#' equivalent to \code{var.type = NULL} in \code{flashier::flash}.
#' @param nsamp Number of samples to use to calculate lfsr.
#'
#' @return A matrix with the lfsr for each element of the data matrix.
#' @export
#'
#' @examples
flash_lfsr_LF_E <- function(fl, fit_E, nsamp = 5000) {

  n <- nrow(fl$flash.fit$Y)
  p <- ncol(fl$flash.fit$Y)

  fx_samp <- flash_sample_fx(fl, nsamp)

  if (fit_E) {

    cat("Generating Posterior Samples for lfsr...\n")
    pb <- txtProgressBar(min = 0, max = nsamp, initial = 0, style = 3)
    # sample from posterior of E | LF'
    for (i in 1:nsamp) {

      setTxtProgressBar(pb, i)
      resid <- fl$flash.fit$Y - fx_samp[[i]]
      E_mean_mat <- resid * (1 / fl$flash.fit$given.S2)
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
    close(pb)

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
mashy_flash <- function(mash_data, ...) {

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

  fl[['LF_E.lfsr']] <- flash_lfsr_LF_E(fl, fit_E)
  fl[['E.pm']] <- flash_get_E_pm(fl, fit_E)
  return(fl)

}

#' Fit Mash and Flash models
#'
#' @param mash_data Output from \code{mashr::mash_set_data}.
#' @param mash_params List of parameters to be passed to
#' \code{mashr::mash}
#' @param flash_params List of parameters to be passed to
#' \code{flashier::flash}.
#'
#' @return
#' @export
#'
#' @examples
fit_multivariate_models <- function(mash_data, mash_params, flash_params) {

  fit_list <- list()
  fit_list['flash'] <- do.call(mashy_flash, flash_params)
  fit_list['mash'] <- do.call(mashr::mash, mash_params)
  return(fit_list)

}

