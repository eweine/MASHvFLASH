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
    var.type = NULL,
    ...
  )

  fl[['lfsr']] <- flash_lfsr(fl)
  return(fl)

}

#' Get Local False Sign Rate (lfsr) for Flash Model
#'
#' @param fl A \code{flash} object.
#'
#' @return A matrix with the lfsr for each element of the data matrix.
#' @export
#'
#' @examples
flash_lfsr <- function(fl, nsamp = 1000) {

  n <- nrow(fl$flash.fit$Y)
  p <- ncol(fl$flash.fit$Y)

  fx_samp <- flash_sample_fx(fl, nsamp)
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

#' Sample from the posterior predictive distribution of a flash object.
#'
#' @param fl flash object
#' @param mash_data output of \code{mashr::mash_set_data}
#' @param nsamp number of samples to draw
#'
#' @return list of sampels
#' @export
#'
#' Note that this sampler is only valid if the error matrix E is not estimated
#' in the flash object. Future versions of the code could take that into account
#' as well.
#'
#' @examples
flash_posterior_pred_sample <- function(
  fl,
  mash_data,
  nsamp
) {

  fx_samp <- flash_sample_fx(fl, nsamp)
  noise_samp <- replicate(
    n = nsamp,
    expr = apply(mash_data$Shat, MARGIN = c(1, 2), FUN = rnorm, n = 1, mean = 0),
    simplify = FALSE
  )

  post_pred_samp <- mapply("+", fx_samp, noise_samp, SIMPLIFY = FALSE)
  return(post_pred_samp)

}

#' Compute a posterior predictive check on a sample from a multivariate model.
#'
#' @param data data matrix for which the model was fit
#' @param sample posterior predictive sample from the model
#' @param distance_norm norm to compute distance between matrices. See
#' documentation of \code{norm} for more details.
#'
#' @return p-value of posterior predictive test
#' @export
#'
#' @examples
matrix_posterior_predictive_test <- function(
  data,
  sample,
  distance_norm = c("O", "I", "F", "M", "2")
) {

  distance_norm <- match.arg(distance_norm)

  nsamp <- length(sample)
  samp_distance_vec <- numeric(as.integer(nsamp * (nsamp - 1) / 2))
  data_distance_vec <- numeric(nsamp)
  counter <- 0

  for (i in 1:(nsamp - 1)) {

    for (j in (i+1):nsamp) {

      samp_norm <- norm(sample[[i]] - sample[[j]], type = distance_norm)
      counter <- counter + 1
      samp_distance_vec[counter] <- samp_norm

    }

    data_distance_vec[i] <- norm(sample[[i]] - data, type = distance_norm)

  }

  data_distance_vec[nsamp] <- norm(sample[[nsamp]] - data, type = distance_norm)

  test_out <- wilcox.test(
    x = data_distance_vec,
    y = samp_distance_vec,
    alternative = "greater"
  )

  return(test_out$p.val)

}


# Generate 10 x 500 matrix with rank 3.
n <- 25
p <- 1000
k <- 10
A <- matrix(
  data = rnorm(n = n * k, sd = 2), nrow = n, ncol = k
)

cov <- A %*% t(A) / (n - 1)

X <- t(MASS::mvrnorm(n = p, mu = rep(0, n), Sigma = cov))
Y <- X + matrix(rnorm(n * p, 0, 1), n, p)

# set mash data
mash_data <- mashr::mash_set_data(Y)
fl <- mashy_flash(mash_data, backfit = TRUE, greedy.Kmax = 5)

fl_pp_samp <- flash_posterior_pred_sample(fl, mash_data, 100)

p_val <- matrix_posterior_predictive_test(Y, fl_pp_samp)

U <- mashr::cov_canonical(mash_data)

m <- mashr::mash(mash_data, U, outputlevel = 4)

