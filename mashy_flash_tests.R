# here, want to write testing utilities for functions in mashy flash
source("./simulation.R")

calculate_lfsr_calibration <- function(X_true, X_fitted, lfsr) {

  is_sign_false <- sign(X_true) != sign(X_fitted)

  lfsr_vec <- c()
  cutoff_seq <- seq(from = .05, to = .45, by = .05)

  for (cutoff in cutoff_seq) {

    lfsr_cutoff <- (lfsr <= cutoff)
    active_signs <- is_sign_false[lfsr_cutoff]

    if (length(active_signs) == 0) {

      est_lfsr <- NA

    } else {

      num_false_signs <- sum(active_signs)
      est_lfsr <- num_false_signs / length(active_signs)

    }

    lfsr_vec <- c(lfsr_vec, est_lfsr)

  }

  lfsr_est_df <- data.frame(
    cutoff = cutoff_seq,
    lfsr = lfsr_vec
  )

  return(lfsr_est_df)

}

# Generate 25 x 1000 matrix with rank 10.
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
fl <- flashier::flash(Y)
fl <- mashy_flash(
  mash_data, backfit = TRUE, greedy.Kmax = 10, var.type = NULL
)

lfsr_calib <- calculate_lfsr_calibration(
  X_true = X,
  X_fitted = fitted(fl),
  lfsr = fl$LF_E.lfsr
)

# now, a second test with E added
distr_list <- list(
  distr::UnivarMixingDistribution(
    distr::Norm(mean = 0, sd = 1),
    distr::Dirac(),
    mixCoeff = c(.2, .8)
  )
)

S_mat_2 <- matrix(data = 2, nrow = n, ncol = p)

S_mat_1 <- matrix(data = 1, nrow = n, ncol = p)

flash_sim_data <- generate_data_from_flash_model(
  n = 25,
  p = 1000,
  k = 5,
  factor_distr_list = distr_list,
  factor_distr_probs = c(1),
  loading_distr_list = distr_list,
  loading_distr_probs = c(1),
  S_1 = S_mat_1,
  S_2 = S_mat_2
)

mash_data <- mashr::mash_set_data(Bhat = flash_sim_data$Y, Shat = S_mat_2)
fl <- mashy_flash(
  mash_data, backfit = TRUE, greedy.Kmax = 10
)

lfsr_calib <- calculate_lfsr_calibration(
  X_true = flash_sim_data$X,
  X_fitted = fitted(fl) + fl$E.pm,
  lfsr = fl$LF_E.lfsr
)
