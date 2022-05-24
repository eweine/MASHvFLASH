library(mashr)
library(flashier)

## SIMULATION FUNCTIONS -------------------------------------------------

# n is number of conditions, p is number of genes

# Noise is i.i.d. N(0, 1)
get_E <- function(n, p, sd = 1) {
  matrix(rnorm(n * p, 0, sd), n, p)
}


# Simulate from null model ----------------------------------------------

null_sim <- function(n, p, seed = NULL) {
  set.seed(seed)
  Y <- get_E(n, p)
  true_Y <- matrix(0, n, p)
  
  list(Y = Y, true_Y = true_Y)
}


# Simulate from MASH model ----------------------------------------------

# Sigma is list of covariance matrices
# pi[j] is probability that effect j has covariance Sigma[[j]]
# s is sparsity (percentage of null effects)
mash_sim <- function(n, p, Sigma, pi = NULL, s = 0.8, seed = NULL) {
  set.seed(NULL)
  if (is.null(pi)) {
    pi = rep(1, length(Sigma)) # default to uniform distribution
  }
  assertthat::are_equal(length(pi), length(Sigma))
  for (j in length(Sigma)) {
    assertthat::are_equal(dim(Sigma[j]), c(n, n))
  }
  
  pi <- pi / sum(pi) # normalize pi to sum to one
  which_sigma <- sample(1:length(pi), p, replace=TRUE, prob=pi)
  nonnull_fx <- sample(1:p, floor((1 - s)*p), replace=FALSE)
  
  X <- matrix(0, n, p)
  for (j in nonnull_fx) {
    X[, j] <- MASS::mvrnorm(1, rep(0, n), Sigma[[which_sigma[j]]])
  }
  Y <- X + get_E(n, p)
  list(Y = Y, true_Y = X)
}


# Simulate from FLASH model ---------------------------------------------

# fs is sparsity of factors (percentage of null effects)
# fvar is variance of effects (generated from normal distribution)
# ls is sparsity of loadings
# lvar is variance of loadings
# UVvar is variance of dense rank-one matrix included to mimic something
#   like unwanted variation (set it to 0 to ignore it)
flash_sim <- function(n, p, k, fs, fvar, ls, lvar, UVvar = 0, seed = NULL) {
  set.seed(seed)
  
  nonnull_ll <- matrix(sample(c(0, 1), n*k, TRUE, c(ls, 1 - ls)), n, k)
  LL <- nonnull_ll * matrix(rnorm(n*k, 0, sqrt(lvar)), nrow=n, ncol=k)
  
  nonnull_ff <- matrix(sample(c(0, 1), p, TRUE, c(fs, 1 - fs)),
                       nrow=k, ncol=p, byrow=TRUE)
  FF <- nonnull_ff * matrix(rnorm(k*p, 0, sqrt(fvar)), nrow=k, ncol=p)
  
  X <- LL %*% FF
  Y <- X + get_E(n, p)
  # add unwanted variation
  Y <- Y + outer(rnorm(n, 0, sqrt(UVvar)), rnorm(p, 0, sqrt(UVvar)))
  list(Y = Y, true_Y = X)
}


## SIMULATIONS ----------------------------------------------------------

# Functions to generate six types of datasets. One is null; three are
# from the MASH model; two are from the FLASH model.

sim_fns <- function(n, p, s,
                    indvar, shvar, uniqvar,
                    r1var, r5var) {
  
  # 1. Everything is null
  sim_null <- function(){ null_sim(n, p) }
  
  Sigma <- list()
  
  # 2. Effects are independent across conditions
  Sigma[[1]] <- diag(rep(indvar, n))
  Sigma_ind <- Sigma
  sim_ind <- function(){ mash_sim(n, p, Sigma_ind, s=s) }
  
  # 3. Effects are either independent or shared
  Sigma[[2]] <- matrix(shvar, n, n)
  Sigma_indsh <- Sigma
  sim_indsh <- function(){ mash_sim(n, p, Sigma_indsh, s=s) }
  
  # 4. Effects are independent, shared, or unique to a single condition
  for (j in 1:n) {
    Sigma[[2 + j]] <- matrix(0, n, n)
    Sigma[[2 + j]][j, j] <- uniqvar
  }
  pi <- c(n, n, rep(1, n))
  sim_mash <- function(){ mash_sim(n, p, Sigma, pi, s=s) }
  
  # 5. Rank one model
  sim_rank1 <- function(){ flash_sim(n, p, 1, s, r1var, 0.2, 1) }
  
  # 6. Rank 5 model
  sim_rank5 <- function(){ flash_sim(n, p, 5, s, r5var, 0.8, 1) }
  
  # 7. Rank 3 model with unwanted variation
  # sim_UV <- function(){ flash_sim(n, p, 3, s, r3var, 0.5, 1, UVvar) }
  
  #c(sim_null, sim_ind, sim_indsh, sim_mash, sim_rank1, sim_rank5)
  c(sim_ind, sim_indsh, sim_mash, sim_rank1, sim_rank5)
}

sim_names <- c(#"Null simulation",
               "All independent effects",
               "Independent and shared",
               "Independent, shared, and unique",
               "Rank 1 FLASH model",
               "Rank 5 FLASH model")

# Fit using FLASH -------------------------------------------------------
fit_flash_zero <- function(data, Kmax, ebnm_fn=ebnm::ebnm_point_normal,
                           init_fn=init.fn.default, greedy=TRUE, backfit=TRUE,
                           warmstart=TRUE) {
  # if (is.matrix(data)) {
  #   data <- flash_set_data(data, S = 1)
  # }
  n <- nrow(data)
  
  flash_obj <- flash.init(
    data = data,
    S = 1,
    var.type = 0
  )
  
  t0 <- Sys.time()
  if (greedy) {
    res <- flash.add.greedy(flash_obj, Kmax, 
                            ebnm.fn=ebnm_fn, init.fn=init_fn,
                            warmstart=warmstart)
    fl <- res$f
  } else {
    fl <- flash_add_factors_from_data(data, Kmax, init_fn=init_fn)
  }
  t1 <- Sys.time()
  if (backfit) {
    res <- flash.backfit(flash_obj)
    fl <- res$f
  }
  t2 <- Sys.time()
  
  t.greedy <- t1 - t0
  t.backfit <- t2 - t1
  
  list(f = fl, t.greedy = t.greedy, t.backfit = t.backfit)
}


fit_flash_OHL <- function(data, Kmax, ebnm_fn=ebnm::ebnm_point_normal,
                          init_fn=init.fn.default, greedy=TRUE, backfit=TRUE,
                          warmstart=TRUE) {
  # if (is.matrix(data)) {
  #   data <- flash_set_data(data, S = 1)
  # }
  n <- nrow(data)
  p <- ncol(data)
  canonical <- cbind(rep(1, n), diag(rep(1, n)))
  k <- ncol(canonical)
  
  zero.res <- fit_flash_zero(data, Kmax, ebnm_fn, init_fn,
                             greedy, backfit=FALSE, warmstart)
  
  # add fixed factors with random loadings
  # I should perhaps add a better initialization of the loadings
  fl <- flash.init.factors(
    flash = zero.res$f,
    init = list(
      canonical,
      matrix(
        data = rnorm(n = p * k),
        nrow = p, 
        ncol = k
        )
    )
  )
  
  # NOTE: I think there needs to be logic here that figures out the dims
  # of the current factors and then figures out what to do
  # this might be easier to see once I've added the factors above
  fl <- flash.fix.factors(
    flash = fl,
    kset = c((fl$n.factors - k + 1):k),
    mode = 1
  )
  
  #fl <- flash_add_fixed_l(data, canonical, zero.res$f, init_fn=init_fn)
  
  t0 <- Sys.time()
  if (backfit) {
    
    res <- flash_backfit(fl)
    fl <- res$f
    
  } else {
    
    res <- flash.backfit(fl, kset=c((fl$n.factors - k + 1):k))
    fl <- res$f
    
  }
  t1 <- Sys.time()
  
  t.backfit <- Sys.time() - t0
  
  list(f = fl, t.greedy = zero.res$t.greedy, t.backfit = t.backfit)
}


fit_flash_OHF <- function(data, Kmax, ebnm_fn=ebnm::ebnm_point_normal,
                          init_fn=init.fn.default, greedy=TRUE, backfit=TRUE,
                          warmstart=TRUE) {
  # if (is.matrix(data)) {
  #   data <- flash_set_data(data, S = 1)
  # }
  n <- nrow(data)
  p <- ncol(data)
  canonical <- cbind(rep(1, n), diag(rep(1, n)))
  k <- ncol(canonical)
  
  t0 <- Sys.time()
  
  # I think here I need to initialize the factors
  
  flash_obj <- flash.init(
    data = data,
    S = 1,
    var.type = 0
  )
  
  fl <- flash.init.factors(
    flash = flash_obj,
    init = list(
      canonical,
      matrix(
        data = rnorm(n = p * k),
        nrow = p, 
        ncol = k
      )
    )
  )
  
  fl <- flash.fix.factors(
    flash = fl,
    kset = c((fl$n.factors - k + 1):k),
    mode = 1
  )
  
  res <- flash.backfit(fl)
  fl <- res$f
  t1 <- Sys.time()
  
  if (greedy) {
    res <- flash.add.greedy(fl, Kmax,
                            ebnm.fn=ebnm_fn, init.fn=init_fn,
                            warmstart=warmstart)
    fl <- res$f
  } else {
    fl <- flash_add_factors_from_data(data, Kmax, fl, init_fn=init_fn)
  }
  t2 <- Sys.time()
  K = fl$n.factors
  if (backfit && K > ncol(canonical)) {
    res <- flash.backfit(fl,
                         kset=(ncol(canonical) + 1):K
                         )
    fl <- res$f
  }
  t3 <- Sys.time()
  
  t.greedy <- t2 - t1
  t.backfit <- (t1 - t0) + (t3 - t2)
  
  list(f = fl, t.greedy = t.greedy, t.backfit = t.backfit)
}


# Fit using MASH -------------------------------------------------------
fit_mash <- function(data) {
  if (is.matrix(data)) {
    data <- mash_set_data(t(data))
  }
  timing <- list()
  
  # time to create canonical matrices is negligible
  U = cov_canonical(data)
  
  t0 <- Sys.time()
  m.1by1 <- mash_1by1(data)
  lvl <- 0.05
  strong <- get_significant_results(m.1by1, lvl)
  while (length(strong) < 5 && lvl < 0.5) {
    lvl <- lvl + 0.05
    strong <- get_significant_results(m.1by1, lvl)
  }
  if (length(strong) >= 5) {
    U.pca <- cov_pca(data, 5, strong)
    U.ed <- cov_ed(data, U.pca, strong)
    U <- c(U, U.ed)
    t.ed <- Sys.time() - t0
  } else {
    t.ed <- as.difftime(0, units="secs")
  }
  
  t0 <- Sys.time()
  m <- mash(data, U)
  t.mash <- Sys.time() - t0
  
  list(m = m, t.ed = t.ed, t.mash = t.mash)
}

flash_diagnostics <- function(fl, Y, true_Y, nsamp) {
  MSE <- flash_mse(fl, true_Y)
  
  # Sample from FLASH fit to estimate CI coverage and TPR vs. FPR
  # NOTE: I don't see a sampler object in the fitted flash object.
  fl_sampler <- flashr::flash_sampler(Y, fl, fixed="loadings")
  fl_samp <- fl_sampler(nsamp)
  
  CI <- flash_ci(fl_samp, true_Y)
  ROC <- flash_roc(fl, fl_samp, true_Y)
  
  list(MSE = MSE, CI = CI, TP = ROC$TP, FP = ROC$FP,
       n_nulls = ROC$n_nulls, n_nonnulls = ROC$n_nonnulls)
}

mash_diagnostics <- function(m, true_Y) {
  MSE <- mash_mse(m, true_Y)
  CI <- mash_ci(m, true_Y)
  ROC <- mash_roc(m, true_Y)
  
  list(MSE = MSE, CI = CI, TP = ROC$TP, FP = ROC$FP,
       n_nulls = ROC$n_nulls, n_nonnulls = ROC$n_nonnulls)
}


# MSE of posterior means (FLASH) ----------------------------------------
flash_mse <- function(fl, true_Y) {
  mean((fitted(fl) - true_Y)^2)
}

# MSE for MASH ----------------------------------------------------------
mash_mse <- function(m, true_Y) {
  mean((get_pm(m) - t(true_Y))^2)
}


# 95% CI coverage for FLASH ---------------------------------------------
flash_ci <- function(fl_samp, true_Y) {
  n <- nrow(true_Y)
  p <- ncol(true_Y)
  nsamp <- length(fl_samp)
  
  flat_samp <- matrix(0, nrow=n*p, ncol=nsamp)
  for (i in 1:nsamp) {
    flat_samp[, i] <- as.vector(fl_samp[[i]])
  }
  CI <- t(apply(flat_samp, 1, function(x) {quantile(x, c(0.025, 0.975))}))
  mean((as.vector(true_Y) >= CI[, 1]) & (as.vector(true_Y) <= CI[, 2]))
}

# 95% CI coverage for MASH ----------------------------------------------
mash_ci <- function(m, true_Y) {
  Y <- t(true_Y)
  mean((Y > get_pm(m) - 1.96 * get_psd(m))
       & (Y < get_pm(m) + 1.96 * get_psd(m)))
}


# LFSR for FLASH --------------------------------------------------------
flash_lfsr <- function(fl_samp) {
  nsamp <- length(fl_samp)
  n <- nrow(fl_samp[[1]])
  p <- ncol(fl_samp[[1]])
  
  pp <- matrix(0, nrow=n, ncol=p)
  pn <- matrix(0, nrow=n, ncol=p)
  for (i in 1:nsamp) {
    pp <- pp + (fl_samp[[i]] > 0)
    pn <- pn + (fl_samp[[i]] < 0)
  }
  1 - pmax(pp, pn) / nsamp
}


# Quantities for plotting ROC curves -----------------------------------
flash_roc <- function(fl, fl_samp, true_Y, step=0.01) {
  roc_data(flash_get_fitted_values(fl), true_Y, flash_lfsr(fl_samp), step)
}

mash_roc <- function(m, true_Y, step=0.01) {
  roc_data(get_pm(m), t(true_Y), get_lfsr(m), step)
}

roc_data <- function(pm, true_Y, lfsr, step) {
  correct_sign <- pm * true_Y > 0
  is_null <- true_Y == 0
  n_nulls <- sum(is_null)
  n_nonnulls <- length(true_Y) - n_nulls
  
  ts <- seq(0, 1, by=step)
  tp <- rep(0, length(ts))
  fp <- rep(0, length(ts))
  
  for (t in 1:length(ts)) {
    signif <- lfsr <= ts[t]
    tp[t] <- sum(signif & correct_sign)
    fp[t] <- sum(signif & is_null)
  }
  
  list(ts = ts, TP = tp, FP = fp, n_nulls = n_nulls, n_nonnulls = n_nonnulls)
}

# Functions for running simulations and combining results ---------------

run_sims <- function(sim_fn, nsims, plot_title, fpath, Kmax=50,
                     backfit=FALSE) {
  if (nsims == 1) {
    res = run_one_sim(sim_fn, Kmax, backfit)
  } else {
    res = run_many_sims(sim_fn, nsims, Kmax, backfit)
  }
  
  saveRDS(output_res_mat(res$res), paste0(fpath, "res.rds"))
  if (!(plot_title == "Null simulation")) {
    png(paste0(fpath, "ROC.png"))
    plot_ROC(res$res, plot_title)
    dev.off()
  }
  png(paste0(fpath, "time.png"))
  plot_timing(res$timing)
  dev.off()
  
  return(res)
}

run_one_sim <- function(sim_fn, Kmax, nsamp=200, backfit) {
  data <- do.call(sim_fn, list())
  
  # If there are no strong signals, trying to run ED throws an error, so
  #   we need to do some error handling to fit the MASH object
  mfit <- fit_mash(data$Y)
  
  flfits <- list()
  flfits$Zero <- fit_flash_zero(data$Y, Kmax, backfit=backfit)
  flfits$OHL <- fit_flash_OHL(data$Y, Kmax, backfit=backfit)
  flfits$OHF <- fit_flash_OHF(data$Y, Kmax, backfit=backfit)
  
  message("Running MASH diagnostics")
  mres <- mash_diagnostics(mfit$m, data$true_Y)
  
  message("Running FLASH diagnostics")
  methods <- names(flfits)
  flres <- list()
  for (method in methods) {
    flres[[method]] <- flash_diagnostics(flfits[[method]]$f, data$Y,
                                         data$true_Y, nsamp)
  }
  
  # Combine results:
  timing = lapply(flfits,
                  function(method) {
                    list(ed.or.greedy = method$t.greedy,
                         mash.or.backfit = method$t.backfit)
                  })
  timing$MASH <- list(ed.or.greedy = mfit$t.ed,
                      mash.or.backfit = mfit$t.mash)
  
  all_res <- flres
  all_res$MASH <- mres
  
  list(timing = timing, res = all_res)
}

run_many_sims <- function(sim_fn, nsims, Kmax, backfit) {
  res <- list()
  combined_res <- list()
  
  for (i in 1:nsims) {
    message(paste("  Simulating dataset", i))
    res[[i]] <- run_one_sim(sim_fn, Kmax, backfit=backfit)
  }
  elems <- names(res[[1]])
  for (elem in elems) {
    combined_res[[elem]] <- list()
    subelems <- names(res[[1]][[elem]])
    for (subelem in subelems) {
      combined_res[[elem]][[subelem]] <- list()
      subsubelems <- names(res[[1]][[elem]][[subelem]])
      for (subsubelem in subsubelems) {
        tmp <- lapply(res, function(x) {x[[elem]][[subelem]][[subsubelem]]})
        combined_res[[elem]][[subelem]][[subsubelem]] <- Reduce(`+`, tmp) / nsims
      }
    }
  }
  return(combined_res)
}


# Plotting and output functions -----------------------------------------

plot_timing <- function(timing, units="secs") {
  data <- sapply(timing, unlist)
  data <- data[, c(4, 1, 2, 3)]
  barplot(colSums(data), axes=T,
          main=paste("Average time to fit in", units),
          names.arg = colnames(data),
          #legend.text = c("ED/Greedy", "MASH/Backfit"),
          ylim = c(0, max(colSums(data))*2),
          cex.names = 0.8)
  # (increasing ylim is easiest way to deal with legend getting in way)
}

plot_ROC <- function(res, main="ROC curve") {
  # Number of nulls and nonnulls are identical across methods
  n_nonnulls <- res[[1]]$n_nonnulls
  n_nulls <- res[[1]]$n_nulls
  
  m_y <- res$MASH$TP / n_nonnulls
  m_x <- res$MASH$FP / n_nulls
  plot(m_x, m_y, xlim=c(0, 1), ylim=c(0, 1), type='l',
       xlab='FPR', ylab='TPR', main=main, col="orange")
  
  idx <- 1:3
  colors <- c("skyblue", "seagreen", "yellow2")
  for (i in idx) {
    y <- res[[i]]$TP / n_nonnulls
    x <- res[[i]]$FP / n_nulls
    lines(x, y, col=colors[i])
  }
  legend("bottomright", c("MASH", names(res)[idx]), lty=1,
         col=c("orange", colors[idx]), cex=0.8)
}

output_res_mat <- function(res) {
  mat <- rbind(sapply(res, function(method) {method$MSE}),
               sapply(res, function(method) {method$CI}))
  
  rownames(mat) = c("MSE", "95% CI cov")
  return(mat)
}


set.seed(1)

# s is sparsity (prop null genes); var parameters define sizes of effects
all_sims <- sim_fns(n=25, p=1000, s=0.8,
                    indvar=4, shvar=2, uniqvar=100,
                    r1var=4, r5var=1)

for (i in 1:length(all_sims)) {
  message(paste("  Beginning simulation #", i, sep=""))
  res <- run_sims(all_sims[[i]], nsims=10, plot_title=sim_names[[i]],
                  fpath = paste0("./output/MASHvFLASHsims/greedy/sim", i),
                  Kmax=25, backfit=FALSE)
}

for (i in 1:length(all_sims)) {
  message(paste("  Beginning simulation #", i, sep=""))
  res <- run_sims(all_sims[[i]], nsims=10, plot_title=sim_names[[i]],
                  fpath = paste0("./output/MASHvFLASHsims/backfit/sim", i),
                  Kmax=25, backfit=TRUE)
}



