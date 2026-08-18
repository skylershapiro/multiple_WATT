# =============================================================================
# Weighted Average Treatment Effect on Treated Under Multple Treatments
# =============================================================================
# Description:  Provides methods for data generation, effect estimates, and 
#               bootstrapped variance estimation. Includes methods adapted from
#               https://github.com/yiliu1998/ATTweights. 
#
# Author:       Skyler Shapiro
# Updated:      2026-07-14
# =============================================================================
library(nnet)
library(MASS)

# -----------------------------------------------------------------------------
# description: Pairwise WATT estimators for multiple-treatments; uses generalized propensity. 
# params:      y, z, X, reference_group, treatment_group, method, alpha, epsilon
# return:      double
# -----------------------------------------------------------------------------

sigmoid <- function(x, k) { 1 / (1 + exp(-k * x)) }

WATT <- function(y, z, X, reference_group = 1, treatment_group = 2, method="att", alpha=0, true_e = NULL){
  z <- factor(z, levels = 1:5)
  
  e_full <- if (is.null(true_e)) { # Optionally use true propensity matrix to calculate truth values for given estimand
    predict(nnet::multinom(z ~ X, trace = FALSE), type = "probs")
  } else {
    true_e
  }
  
  # propensity for reference and treatment
  e_j <- e_full[, reference_group]
  e_i <- e_full[, treatment_group] 
  
  # group indicators
  d_lj <- z == reference_group
  d_li <- z == treatment_group
  
  t = 1000
  
  if(method=="att"){h <- 1}
  if(method=="overlap"){h <- 1 / rowSums(1 / e_full)}
  if(method=="trimming"){ # refit per crump et all and li^2
    h <- d_lj | ( e_i >= alpha & e_i <= 1 )
    # X <- X[compare_inds,]
    X <- X[h,]
    y <- y[h]
    z <- z[h]
    
    e_full <- predict(nnet::multinom(z ~ X, trace = FALSE), type = "probs")
  
    e_j <- e_full[, reference_group]
    e_i <- e_full[, treatment_group] 
    d_li <- d_li[h]
    d_lj <- d_lj[h]
    h <- 1
    }
  #if(method=="trimming"){h <- sigmoid(e_j - alpha, t) - sigmoid(e_j - (1 - alpha), t)}
  if(method=="influence_trimming"){
    h <- (d_li * e_j * e_i^(-1) / sum(d_li * e_j * e_i^(-1))) <= alpha
    X <- X[h,]
    y <- y[h]
    z <- z[h]
  
    e_full <- predict(nnet::multinom(z ~ X, trace = FALSE), type = "probs")
  
  e_j <- e_full[, reference_group]
  e_i <- e_full[, treatment_group] 
  d_li <- d_li[h]
  d_lj <- d_lj[h]
  h <- 1
}
  if(method=="truncation"){
    numerator_sum <- rowSums(e_full) - e_i - e_j
    h <- I(e_i >= alpha & e_i <= 1) + I(e_i < alpha)*(1 - alpha - numerator_sum)/alpha 
  }
  if(method=="matching_weights"){h <- apply(e_full,1,min)}
  if(method=="shannon_entropy"){h <- -1*(e_j*log(e_j) + e_i*log(e_i))}
  if(method=="beta_weights"){h <- (e_j*e_i)^(alpha-1)}
  
  t1 <- sum(y * d_lj) / sum(d_lj)
  t2 <- sum(e_j * e_i^(-1) * h * y * d_li) / sum(e_j * e_i^(-1) *  h * d_li)
  tau <- t1 - t2
  
  return(tau)
}

# -----------------------------------------------------------------------------
# purpose:  Generates covariates X with sample size N.
# params:   N
# return:   matrix
# -----------------------------------------------------------------------------
make_data <- function(N = 6000) {
  mu_mvn    <- c(2, 1, 1)
  Sigma_mvn <- matrix(c( 1.0,  0.5,  0.3,
                         0.5,  1.0,  0.2,
                         0.3,  0.2,  1.0), nrow = 3)
  
  MVN <- mvrnorm(N, mu = mu_mvn, Sigma = Sigma_mvn)
  X1  <- MVN[, 1]
  X2  <- MVN[, 2]
  X3  <- MVN[, 3]
  X4 <- runif(N, -3, 3)
  X5 <- rchisq(N, df = 1)
  X6 <- rbinom(N,1,0.5)
  dat <-  cbind(X1, X2, X3, X4, X5, X6)
  return(dat)
}

# -----------------------------------------------------------------------------
# purpose: Generate treatment Z under defined overlap setting.
# params:  X, N, kappas, alphas, report_ps
# return:  numeric OR matrix
# -----------------------------------------------------------------------------
make_treatment <- function(X, N = nrow(X),kappas = c(0, 0.2, 0.2, 0.2, 0.2), alphas = NULL, report_ps  = FALSE) {
  
  # BETA_DIRS <- matrix(
  #   c(c( 0,  0,  0,  0,  0,  0),   # Grp 1: Ref
  #     c( 1, 1, 1.5,  1,  1,  1),   
  #     c( 1, 1, 1,  1,  1,  -2),  
  #     c( 1, 1, 1,  1,  1,  3),  
  #     c( 1, 1, 1,  -2,  1,  1)), nrow = 6, ncol = 5)

  BETA_DIRS <- matrix(
    c(c(  0,   0,   0,   0,   0,   0),   # Grp 1: Ref
      c(0.3, 0.3, 3.0,   0,   0,   0),   # Grp 2: driven by X3
      c(0.3, 0.3,   0, 3.0,   0,   0),   # Grp 3: driven by X4
      c(0.3, 0.3,   0,   0, 3.0,   0),   # Grp 4: driven by X5
      c(0.3, 0.3,   0,   0,   0, 3.0)),  # Grp 5: driven by X6
    nrow = 6, ncol = 5)
  
  B <- kappas * t(BETA_DIRS)
  if(is.null(alphas)) alphas <- c(2.5,0.05,0.15,0.2, 0.1)
  eta <- X %*% t(B) + matrix(alphas, N, 5, byrow = TRUE)
  ps_matrix <- exp(eta) / rowSums(exp(eta))
  
  Z <- apply(ps_matrix, 1, function(probs) sample(1:5, size = 1, prob = probs))
  
  if (report_ps) return(cbind(Z, ps_matrix))
  return(Z)
}

# -----------------------------------------------------------------------------
# purpose: Generates true potential outcomes Y1-Y5 and observed outcome Y_obs.
# params:  X, Z, N
# return:  matrix
# -----------------------------------------------------------------------------
make_outcome <- function(X, Z, N = 6000) {
  gam1 <- c(-1.5, 1, 1, 1, 1, 1, 1)
  gam2 <- c(-4, 2, 3, 1, 2, 2, 2)
  gam3 <- c(4, 3, 1, 2, -1, -1, -4)
  gam4 <- c(1, 4, 1, 2, -1, -1, -3)
  gam5 <- c(3.5, 5, 1, 2, -1, -1, -2)
  GAMMA <- matrix(c(gam1, gam2, gam3, gam4, gam5),
                  byrow = T,
                  ncol = 7)
  noise <- rnorm(N)
  Y <- cbind(1, X) %*% t(GAMMA) + noise
  Y_obs <- Y[cbind(1:N, Z)]
  Ys <- cbind(Y, Y_obs)
  return(Ys)
}

# -----------------------------------------------------------------------------
# purpose: Wrapper for `WATT()`. Calculates pairwise effect estimates for single
#          simulation setting against reference group.
# params:  Y_obs, Z, X, method
# return:  named list[double] 
# -----------------------------------------------------------------------------
calculate_watts <- function(Y_obs, Z, X, method = "att", alpha = 0, true_e = NULL){
  watts <- rep(NA, 4)
  for (i in 2:5){ 
    watts[i-1] <- WATT(Y_obs,Z,X, method = method, reference_group = 1, treatment_group = i, alpha = alpha, true_e = true_e) # 2v1 3v1 4v1 5v1
  }
  return(watts)
}

# -----------------------------------------------------------------------------
# purpose: Compute true ATT in superpopulation under defined data/overlap setting.
# params:  kappas, alphas, N, n_reps, seed
# return:  named list[double]
# -----------------------------------------------------------------------------

# true_taus <- compute_true_tau(kappas = settings[[s]]$kappas,
#                              alphas = rep(0, 5), tasks = my_tasks)

compute_true_tau <- function(kappas = c(0, 0.2, 0.2, 0.2, 0.2), tasks , alphas = NULL,
                             N = 100000, epsilon = 1e-3, seed = 42) {
  
  
  true_taus <- vector("list", length(tasks))
  names(true_taus) <- names(tasks)
  
 
  runs <- vector("list", length(tasks))
  curr_tau <- vector("list", length(tasks))
  converged <- rep(FALSE, length(tasks))
  r <- 0
  
  while (!all(converged)) {
    r <- r + 1
    set.seed(seed + r)
    X  <- make_data(N)
    ps <- make_treatment(X, N, kappas = kappas, alphas = alphas, report_ps = TRUE)
    Z  <- ps[, 1]
    e  <- ps[, -1]
    Y  <- make_outcome(X, Z, N)
    for (j in seq_along(tasks)) {
      t <- tasks[[j]]
    val <- sapply(2:5, function(i)
        WATT(Y[, 6], Z, X, reference_group = 1, treatment_group = i,
                        method = t$method, alpha = t$alpha, true_e = e))
   
    runs[[j]] <- rbind(runs[[j]], val)
    prev_tau <- curr_tau[[j]]
    curr_tau[[j]] <- colMeans(runs[[j]])
    
    if (!is.null(prev_tau) && !converged[j]) {
      diff <- max(abs(prev_tau - curr_tau[[j]]))
      cat(sprintf("%s - Run %d, max diff: %g\n", names(tasks)[j], r, diff))
      converged[j] <- diff < epsilon
    }
    }
  }
  for (j in seq_along(tasks)) {
    cat(sprintf("%s converged after %d runs (epsilon = %g)\n", names(tasks)[j], nrow(runs[[j]]), epsilon))
    true_taus[[j]] <- setNames(curr_tau[[j]], paste0("ATT_", 2:5, "_vs_1"))
  }
  
  true_taus
}


# -----------------------------------------------------------------------------
# purpose: Simulation function for WATT method comparison in single setting.
#          Uses n_boot bootstrap replicates for SE estimates.
# params:  N, seed, kappas, alphas
# return:  named list[double]
# -----------------------------------------------------------------------------
single_sim <- function(N = 6000, seed = 1, kappas = c(0, 0.2, 0.2, 0.2, 0.2),
                       tasks = NULL, alphas = NULL, n_boot = NULL) { 
  set.seed(seed)
  X <- make_data(N)
  Z <- make_treatment(X, N, kappas = kappas, alphas = alphas)
  Y <- make_outcome(X, Z, N)
  
  point_est <- lapply(tasks, function(t) {
    calculate_watts(Y[, 6], Z, X, method = t$method, alpha = t$alpha)
  })
  
  if (is.null(n_boot)) { n_boot <- 1 }
  
  # bootstrap implementation: sample 1:N with replacement `n_boot` times 
  boot_reps <- lapply(seq_len(n_boot), function(b) {
    idx <- sample(N, N, replace = TRUE)
    lapply(tasks, function(t) {
      calculate_watts(Y[idx, 6], Z[idx], X[idx, ],
                      method = t$method, alpha = t$alpha)
    })
  })
  
  # boot_se[[task_name]] = vector of length 4 (one SE per comparison)
  boot_se <- lapply(names(tasks), function(nm) {
    boot_mat <- do.call(rbind, lapply(boot_reps, `[[`, nm))
    apply(boot_mat, 2, sd)
  })
  
  names(boot_se) <- names(tasks)
  list(point_est = point_est, boot_se = boot_se)
}

# -----------------------------------------------------------------------------
# purpose: Density plot of generalized propensities for all groups
# params:  N, seed, kappas, alphas
# return:  NONE
# -----------------------------------------------------------------------------
plot_overlap <- function(kappas = c(0, 0.2, 0.2, 0.2, 0.2),
                         alphas = NULL,
                         N      = 5000,
                         seed   = 1) {
  set.seed(seed)
  X      <- make_data(N)
  ps_obj <- make_treatment(X, N, kappas = kappas, alphas = alphas, report_ps = TRUE)
  
  Z      <- ps_obj[, 1]
  ps_mat <- ps_obj[, -1]
  
  old_par <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  on.exit(par(old_par))
  
  group_colors <- c("#E69F00", "#56B4E9", "#CC79A7", "#D55E00", "#009E73")
  
  for (j in 1:5) {
    # Pre-calculate densities to find the global max Y for this subplot
    dens_list <- lapply(1:5, function(g) {
      scores <- ps_mat[Z == g, j]
      if (length(scores) > 1) density(scores, from = 0, to = 1) else NULL
    })
    
    # Find max Y across all valid density objects, default to 5 if empty
    max_y <- max(sapply(dens_list, function(d) if (!is.null(d)) max(d$y) else 0))
    if (max_y == 0) max_y <- 5
    
    # Initialize plot with the dynamic ylim (adding 10% buffer for legend space)
    plot(NULL, xlim = c(0, 1), ylim = c(0, max_y * 1.1), 
         xlab = paste0("e", j, "(X)"), ylab = "Density",
         main = paste("Distribution of e", j), bty = "l")
    
    # Draw the stored density lines
    for (g in 1:5) {
      if (!is.null(dens_list[[g]])) {
        lines(dens_list[[g]], col = group_colors[g], lwd = 1.5)
      }
    }
    
    if (j == 1) {
      legend("topright", legend = paste("Group", 1:5), 
             col = group_colors, lwd = 1.5, bty = "n", cex = 0.8)
    }
  }
}

# -----------------------------------------------------------------------------
# purpose: Calculates effective sample size
# params:  Z, e, tasks
# return:  list
# -----------------------------------------------------------------------------
calculate_ess <- function(Z, e, tasks, reference_group = 1) {
  sapply(tasks, function(t) {
    sapply(2:5, function(i) {
      compare_inds <- Z %in% c(reference_group, i)
      e_sub <- e[compare_inds, ]
      z_sub <- Z[compare_inds]
      
      e_j <- e_sub[, reference_group]
      e_i <- e_sub[, i]
      d_li <- z_sub == i
      
      if (t$method == "att")                h <- rep(1, sum(compare_inds))
      if (t$method == "overlap")            h <- 1 / rowSums(1 / e_sub)
      if (t$method == "matching_weights")   h <- apply(e_sub, 1, min)
      if (t$method == "trimming")           h <- as.numeric(e_j >= t$alpha)
      if (t$method == "truncation")         h <- pmax(e_j, t$alpha) / e_j
      if (t$method == "influence_trimming") h <- (e_j * e_i^(-1) / sum(d_li * e_j * e_i^(-1))) <= t$alpha
      
      weights <- e_j * e_i^(-1) * h
      weights <- weights[d_li]  # only treated units
      return(sum(weights)^2 / sum(weights^2))
    })
  })
}

# -----------------------------------------------------------------------------
# purpose: Calculates bias, ARbias, RMSE from bootstrapped simulation run
# params:  res, true_tau
# return:  list
# -----------------------------------------------------------------------------
calculate_metrics <- function(res, true_taus, se_source = "boot_se") {
  do.call(rbind, lapply(names(res[[1]]$point_est), function(nm) {
    est_mat  <- do.call(rbind, lapply(res, function(r) r$point_est[[nm]]))
    se_mat   <- do.call(rbind, lapply(res, function(r) {
      if (se_source == "sandwich_se") {
        sapply(r$sandwich_se[[nm]], `[[`, "se")
      } else {
        r$boot_se[[nm]]
      }
    }))
    true_tau <- true_taus[[nm]]
    bias_mat <- sweep(est_mat, 2, true_tau, "-")
    lower    <- est_mat - 1.959964 * se_mat
    upper    <- est_mat + 1.959964 * se_mat
    covered  <- sweep(lower, 2, true_tau, "<=") & sweep(upper, 2, true_tau, ">=")
    avg_estimated_se <- colMeans(se_mat, na.rm = TRUE)
    
    data.frame(
      Method     = nm,
      Treatment  = names(true_tau),
      Bias       = round(colMeans(bias_mat), 4),
      ARBias_pct = round(abs(colMeans(bias_mat) / true_tau) * 100, 2),
      Avg_Estimated_SE = round(avg_estimated_se, 4),
      RMSE       = round(sqrt(colMeans(bias_mat^2)), 4),
      CP         = round(colMeans(covered, na.rm = T), 4),
      N_NA_SE      = colSums(is.na(se_mat)),
      row.names  = NULL
    )
  }))
}