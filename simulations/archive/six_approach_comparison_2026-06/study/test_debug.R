#!/usr/bin/env Rscript

# Debug script - test approach 1 with simple DGP

cat("Loading packages from source...\n")
suppressPackageStartupMessages({
  devtools::load_all("../../../optimaltrees", quiet = TRUE)
  devtools::load_all("../..", quiet = TRUE)
})

cat("Loading simulation code...\n")
source("code/dgps.R")
source("code/estimators.R")

# Test approach 1 with DGP 1
set.seed(12345)
data <- generate_dgp_simple(n = 500)

cat("\nData summary:\n")
cat(sprintf("  n = %d\n", nrow(data$X)))
cat(sprintf("  n_treated = %d\n", sum(data$A)))
cat(sprintf("  n_control = %d\n", sum(1 - data$A)))
cat(sprintf("  Y range: [%.3f, %.3f]\n", min(data$Y), max(data$Y)))
cat(sprintf("  Y mean: %.3f\n", mean(data$Y)))
cat(sprintf("  True ATT: %.3f\n", data$true_att))

# Fit propensity tree
cat("\nFitting propensity tree...\n")
e_tree <- optimaltrees::fit_tree(
  X = data$X,
  y = data$A,
  loss_function = "log_loss",
  regularization = 0.1,
  verbose = TRUE
)

cat("\nPropensity predictions:\n")
e_hat <- predict(e_tree, data$X)
cat(sprintf("  e_hat range: [%.3f, %.3f]\n", min(e_hat), max(e_hat)))
cat(sprintf("  e_hat mean: %.3f\n", mean(e_hat)))
cat(sprintf("  Any NA: %s\n", any(is.na(e_hat))))
cat(sprintf("  Any NaN: %s\n", any(is.nan(e_hat))))
cat(sprintf("  Any <0: %s\n", any(e_hat < 0)))
cat(sprintf("  Any >1: %s\n", any(e_hat > 1)))

# Fit outcome tree
cat("\nFitting outcome tree (controls only)...\n")
control_idx <- data$A == 0
m0_tree <- optimaltrees::fit_tree(
  X = data$X[control_idx, , drop = FALSE],
  y = data$Y[control_idx],
  loss_function = "log_loss",
  regularization = 0.1,
  verbose = TRUE
)

cat("\nOutcome predictions:\n")
m0_hat <- predict(m0_tree, data$X)
cat(sprintf("  m0_hat range: [%.3f, %.3f]\n", min(m0_hat), max(m0_hat)))
cat(sprintf("  m0_hat mean: %.3f\n", mean(m0_hat)))
cat(sprintf("  Any NA: %s\n", any(is.na(m0_hat))))
cat(sprintf("  Any NaN: %s\n", any(is.nan(m0_hat))))

# Compute EIF
cat("\nComputing EIF...\n")
n <- length(data$Y)
pi_hat <- mean(data$A)
cat(sprintf("  pi_hat = %.3f\n", pi_hat))

psi <- (data$A / pi_hat) * (data$Y - m0_hat) +
       ((1 - data$A) * e_hat) / (pi_hat * (1 - e_hat)) * (data$Y - m0_hat)

cat(sprintf("  psi range: [%.3f, %.3f]\n", min(psi, na.rm = TRUE), max(psi, na.rm = TRUE)))
cat(sprintf("  psi mean: %.3f\n", mean(psi, na.rm = TRUE)))
cat(sprintf("  Any NA: %s\n", any(is.na(psi))))
cat(sprintf("  Any NaN: %s\n", any(is.nan(psi))))
cat(sprintf("  Any Inf: %s\n", any(is.infinite(psi))))

theta_hat <- mean(psi, na.rm = TRUE)
cat(sprintf("\ntheta_hat = %.4f\n", theta_hat))
