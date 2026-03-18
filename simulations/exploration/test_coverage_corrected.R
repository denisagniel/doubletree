# test_coverage_corrected.R
# Coverage test with correctly specified DGPs

suppressMessages({
  devtools::load_all("/Users/dagniel/OneDrive - RAND Corporation/rprojects/global-scholars/optimaltrees")
  devtools::load_all()
  source("simulations/dgps_att_correct.R")
})

cat("===== Coverage Test with Corrected DGPs =====\n\n")

# Parameters
n <- 400
tau <- 0.10  # Using 0.10 for better signal-to-noise
K <- 5
n_reps <- 100

lambda_fixed <- log(n) / n

cat(sprintf("DGP: binary_att (4 features, controlled ATT)\n"))
cat(sprintf("Sample size: n=%d\n", n))
cat(sprintf("True ATT: τ=%.2f\n", tau))
cat(sprintf("Fixed λ: log(n)/n = %.4f\n", lambda_fixed))
cat(sprintf("DML folds: K=%d\n", K))
cat(sprintf("Replications: %d\n\n", n_reps))

# Run simulation
results <- vector("list", n_reps)
cat("Running replications")

for (i in 1:n_reps) {
  if (i %% 10 == 0) cat(".")

  d <- generate_dgp_binary_att(n, tau = tau, seed = 1000 + i)

  capture.output({
    results[[i]] <- tryCatch({
      doubletree::estimate_att(
        d$X, d$A, d$Y,
        K = K,
        regularization = lambda_fixed,
        cv_regularization = FALSE,
        verbose = FALSE
      )
    }, error = function(e) {
      list(theta = NA, sigma = NA, ci = c(NA, NA), error = e$message)
    })
  }, file = tempfile())
}

cat(" done!\n\n")

# Extract results
theta_hats <- sapply(results, function(r) r$theta)
sigma_hats <- sapply(results, function(r) r$sigma)
cis_lower <- sapply(results, function(r) r$ci[1])
cis_upper <- sapply(results, function(r) r$ci[2])

# Remove failed replications
valid <- !is.na(theta_hats)
n_valid <- sum(valid)
n_failed <- sum(!valid)

if (n_failed > 0) {
  cat(sprintf("Warning: %d replications failed\n\n", n_failed))
}

theta_hats <- theta_hats[valid]
sigma_hats <- sigma_hats[valid]
cis_lower <- cis_lower[valid]
cis_upper <- cis_upper[valid]

# Compute statistics
bias <- mean(theta_hats) - tau
coverage <- mean(cis_lower <= tau & cis_upper >= tau)
mean_width <- mean(cis_upper - cis_lower)
rmse <- sqrt(mean((theta_hats - tau)^2))

# Display results
cat("===== RESULTS =====\n")
cat(sprintf("Valid replications: %d/%d\n", n_valid, n_reps))
cat(sprintf("\nEstimation:\n"))
cat(sprintf("  True ATT:      %.4f\n", tau))
cat(sprintf("  Mean estimate: %.4f\n", mean(theta_hats)))
cat(sprintf("  Bias:          %.4f (%.1f%% of truth)\n", bias, 100 * bias / tau))
cat(sprintf("  RMSE:          %.4f\n", rmse))
cat(sprintf("\nInference:\n"))
cat(sprintf("  Coverage:      %.1f%% (target: 95%%)\n", coverage * 100))
cat(sprintf("  Mean CI width: %.4f\n", mean_width))
cat(sprintf("  Mean SE:       %.4f\n", mean(sigma_hats)))

# Evaluation
cat("\n===== EVALUATION =====\n")
if (coverage >= 0.92 && coverage <= 0.98 && abs(bias) < 0.02) {
  cat("✓✓ SUCCESS! Valid inference achieved.\n")
  cat("\nKey findings:\n")
  cat(sprintf("  - Coverage: %.1f%% (within 92-98%%)\n", coverage * 100))
  cat(sprintf("  - Bias: %.1f%% of truth (acceptable)\n", 100 * abs(bias) / tau))
  cat("\nNext steps:\n")
  cat("  1. Test CV-based regularization\n")
  cat("  2. Test continuous-feature DGP\n")
  cat("  3. Document findings\n")
} else if (coverage >= 0.90) {
  cat("⚠ PARTIAL SUCCESS\n")
  cat(sprintf("  Coverage %.1f%% is acceptable but not ideal\n", coverage * 100))
} else {
  cat("✗ FAILED\n")
  cat("  Coverage or bias outside acceptable range\n")
}
