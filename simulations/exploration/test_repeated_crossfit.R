# test_repeated_crossfit.R
# Test if repeated cross-fitting fixes the coverage problem

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools required")
devtools::load_all()
if (!requireNamespace("optimaltrees", quietly = TRUE)) stop("optimaltrees required")

source("simulations/run_simulations_extended.R", local = TRUE)

message("=== Testing Repeated Cross-Fitting with CV Regularization ===\n")
message("Comparing single vs repeated cross-fitting on coverage\n")
message("Using CV to select regularization parameter\n")

# Test parameters
n <- 400  # Use n=400 (moderate sample size)
tau <- 0.15
n_reps <- 50  # Fewer reps for speed
K <- 5

# Test DGP1 (smooth, binary features)
dgp_fn <- generate_data_dgp1

results_single <- list()
results_repeated <- list()

message("Running ", n_reps, " replications...")
message("  Single cross-fit (current approach)")
for (i in seq_len(n_reps)) {
  if (i %% 10 == 0) message("    Rep ", i, "/", n_reps)

  d <- dgp_fn(n, tau, seed = i * 100)

  results_single[[i]] <- estimate_att(
    d$X, d$A, d$Y, K = K,
    use_rashomon = FALSE,
    cv_regularization = TRUE,
    cv_K = 5,
    verbose = FALSE,
    seed = i * 100
  )
}

message("\n  Repeated cross-fit (M=10 splits)")
for (i in seq_len(n_reps)) {
  if (i %% 10 == 0) message("    Rep ", i, "/", n_reps)

  d <- dgp_fn(n, tau, seed = i * 100)

  results_repeated[[i]] <- dml_att_repeated(
    d$X, d$A, d$Y, K = K,
    use_rashomon = FALSE,
    cv_regularization = TRUE,
    cv_K = 5,
    verbose = FALSE,
    seed = i * 100,
    n_splits = 10,
    aggregation = "median"
  )
}

# Compute metrics
compute_coverage <- function(results, true_att, n) {
  theta_hats <- vapply(results, function(r) r$theta, numeric(1))
  ci_low <- vapply(results, function(r) r$ci_95[1L], numeric(1))
  ci_high <- vapply(results, function(r) r$ci_95[2L], numeric(1))
  sigma_hats <- vapply(results, function(r) r$sigma, numeric(1))

  bias <- mean(theta_hats) - true_att
  empirical_se <- sd(theta_hats)
  # sigma is on sqrt(n) scale, convert to theta scale
  mean_estimated_se <- mean(sigma_hats / sqrt(n))
  coverage <- mean(ci_low <= true_att & true_att <= ci_high)
  mean_ci_width <- mean(ci_high - ci_low)

  data.frame(
    bias = bias,
    empirical_se = empirical_se,
    mean_estimated_se = mean_estimated_se,
    se_ratio = mean_estimated_se / empirical_se,
    coverage = coverage,
    mean_ci_width = mean_ci_width
  )
}

metrics_single <- compute_coverage(results_single, tau, n)
metrics_repeated <- compute_coverage(results_repeated, tau, n)

message("\n=== Results ===\n")
message("Single Cross-Fit (K=5, M=1):")
message(sprintf("  Empirical SE:    %.4f", metrics_single$empirical_se))
message(sprintf("  Mean Estimated SE: %.4f", metrics_single$mean_estimated_se))
message(sprintf("  SE Ratio:        %.3f (should be ~1.0)", metrics_single$se_ratio))
message(sprintf("  Coverage:        %.3f (should be ~0.95)", metrics_single$coverage))
message(sprintf("  Mean CI Width:   %.4f", metrics_single$mean_ci_width))

message("\nRepeated Cross-Fit (K=5, M=10):")
message(sprintf("  Empirical SE:    %.4f", metrics_repeated$empirical_se))
message(sprintf("  Mean Estimated SE: %.4f", metrics_repeated$mean_estimated_se))
message(sprintf("  SE Ratio:        %.3f (should be ~1.0)", metrics_repeated$se_ratio))
message(sprintf("  Coverage:        %.3f (should be ~0.95)", metrics_repeated$coverage))
message(sprintf("  Mean CI Width:   %.4f", metrics_repeated$mean_ci_width))

# Variance decomposition for repeated
if (n_reps > 0 && length(results_repeated) > 0) {
  within_frac <- mean(vapply(results_repeated, function(r) r$within_var_frac, numeric(1)))
  between_frac <- mean(vapply(results_repeated, function(r) r$between_var_frac, numeric(1)))

  message("\nVariance Decomposition (Repeated):")
  message(sprintf("  Within-fold:     %.1f%%", within_frac * 100))
  message(sprintf("  Between-splits:  %.1f%%", between_frac * 100))
}

message("\n=== Conclusion ===")
if (metrics_repeated$coverage > metrics_single$coverage) {
  improvement <- (metrics_repeated$coverage - metrics_single$coverage) * 100
  message(sprintf("Coverage improved by %.1f percentage points", improvement))
  if (metrics_repeated$coverage >= 0.94) {
    message("✓ Repeated cross-fitting FIXES the coverage problem!")
  } else {
    message("⚠ Coverage improved but still below nominal level")
  }
} else {
  message("✗ Repeated cross-fitting did not improve coverage")
}

message("\nIf successful, rerun simulations with n_splits=10 to get proper coverage.")
