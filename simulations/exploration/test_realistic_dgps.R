# test_realistic_dgps.R
# Quick test: Do realistic DGPs make CV informative?

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools required")
devtools::load_all()
if (!requireNamespace("optimaltrees", quietly = TRUE)) stop("optimaltrees required")

source("simulations/dgps_realistic.R")

message("=== Testing Realistic DGPs for CV Regularization ===\n")

n <- 400
tau <- 0.15
K <- 5
n_reps <- 10  # Quick test

dgps <- list(
  simple = generate_dgp_simple,
  moderate = generate_dgp_moderate,
  complex = generate_dgp_complex,
  interactions = generate_dgp_interactions
)

test_dgp <- function(dgp_fn, dgp_name, n, tau, n_reps) {
  message("\n=== ", dgp_name, " ===")

  # Step 1: Check CV informativeness
  d <- dgp_fn(n, tau, seed = 123)
  message(sprintf("  Features: %d (patterns: %d)", ncol(d$X), 2^ncol(d$X)))

  cv_e <- optimaltrees::cv_regularization(d$X, d$A, loss_function = "log_loss",
                                       K = 5, refit = FALSE, verbose = FALSE)
  cv_range <- max(cv_e$cv_loss) - min(cv_e$cv_loss)
  cv_rel_range <- cv_range / mean(cv_e$cv_loss)

  message(sprintf("  CV loss range: %.4f (%.1f%% of mean)", cv_range, cv_rel_range * 100))
  message(sprintf("  Selected lambda: %.5f vs Fixed: %.5f",
                  cv_e$best_lambda, log(n)/n))

  is_informative <- cv_rel_range > 0.01  # >1% variation

  if (!is_informative) {
    message("  ✗ CV NOT informative (flat loss)")
    return(NULL)
  }

  message("  ✓ CV is informative")

  # Step 2: Quick coverage test
  message(sprintf("  Running %d reps...", n_reps))

  results_fixed <- list()
  results_cv <- list()

  for (i in seq_len(n_reps)) {
    if (i %% 5 == 0) message(sprintf("    Rep %d/%d", i, n_reps))
    d <- dgp_fn(n, tau, seed = i * 100)

    # Fixed regularization
    results_fixed[[i]] <- estimate_att(
      d$X, d$A, d$Y, K = K,
      use_rashomon = FALSE,
      regularization = log(n)/n,
      verbose = FALSE,
      seed = i * 100
    )

    # CV regularization
    results_cv[[i]] <- estimate_att(
      d$X, d$A, d$Y, K = K,
      use_rashomon = FALSE,
      cv_regularization = TRUE,
      cv_K = 5,
      verbose = FALSE,
      seed = i * 100
    )
  }

  # Compute coverage
  compute_coverage <- function(results, true_att) {
    ci_low <- sapply(results, function(r) r$ci_95[1])
    ci_high <- sapply(results, function(r) r$ci_95[2])
    mean(ci_low <= true_att & true_att <= ci_high)
  }

  cov_fixed <- compute_coverage(results_fixed, tau)
  cov_cv <- compute_coverage(results_cv, tau)

  message(sprintf("  Coverage - Fixed: %.2f, CV: %.2f", cov_fixed, cov_cv))

  # Check if CV improved coverage
  improvement <- cov_cv - cov_fixed
  if (improvement >= 0.05) {
    message(sprintf("  ✓ CV IMPROVED coverage by %.1f pp", improvement * 100))
  } else if (improvement >= -0.05) {
    message(sprintf("  = CV similar to fixed (%.1f pp diff)", improvement * 100))
  } else {
    message(sprintf("  ✗ CV WORSE than fixed (%.1f pp drop)", abs(improvement) * 100))
  }

  list(
    dgp = dgp_name,
    n_features = ncol(d$X),
    cv_informative = TRUE,
    cv_rel_range = cv_rel_range,
    coverage_fixed = cov_fixed,
    coverage_cv = cov_cv,
    improvement = improvement
  )
}

# Run tests
results <- lapply(names(dgps), function(name) {
  test_dgp(dgps[[name]], name, n, tau, n_reps)
})
results <- Filter(Negate(is.null), results)  # Remove NULLs (non-informative)

message("\n\n=== SUMMARY ===\n")

if (length(results) == 0) {
  message("✗ None of the DGPs produced informative CV!")
  message("Need to redesign DGPs or reconsider approach.")
} else {
  for (r in results) {
    message(sprintf("%s (%d features):", r$dgp, r$n_features))
    message(sprintf("  CV informative: %.1f%% variation", r$cv_rel_range * 100))
    message(sprintf("  Coverage: Fixed=%.2f, CV=%.2f (%.1f pp %s)",
                    r$coverage_fixed, r$coverage_cv,
                    abs(r$improvement) * 100,
                    if (r$improvement > 0) "improvement" else "drop"))
    message("")
  }

  # Find best
  best_idx <- which.max(sapply(results, function(r) r$coverage_cv))
  best <- results[[best_idx]]

  message("\n=== RECOMMENDATION ===")
  message(sprintf("Best DGP: %s", best$dgp))
  message(sprintf("  %d features (%d patterns)", best$n_features, 2^best$n_features))
  message(sprintf("  CV coverage: %.2f", best$coverage_cv))

  if (best$coverage_cv >= 0.92) {
    message("  ✓ Excellent coverage - use this DGP!")
  } else if (best$coverage_cv >= 0.85) {
    message("  ✓ Good coverage - acceptable for tree-based DML")
  } else {
    message("  ⚠ Coverage still below 85% - may need larger n or different approach")
  }
}
