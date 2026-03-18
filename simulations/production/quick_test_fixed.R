# Quick test: Verify fixes work
# Tests model_limit fix + threshold encoding with 10 replications

library(optimaltrees)
library(dplyr)

cat("\n=== Quick Test: model_limit Fix + Threshold Encoding ===\n\n")

# Load dmltree
devtools::load_all("../..")

# Source DGPs with continuous features
source("dgps/dgps_beta_continuous.R")

# Small test grid
N_VALUES <- c(800)
TAU <- 0.10
K_FOLDS <- 5
N_REPS <- 10  # Quick test
SEED_OFFSET <- 30000

DGPS <- list(
  beta_high = generate_dgp_beta3_continuous  # Just test one regime for speed
)

cat("Test configuration:\n")
cat("  Sample size: n =", N_VALUES, "\n")
cat("  Replications:", N_REPS, "\n")
cat("  DGP: β=3 (cubic, valid regime)\n")
cat("  Expected: 100% convergence, ~95% coverage\n\n")

# Run simulations
results_list <- vector("list", N_REPS)

cat("Running simulations...\n")
start_time <- Sys.time()

for (rep in 1:N_REPS) {
  cat(sprintf("  Rep %d/%d\r", rep, N_REPS))

  seed <- SEED_OFFSET + rep
  d <- generate_dgp_beta3_continuous(n = N_VALUES, tau = TAU, seed = seed)

  result <- tryCatch({
    fit <- estimate_att(
      X = d$X, A = d$A, Y = d$Y,
      K = K_FOLDS,
      outcome_type = "binary",
      regularization = log(N_VALUES) / N_VALUES,
      cv_regularization = FALSE,
      use_rashomon = FALSE,
      discretize_method = "quantiles",
      discretize_bins = "adaptive",
      verbose = FALSE
    )

    # Check for constant predictions (the old bug)
    e_vals <- sapply(fit$nuisance_fits, function(f) {
      pred <- predict(f$e_model, d$X, type = "prob")
      if (is.matrix(pred)) pred[, 2] else rep(0.5, nrow(d$X))
    })

    pred_sd <- sd(as.vector(e_vals))

    list(
      theta = fit$theta,
      sigma = fit$sigma,
      ci_lower = fit$ci_95[1],
      ci_upper = fit$ci_95[2],
      converged = TRUE,
      pred_sd = pred_sd,  # Should be > 0 (not constant)
      error = NA
    )
  }, error = function(e) {
    error_msg <- conditionMessage(e)
    list(
      theta = NA,
      sigma = NA,
      ci_lower = NA,
      ci_upper = NA,
      converged = FALSE,
      pred_sd = NA,
      error = error_msg
    )
  })

  results_list[[rep]] <- data.frame(
    rep = rep,
    true_att = d$true_att,
    theta = result$theta,
    sigma = result$sigma,
    ci_lower = result$ci_lower,
    ci_upper = result$ci_upper,
    converged = result$converged,
    pred_sd = result$pred_sd,
    error = result$error,
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, results_list)

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat(sprintf("\n\nTest complete in %.1f seconds (%.1f sec/rep)\n\n",
            elapsed, elapsed/N_REPS))

# Check for model_limit errors
model_limit_errors <- sum(grepl("model.limit|model_limit", results$error, ignore.case = TRUE), na.rm = TRUE)

# Check for constant predictions (old predict bug)
constant_predictions <- sum(results$pred_sd < 1e-6, na.rm = TRUE)

# Calculate statistics
converged <- results[results$converged, ]
n_converged <- nrow(converged)
convergence_rate <- 100 * n_converged / N_REPS

if (n_converged > 0) {
  coverage <- 100 * mean(converged$ci_lower <= converged$true_att &
                         converged$ci_upper >= converged$true_att)
  bias <- mean(converged$theta - converged$true_att)
  rmse <- sqrt(mean((converged$theta - converged$true_att)^2))
  mean_se <- mean(converged$sigma)
} else {
  coverage <- NA
  bias <- NA
  rmse <- NA
  mean_se <- NA
}

# Print results
cat("=== TEST RESULTS ===\n\n")

cat("Infrastructure checks:\n")
cat(sprintf("  ✓ model_limit errors: %d (expect 0)\n", model_limit_errors))
cat(sprintf("  ✓ Constant predictions: %d (expect 0)\n", constant_predictions))
cat(sprintf("  ✓ Convergence rate: %.0f%% (expect 100%%)\n", convergence_rate))
cat("\n")

if (n_converged > 0) {
  cat("Statistical performance:\n")
  cat(sprintf("  Coverage: %.0f%% (expect ~95%%)\n", coverage))
  cat(sprintf("  Bias: %.3f (expect ~0)\n", bias))
  cat(sprintf("  RMSE: %.3f\n", rmse))
  cat(sprintf("  Mean SE: %.3f\n", mean_se))
  cat(sprintf("  Prediction SD: %.3f (expect > 0)\n", mean(converged$pred_sd, na.rm = TRUE)))
  cat("\n")

  # Overall assessment
  cat("OVERALL ASSESSMENT:\n")

  if (model_limit_errors > 0) {
    cat("  ✗ FAIL: model_limit errors detected\n")
  } else if (constant_predictions > 0) {
    cat("  ✗ FAIL: Constant predictions (predict bug not fixed)\n")
  } else if (convergence_rate < 90) {
    cat("  ✗ FAIL: Low convergence rate\n")
  } else if (coverage < 80 || coverage > 100) {
    cat("  ⚠ WARNING: Coverage outside acceptable range\n")
  } else {
    cat("  ✓ PASS: All checks successful!\n")
    cat("    - No model_limit errors\n")
    cat("    - Predictions vary correctly\n")
    cat("    - Convergence excellent\n")
    cat("    - Coverage reasonable\n")
  }
} else {
  cat("✗ FAIL: No successful replications\n")
  cat("\nErrors:\n")
  print(results[!results$converged, c("rep", "error")])
}

cat("\n")

# Save results
output_file <- sprintf("results/quick_test_fixed_%s.rds", Sys.Date())
dir.create("results", showWarnings = FALSE)
saveRDS(results, output_file)
cat(sprintf("Results saved to: %s\n", output_file))
