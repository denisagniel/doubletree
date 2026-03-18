# Quick test: Verify fixes work (simple version)
# Avoids package loading issues by sourcing directly

library(optimaltrees)
library(dplyr)

cat("\n=== Quick Test: model_limit Fix + Threshold Encoding ===\n\n")

# Source dmltree functions directly
source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/nuisance_trees.R")
source("../../R/estimate_att.R")

# Source DGPs
source("dgps/dgps_beta_continuous.R")

# Test configuration
N <- 800
TAU <- 0.10
K_FOLDS <- 5
N_REPS <- 10
SEED_OFFSET <- 30000

cat("Test configuration:\n")
cat("  Sample size: n =", N, "\n")
cat("  Replications:", N_REPS, "\n")
cat("  DGP: β=3 (cubic polynomials, valid regime)\n")
cat("  Expected: 100% convergence, ~95% coverage\n\n")

# Run simulations
results_list <- vector("list", N_REPS)

cat("Running simulations...\n")
start_time <- Sys.time()

for (rep in 1:N_REPS) {
  cat(sprintf("  Rep %d/%d...", rep, N_REPS))

  seed <- SEED_OFFSET + rep
  d <- generate_dgp_beta_high(n = N, tau = TAU, seed = seed)

  result <- tryCatch({
    fit <- estimate_att(
      X = d$X, A = d$A, Y = d$Y,
      K = K_FOLDS,
      outcome_type = "binary",
      regularization = log(N) / N,
      cv_regularization = FALSE,
      use_rashomon = FALSE,
      discretize_method = "quantiles",
      discretize_bins = "adaptive",
      verbose = FALSE
    )

    # Extract predictions to check they vary
    fold_1_model <- fit$nuisance_fits[[1]]
    test_pred <- tryCatch({
      pred <- predict(fold_1_model$e_model, d$X[1:100, ], type = "prob")
      if (is.matrix(pred)) sd(pred[, 2]) else 0
    }, error = function(e) 0)

    cat(" OK\n")

    list(
      theta = fit$theta,
      sigma = fit$sigma,
      ci_lower = fit$ci_95[1],
      ci_upper = fit$ci_95[2],
      converged = TRUE,
      pred_sd = test_pred,
      error = NA
    )
  }, error = function(e) {
    error_msg <- conditionMessage(e)
    cat(sprintf(" FAILED: %s\n", substr(error_msg, 1, 60)))

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
    error = as.character(result$error),
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, results_list)

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat(sprintf("\nTest complete in %.1f seconds (%.1f sec/rep)\n\n",
            elapsed, elapsed/N_REPS))

# Check for errors
model_limit_errors <- sum(grepl("model.limit|model_limit", results$error, ignore.case = TRUE), na.rm = TRUE)
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
  mean_pred_sd <- mean(converged$pred_sd, na.rm = TRUE)
} else {
  coverage <- bias <- rmse <- mean_se <- mean_pred_sd <- NA
}

# Print results
cat("=== TEST RESULTS ===\n\n")

cat("Infrastructure checks:\n")
cat(sprintf("  model_limit errors: %d (expect 0) %s\n",
            model_limit_errors,
            ifelse(model_limit_errors == 0, "✓", "✗")))
cat(sprintf("  Constant predictions: %d (expect 0) %s\n",
            constant_predictions,
            ifelse(constant_predictions == 0, "✓", "✗")))
cat(sprintf("  Convergence rate: %.0f%% (expect 100%%) %s\n",
            convergence_rate,
            ifelse(convergence_rate == 100, "✓", "⚠")))
cat("\n")

if (n_converged > 0) {
  cat("Statistical performance:\n")
  cat(sprintf("  Coverage: %.0f%% (expect ~95%%) %s\n",
              coverage,
              ifelse(coverage >= 80 && coverage <= 100, "✓", "⚠")))
  cat(sprintf("  Bias: %.3f (expect ~0)\n", bias))
  cat(sprintf("  RMSE: %.3f\n", rmse))
  cat(sprintf("  Mean SE: %.3f\n", mean_se))
  cat(sprintf("  Prediction SD: %.3f (expect > 0) %s\n",
              mean_pred_sd,
              ifelse(mean_pred_sd > 0.01, "✓", "✗")))
  cat("\n")

  # Overall assessment
  cat("OVERALL ASSESSMENT:\n")

  all_pass <- TRUE

  if (model_limit_errors > 0) {
    cat("  ✗ FAIL: model_limit errors detected\n")
    all_pass <- FALSE
  }

  if (constant_predictions > 0) {
    cat("  ✗ FAIL: Constant predictions (predict bug not fixed)\n")
    all_pass <- FALSE
  }

  if (convergence_rate < 90) {
    cat("  ✗ FAIL: Low convergence rate\n")
    all_pass <- FALSE
  }

  if (coverage < 80 || coverage > 100) {
    cat("  ⚠ WARNING: Coverage outside acceptable range (80-100%)\n")
    if (coverage >= 70) {
      cat("    (Still reasonable for small test, may improve with more reps)\n")
    }
  }

  if (all_pass) {
    cat("  ✓✓ PASS: All infrastructure checks successful!\n")
    cat("    - No model_limit errors\n")
    cat("    - Predictions vary correctly (not constant)\n")
    cat("    - Convergence excellent\n")
    if (coverage >= 80 && coverage <= 100) {
      cat("    - Coverage good\n")
    }
    cat("\n  Ready for full simulation!\n")
  }
} else {
  cat("✗ FAIL: No successful replications\n\n")
  cat("Errors encountered:\n")
  failed <- results[!results$converged, ]
  for (i in 1:nrow(failed)) {
    cat(sprintf("  Rep %d: %s\n", failed$rep[i],
                substr(failed$error[i], 1, 100)))
  }
}

cat("\n")

# Save results
dir.create("results", showWarnings = FALSE, recursive = TRUE)
output_file <- sprintf("results/quick_test_%s.rds", format(Sys.time(), "%Y%m%d_%H%M"))
saveRDS(results, output_file)
cat(sprintf("Results saved to: %s\n", output_file))
