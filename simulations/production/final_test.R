# Final test with all fixes applied
# C++ fix: model_limit=0 now works correctly
# R fix: threshold encoding
# Error checks: use n_trees instead of tree_json

library(optimaltrees)
library(dplyr)

cat("\n=== FINAL TEST: All Fixes Applied ===\n\n")

# Source doubletree
source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/nuisance_trees.R")
source("../../R/estimate_att.R")

# Source DGPs
source("dgps/dgps_beta_continuous.R")

N <- 800
TAU <- 0.10
K_FOLDS <- 5
N_REPS <- 5  # Just 5 reps for speed
SEED_OFFSET <- 40000

cat("Configuration:\n")
cat("  n =", N, "\n")
cat("  Replications:", N_REPS, "\n")
cat("  DGP: β=3 (valid regime)\n\n")

results_list <- vector("list", N_REPS)
start_time <- Sys.time()

for (rep in 1:N_REPS) {
  cat(sprintf("Rep %d/%d...", rep, N_REPS))

  d <- generate_dgp_beta_high(n = N, tau = TAU, seed = SEED_OFFSET + rep)

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

    cat(" OK (theta=", round(fit$theta, 3), ")\n", sep="")

    list(
      theta = fit$theta,
      sigma = fit$sigma,
      ci_lower = fit$ci_95[1],
      ci_upper = fit$ci_95[2],
      converged = TRUE,
      error = NA
    )
  }, error = function(e) {
    cat(" FAILED\n")
    list(
      theta = NA, sigma = NA, ci_lower = NA, ci_upper = NA,
      converged = FALSE, error = conditionMessage(e)
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
    error = as.character(result$error),
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, results_list)
elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat(sprintf("\nCompleted in %.1f seconds (%.1f sec/rep)\n\n", elapsed, elapsed/N_REPS))

# Results
converged <- results[results$converged, ]
n_conv <- nrow(converged)

cat("=== RESULTS ===\n\n")
cat(sprintf("Convergence: %d/%d (%.0f%%)\n", n_conv, N_REPS, 100*n_conv/N_REPS))

if (n_conv > 0) {
  coverage <- 100 * mean(converged$ci_lower <= converged$true_att &
                         converged$ci_upper >= converged$true_att)
  bias <- mean(converged$theta - converged$true_att)
  rmse <- sqrt(mean((converged$theta - converged$true_att)^2))

  cat(sprintf("Coverage: %.0f%% (expect ~95%%)\n", coverage))
  cat(sprintf("Bias: %.4f\n", bias))
  cat(sprintf("RMSE: %.4f\n\n", rmse))

  if (n_conv == N_REPS && coverage >= 80 && coverage <= 100) {
    cat("✓✓ ALL TESTS PASS!\n")
    cat("   - 100% convergence\n")
    cat("   - Coverage in acceptable range\n")
    cat("   - Ready for full simulations!\n")
  } else if (n_conv >= 0.9 * N_REPS) {
    cat("✓ MOSTLY PASSING\n")
    cat("   - Good convergence rate\n")
    if (coverage < 80 || coverage > 100) {
      cat("   ⚠ Coverage needs investigation\n")
    }
  } else {
    cat("⚠ NEEDS ATTENTION\n")
    cat("   - Low convergence rate\n")
  }
} else {
  cat("✗ NO SUCCESSES\n\n")
  cat("Errors:\n")
  for (i in 1:nrow(results)) {
    if (!results$converged[i]) {
      cat(sprintf("  Rep %d: %s\n", i, substr(results$error[i], 1, 80)))
    }
  }
}

# Save
dir.create("results", showWarnings = FALSE)
output_file <- sprintf("results/final_test_%s.rds", format(Sys.time(), "%Y%m%d_%H%M"))
saveRDS(results, output_file)
cat(sprintf("\nResults saved: %s\n", output_file))
