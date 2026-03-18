# Micro Test: Verify Predict Bug Fix (Simplified)
# 3 reps, ~2 minutes, ~300MB memory
# Uses accessible DGP file to test predictions vary (not constant 0.5)

library(optimaltrees)
library(dplyr)

cat("\n=== MICRO TEST: Verify Predict Bug Fix ===\n\n")

# Source dmltree
source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/nuisance_trees.R")
source("../../R/estimate_att.R")

# Source accessible DGPs
source("dgps/dgps_smooth.R")

N <- 800
TAU <- 0.10
K_FOLDS <- 5
N_REPS <- 3
SEED_OFFSET <- 60000

cat("Configuration:\n")
cat("  n =", N, "\n")
cat("  Replications:", N_REPS, "\n")
cat("  DGP: Binary (smooth nuisances)\n")
cat("  Memory check: ON\n\n")

# Initial memory
mem_start <- gc()
mem_start_mb <- sum(mem_start[, "used"]) * 0.001
cat(sprintf("Starting memory: %.0f MB\n\n", mem_start_mb))

results_list <- vector("list", N_REPS)
start_time <- Sys.time()

for (rep in 1:N_REPS) {
  cat(sprintf("Rep %d/%d...", rep, N_REPS))

  # Use simple binary DGP
  d <- generate_dgp_binary_att(n = N, tau = TAU, seed = SEED_OFFSET + rep)

  result <- tryCatch({
    fit <- estimate_att(
      X = d$X, A = d$A, Y = d$Y,
      K = K_FOLDS,
      outcome_type = "binary",
      regularization = log(N) / N,
      cv_regularization = FALSE,
      use_rashomon = FALSE,
      verbose = FALSE
    )

    cat(sprintf(" OK (theta=%.3f)\n", fit$theta))

    list(
      theta = fit$theta,
      sigma = fit$sigma,
      ci_lower = fit$ci_95[1],
      ci_upper = fit$ci_95[2],
      converged = TRUE,
      error = NA,
      # Store predictions to check they vary
      propensity_sd = if(!is.null(fit$nuisance_fits$propensity)) sd(fit$nuisance_fits$propensity) else NA,
      outcome_sd = if(!is.null(fit$nuisance_fits$outcome_control)) sd(fit$nuisance_fits$outcome_control) else NA
    )
  }, error = function(e) {
    cat(" FAILED\n")
    cat(sprintf("  Error: %s\n", substr(conditionMessage(e), 1, 100)))
    list(
      theta = NA, sigma = NA, ci_lower = NA, ci_upper = NA,
      converged = FALSE, error = conditionMessage(e),
      propensity_sd = NA, outcome_sd = NA
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
    propensity_sd = result$propensity_sd,
    outcome_sd = result$outcome_sd,
    stringsAsFactors = FALSE
  )

  # Force garbage collection
  if (rep %% 2 == 0) gc(verbose = FALSE)
}

results <- do.call(rbind, results_list)
elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

# Final memory
mem_end <- gc()
mem_end_mb <- sum(mem_end[, "used"]) * 0.001
mem_used <- mem_end_mb - mem_start_mb

cat(sprintf("\nCompleted in %.1f seconds (%.1f sec/rep)\n", elapsed, elapsed/N_REPS))
cat(sprintf("Memory used: %.0f MB (%.0f MB/rep)\n\n", mem_used, mem_used/N_REPS))

# Results
converged <- results[results$converged, ]
n_conv <- nrow(converged)

cat("=== RESULTS ===\n\n")
cat(sprintf("Convergence: %d/%d (%.0f%%)\n\n", n_conv, N_REPS, 100*n_conv/N_REPS))

if (n_conv > 0) {
  # Check 1: Predictions vary (not constant 0.5)
  cat("1. PREDICT BUG CHECK:\n")
  prop_sd <- mean(converged$propensity_sd, na.rm=TRUE)
  out_sd <- mean(converged$outcome_sd, na.rm=TRUE)

  cat(sprintf("   Propensity SD: %.4f (expect > 0.01)\n", prop_sd))
  cat(sprintf("   Outcome SD:    %.4f (expect > 0.01)\n", out_sd))

  predict_ok <- (!is.na(prop_sd) && prop_sd > 0.01) || (!is.na(out_sd) && out_sd > 0.01)

  if (predict_ok) {
    cat("   Status: PASS - Predictions vary correctly\n\n")
  } else {
    cat("   Status: FAIL - Predictions still constant!\n\n")
  }

  # Check 2: Coverage
  coverage <- 100 * mean(converged$ci_lower <= converged$true_att &
                         converged$ci_upper >= converged$true_att)
  bias <- mean(converged$theta - converged$true_att)
  rmse <- sqrt(mean((converged$theta - converged$true_att)^2))

  cat("2. STATISTICAL CHECKS:\n")
  cat(sprintf("   Coverage: %.0f%% (expect ~95%%, but only %d reps)\n", coverage, n_conv))
  cat(sprintf("   Bias:     %.4f\n", bias))
  cat(sprintf("   RMSE:     %.4f\n\n", rmse))

  # Overall verdict
  cat("=== VERDICT ===\n\n")

  if (n_conv == N_REPS && predict_ok) {
    cat("PASS: All tests successful!\n")
    cat("  - 100% convergence\n")
    cat("  - Predictions varying correctly\n")
    cat(sprintf("  - Memory usage reasonable (%.0f MB/rep)\n", mem_used/N_REPS))
    cat("\nReady to proceed with larger simulations.\n")
  } else if (predict_ok && n_conv >= 2) {
    cat("MOSTLY PASS:\n")
    cat("  - Predictions working correctly\n")
    if (n_conv < N_REPS) {
      cat(sprintf("  - Some convergence failures (%d/%d)\n", N_REPS-n_conv, N_REPS))
    }
  } else if (!predict_ok) {
    cat("FAIL: Predict bug still present!\n")
    cat("  - Predictions are constant (SD < 0.01)\n")
    cat("  - Need to investigate predict.R fix\n")
  } else {
    cat("FAIL: Too many convergence failures\n")
  }

} else {
  cat("FAIL: NO SUCCESSFUL RUNS\n\n")
  cat("Errors:\n")
  for (i in 1:nrow(results)) {
    if (!results$converged[i]) {
      err_msg <- substr(results$error[i], 1, 150)
      cat(sprintf("  Rep %d: %s\n", i, err_msg))
    }
  }
}

cat("\n")
