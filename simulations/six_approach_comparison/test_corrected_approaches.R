# Test Corrected Approaches (4 and 6) - Averaging Implementation
# Created: 2026-05-20
#
# Quick local test to verify approaches 4 and 6 work with averaging

library(optimaltrees)

# Load development version of doubletree (not installed version)
# This ensures we test the latest fixes
devtools::load_all("/Users/dagniel/RAND/rprojects/global-scholars/doubletree", quiet = TRUE)

# Source simulation code
setwd("/Users/dagniel/RAND/rprojects/global-scholars/doubletree/simulations/six_approach_comparison")
source("code/dgps.R")
source("code/estimators.R")
# Note: tree_averaging functions now in doubletree package (estimate_att_*_averaged)

cat("=====================================\n")
cat("Testing Corrected Approaches 4 and 6\n")
cat("=====================================\n\n")

# Set seed for reproducibility
set.seed(12345)

# Generate small test data
cat("Generating test data (n=500)...\n")
data <- generate_dgp_simple(n = 500)
cat("  X: ", nrow(data$X), " x ", ncol(data$X), "\n")
cat("  A: ", sum(data$A == 1), " treated, ", sum(data$A == 0), " control\n")
cat("  Y: ", sum(data$Y == 1), " / ", length(data$Y), " (",
    round(100*mean(data$Y), 1), "%)\n")
cat("  True ATT: ", data$true_att, "\n\n")

# Test all 6 approaches
approaches <- list(
  list(num = 1, name = "Full-sample", fun = estimate_att_fullsample),
  list(num = 2, name = "Cross-fit", fun = estimate_att_crossfit),
  list(num = 3, name = "Doubletree", fun = estimate_att_doubletree),
  list(num = 4, name = "Doubletree averaged", fun = estimate_att_doubletree_averaged),
  list(num = 5, name = "M-split", fun = estimate_att_msplit),
  list(num = 6, name = "M-split averaged", fun = estimate_att_msplit_averaged)
)

results <- list()

for (app in approaches) {
  cat("-----------------------------------\n")
  cat(sprintf("Approach %d: %s\n", app$num, app$name))
  cat("-----------------------------------\n")

  result <- tryCatch({
    start_time <- Sys.time()

    # Run estimator
    if (app$num %in% c(5, 6)) {
      # M-split approaches (smaller M for testing)
      res <- app$fun(data$X, data$A, data$Y, M = 5, K = 3, regularization = 0.1)
    } else if (app$num %in% c(1)) {
      # Full-sample (no K parameter)
      res <- app$fun(data$X, data$A, data$Y, regularization = 0.1)
    } else {
      # Standard cross-fit approaches
      res <- app$fun(data$X, data$A, data$Y, K = 3, regularization = 0.1)
    }

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

    # Check for errors
    if (!is.null(res$error)) {
      cat("  ERROR: ", res$error, "\n\n")
      list(success = FALSE, error = res$error, time = elapsed)
    } else {
      # Success - print results
      cat(sprintf("  Theta:  %.4f (true: %.4f)\n", res$theta, data$true_att))
      cat(sprintf("  SE:     %.4f\n", res$se))
      cat(sprintf("  95%% CI: [%.4f, %.4f]\n",
                  res$theta - 1.96*res$se, res$theta + 1.96*res$se))
      cat(sprintf("  Bias:   %.4f (%.1f%% of SE)\n",
                  res$theta - data$true_att,
                  100*(res$theta - data$true_att)/res$se))
      cat(sprintf("  Time:   %.2f seconds\n", elapsed))

      # Additional diagnostics for approaches 4 and 6
      if (app$num == 4 && !is.null(res$averaged_trees)) {
        cat("  Averaged trees: propensity + outcome (from K=3 folds)\n")
      }
      if (app$num == 6 && !is.null(res$n_trees_averaged)) {
        cat(sprintf("  Trees averaged: e=%d, m0=%d (from M×K=5×3=15)\n",
                    res$n_trees_averaged$e, res$n_trees_averaged$m0))
      }

      cat("\n")
      list(success = TRUE, theta = res$theta, se = res$se,
           bias = res$theta - data$true_att, time = elapsed)
    }
  }, error = function(e) {
    cat("  EXCEPTION: ", e$message, "\n\n")
    list(success = FALSE, error = e$message, time = NA)
  })

  results[[app$num]] <- result
}

# Summary
cat("=====================================\n")
cat("Summary\n")
cat("=====================================\n\n")

cat(sprintf("%-25s %10s %10s %10s %8s\n", "Approach", "Theta", "SE", "Bias", "Time(s)"))
cat(sprintf("%s\n", paste(rep("-", 65), collapse = "")))

for (app in approaches) {
  res <- results[[app$num]]
  if (res$success) {
    cat(sprintf("%-25s %10.4f %10.4f %10.4f %8.2f\n",
                app$name, res$theta, res$se, res$bias, res$time))
  } else {
    cat(sprintf("%-25s %10s %10s %10s %8s\n",
                app$name, "FAILED", "-", "-", "-"))
  }
}

cat("\n")

# Check success
n_success <- sum(sapply(results, function(r) r$success))
cat(sprintf("Success: %d / %d approaches\n", n_success, length(approaches)))

if (n_success == length(approaches)) {
  cat("\n✓ All approaches working!\n\n")

  # Compare approaches 4 and 6 to baseline
  if (results[[2]]$success && results[[4]]$success) {
    theta_diff <- results[[4]]$theta - results[[2]]$theta
    cat(sprintf("Approach 4 vs 2 (baseline): Δtheta = %.4f (%.1f%% difference)\n",
                theta_diff, 100*abs(theta_diff)/results[[2]]$theta))
  }

  if (results[[5]]$success && results[[6]]$success) {
    theta_diff <- results[[6]]$theta - results[[5]]$theta
    cat(sprintf("Approach 6 vs 5 (baseline): Δtheta = %.4f (%.1f%% difference)\n",
                theta_diff, 100*abs(theta_diff)/results[[5]]$theta))
  }
} else {
  cat("\n✗ Some approaches failed. See details above.\n")
}

cat("\n")
cat("Test complete.\n")
