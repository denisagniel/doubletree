#!/usr/bin/env Rscript
# Test all 6 approaches after lambda fix
# Quick sanity check before cluster deployment

library(optimaltrees)
library(doubletree)

# Source simulation code
setwd("doubletree/simulations/six_approach_comparison")
source("code/dgps.R")
source("code/estimators.R")
source("code/metrics.R")

cat("================================================================================\n")
cat("TESTING ALL 6 APPROACHES AFTER LAMBDA FIX\n")
cat("================================================================================\n\n")

cat("Test configuration:\n")
cat("  DGP: Complex (the one that showed bias)\n")
cat("  Sample size: n=500\n")
cat("  True ATT: 0.15\n")
cat("  Replications: 1 per approach (quick test)\n\n")

# Generate test data
set.seed(42)
data <- generate_dgp_complex(500)

approaches <- list(
  list(num = 1, name = "Full-sample", fun = estimate_att_fullsample),
  list(num = 2, name = "Cross-fit separate", fun = estimate_att_crossfit),
  list(num = 3, name = "Doubletree", fun = estimate_att_doubletree),
  list(num = 4, name = "Doubletree averaged", fun = estimate_att_doubletree_averaged),
  list(num = 5, name = "M-split", fun = estimate_att_msplit),
  list(num = 6, name = "M-split averaged", fun = estimate_att_msplit_averaged)
)

results <- list()
cat("================================================================================\n")
cat("RUNNING APPROACHES\n")
cat("================================================================================\n\n")

for (i in seq_along(approaches)) {
  approach <- approaches[[i]]

  cat(sprintf("Approach %d: %s\n", approach$num, approach$name))
  cat(strrep("-", 80), "\n")

  start_time <- Sys.time()

  result <- tryCatch({
    approach$fun(X = data$X, A = data$A, Y = data$Y)
  }, error = function(e) {
    list(theta = NA_real_, se = NA_real_, error = conditionMessage(e))
  })

  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")

  if (!is.null(result$error) && !is.na(result$error)) {
    cat(sprintf("  ✗ FAILED\n"))
    cat(sprintf("  Error: %s\n", result$error))
    cat(sprintf("  Time: %.1f sec\n\n", elapsed))

    results[[i]] <- list(
      approach = approach$num,
      name = approach$name,
      success = FALSE,
      theta = NA,
      se = NA,
      bias = NA,
      time = elapsed,
      error = result$error
    )
  } else if (is.na(result$theta) || is.na(result$se)) {
    cat(sprintf("  ✗ FAILED (returned NA)\n"))
    cat(sprintf("  Theta: %s, SE: %s\n", result$theta, result$se))
    cat(sprintf("  Time: %.1f sec\n\n", elapsed))

    results[[i]] <- list(
      approach = approach$num,
      name = approach$name,
      success = FALSE,
      theta = result$theta,
      se = result$se,
      bias = NA,
      time = elapsed,
      error = "Returned NA"
    )
  } else {
    bias <- result$theta - data$true_att

    cat(sprintf("  ✓ SUCCESS\n"))
    cat(sprintf("  Estimate: %.4f\n", result$theta))
    cat(sprintf("  True ATT: %.4f\n", data$true_att))
    cat(sprintf("  Bias:     %.4f\n", bias))
    cat(sprintf("  SE:       %.4f\n", result$se))
    cat(sprintf("  Time:     %.1f sec\n\n", elapsed))

    results[[i]] <- list(
      approach = approach$num,
      name = approach$name,
      success = TRUE,
      theta = result$theta,
      se = result$se,
      bias = bias,
      time = elapsed,
      error = NA
    )
  }
}

cat("================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n\n")

n_success <- sum(sapply(results, function(r) r$success))
n_total <- length(results)

cat(sprintf("Success rate: %d / %d (%.1f%%)\n\n", n_success, n_total, 100 * n_success / n_total))

if (n_success > 0) {
  cat("Successful approaches:\n")
  cat(sprintf("%-5s %-25s %10s %10s %10s %10s\n",
              "Num", "Name", "Estimate", "Bias", "SE", "Time (s)"))
  cat(strrep("-", 80), "\n")

  for (r in results) {
    if (r$success) {
      cat(sprintf("%-5d %-25s %10.4f %10.4f %10.4f %10.1f\n",
                  r$approach, r$name, r$theta, r$bias, r$se, r$time))
    }
  }
  cat("\n")

  # Summary stats
  successful_results <- results[sapply(results, function(r) r$success)]
  biases <- sapply(successful_results, function(r) r$bias)

  cat(sprintf("Bias range: [%.4f, %.4f]\n", min(biases), max(biases)))
  cat(sprintf("Mean bias:  %.4f\n", mean(biases)))
  cat(sprintf("Max |bias|: %.4f\n\n", max(abs(biases))))

  if (max(abs(biases)) < 0.01) {
    cat("✓✓✓ EXCELLENT: All biases < 0.01 (nearly unbiased)\n")
  } else if (max(abs(biases)) < 0.015) {
    cat("✓✓ VERY GOOD: All biases < 0.015\n")
  } else if (max(abs(biases)) < 0.020) {
    cat("✓ GOOD: All biases < 0.020 (improved from old ~-0.020)\n")
  } else {
    cat("⚠ WARNING: Some biases still >= 0.020\n")
  }
}

if (n_success < n_total) {
  cat("\nFailed approaches:\n")
  cat(strrep("-", 80), "\n")

  for (r in results) {
    if (!r$success) {
      cat(sprintf("Approach %d (%s):\n", r$approach, r$name))
      cat(sprintf("  Error: %s\n", r$error))
    }
  }
}

cat("\n")
if (n_success == n_total) {
  cat("✓ ALL APPROACHES WORKING - Ready for cluster deployment\n")
} else if (n_success >= 4) {
  cat("⚠ MOST APPROACHES WORKING - Review failures before deployment\n")
} else {
  cat("✗ MULTIPLE FAILURES - Fix issues before deployment\n")
}

cat("\n")
