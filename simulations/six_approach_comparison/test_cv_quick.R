#!/usr/bin/env Rscript

# Quick test of CV default for doubletree approach
# Tests 2 replications with simple and complex DGPs

cat("Loading packages from source...\n")
suppressPackageStartupMessages({
  devtools::load_all("../../../optimaltrees", quiet = TRUE)
  devtools::load_all("../..", quiet = TRUE)  # doubletree
})

cat("Loading simulation code...\n")
source("code/dgps.R")
source("code/estimators.R")
source("code/metrics.R")

cat("\n===========================================\n")
cat("Testing Doubletree with CV Default\n")
cat("===========================================\n\n")

# Test with 2 DGPs and 2 replications each
dgps <- list(
  list(id = 1, name = "simple", fn = generate_dgp_simple),
  list(id = 3, name = "complex", fn = generate_dgp_complex)
)

results <- list()

for (dgp_info in dgps) {
  cat(sprintf("\n--- DGP %d (%s) ---\n", dgp_info$id, dgp_info$name))

  for (rep in 1:2) {
    cat(sprintf("\n  Replication %d:\n", rep))

    set.seed(1000 + dgp_info$id * 100 + rep)
    data <- dgp_info$fn(n = 500)

    start_time <- Sys.time()
    result <- tryCatch({
      estimate_att_doubletree(X = data$X, A = data$A, Y = data$Y, K = 5)
    }, error = function(e) {
      list(theta = NA, se = NA, error = as.character(e))
    })
    elapsed <- as.numeric(Sys.time() - start_time, units = "secs")

    if (!is.null(result$error)) {
      cat(sprintf("    ERROR: %s\n", result$error))
      status <- "FAILED"
    } else {
      bias <- result$theta - data$true_att
      z_score <- abs(bias / result$se)
      covered <- z_score <= 1.96

      cat(sprintf("    theta_hat = %.4f (true = %.4f)\n", result$theta, data$true_att))
      cat(sprintf("    se = %.4f\n", result$se))
      cat(sprintf("    bias = %.4f\n", bias))
      cat(sprintf("    z-score = %.2f\n", z_score))
      cat(sprintf("    covered = %s\n", covered))
      cat(sprintf("    time = %.1f sec\n", elapsed))
      status <- "SUCCESS"
    }

    results[[length(results) + 1]] <- list(
      dgp = dgp_info$name,
      rep = rep,
      theta = result$theta,
      se = result$se,
      true_att = data$true_att,
      elapsed = elapsed,
      status = status
    )
  }
}

cat("\n\n===========================================\n")
cat("Summary\n")
cat("===========================================\n\n")

results_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(
    dgp = x$dgp,
    rep = x$rep,
    theta = x$theta,
    se = x$se,
    true_att = x$true_att,
    elapsed = x$elapsed,
    status = x$status
  )
}))

successes <- sum(results_df$status == "SUCCESS")
failures <- sum(results_df$status == "FAILED")

cat(sprintf("Total tests: %d\n", nrow(results_df)))
cat(sprintf("Successful: %d\n", successes))
cat(sprintf("Failed: %d\n", failures))

if (successes > 0) {
  cat(sprintf("\nMean elapsed time: %.1f sec\n", mean(results_df$elapsed[results_df$status == "SUCCESS"], na.rm = TRUE)))
}

cat("\nDetailed results:\n")
print(results_df)

if (failures == 0) {
  cat("\n✓ All tests passed! CV default is working correctly.\n")
} else {
  cat("\n✗ Some tests failed. Review errors above.\n")
}
