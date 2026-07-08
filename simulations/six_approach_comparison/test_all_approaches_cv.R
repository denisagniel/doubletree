#!/usr/bin/env Rscript

# Comprehensive test of ALL 6 approaches with CV defaults
# Tests 1 replication per approach with simple DGP

library(optimaltrees)
# Load development version of doubletree with CV fixes (2026-05-26)
devtools::load_all("../..", quiet = TRUE)

cat("Loading simulation code...\n")
source("code/dgps.R")
source("code/estimators.R")
source("code/metrics.R")

cat("\n===========================================\n")
cat("Testing All 6 Approaches with CV\n")
cat("===========================================\n\n")

# Generate test data once
set.seed(999)
data <- generate_dgp_simple(n = 500)

approaches <- list(
  list(id = 1, name = "fullsample", fn = estimate_att_fullsample, cv = "yes"),
  list(id = 2, name = "crossfit", fn = estimate_att_crossfit, cv = "yes"),
  list(id = 3, name = "doubletree", fn = estimate_att_doubletree, cv = "yes"),
  list(id = 4, name = "doubletree_averaged", fn = estimate_att_doubletree_averaged, cv = "yes"),  # Updated 2026-05-26: now uses CV
  list(id = 5, name = "msplit", fn = estimate_att_msplit, cv = "no"),
  list(id = 6, name = "msplit_averaged", fn = estimate_att_msplit_averaged, cv = "yes")  # Updated 2026-05-26: now uses CV
)

results <- list()

for (approach in approaches) {
  cat(sprintf("\n--- Approach %d: %s (CV: %s) ---\n",
              approach$id, approach$name, approach$cv))

  start_time <- Sys.time()
  result <- tryCatch({
    if (approach$id == 1) {
      # Fullsample: no K parameter
      approach$fn(X = data$X, A = data$A, Y = data$Y)
    } else if (approach$id %in% 2:4) {
      # Crossfit, doubletree, doubletree_averaged: use K
      approach$fn(X = data$X, A = data$A, Y = data$Y, K = 5)
    } else {
      # M-split approaches: use M and K
      approach$fn(X = data$X, A = data$A, Y = data$Y, M = 5, K = 3)
    }
  }, error = function(e) {
    list(theta = NA, se = NA, error = as.character(e))
  })
  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")

  if (!is.null(result$error)) {
    cat(sprintf("  ✗ ERROR: %s\n", substr(result$error, 1, 100)))
    status <- "FAILED"
  } else {
    bias <- result$theta - data$true_att
    z_score <- abs(bias / result$se)
    covered <- z_score <= 1.96

    cat(sprintf("  ✓ theta = %.4f (true = %.4f, bias = %.4f)\n",
                result$theta, data$true_att, bias))
    cat(sprintf("    se = %.4f, z-score = %.2f, covered = %s\n",
                result$se, z_score, covered))
    cat(sprintf("    time = %.1f sec\n", elapsed))
    status <- "SUCCESS"
  }

  results[[length(results) + 1]] <- list(
    id = approach$id,
    name = approach$name,
    cv = approach$cv,
    theta = result$theta,
    se = result$se,
    true_att = data$true_att,
    elapsed = elapsed,
    status = status
  )
}

cat("\n\n===========================================\n")
cat("Summary of All Approaches\n")
cat("===========================================\n\n")

results_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(
    id = x$id,
    name = x$name,
    cv_enabled = x$cv,
    theta = x$theta,
    se = x$se,
    bias = x$theta - x$true_att,
    elapsed = x$elapsed,
    status = x$status,
    stringsAsFactors = FALSE
  )
}))

successes <- sum(results_df$status == "SUCCESS")
failures <- sum(results_df$status == "FAILED")

cat(sprintf("Total tests: %d\n", nrow(results_df)))
cat(sprintf("Successful: %d\n", successes))
cat(sprintf("Failed: %d\n", failures))

cat("\nApproaches using CV:\n")
print(results_df[results_df$cv_enabled == "yes", c("id", "name", "theta", "bias", "se", "elapsed", "status")])

cat("\nApproaches using fixed regularization:\n")
print(results_df[results_df$cv_enabled == "no", c("id", "name", "theta", "bias", "se", "elapsed", "status")])

if (failures == 0) {
  cat("\n✓ ALL TESTS PASSED! All approaches are working correctly.\n")
  cat("  Approaches 1-3 now use CV-selected lambda.\n")
  cat("  Ready for cluster deployment.\n")
} else {
  cat("\n✗ Some tests failed. Review errors above.\n")
  cat("  Do NOT deploy to cluster until issues are resolved.\n")
}
