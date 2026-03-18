# test_continuous_cv.R
# Phase 2: Test that continuous DGPs produce varying predictions and informative CV

# Load packages
devtools::load_all("/Users/dagniel/OneDrive - RAND Corporation/rprojects/global-scholars/optimaltrees")
devtools::load_all()
source("simulations/dgps_continuous.R")

cat("===== Phase 2: Testing Continuous DGPs with CV =====\n\n")

# Test parameters
n <- 400
tau <- 0.15
seed <- 123

# Test each DGP
dgps <- list(
  sparse = generate_dgp_continuous_sparse,
  moderate = generate_dgp_continuous_moderate,
  complex = generate_dgp_continuous_complex,
  rct_like = generate_dgp_continuous_rct_like
)

results <- list()

for (dgp_name in names(dgps)) {
  cat(sprintf("Testing DGP: %s\n", dgp_name))
  cat(rep("-", 50), "\n", sep = "")

  # Generate data
  d <- dgps[[dgp_name]](n, tau = tau, seed = seed)

  cat(sprintf("  Generated n=%d with %d continuous features\n", n, ncol(d$X)))

  # Note: fit_tree() will automatically discretize continuous features
  # using discretize_bins = "adaptive" (default creates ~log(n)/3 bins per feature)

  # Test 1: Single tree predictions
  cat("\n  Test 1: Tree predictions (should vary)...\n")
  tree <- tryCatch({
    optimaltrees::fit_tree(
      d$X, d$A,
      loss_function = "log_loss",
      regularization = 0.01,
      discretize_method = "quantiles",
      discretize_bins = "adaptive",
      verbose = FALSE
    )
  }, error = function(e) {
    cat(sprintf("    ✗ FAILED: %s\n", e$message))
    return(NULL)
  })

  if (!is.null(tree)) {
    pred <- predict(tree, d$X, type = "prob")
    pred_sd <- sd(pred[, 2])
    cat(sprintf("    Predictions SD: %.4f\n", pred_sd))

    if (pred_sd > 0.02) {
      cat("    ✓ SUCCESS: Predictions vary (SD > 0.02)\n")
      results[[dgp_name]]$pred_sd <- pred_sd
      results[[dgp_name]]$pred_test <- "pass"
    } else {
      cat("    ✗ FAILED: Predictions too constant (SD ≤ 0.02)\n")
      results[[dgp_name]]$pred_sd <- pred_sd
      results[[dgp_name]]$pred_test <- "fail"
    }
  } else {
    results[[dgp_name]]$pred_test <- "error"
  }

  # Test 2: CV informativeness
  cat("\n  Test 2: CV regularization (should be informative)...\n")
  cv <- tryCatch({
    optimaltrees::cv_regularization(
      d$X, d$A,
      loss_function = "log_loss",
      K = 5,
      refit = FALSE,
      discretize_method = "quantiles",
      discretize_bins = "adaptive",
      verbose = FALSE
    )
  }, error = function(e) {
    cat(sprintf("    ✗ FAILED: %s\n", e$message))
    return(NULL)
  })

  if (!is.null(cv)) {
    range_cv <- max(cv$cv_loss) - min(cv$cv_loss)
    rel_range <- range_cv / mean(cv$cv_loss)

    cat(sprintf("    CV loss range: %.6f\n", range_cv))
    cat(sprintf("    Relative range: %.2f%% of mean\n", rel_range * 100))
    cat(sprintf("    Min λ: %.4f (CV loss: %.4f)\n",
                cv$regularization[which.min(cv$cv_loss)],
                min(cv$cv_loss)))
    cat(sprintf("    Max λ: %.4f (CV loss: %.4f)\n",
                cv$regularization[which.max(cv$cv_loss)],
                max(cv$cv_loss)))

    if (rel_range > 0.01) {
      cat("    ✓ SUCCESS: CV informative (range > 1% of mean)\n")
      results[[dgp_name]]$cv_range <- rel_range
      results[[dgp_name]]$cv_test <- "pass"
    } else {
      cat("    ✗ WARNING: CV weakly informative (range ≤ 1% of mean)\n")
      results[[dgp_name]]$cv_range <- rel_range
      results[[dgp_name]]$cv_test <- "warn"
    }
  } else {
    results[[dgp_name]]$cv_test <- "error"
  }

  cat("\n")
}

# Summary
cat("\n")
cat("===== SUMMARY =====\n")
cat(sprintf("%-15s | %10s | %10s | %s\n", "DGP", "Pred Test", "CV Test", "Status"))
cat(rep("=", 70), "\n", sep = "")

all_pass <- TRUE
for (dgp_name in names(results)) {
  pred_status <- results[[dgp_name]]$pred_test
  cv_status <- results[[dgp_name]]$cv_test

  overall <- if (pred_status == "pass" && cv_status %in% c("pass", "warn")) {
    "✓ READY"
  } else {
    all_pass <- FALSE
    "✗ FAIL"
  }

  cat(sprintf("%-15s | %10s | %10s | %s\n",
              dgp_name, pred_status, cv_status, overall))
}

cat("\n")
if (all_pass) {
  cat("✓✓ ALL TESTS PASSED! Continuous DGPs are ready for Phase 3 (coverage testing).\n")
  cat("\nNext step: Run simulations/test_continuous_coverage.R\n")
} else {
  cat("✗✗ SOME TESTS FAILED. Review output above.\n")
}
