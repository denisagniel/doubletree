# test_continuous_quick.R
# Quick test of one continuous DGP

# Load packages quietly
suppressMessages({
  devtools::load_all("/Users/dagniel/OneDrive - RAND Corporation/rprojects/global-scholars/optimaltrees")
  devtools::load_all()
  source("simulations/dgps_continuous.R")
})

cat("Testing continuous DGP (sparse)...\n\n")

# Generate data
n <- 400
d <- generate_dgp_continuous_sparse(n, tau = 0.15, seed = 123)

cat(sprintf("Generated n=%d with %d continuous features\n", n, ncol(d$X)))

# Test 1: Single tree
cat("\nTest 1: Fitting tree with adaptive discretization...\n")

# Capture all output from GOSDT (including C++ stdout)
capture.output({
  tree <- optimaltrees::fit_tree(
    d$X, d$A,
    loss_function = "log_loss",
    regularization = 0.01,
    discretize_method = "quantiles",
    discretize_bins = "adaptive",
    verbose = FALSE
  )
}, file = tempfile())

# Check predictions
pred <- predict(tree, d$X, type = "prob")
pred_sd <- sd(pred[, 2])

cat(sprintf("  Predictions SD: %.4f\n", pred_sd))
if (pred_sd > 0.02) {
  cat("  ✓ Predictions vary (SD > 0.02)\n")
} else {
  cat("  ✗ Predictions too constant\n")
}

# Test 2: CV
cat("\nTest 2: CV regularization...\n")

# Capture all output from GOSDT
capture.output({
  cv <- optimaltrees::cv_regularization(
    d$X, d$A,
    loss_function = "log_loss",
    K = 5,
    refit = FALSE,
    discretize_method = "quantiles",
    discretize_bins = "adaptive",
    verbose = FALSE
  )
}, file = tempfile())

range_cv <- max(cv$cv_loss) - min(cv$cv_loss)
rel_range <- range_cv / mean(cv$cv_loss)

cat(sprintf("  CV range: %.2f%% of mean\n", rel_range * 100))
cat(sprintf("  Best λ: %.4f (CV loss: %.4f)\n",
            cv$regularization[which.min(cv$cv_loss)],
            min(cv$cv_loss)))

if (rel_range > 0.01) {
  cat("  ✓ CV informative (range > 1%)\n")
} else {
  cat("  ⚠ CV weakly informative\n")
}

# Summary
cat("\n")
if (pred_sd > 0.02 && rel_range > 0.005) {
  cat("✓✓ SUCCESS! Continuous DGPs work.\n")
  cat("Ready for Phase 3: Coverage testing\n")
} else {
  cat("✗ Issues detected. Review above.\n")
}
