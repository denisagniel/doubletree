# test_binary_5feat.R
# Test if 5 binary features (32 patterns) give better CV than 4 features

suppressMessages({
  devtools::load_all("/Users/dagniel/OneDrive - RAND Corporation/rprojects/global-scholars/optimaltrees")
  devtools::load_all()
  source("simulations/dgps_realistic.R")
})

cat("Testing 5 binary features (32 patterns)...\n\n")

# Generate moderate DGP (5 binary features)
n <- 400
d <- generate_dgp_moderate(n, tau = 0.15, seed = 123)

cat(sprintf("Generated n=%d with %d binary features\n", n, ncol(d$X)))
cat(sprintf("Covariate patterns: 2^%d = %d\n", ncol(d$X), 2^ncol(d$X)))

# Test 1: Tree predictions
cat("\nTest 1: Fitting tree...\n")
capture.output({
  tree <- optimaltrees::fit_tree(
    d$X, d$A,
    loss_function = "log_loss",
    regularization = 0.01,
    verbose = FALSE
  )
}, file = tempfile())

pred <- predict(tree, d$X, type = "prob")
pred_sd <- sd(pred[, 2])
cat(sprintf("  Predictions SD: %.4f\n", pred_sd))

# Test 2: CV
cat("\nTest 2: CV regularization...\n")
capture.output({
  cv <- optimaltrees::cv_regularization(
    d$X, d$A,
    loss_function = "log_loss",
    K = 5,
    refit = FALSE,
    verbose = FALSE
  )
}, file = tempfile())

range_cv <- max(cv$cv_loss) - min(cv$cv_loss)
rel_range <- range_cv / mean(cv$cv_loss)

cat(sprintf("  CV range: %.2f%% of mean\n", rel_range * 100))
cat(sprintf("  Best λ: %.4f\n", cv$regularization[which.min(cv$cv_loss)]))

# Summary
cat("\n")
if (pred_sd > 0.02 && rel_range > 0.01) {
  cat("✓ SUCCESS: 5 binary features work well\n")
} else if (pred_sd > 0.02) {
  cat("⚠ PARTIAL: Predictions vary but CV not very informative\n")
} else {
  cat("✗ FAILED\n")
}
