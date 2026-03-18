# test_8_features.R
# Quick test: Does 8 features make CV informative?

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools required")
devtools::load_all()
if (!requireNamespace("optimaltrees", quietly = TRUE)) stop("optimaltrees required")

source("simulations/dgps_final.R")

message("=== Testing 8-feature DGPs ===\n")

n <- 400
d1 <- generate_dgp1_final(n, tau = 0.15, seed = 123)
d2 <- generate_dgp2_final(n, tau = 0.15, seed = 123)

message("DGP1 (sparse): 8 features, 256 patterns\n")

# Test CV on propensity
cv1 <- optimaltrees::cv_regularization(d1$X, d1$A, loss_function = "log_loss",
                                    K = 5, refit = FALSE, verbose = FALSE)
range1 <- max(cv1$cv_loss) - min(cv1$cv_loss)
rel_range1 <- range1 / mean(cv1$cv_loss)

message("Propensity model:")
message(sprintf("  Lambda grid: %s", paste(sprintf("%.5f", cv1$lambda_grid), collapse=", ")))
message(sprintf("  CV losses:   %s", paste(sprintf("%.6f", cv1$cv_loss), collapse=", ")))
message(sprintf("  Range: %.6f (%.2f%% of mean)", range1, rel_range1 * 100))

if (rel_range1 > 0.01) {
  message("  ✓ CV IS INFORMATIVE!\n")
} else {
  message("  ✗ CV still flat\n")
}

message("DGP2 (moderate): 8 features, 256 patterns\n")

cv2 <- optimaltrees::cv_regularization(d2$X, d2$A, loss_function = "log_loss",
                                    K = 5, refit = FALSE, verbose = FALSE)
range2 <- max(cv2$cv_loss) - min(cv2$cv_loss)
rel_range2 <- range2 / mean(cv2$cv_loss)

message("Propensity model:")
message(sprintf("  Lambda grid: %s", paste(sprintf("%.5f", cv2$lambda_grid), collapse=", ")))
message(sprintf("  CV losses:   %s", paste(sprintf("%.6f", cv2$cv_loss), collapse=", ")))
message(sprintf("  Range: %.6f (%.2f%% of mean)", range2, rel_range2 * 100))

if (rel_range2 > 0.01) {
  message("  ✓ CV IS INFORMATIVE!\n")
} else {
  message("  ✗ CV still flat\n")
}

# Check tree complexity at different lambdas
message("\nTree complexity comparison (DGP1):")
for (lam in c(0.005, 0.015, 0.05)) {
  tree <- optimaltrees::fit_tree(d1$X, d1$A, loss_function = "log_loss",
                               regularization = lam, verbose = FALSE)
  tree_str <- as.character(tree)
  n_leaves <- length(gregexpr("prediction", tree_str)[[1]])
  if (n_leaves < 1) n_leaves <- 1
  message(sprintf("  lambda = %.3f: ~%d leaves", lam, n_leaves))
}

message("\n=== Conclusion ===")
if (rel_range1 > 0.01 || rel_range2 > 0.01) {
  message("✓ 8 features creates sufficient feature space!")
  message("CV regularization can now distinguish between lambda values.")
  message("Ready to run full simulations with these DGPs.")
} else {
  message("✗ Even 8 features insufficient - may need 10+ or different approach")
}
