# test_8_features_fixed.R
# Test if 8 features now produces informative CV after the fix

devtools::load_all("/Users/dagniel/OneDrive - RAND Corporation/rprojects/global-scholars/optimaltrees")
devtools::load_all()

source("simulations/dgps_final.R")

message("=== Testing 8-Feature DGPs After Fix ===\n")

n <- 400
d1 <- generate_dgp1_final(n, tau = 0.15, seed = 123)

message("DGP1 (sparse): 8 features, 256 patterns\n")

# Test single tree first
message("Testing single tree fit...")
tree <- optimaltrees::fit_tree(d1$X, d1$A, loss_function = "log_loss",
                            regularization = 0.01, verbose = FALSE)
pred <- predict(tree, d1$X, type = "prob")
message(sprintf("  Predictions SD: %.6f", sd(pred[,2])))

if (sd(pred[,2]) > 0.01) {
  message("  ✓ Tree learning successfully\n")

  # Now test CV
  message("Testing CV...")
  cv1 <- optimaltrees::cv_regularization(d1$X, d1$A, loss_function = "log_loss",
                                      K = 5, refit = FALSE, verbose = FALSE)
  range1 <- max(cv1$cv_loss) - min(cv1$cv_loss)
  rel_range1 <- range1 / mean(cv1$cv_loss)

  message(sprintf("  Lambda grid: %s", paste(sprintf("%.5f", cv1$lambda_grid), collapse=", ")))
  message(sprintf("  CV losses:   %s", paste(sprintf("%.6f", cv1$cv_loss), collapse=", ")))
  message(sprintf("  Range: %.6f (%.2f%% of mean)", range1, rel_range1 * 100))

  if (rel_range1 > 0.01) {
    message("\n✓✓✓ CV IS INFORMATIVE with 8 features!")
    message("The bug fix has solved the problem. Ready to run full simulations.")
  } else {
    message("\n? CV still relatively flat (", round(rel_range1 * 100, 2), "%)")
    message("May need continuous features or different DGP design")
  }
} else {
  message("  ✗ Tree still not learning - something else is wrong")
}
