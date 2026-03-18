# test_fix.R
# Test that the environment fix works

# Reload optimaltrees first (where the fix was made)
message("Reloading optimaltrees...")
devtools::load_all("/Users/dagniel/OneDrive - RAND Corporation/rprojects/global-scholars/optimaltrees")

message("Reloading dmltree...")
devtools::load_all()

source("simulations/dgps_realistic.R")

message("\n=== Testing Fix with 4 Features ===\n")

n <- 400
d <- generate_dgp_simple(n, tau = 0.15, seed = 123)

message("Fitting tree...")
tree <- optimaltrees::fit_tree(
  d$X, d$A,
  loss_function = "log_loss",
  regularization = 0.01,
  verbose = FALSE
)

pred <- predict(tree, d$X, type = "prob")
if (is.matrix(pred)) {
  e_pred <- pred[, 2]
  message(sprintf("Predictions: Min=%.3f, Max=%.3f, Mean=%.3f, SD=%.6f",
                  min(e_pred), max(e_pred), mean(e_pred), sd(e_pred)))

  if (sd(e_pred) > 0.01) {
    message("\n✓✓ SUCCESS! Tree is now learning!")

    # Test CV
    message("\nTesting CV with fixed predictions...")
    cv_result <- optimaltrees::cv_regularization(
      d$X, d$A,
      loss_function = "log_loss",
      K = 5,
      refit = FALSE,
      verbose = FALSE
    )

    range_loss <- max(cv_result$cv_loss) - min(cv_result$cv_loss)
    rel_range <- range_loss / mean(cv_result$cv_loss)

    message(sprintf("Lambda grid: %s", paste(sprintf("%.5f", cv_result$lambda_grid), collapse=", ")))
    message(sprintf("CV losses:   %s", paste(sprintf("%.6f", cv_result$cv_loss), collapse=", ")))
    message(sprintf("Range: %.6f (%.2f%% of mean)", range_loss, rel_range * 100))

    if (rel_range > 0.01) {
      message("\n✓✓✓ CV IS NOW INFORMATIVE!")
      message("The fix worked! CV can now distinguish between lambda values.")
    } else {
      message("\n? CV still flat (may need more features)")
    }
  } else {
    message("\n✗ STILL BROKEN - predictions constant")
  }
}
