# test_low_regularization.R
# Test if very low regularization allows model generation

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools required")
devtools::load_all()

source("simulations/dgps_final.R")

message("=== Testing Very Low Regularization ===\n")

n <- 400
d <- generate_dgp1_final(n, tau = 0.15, seed = 123)

# Try progressively lower regularization
for (lam in c(0.01, 0.005, 0.001, 0.0001)) {
  message(sprintf("\nTesting lambda = %.4f...", lam))

  tree <- optimaltrees::fit_tree(
    d$X, d$A,
    loss_function = "log_loss",
    regularization = lam,
    verbose = FALSE
  )

  pred <- predict(tree, d$X, type = "prob")
  if (is.matrix(pred)) {
    e_pred <- pred[, 2]
    sd_pred <- sd(e_pred)
    message(sprintf("  Predictions SD: %.6f", sd_pred))

    if (sd_pred > 0.01) {
      message(sprintf("  ✓ SUCCESS at lambda = %.4f!", lam))

      # Now test CV with very low lambda grid
      message("\n=== Testing CV with low lambda grid ===")
      low_grid <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)

      cv_result <- optimaltrees::cv_regularization(
        d$X, d$A,
        loss_function = "log_loss",
        lambda_grid = low_grid,
        K = 5,
        refit = FALSE,
        verbose = FALSE
      )

      range_loss <- max(cv_result$cv_loss) - min(cv_result$cv_loss)
      rel_range <- range_loss / mean(cv_result$cv_loss)

      message(sprintf("Lambda grid: %s", paste(sprintf("%.4f", cv_result$lambda_grid), collapse=", ")))
      message(sprintf("CV losses:   %s", paste(sprintf("%.6f", cv_result$cv_loss), collapse=", ")))
      message(sprintf("Range: %.6f (%.2f%% of mean)", range_loss, rel_range * 100))

      if (rel_range > 0.01) {
        message("\n✓✓ CV IS INFORMATIVE with low regularization!")
      } else {
        message("\n✗ CV still flat")
      }

      break  # Success, stop testing
    } else {
      message(sprintf("  ✗ Still constant at lambda = %.4f", lam))
    }
  }
}
