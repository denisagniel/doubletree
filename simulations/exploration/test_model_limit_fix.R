# test_model_limit_fix.R
# Test if higher model_limit allows 8-feature trees to fit

devtools::load_all("/Users/dagniel/OneDrive - RAND Corporation/rprojects/global-scholars/optimaltrees")
devtools::load_all()

source("simulations/dgps_final.R")

message("=== Testing 8 Features with Higher model_limit ===\n")

n <- 400
d <- generate_dgp1_final(n, tau = 0.15, seed = 123)

# Try progressively higher model_limit values
for (limit in c(50000, 100000, 500000)) {
  message(sprintf("Testing model_limit = %d...", limit))

  tree <- optimaltrees::fit_tree(
    d$X, d$A,
    loss_function = "log_loss",
    regularization = 0.01,
    model_limit = limit,
    verbose = FALSE
  )

  pred <- predict(tree, d$X, type = "prob")
  sd_pred <- sd(pred[, 2])

  message(sprintf("  Predictions SD: %.6f", sd_pred))

  if (sd_pred > 0.01) {
    message(sprintf("  ✓ SUCCESS with model_limit = %d!\n", limit))

    # Test CV with this limit
    message("Testing CV...")
    cv <- optimaltrees::cv_regularization(
      d$X, d$A,
      loss_function = "log_loss",
      K = 5,
      refit = FALSE,
      verbose = FALSE,
      model_limit = limit
    )

    range_loss <- max(cv$cv_loss) - min(cv$cv_loss)
    rel_range <- range_loss / mean(cv$cv_loss)

    message(sprintf("  CV losses: %s", paste(sprintf("%.6f", cv$cv_loss), collapse=", ")))
    message(sprintf("  Range: %.6f (%.2f%% of mean)", range_loss, rel_range * 100))

    if (rel_range > 0.01) {
      message(sprintf("\n✓✓✓ CV IS INFORMATIVE with model_limit = %d!", limit))
      message("Solution: Pass model_limit to fit_tree() and cv_regularization()")
    } else {
      message("\n? CV still flat even with higher model_limit")
    }

    break  # Found working limit, stop testing
  } else {
    message(sprintf("  ✗ Still fails with model_limit = %d\n", limit))
  }
}
