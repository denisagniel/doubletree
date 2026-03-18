# test_4_features.R
# Test if problem only occurs with 8 features

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools required")
devtools::load_all()

source("simulations/dgps_realistic.R")  # Has 4-feature DGPs

message("=== Testing with 4 Features ===\n")

n <- 400
d <- generate_dgp_simple(n, tau = 0.15, seed = 123)

message("DGP: 4 binary features = 16 covariate patterns\n")

# Test single tree fit
message("Testing tree fit with lambda = 0.01...")
tree <- optimaltrees::fit_tree(
  d$X, d$A,
  loss_function = "log_loss",
  regularization = 0.01,
  verbose = FALSE
)

pred <- predict(tree, d$X, type = "prob")
if (is.matrix(pred)) {
  e_pred <- pred[, 2]
  message(sprintf("Predictions: Min=%.3f, Max=%.3f, SD=%.6f",
                  min(e_pred), max(e_pred), sd(e_pred)))

  if (sd(e_pred) > 0.01) {
    message("✓ Tree IS learning with 4 features!")

    # Test CV
    message("\nTesting CV...")
    cv_result <- optimaltrees::cv_regularization(
      d$X, d$A,
      loss_function = "log_loss",
      K = 5,
      refit = FALSE,
      verbose = FALSE
    )

    range_loss <- max(cv_result$cv_loss) - min(cv_result$cv_loss)
    rel_range <- range_loss / mean(cv_result$cv_loss)

    message(sprintf("CV losses: %s", paste(sprintf("%.6f", cv_result$cv_loss), collapse=", ")))
    message(sprintf("Range: %.6f (%.2f%% of mean)", range_loss, rel_range * 100))

    if (rel_range > 0.01) {
      message("✓ CV is informative with 4 features")
    } else {
      message("✗ CV still flat with 4 features")
    }
  } else {
    message("✗ Tree NOT learning even with 4 features - deeper problem!")
  }
}
