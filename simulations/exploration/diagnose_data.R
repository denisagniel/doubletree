# diagnose_data.R
# Check if DGPs are generating valid data

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools required")
devtools::load_all()

source("simulations/dgps_final.R")

message("=== Data Quality Check ===\n")

n <- 400
d <- generate_dgp1_final(n, tau = 0.15, seed = 123)

message("Sample size: ", n)
message("Features: ", ncol(d$X))
message("\nFeature summaries:")
print(summary(d$X))

message("\nTreatment:")
message(sprintf("  A=0: %d (%.1f%%)", sum(d$A == 0), 100 * mean(d$A == 0)))
message(sprintf("  A=1: %d (%.1f%%)", sum(d$A == 1), 100 * mean(d$A == 1)))

message("\nOutcome:")
message(sprintf("  Y=0: %d (%.1f%%)", sum(d$Y == 0), 100 * mean(d$Y == 0)))
message(sprintf("  Y=1: %d (%.1f%%)", sum(d$Y == 1), 100 * mean(d$Y == 1)))

message("\nTrue propensity distribution:")
message(sprintf("  Min: %.3f, Max: %.3f, Mean: %.3f",
                min(d$true_e), max(d$true_e), mean(d$true_e)))

# Try fitting a single tree
message("\n=== Testing Tree Fitting ===\n")
message("Fitting propensity tree with lambda = 0.01...")

tree <- optimaltrees::fit_tree(
  d$X, d$A,
  loss_function = "log_loss",
  regularization = 0.01,
  verbose = TRUE
)

message("\nTree fitted. Checking predictions...")
pred <- predict(tree, d$X, type = "prob")
if (is.matrix(pred)) {
  e_pred <- pred[, 2]
  message(sprintf("Predictions: Min=%.3f, Max=%.3f, Mean=%.3f, SD=%.3f",
                  min(e_pred), max(e_pred), mean(e_pred), sd(e_pred)))

  if (sd(e_pred) < 0.01) {
    message("✗ PROBLEM: Predictions are nearly constant!")
    message("Tree is not learning from the data.")
  } else {
    message("✓ Predictions vary - tree is learning")
  }
} else {
  message("✗ PROBLEM: Unexpected prediction format")
}
