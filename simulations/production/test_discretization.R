#' Test Manual Discretization + TreeFARMS

library(dplyr)
devtools::load_all("../../../optimaltrees")

source("../../R/discretize_adaptive.R")
source("dgps/dgps_beta_continuous.R")

cat("Testing manual discretization + optimaltrees...\n\n")

# Generate test data
set.seed(123)
d <- generate_dgp_beta_high(n = 800, tau = 0.10, seed = 123)

cat("Original continuous X:\n")
cat("  X1 unique values:", length(unique(d$X$X1)), "\n")
cat("  X2 unique values:", length(unique(d$X$X2)), "\n\n")

# Discretize
disc_result <- discretize_adaptive(d$X, n_bins = "adaptive", method = "quantiles")

cat("After discretization:\n")
cat("  Number of bins (b_n):", disc_result$n_bins, "\n")
cat("  X1_disc unique values:", length(unique(disc_result$X_discrete$X1)), "\n")
cat("  X2_disc unique values:", length(unique(disc_result$X_discrete$X2)), "\n")
cat("  X1_disc distribution:", table(disc_result$X_discrete$X1), "\n")
cat("  X2_disc distribution:", table(disc_result$X_discrete$X2), "\n\n")

# Fit tree on discretized features
cat("Fitting tree on discretized features...\n")
tree_model <- optimaltrees::fit_tree(
  disc_result$X_discrete,
  d$A,
  loss_function = "log_loss",
  regularization = log(800) / 800,
  verbose = FALSE
)

# Predict
pred_prob <- predict(tree_model, disc_result$X_discrete, type = "prob")
pred_class1 <- if (is.matrix(pred_prob)) pred_prob[, 2] else rep(0.5, nrow(d$X))

cat("\nPredictions:\n")
cat("  Mean:", mean(pred_class1), "\n")
cat("  SD:", sd(pred_class1), "\n")
cat("  Range: [", min(pred_class1), ",", max(pred_class1), "]\n")
cat("  Unique:", length(unique(pred_class1)), "\n\n")

if (sd(pred_class1) > 0.01) {
  cat("✓ Manual discretization WORKS!\n")
  cat("  Correlation with true e(X):", cor(pred_class1, d$true_e), "\n")
} else {
  cat("❌ Still producing constant predictions\n")
  cat("\nDebugging info:\n")
  cat("  Tree model class:", class(tree_model), "\n")
  cat("  Tree has model?:", !is.null(tree_model$model), "\n")
  if (!is.null(tree_model$model)) {
    cat("  Model has tree_json?:", !is.null(tree_model$model$tree_json), "\n")
  }
}
