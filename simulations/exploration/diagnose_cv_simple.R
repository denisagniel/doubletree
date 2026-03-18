# diagnose_cv_simple.R
# Deep dive: Why is CV loss flat even with 4 features?

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools required")
devtools::load_all()
if (!requireNamespace("optimaltrees", quietly = TRUE)) stop("optimaltrees required")

source("simulations/dgps_realistic.R")

n <- 400
d <- generate_dgp_simple(n, tau = 0.15, seed = 123)

message("=== Diagnosing CV with 4 features ===\n")
message(sprintf("n = %d, features = %d, patterns = %d\n", n, ncol(d$X), 2^ncol(d$X)))

# Run CV with default grid
message("Testing default CV grid:")
cv_result <- optimaltrees::cv_regularization(
  d$X, d$A,
  loss_function = "log_loss",
  K = 5,
  refit = FALSE,
  verbose = FALSE
)

message(sprintf("Lambda grid: %s", paste(sprintf("%.5f", cv_result$lambda_grid), collapse=", ")))
message(sprintf("CV losses:   %s", paste(sprintf("%.6f", cv_result$cv_loss), collapse=", ")))
message(sprintf("Range: %.6f (%.2f%% of mean)\n",
                max(cv_result$cv_loss) - min(cv_result$cv_loss),
                100 * (max(cv_result$cv_loss) - min(cv_result$cv_loss)) / mean(cv_result$cv_loss)))

# Try a MUCH wider grid
message("Testing wide lambda grid (0.001 to 0.5):")
wide_grid <- c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
cv_wide <- optimaltrees::cv_regularization(
  d$X, d$A,
  loss_function = "log_loss",
  lambda_grid = wide_grid,
  K = 5,
  refit = FALSE,
  verbose = FALSE
)

message(sprintf("Lambda grid: %s", paste(sprintf("%.3f", cv_wide$lambda_grid), collapse=", ")))
message(sprintf("CV losses:   %s", paste(sprintf("%.6f", cv_wide$cv_loss), collapse=", ")))
message(sprintf("Range: %.6f (%.2f%% of mean)\n",
                max(cv_wide$cv_loss) - min(cv_wide$cv_loss),
                100 * (max(cv_wide$cv_loss) - min(cv_wide$cv_loss)) / mean(cv_wide$cv_loss)))

# Fit trees with different lambdas to see tree sizes
message("Tree complexity at different lambda values:")
for (lam in c(0.001, 0.01, 0.05, 0.2)) {
  tree <- optimaltrees::fit_tree(d$X, d$A, loss_function = "log_loss",
                               regularization = lam, verbose = FALSE)
  # Count leaves (prediction nodes)
  tree_str <- as.character(tree)
  n_leaves <- length(gregexpr("prediction", tree_str)[[1]])
  if (n_leaves < 1) n_leaves <- 1
  message(sprintf("  lambda = %.3f: ~%d leaves", lam, n_leaves))
}

message("\n=== Hypothesis ===")
message("If CV loss is flat AND all trees have similar complexity,")
message("then the feature space is saturated - all lambda values")
message("produce functionally equivalent trees.")
message("\nSolution: Need MORE features or CONTINUOUS features")
message("to create sufficient feature space for regularization to matter.")
