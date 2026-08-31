#!/usr/bin/env Rscript
# Test: Does optimaltrees actually find multiple trees in Rashomon set?

library(optimaltrees)
source("code/dgps.R")

set.seed(12345)
n <- 500
data <- generate_dgp_complex(n)

# Take first fold's training data
fold_indices <- sample(rep(1:5, length.out = n))
train_idx <- which(fold_indices != 1)
X_train <- data$X[train_idx, ]
A_train <- data$A[train_idx]

cat("==============================================\n")
cat("Testing Rashomon Set Computation\n")
cat("==============================================\n\n")

cat(sprintf("Training data: n=%d\n", length(train_idx)))
cat("Fitting with epsilon=0.11 (should find multiple trees if landscape is flat)\n\n")

# Fit with Rashomon set computation
model <- optimaltrees::fit_rashomon(
  X_train, A_train,
  loss_function = "log_loss",
  regularization = 0.005,  # Similar to what CV selected
  rashomon_bound_multiplier = 0.11,  # c=1 for n=500
  verbose = TRUE
)

cat("\n==============================================\n")
cat("Results\n")
cat("==============================================\n")

n_trees <- model@n_trees
cat(sprintf("Number of trees in Rashomon set: %d\n\n", n_trees))

if (n_trees == 1) {
  cat("⚠️  WARNING: Only 1 tree found!\n")
  cat("This suggests either:\n")
  cat("  1. treefarms isn't enumerating alternative trees\n")
  cat("  2. The optimal tree is uniquely optimal (no alternatives within epsilon)\n\n")

  cat("Let's try with a MUCH larger epsilon:\n\n")

  model2 <- optimaltrees::fit_rashomon(
    X_train, A_train,
    loss_function = "log_loss",
    regularization = 0.005,
    rashomon_bound_multiplier = 1.0,  # 10× larger!
    verbose = TRUE
  )

  n_trees2 <- model2@n_trees
  cat(sprintf("\nWith epsilon=1.0: %d trees\n", n_trees2))

  if (n_trees2 == 1) {
    cat("\n❌ PROBLEM: Even at epsilon=1.0, only 1 tree!\n")
    cat("This strongly suggests treefarms is NOT enumerating alternatives.\n")
  } else {
    cat(sprintf("\n✓ With larger epsilon, found %d trees\n", n_trees2))
  }
}

cat("\n==============================================\n")
