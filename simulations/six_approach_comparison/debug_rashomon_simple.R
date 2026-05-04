#!/usr/bin/env Rscript
#
# Simplified Rashomon Debug
# Just show key structure info for each fold
#
# Created: 2026-05-04

suppressPackageStartupMessages({
  devtools::load_all("../../../optimaltrees", quiet = TRUE)
  devtools::load_all("../..", quiet = TRUE)
})

source("code/dgps.R")

cat("===================================\n")
cat("Simplified Rashomon Structure Debug\n")
cat("===================================\n\n")

set.seed(123)
n <- 1000
data <- generate_dgp_simple(n = n)

K <- 5
folds <- sample(rep(1:K, length.out = n))

cat("--- Propensity Trees ---\n\n")
for (k in 1:K) {
  train_idx <- folds != k

  tree <- optimaltrees::fit_tree(
    X = data$X[train_idx, , drop = FALSE],
    y = data$A[train_idx],
    loss_function = "log_loss",
    regularization = 0.1,
    verbose = FALSE
  )

  struct <- optimaltrees::extract_tree_structure(tree)

  # Get split features
  split_features <- sapply(struct@splits, function(s) s$feature_name)

  cat(sprintf("Fold %d: depth=%d, leaves=%d, splits=[%s]\n",
              k, struct@max_depth, struct@n_leaves,
              paste(split_features, collapse=", ")))
}

cat("\n--- Outcome Trees ---\n\n")
control_idx <- data$A == 0

for (k in 1:K) {
  train_idx <- (folds != k) & control_idx

  tree <- optimaltrees::fit_tree(
    X = data$X[train_idx, , drop = FALSE],
    y = data$Y[train_idx],
    loss_function = "log_loss",
    regularization = 0.1,
    verbose = FALSE
  )

  struct <- optimaltrees::extract_tree_structure(tree)

  # Get split features
  split_features <- sapply(struct@splits, function(s) s$feature_name)

  cat(sprintf("Fold %d: depth=%d, leaves=%d, splits=[%s]\n",
              k, struct@max_depth, struct@n_leaves,
              paste(split_features, collapse=", ")))
}

cat("\n")
