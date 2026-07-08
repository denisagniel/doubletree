#!/usr/bin/env Rscript
# Inspect the actual tree structures from each fold to see why they differ

library(optimaltrees)
library(doubletree)

source("code/dgps.R")

cat("==============================================\n")
cat("Inspecting Fold-Specific Tree Structures\n")
cat("==============================================\n\n")

set.seed(12345)
n <- 500
data <- generate_dgp_complex(n)

# Create folds
K <- 5
fold_indices <- sample(rep(1:K, length.out = n))

cat(sprintf("Data: n=%d, DGP=complex\n", n))
cat(sprintf("Folds: K=%d\n\n", K))

# Helper function to count leaves
count_leaves <- function(node) {
  if (is.null(node$feature)) {
    return(1)  # Leaf
  } else {
    left <- if (!is.null(node$true)) count_leaves(node$true) else 0
    right <- if (!is.null(node$false)) count_leaves(node$false) else 0
    return(left + right)
  }
}

# Fit propensity model for each fold
cat("==============================================\n")
cat("Propensity Trees (predicting A from X)\n")
cat("==============================================\n\n")

for (k in 1:K) {
  train_idx <- which(fold_indices != k)
  X_train <- data$X[train_idx, ]
  A_train <- data$A[train_idx]

  cat(sprintf("--- Fold %d (n_train=%d) ---\n", k, length(train_idx)))

  # Use CV to select lambda (same as the simulation)
  cv_result <- optimaltrees::cv_regularization_adaptive(
    X_train, A_train,
    loss_function = "log_loss",
    K = 5,
    max_iterations = 10,
    refit = FALSE,
    verbose = FALSE
  )

  lambda_selected <- cv_result$best_lambda
  cat(sprintf("CV selected lambda: %.5f\n", lambda_selected))

  # Fit tree with selected lambda
  model <- optimaltrees::fit_tree(
    X_train, A_train,
    loss_function = "log_loss",
    regularization = lambda_selected
  )

  # Extract tree structure
  tree <- model@trees[[1]]

  # Print tree details
  n_leaves <- count_leaves(tree)
  cat(sprintf("Tree has %d leaves\n", n_leaves))

  # Get splits used
  get_splits <- function(node, splits = c()) {
    if (!is.null(node$feature)) {
      splits <- c(splits, node$feature)
      if (!is.null(node$true)) splits <- get_splits(node$true, splits)
      if (!is.null(node$false)) splits <- get_splits(node$false, splits)
    }
    return(splits)
  }

  splits_used <- unique(get_splits(tree))
  feat_names <- colnames(data$X)[splits_used + 1]
  cat(sprintf("Features used in splits: %s\n", paste(feat_names, collapse = ", ")))

  # Only print full tree if it's small
  if (n_leaves <= 10) {
    # Print splits
    print_tree_structure <- function(node, depth = 0, prefix = "") {
    indent <- paste(rep("  ", depth), collapse = "")

    if (!is.null(node$feature)) {
      # Internal node
      feat_name <- colnames(data$X)[node$feature + 1]  # 0-indexed
      cat(sprintf("%s%sSplit on %s (feature %d)\n",
                  indent, prefix, feat_name, node$feature))
      if (!is.null(node$true)) {
        print_tree_structure(node$true, depth + 1, "├─ TRUE:  ")
      }
      if (!is.null(node$false)) {
        print_tree_structure(node$false, depth + 1, "└─ FALSE: ")
      }
    } else {
      # Leaf node
      if (!is.null(node$prediction)) {
        n_obs <- if (!is.null(node$n)) node$n else "?"
        cat(sprintf("%s%sLeaf: P(A=1) = %.3f (n=%s)\n",
                    indent, prefix, node$prediction[2], as.character(n_obs)))
      }
    }
  }

    cat("\nTree structure:\n")
    print_tree_structure(tree)
  } else {
    cat(sprintf("\n(Tree too large to display - %d leaves)\n", n_leaves))
  }
  cat("\n")
}

cat("\n==============================================\n")
cat("Summary\n")
cat("==============================================\n\n")

cat("Question: Why do folds have different tree structures?\n\n")

cat("Possible reasons:\n")
cat("1. Small training sets (n_train ≈ 400) → high variance\n")
cat("2. Binary features → limited smoothing, sharp decisions\n")
cat("3. Complex DGP → multiple \"good\" ways to partition\n")
cat("4. CV selects different lambda per fold → different complexity\n\n")
