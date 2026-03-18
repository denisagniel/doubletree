# test_traverse_debug.R
# Create a debug version of get_probabilities_from_tree with verbose output

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools required")
devtools::load_all()

source("simulations/dgps_realistic.R")

n <- 10  # Very small for debugging
d <- generate_dgp_simple(n, tau = 0.15, seed = 123)

tree <- optimaltrees::fit_tree(
  d$X, d$A,
  loss_function = "log_loss",
  regularization = 0.01,
  verbose = FALSE
)

tree_json <- tree$model$tree_json
X <- d$X

message("=== Debug version of get_probabilities_from_tree ===\n")

# Copy of the function with debug output
get_probabilities_debug <- function(tree_json, X) {
  n_samples <- nrow(X)
  probabilities <- matrix(0.5, nrow = n_samples, ncol = 2)
  max_depth <- 100L
  update_count <- 0

  leaf_probs_from_node <- function(node) {
    if (is.null(node) || !is.list(node)) return(c(0.5, 0.5))
    if (is.null(node$prediction)) return(c(0.5, 0.5))
    if (!is.null(node$probabilities) && length(node$probabilities) >= 2) {
      probs <- as.numeric(node$probabilities)
      if (length(probs) == 2 && all(is.finite(probs)) && all(probs >= 0)) {
        prob_sum <- sum(probs)
        if (prob_sum > 0) probs <- probs / prob_sum else probs <- c(0.5, 0.5)
        return(probs)
      }
    }
    pred <- as.numeric(node$prediction)
    if (pred == 0) c(1.0, 0.0) else c(0.0, 1.0)
  }

  traverse_batch <- function(node, row_indices, depth) {
    if (length(row_indices) == 0) return()
    if (depth > max_depth) {
      message("  WARNING: Max depth exceeded")
      return()
    }
    if (is.null(node) || !is.list(node)) {
      message("  WARNING: NULL or non-list node at depth ", depth)
      return()
    }

    # Check if leaf
    if (!is.null(node$prediction)) {
      probs <- leaf_probs_from_node(node)
      message("  LEAF at depth ", depth, ": rows ", min(row_indices), "-", max(row_indices),
              " (", length(row_indices), " obs), probs = [",
              round(probs[1], 3), ", ", round(probs[2], 3), "]")
      probabilities[row_indices, 1] <<- probs[1]
      probabilities[row_indices, 2] <<- probs[2]
      update_count <<- update_count + length(row_indices)
      return()
    }

    # Split node
    if (!is.null(node$feature)) {
      feature_idx <- as.integer(as.numeric(node$feature) + 1)
      message("  SPLIT at depth ", depth, ": feature ", feature_idx,
              " (", names(X)[feature_idx], "), rows ", min(row_indices), "-", max(row_indices))

      if (feature_idx < 1 || feature_idx > ncol(X)) {
        message("    WARNING: Invalid feature index ", feature_idx)
        return()
      }

      feature_vals <- X[row_indices, feature_idx, drop = TRUE]
      go_true <- (feature_vals == 1 | feature_vals == TRUE)
      true_idx <- row_indices[go_true]
      false_idx <- row_indices[!go_true]

      message("    → TRUE: ", length(true_idx), " rows, FALSE: ", length(false_idx), " rows")

      if (length(true_idx) > 0 && !is.null(node$true) && is.list(node$true)) {
        traverse_batch(node$true, true_idx, depth + 1)
      }
      if (length(false_idx) > 0 && !is.null(node$false) && is.list(node$false)) {
        traverse_batch(node$false, false_idx, depth + 1)
      }
    }
  }

  traverse_batch(tree_json, seq_len(n_samples), 0L)

  message("\n", update_count, " / ", n_samples, " rows updated")
  colnames(probabilities) <- c("P(class=0)", "P(class=1)")
  return(probabilities)
}

probs <- get_probabilities_debug(tree_json, X)

message("\nFirst 5 predictions:")
print(head(probs, 5))
message("SD: ", sd(probs[,2]))
