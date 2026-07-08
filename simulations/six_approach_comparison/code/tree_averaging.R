# Tree Averaging Utilities for Six-Approach Comparison
# Created: 2026-05-20
#
# Functions for extracting and averaging leaf values across K trees
# with the same structure but different leaf values (from cross-fitting).

#' Extract Leaf Values from a Tree
#'
#' Recursively traverse a tree (nested list structure) and extract
#' all leaf values with their paths from the root.
#'
#' @param tree_node Nested list representing a tree (from refit_structure_on_data)
#' @param path Integer vector tracking current path from root
#' @return Named numeric vector: names = leaf paths (e.g., "0-1", "1-0-1"),
#'   values = probabilities P(Y=1|leaf) for binary outcomes
#' @examples
#' # tree_node structure:
#' # Internal node: list(feature = 0, relation = "==", reference = 1, true = ..., false = ...)
#' # Leaf node: list(prediction = 0 or 1, probabilities = c(p0, p1))
extract_leaf_values <- function(tree_node, path = integer(0)) {
  # Validate input
  if (!is.list(tree_node) || length(tree_node) == 0) {
    stop("tree_node must be a non-empty list", call. = FALSE)
  }

  # Check if this is a leaf (has 'prediction' field)
  if (!is.null(tree_node$prediction)) {
    # This is a leaf
    if (is.null(tree_node$probabilities) || length(tree_node$probabilities) != 2) {
      stop("Leaf node missing 'probabilities' or wrong length. Expected c(p0, p1).", call. = FALSE)
    }

    # Extract P(Y=1) from probabilities = c(P(Y=0), P(Y=1))
    p1 <- tree_node$probabilities[2]

    # Create path string
    if (length(path) == 0) {
      path_str <- "root"
    } else {
      path_str <- paste(path, collapse = "-")
    }

    # Return named vector
    result <- p1
    names(result) <- path_str
    return(result)
  }

  # This is an internal node (has 'feature')
  if (is.null(tree_node$feature)) {
    stop("Internal node must have 'feature' field", call. = FALSE)
  }

  if (is.null(tree_node$true) || is.null(tree_node$false)) {
    stop("Internal node must have 'true' and 'false' children", call. = FALSE)
  }

  # Recursively extract from left (false) and right (true) children
  # Path convention: 0 = left/false, 1 = right/true
  left_values <- extract_leaf_values(tree_node$false, c(path, 0L))
  right_values <- extract_leaf_values(tree_node$true, c(path, 1L))

  # Combine and return
  c(left_values, right_values)
}

#' Average Leaf Values Across K Trees
#'
#' Given K trees with the same structure but different leaf values,
#' compute the average leaf value for each leaf position.
#'
#' @param tree_list List of K trees (each from refit_structure_on_data)
#' @return Named numeric vector of averaged leaf values
#' @examples
#' # tree_list[[1]] has leaves: c("0" = 0.3, "1-0" = 0.6, "1-1" = 0.8)
#' # tree_list[[2]] has leaves: c("0" = 0.4, "1-0" = 0.5, "1-1" = 0.9)
#' # Result: c("0" = 0.35, "1-0" = 0.55, "1-1" = 0.85)
average_leaf_values <- function(tree_list) {
  # Validate input
  if (!is.list(tree_list) || length(tree_list) == 0) {
    stop("tree_list must be a non-empty list of trees", call. = FALSE)
  }

  K <- length(tree_list)

  # Extract leaf values from each tree
  leaf_values_list <- lapply(tree_list, extract_leaf_values)

  # Check all trees have same structure (same leaf paths)
  leaf_paths_list <- lapply(leaf_values_list, names)
  first_paths <- sort(leaf_paths_list[[1]])

  for (k in 2:K) {
    current_paths <- sort(leaf_paths_list[[k]])
    if (!identical(first_paths, current_paths)) {
      stop("Trees have different structures. Cannot average.\n",
           "Tree 1 leaves: ", paste(first_paths, collapse = ", "), "\n",
           "Tree ", k, " leaves: ", paste(current_paths, collapse = ", "),
           call. = FALSE)
    }
  }

  # Convert to matrix: rows = trees, cols = leaves
  leaf_paths <- first_paths
  n_leaves <- length(leaf_paths)
  leaf_matrix <- matrix(NA_real_, nrow = K, ncol = n_leaves)
  colnames(leaf_matrix) <- leaf_paths

  for (k in 1:K) {
    # Match leaf values to standardized order
    leaf_matrix[k, ] <- leaf_values_list[[k]][leaf_paths]
  }

  # Average across K trees (down columns)
  averaged <- colMeans(leaf_matrix, na.rm = FALSE)

  # Validate no NAs in result
  if (any(is.na(averaged))) {
    stop("Averaged leaf values contain NA. This should not happen.", call. = FALSE)
  }

  return(averaged)
}

#' Create Single Tree with Averaged Leaf Values
#'
#' Rebuild a tree structure with averaged leaf values.
#' Uses the structure from the first tree and replaces all leaf values
#' with the averaged values.
#'
#' @param tree_template First tree (used for structure)
#' @param averaged_values Named numeric vector from average_leaf_values()
#' @return Tree (nested list) with averaged leaf values
rebuild_tree_with_averaged_values <- function(tree_template, averaged_values, path = integer(0)) {
  # Validate inputs
  if (!is.list(tree_template) || length(tree_template) == 0) {
    stop("tree_template must be a non-empty list", call. = FALSE)
  }

  if (!is.numeric(averaged_values) || is.null(names(averaged_values))) {
    stop("averaged_values must be a named numeric vector", call. = FALSE)
  }

  # Check if this is a leaf
  if (!is.null(tree_template$prediction)) {
    # Create path string
    if (length(path) == 0) {
      path_str <- "root"
    } else {
      path_str <- paste(path, collapse = "-")
    }

    # Get averaged value for this leaf
    if (!(path_str %in% names(averaged_values))) {
      stop("Leaf path '", path_str, "' not found in averaged_values.\n",
           "Available: ", paste(names(averaged_values), collapse = ", "),
           call. = FALSE)
    }

    p1_avg <- averaged_values[[path_str]]

    # Validate probability
    if (p1_avg < 0 || p1_avg > 1) {
      stop("Averaged probability out of [0,1]: ", p1_avg, " at path ", path_str,
           call. = FALSE)
    }

    # Return leaf node with averaged value
    return(list(
      prediction = as.integer(p1_avg >= 0.5),
      probabilities = c(1 - p1_avg, p1_avg)
    ))
  }

  # This is an internal node - rebuild with same structure
  if (is.null(tree_template$feature)) {
    stop("Internal node must have 'feature' field", call. = FALSE)
  }

  # Recursively rebuild children
  false_child <- rebuild_tree_with_averaged_values(
    tree_template$false, averaged_values, c(path, 0L)
  )
  true_child <- rebuild_tree_with_averaged_values(
    tree_template$true, averaged_values, c(path, 1L)
  )

  # Return internal node with same split, averaged children
  list(
    feature = tree_template$feature,
    relation = tree_template$relation,
    reference = tree_template$reference,
    false = false_child,
    true = true_child
  )
}

#' Average K Trees and Create Single Tree
#'
#' Main function: Takes K trees with same structure, averages leaf values,
#' returns single tree.
#'
#' @param tree_list List of K trees (from fold_refits)
#' @return Single tree with averaged leaf values
#' @export
average_trees <- function(tree_list) {
  # Validate
  if (!is.list(tree_list) || length(tree_list) == 0) {
    stop("tree_list must be a non-empty list", call. = FALSE)
  }

  # Compute averaged leaf values
  averaged_values <- average_leaf_values(tree_list)

  # Rebuild tree using first tree as template
  single_tree <- rebuild_tree_with_averaged_values(tree_list[[1]], averaged_values)

  return(single_tree)
}

#' Predict from Tree (Helper)
#'
#' Predict P(Y=1|X) for new observations using a tree.
#' This is a simplified prediction function for trees from refit_structure_on_data.
#'
#' @param tree Tree (nested list)
#' @param X Data.frame of covariates
#' @return Numeric vector of predictions P(Y=1|X)
predict_from_tree <- function(tree, X) {
  n <- nrow(X)
  predictions <- numeric(n)

  for (i in 1:n) {
    predictions[i] <- predict_single_obs(tree, X[i, , drop = FALSE])
  }

  return(predictions)
}

#' Predict for Single Observation
#'
#' Traverse tree to get prediction for one observation.
#'
#' @param tree_node Current node (nested list)
#' @param x_row Single-row data.frame
#' @return Numeric P(Y=1)
predict_single_obs <- function(tree_node, x_row) {
  # Leaf?
  if (!is.null(tree_node$prediction)) {
    return(tree_node$probabilities[2])  # P(Y=1)
  }

  # Internal node - evaluate split
  feature_idx <- tree_node$feature + 1  # 0-indexed to 1-indexed
  feature_val <- as.numeric(x_row[[feature_idx]])

  # For binary features with relation "=="
  if (tree_node$relation == "==") {
    # reference is typically 1 for binary
    if (abs(feature_val - tree_node$reference) < 1e-10) {
      # Goes right (true)
      return(predict_single_obs(tree_node$true, x_row))
    } else {
      # Goes left (false)
      return(predict_single_obs(tree_node$false, x_row))
    }
  } else if (tree_node$relation == "<=") {
    if (feature_val <= tree_node$reference) {
      return(predict_single_obs(tree_node$true, x_row))
    } else {
      return(predict_single_obs(tree_node$false, x_row))
    }
  } else if (tree_node$relation == ">") {
    if (feature_val > tree_node$reference) {
      return(predict_single_obs(tree_node$true, x_row))
    } else {
      return(predict_single_obs(tree_node$false, x_row))
    }
  } else {
    stop("Unknown relation: ", tree_node$relation, call. = FALSE)
  }
}
