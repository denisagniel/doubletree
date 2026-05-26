# Tree Averaging Utilities
#
# Functions for averaging leaf values across K trees with the same structure
# but different leaf values (from cross-fitting). Used by estimate_att_*_averaged()
# functions to create a single interpretable tree while maintaining cross-fit validity.

#' Validate Tree Structure for Averaging
#'
#' Check that a tree object is in the correct nested-list format for averaging.
#' Recursively validates internal nodes and leaf nodes.
#'
#' @param tree Nested list representing a tree
#' @param tree_name Character name for error messages (default "tree")
#'
#' @return NULL if valid, stops with error if invalid
#'
#' @keywords internal
validate_tree_for_averaging <- function(tree, tree_name = "tree") {
  if (!is.list(tree) || length(tree) == 0) {
    stop("Tree '", tree_name, "' must be a non-empty list, got: ", class(tree)[1], call. = FALSE)
  }

  # Check if S7 object (wrong type)
  if (inherits(tree, "S7_object")) {
    stop("Tree '", tree_name, "' is an S7 object (", class(tree)[1], "). ",
         "Expected nested list with fields like 'feature', 'true', 'false'. ",
         "Hint: Use model@trees[[1]] instead of extract_tree_structure(model).",
         call. = FALSE)
  }

  # Check for required fields
  is_leaf <- !is.null(tree$prediction)
  is_internal <- !is.null(tree$feature)

  if (!is_leaf && !is_internal) {
    stop("Tree '", tree_name, "' has neither 'prediction' (leaf) nor 'feature' (internal node). ",
         "Available fields: ", paste(names(tree), collapse = ", "), call. = FALSE)
  }

  # Validate leaf node
  if (is_leaf) {
    if (is.null(tree$probabilities) || length(tree$probabilities) != 2) {
      stop("Leaf node in '", tree_name, "' missing 'probabilities' or wrong length. ",
           "Expected c(p0, p1), got length: ", length(tree$probabilities), call. = FALSE)
    }
  }

  # Validate internal node
  if (is_internal) {
    if (is.null(tree$true) || is.null(tree$false)) {
      stop("Internal node in '", tree_name, "' missing 'true' or 'false' children.", call. = FALSE)
    }
    # Recursively validate children
    validate_tree_for_averaging(tree$true, paste0(tree_name, "$true"))
    validate_tree_for_averaging(tree$false, paste0(tree_name, "$false"))
  }

  NULL  # Valid
}

#' Extract Leaf Values from a Tree
#'
#' Recursively traverse a tree (nested list structure) and extract
#' all leaf values with their paths from the root.
#'
#' @param tree_node Nested list representing a tree (from refit_structure_on_data).
#'   Internal nodes have: feature, relation, reference, true, false.
#'   Leaf nodes have: prediction, probabilities.
#' @param path Integer vector tracking current path from root (internal use)
#'
#' @return Named numeric vector where names are leaf paths (e.g., "0-1", "1-0-1")
#'   and values are probabilities P(Y=1|leaf) for binary outcomes
#'
#' @keywords internal
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
    # Note: probabilities might be a list from JSON, convert to numeric vector
    probs <- as.numeric(tree_node$probabilities)
    p1 <- probs[2]

    # Create path string
    if (length(path) == 0) {
      path_str <- "root"
    } else {
      path_str <- paste(path, collapse = "-")
    }

    # Return named numeric vector (not list)
    result <- setNames(p1, path_str)
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
#' @param tree_list List of K trees (each from refit_structure_on_data or fold_refits)
#'
#' @return Named numeric vector of averaged leaf values
#'
#' @keywords internal
average_leaf_values <- function(tree_list) {
  # Validate input
  if (!is.list(tree_list) || length(tree_list) == 0) {
    stop("tree_list must be a non-empty list of trees", call. = FALSE)
  }

  K <- length(tree_list)

  # Extract leaf values from each tree with better error handling
  leaf_values_list <- vector("list", K)
  for (i in seq_len(K)) {
    leaf_values_list[[i]] <- tryCatch({
      extract_leaf_values(tree_list[[i]])
    }, error = function(e) {
      stop("Failed to extract leaf values from tree ", i, " of ", K,
           ": ", e$message, call. = FALSE)
    })
  }

  # Check all trees have same structure (same leaf paths)
  leaf_paths_list <- lapply(leaf_values_list, names)

  # Validate all extractions succeeded
  if (length(leaf_paths_list) != K) {
    stop("Expected ", K, " trees but only extracted ", length(leaf_paths_list),
         " leaf path sets", call. = FALSE)
  }

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

#' Rebuild Tree with Averaged Leaf Values
#'
#' Rebuild a tree structure with averaged leaf values.
#' Uses the structure from the first tree and replaces all leaf values
#' with the averaged values.
#'
#' @param tree_template First tree (used for structure)
#' @param averaged_values Named numeric vector from average_leaf_values()
#' @param path Integer vector tracking current path from root (internal use)
#'
#' @return Tree (nested list) with averaged leaf values
#'
#' @keywords internal
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

#' Average Trees with Same Structure
#'
#' Takes K trees with the same structure but different leaf values (from cross-fitting),
#' averages the leaf values, and returns a single tree. This maintains cross-fit validity
#' while producing one interpretable tree.
#'
#' @param tree_list List of K trees (from fold_refits). Each tree should be a nested list
#'   from \code{optimaltrees::refit_structure_on_data()} or similar.
#'
#' @return Single tree (nested list) with averaged leaf values
#'
#' @details
#' This function is used by \code{\link{estimate_att_doubletree_averaged}} and
#' \code{\link{estimate_att_msplit_averaged}} to create a single interpretable tree
#' from K cross-fitted trees.
#'
#' The averaging maintains cross-fit validity because each tree's leaf values were
#' fitted on training data (excluding the test fold), so averaging them preserves
#' the cross-fit property.
#'
#' All trees must have identical structure (same splits). If trees have different
#' structures, the function will error.
#'
#' @seealso \code{\link{estimate_att_doubletree_averaged}}, \code{\link{estimate_att_msplit_averaged}}
#'
#' @export
average_trees <- function(tree_list) {
  # Validate
  if (!is.list(tree_list) || length(tree_list) == 0) {
    stop("tree_list must be a non-empty list", call. = FALSE)
  }

  # Validate each tree structure
  for (i in seq_along(tree_list)) {
    tryCatch({
      validate_tree_for_averaging(tree_list[[i]], tree_name = paste0("tree_", i))
    }, error = function(e) {
      stop("Tree ", i, " in tree_list is invalid: ", e$message, call. = FALSE)
    })
  }

  # Compute averaged leaf values
  averaged_values <- average_leaf_values(tree_list)

  # Rebuild tree using first tree as template
  single_tree <- rebuild_tree_with_averaged_values(tree_list[[1]], averaged_values)

  return(single_tree)
}

#' Predict from Tree
#'
#' Predict P(Y=1|X) for new observations using a tree (nested list structure).
#' This is a helper for trees from \code{average_trees()}.
#'
#' @param tree Tree (nested list) from \code{average_trees()}
#' @param X Data.frame of covariates
#'
#' @return Numeric vector of predictions P(Y=1|X)
#'
#' @keywords internal
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
#'
#' @return Numeric P(Y=1)
#'
#' @keywords internal
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
