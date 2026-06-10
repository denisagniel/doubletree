# Averaged-Tree ATT Estimators
#
# Alternative inference methods that average leaf values across cross-fitted trees
# to produce a single interpretable tree while maintaining cross-fit validity.

# =============================================================================
# Helper Functions for Tree Averaging
# =============================================================================

#' Collect Rashomon Trees at Given Tolerance with CV-Selected Lambda
#'
#' Fits K Rashomon sets (one per fold) at specified tolerance and collects
#' all trees from all folds. Uses CV to select lambda for each fold.
#'
#' @param X Data.frame of covariates
#' @param outcome Numeric vector (A for propensity, Y[A==0] for outcome)
#' @param K Number of folds
#' @param fold_indices List of length K with test indices per fold
#' @param tolerance Numeric. Rashomon bound multiplier (epsilon_n)
#' @param regularization Numeric. Not used (CV selects lambda). Kept for backward compatibility.
#' @param outcome_type Character. "binary" or "continuous"
#' @param verbose Logical
#' @param ... Additional arguments for optimaltrees
#'
#' @return List with:
#'   \item{all_trees}{List of all tree structures (nested lists) from all K folds}
#'   \item{all_structures}{List of all TreeStructure objects corresponding to all_trees}
#'   \item{intersection_trees}{List of trees in intersection (could be empty)}
#'   \item{n_total}{Total trees collected across all K folds}
#'   \item{n_intersecting}{Number of trees in intersection (0 if empty)}
#'   \item{tolerance_used}{The tolerance value used}
#'
#' @keywords internal
collect_rashomon_trees_at_tolerance <- function(X, outcome, K, fold_indices,
                                                tolerance, regularization,
                                                outcome_type, verbose = FALSE, ...) {
  n <- nrow(X)
  all_trees <- list()               # Tree structures (nested lists) for averaging
  all_structures <- list()          # TreeStructure objects for modal selection
  fold_tree_structures <- list()    # For intersection computation

  for (k in 1:K) {
    test_idx <- fold_indices[[k]]
    train_idx <- setdiff(1:n, test_idx)

    X_train <- X[train_idx, , drop = FALSE]
    outcome_train <- outcome[train_idx]

    # Select lambda via adaptive CV on training fold (2026-05-27)
    loss_fn <- if (outcome_type == "binary") "log_loss" else "squared_error"
    cv_result <- tryCatch({
      optimaltrees::cv_regularization_adaptive(
        X = X_train,
        y = outcome_train,
        loss_function = loss_fn,
        K = 5,
        max_iterations = 10,
        refit = FALSE,
        verbose = FALSE
      )
    }, error = function(e) {
      if (verbose) message("Fold ", k, " CV failed: ", e$message)
      return(list(best_lambda = NA_real_))
    })

    # No fallback - CV must succeed
    if (is.na(cv_result$best_lambda)) {
      stop(
        "CV failed for ", nuisance_type, " in fold ", k, " during Rashomon collection.\n",
        "Possible fixes:\n",
        "  1. Check data quality (enough variation in outcome for this fold?)\n",
        "  2. Try different lambda_grid in cv_regularization()\n",
        "  3. Increase K in cv_regularization() for more stable CV\n",
        "  4. Check for numerical issues (NaN, Inf in data)\n",
        "  5. Try different random seed (fold might have unusual split)",
        call. = FALSE
      )
    }

    lambda_k <- cv_result$best_lambda

    # Fit Rashomon set for this fold with selected lambda
    rashomon_model <- tryCatch({
      optimaltrees::fit_rashomon(
        X = X_train,
        y = outcome_train,
        loss_function = loss_fn,
        regularization = lambda_k,
        bound_multiplier = tolerance,
        bound_adder = 0,
        verbose = verbose,
        ...
      )
    }, error = function(e) {
      if (verbose) message("Fold ", k, " Rashomon fit failed: ", e$message)
      return(NULL)
    })

    if (is.null(rashomon_model)) next

    # Get number of trees in Rashomon set
    n_trees <- rashomon_model@n_trees
    if (n_trees == 0) next

    # Extract tree structures for intersection computation and modal selection
    fold_structures <- vector("list", n_trees)
    for (i in 1:n_trees) {
      fold_structures[[i]] <- optimaltrees::extract_tree_structure(rashomon_model, tree_index = i)
    }
    fold_tree_structures[[k]] <- fold_structures
    all_structures <- c(all_structures, fold_structures)  # Collect for modal selection

    # Extract actual tree structures (nested lists) for averaging
    fold_trees <- optimaltrees::get_rashomon_trees(rashomon_model)
    all_trees <- c(all_trees, fold_trees)
  }

  # Compute intersection across K folds using optimaltrees utilities
  intersection_trees <- list()
  n_intersecting <- 0

  if (length(fold_tree_structures) == K && all(sapply(fold_tree_structures, function(x) length(x) > 0))) {
    # Find structures that appear in all K folds
    # Use hash-based comparison for efficiency
    hash_counts <- list()

    for (k in 1:K) {
      fold_hashes <- sapply(fold_tree_structures[[k]], optimaltrees::structure_hash)
      unique_hashes <- unique(fold_hashes)
      for (h in unique_hashes) {
        if (is.null(hash_counts[[h]])) {
          hash_counts[[h]] <- 1
        } else {
          hash_counts[[h]] <- hash_counts[[h]] + 1
        }
      }
    }

    # Find hashes that appear in all K folds
    intersecting_hashes <- names(hash_counts)[sapply(hash_counts, function(x) x == K)]

    if (length(intersecting_hashes) > 0) {
      # Extract one tree per intersecting structure (from first fold)
      for (h in intersecting_hashes) {
        for (k in 1:K) {
          fold_hashes <- sapply(fold_tree_structures[[k]], optimaltrees::structure_hash)
          match_idx <- which(fold_hashes == h)[1]
          if (!is.na(match_idx)) {
            # Map back to original tree using all_structures (which correspond to all_trees)
            all_hashes <- sapply(all_structures, optimaltrees::structure_hash)
            tree_idx <- which(all_hashes == h)[1]
            if (!is.na(tree_idx)) {
              intersection_trees[[length(intersection_trees) + 1]] <- all_trees[[tree_idx]]
            }
            break
          }
        }
      }
      n_intersecting <- length(intersection_trees)
    }
  }

  list(
    all_trees = all_trees,
    all_structures = all_structures,
    intersection_trees = intersection_trees,
    n_total = length(all_trees),
    n_intersecting = n_intersecting,
    tolerance_used = tolerance
  )
}

#' Filter Trees by Structure
#'
#' Given lists of trees and their structures, return only trees that match
#' the target structure.
#'
#' @param tree_list List of trees (nested lists for averaging)
#' @param structure_list List of TreeStructure objects corresponding to tree_list
#' @param target_structure TreeStructure object to match
#'
#' @return List with:
#'   \item{matched_trees}{List of trees matching target structure}
#'   \item{n_matched}{Number of matched trees}
#'
#' @keywords internal
filter_trees_by_structure <- function(tree_list, structure_list, target_structure) {
  if (length(tree_list) == 0 || length(structure_list) == 0) {
    return(list(matched_trees = list(), n_matched = 0))
  }

  if (length(tree_list) != length(structure_list)) {
    stop("tree_list and structure_list must have same length", call. = FALSE)
  }

  target_hash <- optimaltrees::structure_hash(target_structure)

  matched_trees <- list()
  for (i in seq_along(tree_list)) {
    tree_hash <- optimaltrees::structure_hash(structure_list[[i]])

    if (tree_hash == target_hash) {
      matched_trees[[length(matched_trees) + 1]] <- tree_list[[i]]
    }
  }

  list(
    matched_trees = matched_trees,
    n_matched = length(matched_trees)
  )
}

#' Find Trees Through Five-Tier Fallback (Approach 4)
#'
#' Progressively relax strategy to find trees for one nuisance function:
#' Tier 1: Intersection at 0.05
#' Tier 2: Expanded intersection (0.10, 0.15, 0.20)
#' Tier 3: Modal at 0.05
#' Tier 4: Modal expanded (0.10, 0.15, 0.20)
#' Tier 5: Fold-specific fallback
#'
#' @param X Data.frame of covariates
#' @param outcome Numeric vector (A for propensity, Y[A==0] for m0)
#' @param K Number of folds
#' @param fold_indices List of K fold indices
#' @param nuisance_name Character. "propensity" or "outcome" (for messages)
#' @param regularization Numeric
#' @param outcome_type Character. "binary" or "continuous"
#' @param verbose Logical
#' @param ... Additional arguments
#'
#' @return List with:
#'   \item{trees}{List of trees (or NULL if all tiers fail)}
#'   \item{structure}{TreeStructure object (or NULL if all tiers fail)}
#'   \item{tier}{Character. Which tier succeeded}
#'   \item{tolerance}{Numeric or NA}
#'   \item{warning}{Character or NULL (theory warning if tolerance > 0.05)}
#'   \item{n_matched}{Integer (for modal tiers only)}
#'   \item{n_total}{Integer (for modal tiers only)}
#'   \item{modal_frequency}{Numeric (for modal tiers only)}
#'
#' @keywords internal
find_trees_through_tiers <- function(X, outcome, K, fold_indices, nuisance_name,
                                     regularization, outcome_type, verbose, ...) {
  n <- nrow(X)

  # Tier 1: Rashomon intersection at ε = 0.05
  if (verbose) message(sprintf("%s: Tier 1 - Intersection at ε = 0.05", nuisance_name))

  tier1_result <- collect_rashomon_trees_at_tolerance(
    X, outcome, K, fold_indices, tolerance = 0.05, regularization, outcome_type, verbose, ...
  )

  if (tier1_result$n_intersecting > 0) {
    if (verbose) message(sprintf("%s: Tier 1 succeeded (n_intersecting = %d)",
                                  nuisance_name, tier1_result$n_intersecting))
    # Extract structure from first tree in intersection
    first_tree_structure <- tier1_result$all_structures[[1]]

    # For averaging, we need one tree per fold (K trees total)
    # Filter all collected trees to only those matching the intersection structure
    filtered <- filter_trees_by_structure(tier1_result$all_trees,
                                          tier1_result$all_structures,
                                          first_tree_structure)

    if (filtered$n_matched < K) {
      warning("Tier 1 intersection succeeded but only ", filtered$n_matched,
              " of ", K, " folds have trees matching structure. Continuing to Tier 2.",
              call. = FALSE)
    } else {
      return(list(
        trees = filtered$matched_trees[1:K],  # Use first K trees (one per fold)
        structure = first_tree_structure,
        tier = "tier1_0.05",
        tolerance = 0.05,
        warning = NULL,
        n_matched = filtered$n_matched,
        n_total = tier1_result$n_total,
        modal_frequency = filtered$n_matched / tier1_result$n_total
      ))
    }
  }

  # Tier 2: Expanded Rashomon intersection (0.10, 0.15, 0.20)
  expanded_tolerances <- c(0.10, 0.15, 0.20)

  for (eps in expanded_tolerances) {
    if (verbose) message(sprintf("%s: Tier 2 - Intersection at ε = %.2f", nuisance_name, eps))

    tier2_result <- collect_rashomon_trees_at_tolerance(
      X, outcome, K, fold_indices, tolerance = eps, regularization, outcome_type, verbose, ...
    )

    if (tier2_result$n_intersecting > 0) {
      warning_msg <- sprintf(
        "Rashomon tolerance %.2f exceeds theoretical bound (0.05) for n=%d. Inference guarantees may not hold.",
        eps, n
      )
      if (verbose) message(sprintf("%s: Tier 2 succeeded at ε = %.2f (n_intersecting = %d) - WARNING ISSUED",
                                    nuisance_name, eps, tier2_result$n_intersecting))
      # Extract structure from first tree in intersection
      first_tree_structure <- tier2_result$all_structures[[1]]

      # For averaging, we need one tree per fold (K trees total)
      filtered <- filter_trees_by_structure(tier2_result$all_trees,
                                            tier2_result$all_structures,
                                            first_tree_structure)

      if (filtered$n_matched < K) {
        if (verbose) message("Tier 2 intersection succeeded but only ", filtered$n_matched,
                            " of ", K, " folds have trees. Continuing...")
        next  # Continue to next epsilon
      }

      return(list(
        trees = filtered$matched_trees[1:K],  # Use first K trees
        structure = first_tree_structure,
        tier = sprintf("tier2_%.2f", eps),
        tolerance = eps,
        warning = warning_msg,
        n_matched = filtered$n_matched,
        n_total = tier2_result$n_total,
        modal_frequency = filtered$n_matched / tier2_result$n_total
      ))
    }
  }

  # Tier 3: Modal at DEFAULT ε = 0.05
  if (verbose) message(sprintf("%s: Tier 3 - Modal structure at ε = 0.05", nuisance_name))

  tier3_result <- collect_rashomon_trees_at_tolerance(
    X, outcome, K, fold_indices, tolerance = 0.05, regularization, outcome_type, verbose, ...
  )

  if (tier3_result$n_total > 0) {
    # Use pre-extracted TreeStructure objects
    all_structures <- tier3_result$all_structures

    # Find modal structure
    modal_result <- select_structure_modal(all_structures)

    # Filter to trees matching modal
    filtered <- filter_trees_by_structure(tier3_result$all_trees, all_structures, modal_result$structure)

    if (filtered$n_matched >= 3) {
      if (verbose) message(sprintf("%s: Tier 3 succeeded (n_matched = %d / %d = %.1f%%)",
                                    nuisance_name, filtered$n_matched, tier3_result$n_total,
                                    modal_result$frequency * 100))
      return(list(
        trees = filtered$matched_trees,
        structure = modal_result$structure,
        tier = "tier3_0.05",
        tolerance = 0.05,
        warning = NULL,
        n_matched = filtered$n_matched,
        n_total = tier3_result$n_total,
        modal_frequency = modal_result$frequency
      ))
    }
  }

  # Tier 4: Modal with expanded Rashomon (0.10, 0.15, 0.20)
  for (eps in expanded_tolerances) {
    if (verbose) message(sprintf("%s: Tier 4 - Modal structure at ε = %.2f", nuisance_name, eps))

    tier4_result <- collect_rashomon_trees_at_tolerance(
      X, outcome, K, fold_indices, tolerance = eps, regularization, outcome_type, verbose, ...
    )

    if (tier4_result$n_total > 0) {
      # Use pre-extracted TreeStructure objects
      all_structures <- tier4_result$all_structures

      modal_result <- select_structure_modal(all_structures)
      filtered <- filter_trees_by_structure(tier4_result$all_trees, all_structures, modal_result$structure)

      if (filtered$n_matched >= 3) {
        warning_msg <- sprintf(
          "Rashomon tolerance %.2f exceeds theoretical bound (0.05) for n=%d. Inference guarantees may not hold.",
          eps, n
        )
        if (verbose) message(sprintf("%s: Tier 4 succeeded at ε = %.2f (n_matched = %d / %d = %.1f%%) - WARNING ISSUED",
                                      nuisance_name, eps, filtered$n_matched, tier4_result$n_total,
                                      modal_result$frequency * 100))
        return(list(
          trees = filtered$matched_trees,
          structure = modal_result$structure,
          tier = sprintf("tier4_%.2f", eps),
          tolerance = eps,
          warning = warning_msg,
          n_matched = filtered$n_matched,
          n_total = tier4_result$n_total,
          modal_frequency = modal_result$frequency
        ))
      }
    }
  }

  # Tier 5: Fold-specific trees with CV-selected lambda (no fallback - CV must succeed)
  if (verbose) message(sprintf("%s: Tier 5 - Fold-specific trees (no Rashomon)", nuisance_name))

  fold_trees <- list()
  for (k in 1:K) {
    test_idx <- fold_indices[[k]]
    train_idx <- setdiff(1:n, test_idx)

    X_train <- X[train_idx, , drop = FALSE]
    outcome_train <- outcome[train_idx]

    # Select lambda via adaptive CV (2026-05-27)
    loss_fn <- if (outcome_type == "binary") "log_loss" else "squared_error"
    cv_result <- tryCatch({
      optimaltrees::cv_regularization_adaptive(
        X = X_train,
        y = outcome_train,
        loss_function = loss_fn,
        K = 5,
        max_iterations = 10,
        refit = TRUE,
        verbose = FALSE
      )
    }, error = function(e) {
      if (verbose) message("Fold ", k, " CV failed: ", e$message)
      return(list(best_lambda = NA_real_, model = NULL))
    })

    # No fallback - CV must succeed
    if (is.na(cv_result$best_lambda) || is.null(cv_result$model)) {
      stop(
        "CV failed for ", nuisance_type, " in fold ", k, " (Tier 5 fallback).\n",
        "Possible fixes:\n",
        "  1. Check data quality (enough variation in outcome for this fold?)\n",
        "  2. Try different lambda_grid in cv_regularization()\n",
        "  3. Increase K in cv_regularization() for more stable CV\n",
        "  4. Check for numerical issues (NaN, Inf in data)\n",
        "  5. Try different random seed (fold might have unusual split)",
        call. = FALSE
      )
    }

    tree_k <- tryCatch({
      cv_result$model
    }, error = function(e) {
      if (verbose) message("Fold ", k, " tree fit failed: ", e$message)
      return(NULL)
    })

    if (!is.null(tree_k)) {
      fold_trees[[length(fold_trees) + 1]] <- tree_k
    }
  }

  if (length(fold_trees) > 0) {
    # Extract structures from fold-specific trees (these are models)
    fold_structures <- lapply(fold_trees, function(tree_model) {
      optimaltrees::extract_tree_structure(tree_model)
    })

    modal_result <- select_structure_modal(fold_structures)

    # For Tier 5, fold_trees are models, not tree structures
    # We need to convert models to tree structures for averaging
    fold_tree_structures_for_averaging <- lapply(fold_trees, function(model) {
      optimaltrees::get_rashomon_trees(model)[[1]]  # Get first tree from single-tree model
    })

    filtered <- filter_trees_by_structure(fold_tree_structures_for_averaging, fold_structures, modal_result$structure)

    if (filtered$n_matched >= 3) {
      if (verbose) message(sprintf("%s: Tier 5 succeeded (n_matched = %d / %d = %.1f%%)",
                                    nuisance_name, filtered$n_matched, length(fold_trees),
                                    modal_result$frequency * 100))
      return(list(
        trees = filtered$matched_trees,
        structure = modal_result$structure,
        tier = "tier5_fold_specific",
        tolerance = NA_real_,
        warning = NULL,
        n_matched = filtered$n_matched,
        n_total = length(fold_trees),
        modal_frequency = modal_result$frequency
      ))
    }
  }

  # All tiers failed
  if (verbose) message(sprintf("%s: All tiers failed (<3 trees matched modal)", nuisance_name))
  return(list(
    trees = NULL,
    structure = NULL,
    tier = "failed",
    tolerance = NA_real_,
    warning = NULL,
    n_matched = NA_integer_,
    n_total = NA_integer_,
    modal_frequency = NA_real_
  ))
}

# =============================================================================
# Approach 4: Doubletree Averaged (K-fold with Rashomon)
# =============================================================================

#' Extract K Trees from Cross-Fitted Rashomon Object
#'
#' Extract one tree per fold from a cross_fitted_rashomon object.
#' Each tree has the same structure (from Rashomon intersection) but different
#' leaf values (fit on that fold's training data).
#'
#' @param cf_rashomon_obj A cross_fitted_rashomon S7 object from fit_nuisances_rashomon()
#'
#' @return List of K trees (nested lists), one per fold
#'
#' @keywords internal
extract_k_trees_from_rashomon <- function(cf_rashomon_obj) {
  # Validate input
  if (cf_rashomon_obj@n_intersecting == 0) {
    stop("Cannot extract K trees: Rashomon intersection is empty", call. = FALSE)
  }

  K <- cf_rashomon_obj@K
  fold_refits <- cf_rashomon_obj@fold_refits

  # Extract first intersecting structure from each fold
  trees <- vector("list", K)
  for (k in 1:K) {
    # fold_refits[[k]] is a list of refit_result objects (one per intersecting structure)
    # We want the FIRST intersecting structure only
    # refit_result is a nested list tree from refit_structure_on_data()
    if (length(fold_refits[[k]]) == 0) {
      stop("Fold ", k, " has no refit results, but intersection claimed to exist", call. = FALSE)
    }

    trees[[k]] <- fold_refits[[k]][[1]]  # First intersecting structure
  }

  return(trees)
}

#' Estimate ATT with Doubletree and Averaged Leaves (Approach 4)
#'
#' Estimates the Average Treatment Effect on the Treated (ATT) using Rashomon
#' intersection to find K trees with a common structure, then averages their
#' leaf values to create a single interpretable tree. No cross-fitting in final
#' predictions (all observations use the same averaged tree).
#'
#' @inheritParams estimate_att
#' @param K Number of cross-fitting folds. Default: 5
#' @param outcome_type Character: "binary" or "continuous"
#' @param cv_regularization Logical. Use CV to select lambda? Default: TRUE
#' @param cv_K Integer. Number of folds for lambda CV. Default: 5
#' @param regularization Numeric. Fixed lambda (only used if cv_regularization=FALSE). Default: 0.1
#' @param stratified Logical. Stratify folds by treatment? Default: TRUE
#' @param seed Integer. Random seed for reproducibility. Default: NULL
#' @param verbose Logical. Print progress? Default: FALSE
#' @param rashomon_bound_multiplier Numeric. Rashomon bound multiplier (epsilon_n). Default: 0.05
#' @param rashomon_bound_adder Numeric. Additive Rashomon bound. Default: 0
#' @param max_leaves Integer. Maximum number of leaves (optional sieve). Default: NULL
#' @param auto_tune_intersecting Logical. Auto-tune epsilon_n to find intersection? Default: FALSE
#' @param discretize_method Character. "quantiles" or other. Default: "quantiles"
#' @param discretize_bins Integer or "adaptive". Default: "adaptive"
#' @param ... Additional arguments passed to optimaltrees functions
#'
#' @return List with elements:
#'   \item{theta}{ATT point estimate}
#'   \item{sigma}{Standard error}
#'   \item{ci_95}{95\% confidence interval}
#'   \item{score_values}{EIF score values}
#'   \item{e_hat}{Propensity score predictions (from averaged tree)}
#'   \item{m0_hat}{Control outcome predictions (from averaged tree)}
#'   \item{averaged_trees}{List with e and m0 averaged trees (nested lists)}
#'   \item{structures}{List with e and m0 TreeStructure objects}
#'   \item{n}{Sample size}
#'   \item{K}{Number of folds}
#'   \item{n_treated}{Number of treated observations}
#'   \item{outcome_type}{Outcome type ("binary" or "continuous")}
#'   \item{converged}{Logical (TRUE if Rashomon intersection succeeded)}
#'   \item{epsilon_n}{Rashomon bound multiplier used}
#'   \item{n_trees_averaged}{Number of trees averaged (always K)}
#'
#' @details
#' \strong{Algorithm:}
#' \enumerate{
#'   \item Run Rashomon intersection via \code{fit_nuisances_rashomon()} for both e and m0
#'   \item Check intersection succeeded (n_intersecting > 0) for both nuisances
#'   \item Extract K trees with common structure from each nuisance
#'   \item Average leaf values across K trees using \code{average_trees()}
#'   \item Predict for ALL observations with averaged trees (no cross-fitting)
#'   \item Compute ATT via EIF
#' }
#'
#' \strong{Key differences from old implementation:}
#' - Uses Rashomon intersection directly (no tier system)
#' - Fails loudly if intersection is empty (no fallback)
#' - Averages exactly K trees (one per fold, all with same structure)
#' - No cross-fitting in final predictions (all obs use averaged tree)
#'
#' \strong{Error handling:}
#' If Rashomon intersection is empty for either nuisance, the function stops
#' with an informative error suggesting:
#' - Increase rashomon_bound_multiplier
#' - Use estimate_att() with use_rashomon=FALSE
#' - Use estimate_att_msplit_averaged() (more robust to empty intersection)
#'
#' @seealso \code{\link{estimate_att}}, \code{\link{estimate_att_msplit_averaged}},
#'   \code{\link{fit_nuisances_rashomon}}, \code{\link{average_trees}}
#'
#' @export
estimate_att_doubletree_averaged <- function(
  X, A, Y,
  K = 5,
  outcome_type = c("binary", "continuous"),
  cv_regularization = TRUE,
  cv_K = 5,
  regularization = 0.1,
  stratified = TRUE,
  seed = NULL,
  verbose = FALSE,
  rashomon_bound_multiplier = 0.05,
  rashomon_bound_adder = 0,
  max_leaves = NULL,
  auto_tune_intersecting = FALSE,
  discretize_method = "quantiles",
  discretize_bins = "adaptive",
  ...
) {
  # Input validation
  outcome_type <- match.arg(outcome_type)
  n <- nrow(X)
  check_att_data(X, A, Y, outcome_type = outcome_type)
  if (is.matrix(X)) X <- as.data.frame(X)

  if (verbose) {
    message("=== Approach 4: Doubletree Averaged ===")
    message("Using Rashomon intersection + tree averaging")
  }

  # Create fold indices
  fold_indices <- create_folds(n, K, strata = if (stratified) A else NULL, seed = seed)

  # Fit both nuisances with Rashomon intersection
  if (verbose) message("\n--- Fitting nuisances with Rashomon intersection ---")

  nuisance_fits <- tryCatch({
    fit_nuisances_rashomon(
      X = X,
      A = A,
      Y = Y,
      fold_indices = fold_indices,
      outcome_type = outcome_type,
      regularization = regularization,
      cv_regularization = cv_regularization,
      cv_K = cv_K,
      verbose = verbose,
      rashomon_bound_multiplier = rashomon_bound_multiplier,
      rashomon_bound_adder = rashomon_bound_adder,
      max_leaves = max_leaves,
      auto_tune_intersecting = auto_tune_intersecting,
      discretize_method = discretize_method,
      discretize_bins = discretize_bins,
      ...
    )
  }, error = function(e) {
    stop("Failed to fit nuisance models: ", e$message, call. = FALSE)
  })

  # Extract Rashomon objects
  cf_e <- nuisance_fits$cf_e
  cf_m0 <- nuisance_fits$cf_m0

  # Check intersection succeeded for propensity
  if (is.null(cf_e) || cf_e@n_intersecting == 0) {
    stop(
      "Rashomon intersection empty for propensity.\n",
      "Suggestions:\n",
      "  1. Increase rashomon_bound_multiplier (current: ", rashomon_bound_multiplier, ")\n",
      "  2. Use estimate_att() with use_rashomon=FALSE\n",
      "  3. Use estimate_att_msplit_averaged() (more robust to empty intersection)",
      call. = FALSE
    )
  }

  # Check intersection succeeded for outcome
  if (is.null(cf_m0) || cf_m0@n_intersecting == 0) {
    stop(
      "Rashomon intersection empty for outcome.\n",
      "Suggestions:\n",
      "  1. Increase rashomon_bound_multiplier (current: ", rashomon_bound_multiplier, ")\n",
      "  2. Use estimate_att() with use_rashomon=FALSE\n",
      "  3. Use estimate_att_msplit_averaged() (more robust to empty intersection)",
      call. = FALSE
    )
  }

  # Extract K trees with common structure from each fold
  if (verbose) message("\n--- Extracting K trees from each nuisance ---")

  e_trees <- extract_k_trees_from_rashomon(cf_e)
  m0_trees <- extract_k_trees_from_rashomon(cf_m0)

  if (verbose) {
    message(sprintf("Extracted %d propensity trees", length(e_trees)))
    message(sprintf("Extracted %d outcome trees", length(m0_trees)))
  }

  # Average trees
  if (verbose) message("\n--- Averaging leaf values ---")

  # TODO: Update to use actual leaf counts from fold_refits
  # For now, use uniform weights (unweighted averaging)
  uniform_weights_e <- lapply(e_trees, function(tree) {
    leaf_values <- extract_leaf_values(tree)
    setNames(rep(1L, length(leaf_values)), names(leaf_values))
  })

  uniform_weights_m0 <- lapply(m0_trees, function(tree) {
    leaf_values <- extract_leaf_values(tree)
    setNames(rep(1L, length(leaf_values)), names(leaf_values))
  })

  e_averaged <- tryCatch({
    average_trees(e_trees, uniform_weights_e)
  }, error = function(e) {
    stop("Failed to average propensity trees: ", e$message, call. = FALSE)
  })

  m0_averaged <- tryCatch({
    average_trees(m0_trees, uniform_weights_m0)
  }, error = function(e) {
    stop("Failed to average outcome trees: ", e$message, call. = FALSE)
  })

  # Predict for ALL observations (no cross-fitting)
  if (verbose) message("\n--- Generating predictions with averaged trees ---")

  # Trees were built on discretized binary features (X_binary).
  # Must apply the same discretization to X before prediction.
  apply_disc <- get("apply_discretization", envir = asNamespace("optimaltrees"))
  X_for_e_pred <- if (!is.null(cf_e@disc_metadata)) {
    apply_disc(X, cf_e@disc_metadata)
  } else {
    X
  }
  X_for_m0_pred <- if (!is.null(cf_m0@disc_metadata)) {
    apply_disc(X, cf_m0@disc_metadata)
  } else {
    X
  }

  e_hat <- predict_from_tree(e_averaged, X_for_e_pred)
  m0_hat <- predict_from_tree(m0_averaged, X_for_m0_pred)

  # Compute ATT via EIF
  if (verbose) message("\n--- Computing ATT estimate ---")

  pi_hat <- mean(A)
  eta <- list(e = e_hat, m0 = m0_hat, m1 = NULL)

  # Solve for theta
  score_at_zero <- psi_att(Y, A, theta = 0, eta, pi_hat)
  sum_a_over_pi <- sum(A / pi_hat)
  theta_hat <- sum(score_at_zero) / sum_a_over_pi

  # Compute standard error
  score_values <- psi_att(Y, A, theta_hat, eta, pi_hat)
  sigma <- sqrt(mean((score_values - mean(score_values))^2) / n)

  # Confidence interval
  ci_95 <- att_ci(theta_hat, sigma)

  if (verbose) {
    message(sprintf("\n=== Results ==="))
    message(sprintf("ATT estimate: %.4f", theta_hat))
    message(sprintf("Standard error: %.4f", sigma))
    message(sprintf("95%% CI: [%.4f, %.4f]", ci_95[1], ci_95[2]))
  }

  # Return structure
  structure(list(
    theta = theta_hat,
    sigma = sigma,
    ci_95 = ci_95,
    score_values = score_values,
    e_hat = e_hat,
    m0_hat = m0_hat,
    averaged_trees = list(e = e_averaged, m0 = m0_averaged),
    structures = NULL,  # Averaged trees contain full structure information
    n = n,
    K = K,
    n_treated = sum(A),
    outcome_type = outcome_type,
    converged = TRUE,
    epsilon_n = rashomon_bound_multiplier,
    n_trees_averaged = K
  ), class = c("doubletree_att_averaged", "list"))
}


# =============================================================================
# Approach 6: M-Split Averaged
# =============================================================================

#' M-Split ATT with Tree Averaging (Approach 6)
#'
#' @description
#' Estimates ATT using M-split with tree averaging: Find modal structure across
#' M×K trees, refit that structure M×K times, average leaf values across all
#' M×K trees to create single averaged tree.
#'
#' @param X Data.frame or matrix of covariates
#' @param A Treatment vector (0/1)
#' @param Y Outcome vector
#' @param M Number of independent sample splits (default 10)
#' @param K Number of cross-validation folds per split (default 5)
#' @param seed_base Base seed for reproducibility
#' @param verbose Logical
#' @param outcome_type "binary" or "continuous"
#'
#' @return List with class "msplit_att_averaged"
#'
#' @details
#' Algorithm:
#' 1. Stage 1: Find modal structure (same as approach 5)
#' 2. Stage 2: Refit modal structure M×K times AND store cross-fitted fold predictions
#' 3. Stage 3: **Average leaf values** across all M×K trees → 1 averaged tree
#' 4. **Inference** via the averaged tree (predict all n observations)
#'
#' The averaged tree is used for both inference and display. Cross-fitted predictions
#' from Stage 2 are also stored in \code{predictions_all_splits} to enable investigation
#' of alternative inference paths (e.g., rowMeans of the n×M prediction matrix).
#'
#' \strong{Known issue:} Using the averaged tree for inference introduces (K-1)/K
#' in-sample contamination (80% at K=5) because each leaf's averaged value is
#' computed from M×K refits, M(K-1) of which used observation i in training.
#' This causes structural positive bias that does not vanish with M or n.
#' The \code{predictions_all_splits} field enables comparison with cross-fitted inference.
#'
#' @export
estimate_att_msplit_averaged <- function(X, A, Y,
                                            M = 10,
                                            K = 5,
                                            seed_base = NULL,
                                            verbose = TRUE,
                                            outcome_type = "binary") {
  n <- nrow(X)

  if (!is.data.frame(X) && !is.matrix(X)) {
    stop("X must be a data.frame or matrix", call. = FALSE)
  }
  if (length(A) != n || length(Y) != n) {
    stop("length(A) and length(Y) must equal nrow(X)", call. = FALSE)
  }

  if (verbose) {
    cat(sprintf("M-split averaged: M=%d splits, K=%d folds, n=%d\n", M, K, n))
  }

  # ============================================================
  # Stage 1: Structure Selection (same as approach 5)
  # ============================================================
  if (verbose) cat("Stage 1: Selecting modal structures...\n")

  structures_e <- vector("list", M)
  structures_m0 <- vector("list", M)

  for (m in seq_len(M)) {
    seed_m <- if (!is.null(seed_base)) seed_base + m else NULL
    folds_m <- create_folds(n, K, strata = A, seed = seed_m)

    # Use fold 1's training set for structure discovery
    train_idx <- which(folds_m != 1)
    X_train <- X[train_idx, , drop = FALSE]
    A_train <- A[train_idx]
    Y_train <- Y[train_idx]

    # Cap lambda at 15× theory value to prevent stump-producing regularization.
    # Without this cap, cv_regularization_adaptive can select lambda large enough
    # that all splits are pruned, yielding constant (stump) predictions which bias DML.
    n_train <- length(train_idx)
    max_lambda_cap <- (log(n_train) / n_train) * 15

    # Fit propensity with adaptive CV
    cv_e <- optimaltrees::cv_regularization_adaptive(
      X = X_train, y = A_train, loss_function = "log_loss",
      K = 5, max_iterations = 10, refit = TRUE, verbose = FALSE,
      max_lambda = max_lambda_cap,
      discretize_bins = "adaptive",
      discretize_method = "quantiles"
    )

    if (is.na(cv_e$best_lambda)) {
      stop("CV failed for propensity in split ", m, call. = FALSE)
    }

    model_e <- cv_e$model
    structures_e[[m]] <- list(
      structure = optimaltrees::extract_tree_structure(model_e),
      discretization_metadata = model_e@discretization_metadata
    )

    # Fit outcome on controls with CV
    control_idx <- which(A_train == 0)
    X_control <- X_train[control_idx, , drop = FALSE]
    Y_control <- Y_train[control_idx]

    outcome_loss <- if (outcome_type == "binary") "log_loss" else "squared_error"

    cv_m0 <- optimaltrees::cv_regularization_adaptive(
      X = X_control, y = Y_control, loss_function = outcome_loss,
      K = 5, max_iterations = 10, refit = TRUE, verbose = FALSE,
      max_lambda = max_lambda_cap,
      discretize_bins = "adaptive",
      discretize_method = "quantiles"
    )

    if (is.na(cv_m0$best_lambda)) {
      stop("CV failed for outcome in split ", m, call. = FALSE)
    }

    model_m0 <- cv_m0$model
    structures_m0[[m]] <- list(
      structure = optimaltrees::extract_tree_structure(model_m0),
      discretization_metadata = model_m0@discretization_metadata
    )

    if (verbose && m %% max(1, M %/% 10) == 0) {
      cat(sprintf("  Structure discovery: %d/%d splits\n", m, M))
    }
  }

  # Select modal structures, excluding stumps (min_leaves = 2): prevents
  # degenerate constant-prediction modal structures from biasing DML estimates.
  s_star_e <- select_structure_modal(structures_e, min_leaves = 2L)
  s_star_m0 <- select_structure_modal(structures_m0, min_leaves = 2L)

  if (verbose) {
    cat(sprintf("  Modal propensity: %.1f%% agreement\n", s_star_e$frequency * 100))
    cat(sprintf("  Modal outcome: %.1f%% agreement\n", s_star_m0$frequency * 100))
  }

  # ============================================================
  # Stage 2: Refit Modal Structure M×K Times
  # ============================================================
  if (verbose) cat("Stage 2: Refitting modal structure M×K times...\n")

  # Store all M×K trees (as tree structure objects, not models)
  trees_e <- vector("list", M * K)
  trees_m0 <- vector("list", M * K)
  leaf_counts_e <- vector("list", M * K)
  leaf_counts_m0 <- vector("list", M * K)
  tree_idx <- 1

  # Cross-fitted predictions (n × M): for valid inference (Option B).
  # Each column m is filled by the K fold test-set predictions from split m.
  # These are proper cross-fitted predictions: observation i is predicted by
  # a tree trained on data that does NOT include i, fixing the in-sample
  # contamination that biases the averaged-tree inference path.
  predictions_e_xfitted  <- matrix(NA_real_, nrow = n, ncol = M)
  predictions_m0_xfitted <- matrix(NA_real_, nrow = n, ncol = M)

  for (m in seq_len(M)) {
    seed_m <- if (!is.null(seed_base)) seed_base + m else NULL
    folds_m <- create_folds(n, K, strata = A, seed = seed_m)

    for (k in seq_len(K)) {
      test_idx <- which(folds_m == k)
      train_idx <- which(folds_m != k)

      # Refit propensity tree with modal structure
      X_train <- X[train_idx, , drop = FALSE]
      A_train <- A[train_idx]

      refit_result_e <- optimaltrees::refit_tree_structure(
        structure = s_star_e$structure,
        X_new = X_train,
        y_new = A_train,
        loss_function = "log_loss",
        store_training_data = FALSE,
        discretization_metadata = s_star_e$discretization_metadata,
        allow_partial_leaves = TRUE
      )

      # Extract tree structure and leaf counts
      trees_e[[tree_idx]] <- refit_result_e$model@trees[[1]]
      leaf_counts_e[[tree_idx]] <- refit_result_e$n_per_leaf

      # Refit outcome tree with modal structure
      control_idx <- which(A_train == 0)
      Y_train_control <- Y[train_idx][A_train == 0]
      X_train_control <- X_train[control_idx, , drop = FALSE]

      outcome_loss <- if (outcome_type == "binary") "log_loss" else "squared_error"

      refit_result_m0 <- optimaltrees::refit_tree_structure(
        structure = s_star_m0$structure,
        X_new = X_train_control,
        y_new = Y_train_control,
        loss_function = outcome_loss,
        store_training_data = FALSE,
        discretization_metadata = s_star_m0$discretization_metadata,
        allow_partial_leaves = TRUE
      )

      # Extract tree structure and leaf counts
      trees_m0[[tree_idx]] <- refit_result_m0$model@trees[[1]]
      leaf_counts_m0[[tree_idx]] <- refit_result_m0$n_per_leaf

      # Store cross-fitted predictions on test fold for valid inference.
      # predict() uses the discretization_metadata stored on the model object,
      # so X_test (original X) is handled correctly without explicit discretization.
      X_test <- X[test_idx, , drop = FALSE]

      preds_e_mk <- predict(refit_result_e$model, X_test, type = "prob")
      predictions_e_xfitted[test_idx, m] <- preds_e_mk[, 2L]

      if (outcome_loss == "log_loss") {
        preds_m0_mk <- predict(refit_result_m0$model, X_test, type = "prob")
        predictions_m0_xfitted[test_idx, m] <- preds_m0_mk[, 2L]
      } else {
        preds_m0_mk <- predict(refit_result_m0$model, X_test)
        predictions_m0_xfitted[test_idx, m] <- preds_m0_mk
      }

      tree_idx <- tree_idx + 1
    }

    if (verbose && m %% max(1, M %/% 10) == 0) {
      cat(sprintf("  Refitting: %d/%d splits complete\n", m, M))
    }
  }

  # ============================================================
  # Stage 3: Average Leaf Values Across All M×K Trees (Weighted)
  # ============================================================
  if (verbose) cat("Stage 3: Averaging leaf values (weighted by sample size)...\n")

  e_averaged <- average_trees(trees_e, leaf_counts_e)
  m0_averaged <- average_trees(trees_m0, leaf_counts_m0)

  # ============================================================
  # Stage 4: Predict and Compute ATT
  # ============================================================
  if (verbose) cat("Stage 4: Computing ATT with averaged tree...\n")

  # Predict for ALL observations using the averaged tree.
  # Trees use discretized binary features; apply_discretization converts X to
  # the same binary feature space before calling predict_from_tree.
  apply_disc <- get("apply_discretization", envir = asNamespace("optimaltrees"))
  X_for_e_pred <- if (!is.null(s_star_e$discretization_metadata)) {
    apply_disc(X, s_star_e$discretization_metadata)
  } else {
    X
  }
  X_for_m0_pred <- if (!is.null(s_star_m0$discretization_metadata)) {
    apply_disc(X, s_star_m0$discretization_metadata)
  } else {
    X
  }

  e_hat  <- pmax(pmin(predict_from_tree(e_averaged,  X_for_e_pred),  0.99), 0.01)
  m0_hat <- pmax(pmin(predict_from_tree(m0_averaged, X_for_m0_pred), 0.99), 0.01)

  # Compute ATT
  pi_hat <- mean(A)
  eta <- list(e = e_hat, m0 = m0_hat, m1 = NULL)
  score <- psi_att(Y, A, theta = 0, eta = eta, pi_hat = pi_hat)
  theta_hat <- sum(score) / sum(A / pi_hat)

  # Standard error
  score_centered <- score - mean(score)
  sigma <- sqrt(mean(score_centered^2) / n)

  # CI
  ci_95 <- att_ci(theta_hat, sigma)

  # ============================================================
  # Return
  # ============================================================
  structure(list(
    theta = theta_hat,
    sigma = sigma,
    ci_95 = ci_95,
    score_values = score,
    e_hat = e_hat,
    m0_hat = m0_hat,
    averaged_trees = list(e = e_averaged, m0 = m0_averaged),
    predictions_all_splits = list(e = predictions_e_xfitted, m0 = predictions_m0_xfitted),
    structures = list(
      modal = list(e = s_star_e$structure, m0 = s_star_m0$structure)
    ),
    structure_selection = list(
      modal_freq_e = s_star_e$frequency,
      modal_freq_m0 = s_star_m0$frequency
    ),
    M = M,
    K = K,
    n = n,
    n_treated = sum(A),
    outcome_type = outcome_type,
    n_trees_averaged = M * K
  ), class = c("msplit_att_averaged", "list"))
}
