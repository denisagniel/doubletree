# Averaged-Tree ATT Estimators
#
# Alternative inference methods that average leaf values across cross-fitted trees
# to produce a single interpretable tree while maintaining cross-fit validity.

# =============================================================================
# Helper Functions for Robust Fallback
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
#' @param regularization Numeric. Fallback tree complexity penalty if CV fails
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

    # Select lambda via CV on training fold (2026-05-26)
    loss_fn <- if (outcome_type == "binary") "log_loss" else "squared_error"
    cv_result <- tryCatch({
      optimaltrees::cv_regularization(
        X = X_train,
        y = outcome_train,
        loss_function = loss_fn,
        K = 5,
        refit = FALSE,
        verbose = FALSE
      )
    }, error = function(e) {
      if (verbose) message("Fold ", k, " CV failed: ", e$message)
      return(list(best_lambda = NA_real_))
    })

    # Use CV-selected lambda or fallback to fixed
    lambda_k <- if (!is.na(cv_result$best_lambda)) {
      cv_result$best_lambda
    } else {
      if (verbose) message("Fold ", k, " using fallback regularization = ", regularization)
      regularization
    }

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

  # Tier 5: Fold-specific fallback with CV-selected lambda
  if (verbose) message(sprintf("%s: Tier 5 - Fold-specific trees (no Rashomon)", nuisance_name))

  fold_trees <- list()
  for (k in 1:K) {
    test_idx <- fold_indices[[k]]
    train_idx <- setdiff(1:n, test_idx)

    X_train <- X[train_idx, , drop = FALSE]
    outcome_train <- outcome[train_idx]

    # Select lambda via CV (2026-05-26)
    loss_fn <- if (outcome_type == "binary") "log_loss" else "squared_error"
    cv_result <- tryCatch({
      optimaltrees::cv_regularization(
        X = X_train,
        y = outcome_train,
        loss_function = loss_fn,
        K = 5,
        refit = TRUE,
        verbose = FALSE
      )
    }, error = function(e) {
      if (verbose) message("Fold ", k, " CV failed: ", e$message)
      return(list(best_lambda = NA_real_, model = NULL))
    })

    # Use CV-selected model or fit with fallback lambda
    tree_k <- if (!is.na(cv_result$best_lambda) && !is.null(cv_result$model)) {
      cv_result$model
    } else {
      if (verbose) message("Fold ", k, " using fallback regularization = ", regularization)
      tryCatch({
        optimaltrees::fit_tree(
          X = X_train,
          y = outcome_train,
          loss_function = loss_fn,
          regularization = regularization,
          verbose = FALSE,
          ...
        )
      }, error = function(e) {
        if (verbose) message("Fold ", k, " tree fit failed: ", e$message)
        return(NULL)
      })
    }

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

#' Estimate ATT with Doubletree and Averaged Leaves
#'
#' Estimates the Average Treatment Effect on the Treated (ATT) using doubletree
#' (Rashomon intersection) with robust five-tier fallback to find a common tree
#' structure, then averages leaf values across K cross-fitted trees to create a
#' single interpretable tree. Uses CV-selected lambda for all tree fits (2026-05-26).
#'
#' @inheritParams estimate_att
#' @param K Number of cross-fitting folds. Default 5.
#' @param regularization Numeric. Fallback tree complexity penalty if CV fails. Default 0.1.
#'
#' @return List with elements:
#'   \item{theta}{ATT point estimate}
#'   \item{sigma}{Standard error}
#'   \item{ci_95}{95\% confidence interval}
#'   \item{e_hat}{Propensity score predictions (from averaged tree)}
#'   \item{m0_hat}{Control outcome predictions (from averaged tree)}
#'   \item{structures}{List with e and m0 structures}
#'   \item{averaged_trees}{List with e and m0 averaged trees}
#'   \item{n_folds}{Number of folds (K)}
#'   \item{tier_used}{List with e and m0 tier information}
#'   \item{rashomon_tolerance_used}{List with e and m0 tolerance values}
#'   \item{theory_warning}{Character or NULL (if any nuisance used ε > 0.05)}
#'   \item{n_trees_matched_modal}{List with e and m0 (for modal tiers only)}
#'   \item{modal_frequency}{List with e and m0 (for modal tiers only)}
#'
#' @details
#' This estimator combines the interpretability of a single tree with the validity
#' of cross-fitting, using a robust five-tier fallback strategy:
#'
#' \strong{Five-Tier Fallback (Per Nuisance):}
#' \enumerate{
#'   \item Rashomon intersection at ε = 0.05 (theory-compliant)
#'   \item Expanded intersection at ε = 0.10, 0.15, 0.20 (with warning)
#'   \item Modal structure across K Rashomon sets at ε = 0.05 (theory-compliant)
#'   \item Modal structure with expanded ε = 0.10, 0.15, 0.20 (with warning)
#'   \item Fold-specific trees with modal structure selection (no Rashomon)
#' }
#'
#' \strong{Key features:}
#' - Each nuisance (e, m0) progresses through tiers independently
#' - Theory warnings issued when ε > 0.05 used (violates rate condition)
#' - Minimum 3 trees required for averaging (ensures robustness)
#' - All fallback paths are explicit (no silent failures)
#'
#' @seealso \code{\link{estimate_att}}, \code{\link{estimate_att_msplit_averaged}},
#'   \code{\link{average_trees}}
#'
#' @export
estimate_att_doubletree_averaged <- function(X, A, Y, K = 5, regularization = 0.1,
                                             outcome_type = c("binary", "continuous"),
                                             rashomon_bound_multiplier = 0.05,
                                             verbose = FALSE, ...) {
  outcome_type <- match.arg(outcome_type)
  n <- nrow(X)
  check_att_data(X, A, Y, outcome_type = outcome_type)
  if (is.matrix(X)) X <- as.data.frame(X)

  # Create fold indices (stratified by A)
  fold_indices <- create_folds(n, K, strata = A, seed = NULL)

  # Separate control outcomes for m0
  Y_control <- Y[A == 0]
  X_control <- X[A == 0, , drop = FALSE]
  n_control <- length(Y_control)
  fold_indices_control <- lapply(fold_indices, function(idx) {
    idx[idx <= n_control]  # Map to control indices
  })

  # Find trees for propensity (e)
  if (verbose) message("\n=== Fitting propensity (e) ===")
  e_result <- find_trees_through_tiers(
    X = X,
    outcome = A,
    K = K,
    fold_indices = fold_indices,
    nuisance_name = "propensity",
    regularization = regularization,
    outcome_type = "binary",
    verbose = verbose,
    ...
  )

  # Find trees for outcome (m0)
  if (verbose) message("\n=== Fitting outcome (m0) ===")
  m0_result <- find_trees_through_tiers(
    X = X_control,
    outcome = Y_control,
    K = K,
    fold_indices = fold_indices_control,
    nuisance_name = "outcome",
    regularization = regularization,
    outcome_type = outcome_type,
    verbose = verbose,
    ...
  )

  # Check both succeeded
  if (is.null(e_result$trees) || is.null(m0_result$trees)) {
    error_parts <- c()
    if (is.null(e_result$trees)) error_parts <- c(error_parts, "propensity (e) failed: all tiers exhausted")
    if (is.null(m0_result$trees)) error_parts <- c(error_parts, "outcome (m0) failed: all tiers exhausted")

    return(list(
      theta = NA_real_,
      sigma = NA_real_,
      ci_95 = c(NA_real_, NA_real_),
      e_hat = rep(NA_real_, n),
      m0_hat = rep(NA_real_, n),
      error = paste(error_parts, collapse = "; ")
    ))
  }

  # Average trees for each nuisance
  if (verbose) {
    message(sprintf("Averaging %d propensity trees", length(e_result$trees)))
    message(sprintf("Averaging %d outcome trees", length(m0_result$trees)))
  }

  e_averaged <- tryCatch({
    average_trees(e_result$trees)
  }, error = function(e) {
    stop("Failed to average propensity trees: ", e$message,
         "\nNumber of trees: ", length(e_result$trees),
         "\nTier used: ", e_result$tier, call. = FALSE)
  })

  m0_averaged <- tryCatch({
    average_trees(m0_result$trees)
  }, error = function(e) {
    stop("Failed to average outcome trees: ", e$message,
         "\nNumber of trees: ", length(m0_result$trees),
         "\nTier used: ", m0_result$tier, call. = FALSE)
  })

  # Predict using averaged trees
  e_hat <- predict_from_tree(e_averaged, X)
  m0_hat <- predict_from_tree(m0_averaged, X)

  # Compute ATT using the core inference functions
  pi_hat <- mean(A)
  eta <- list(e = e_hat, m0 = m0_hat, m1 = NULL)
  score_at_zero <- psi_att(Y, A, theta = 0, eta, pi_hat)
  sum_a_over_pi <- sum(A / pi_hat)
  theta_hat <- sum(score_at_zero) / sum_a_over_pi

  score_values <- psi_att(Y, A, theta_hat, eta, pi_hat)
  sigma <- att_se(score_values, n)
  ci_95 <- att_ci(theta_hat, sigma, n, level = 0.95)

  # Combine warnings if both nuisances have them
  combined_warning <- NULL
  if (!is.null(e_result$warning) || !is.null(m0_result$warning)) {
    warning_parts <- c()
    if (!is.null(e_result$warning)) warning_parts <- c(warning_parts, paste("e:", e_result$warning))
    if (!is.null(m0_result$warning)) warning_parts <- c(warning_parts, paste("m0:", m0_result$warning))
    combined_warning <- paste(warning_parts, collapse = " | ")
  }

  # Use structures from tier results
  e_structure <- e_result$structure
  m0_structure <- m0_result$structure

  list(
    theta = theta_hat,
    sigma = sigma,
    ci_95 = ci_95,
    score_values = score_values,
    e_hat = e_hat,
    m0_hat = m0_hat,
    structures = list(e = e_structure, m0 = m0_structure),
    averaged_trees = list(e = e_averaged, m0 = m0_averaged),
    n_folds = K,
    n = n,
    converged = TRUE,
    tier_used = list(e = e_result$tier, m0 = m0_result$tier),
    rashomon_tolerance_used = list(e = e_result$tolerance, m0 = m0_result$tolerance),
    theory_warning = combined_warning,
    n_trees_matched_modal = list(e = e_result$n_matched, m0 = m0_result$n_matched),
    n_trees_total = list(e = e_result$n_total, m0 = m0_result$n_total),
    modal_frequency = list(e = e_result$modal_frequency, m0 = m0_result$modal_frequency)
  )
}

# =============================================================================
# Approach 6: M-Split Averaged
# =============================================================================

#' Try M-Split with Fold-Specific Trees (Tier 1)
#'
#' Collect M×K fold-specific trees, find modal structure, filter and average.
#'
#' @keywords internal
try_msplit_fold_specific <- function(X, A, Y, M, K, regularization, outcome_type, verbose, ...) {
  n <- nrow(X)

  # Separate control data for m0
  Y_control <- Y[A == 0]
  X_control <- X[A == 0, , drop = FALSE]

  # Collect M×K fold-specific trees for each nuisance
  all_e_trees <- list()
  all_m0_trees <- list()

  for (m in 1:M) {
    if (verbose && m %% max(1, M %/% 4) == 0) {
      message(sprintf("M-split fold-specific: Split %d/%d", m, M))
    }

    # Create fold indices for this split
    fold_indices <- create_folds(n, K, strata = A, seed = NULL)

    # Fit K trees for propensity with CV-selected lambda (2026-05-26)
    for (k in 1:K) {
      test_idx <- fold_indices[[k]]
      train_idx <- setdiff(1:n, test_idx)

      X_train <- X[train_idx, , drop = FALSE]
      A_train <- A[train_idx]

      # Select lambda via CV
      cv_e <- tryCatch({
        optimaltrees::cv_regularization(
          X = X_train, y = A_train, loss_function = "log_loss",
          K = 5, refit = TRUE, verbose = FALSE
        )
      }, error = function(e) list(best_lambda = NA_real_, model = NULL))

      # Use CV-selected model or fallback
      tree_e <- if (!is.na(cv_e$best_lambda) && !is.null(cv_e$model)) {
        cv_e$model
      } else {
        tryCatch({
          optimaltrees::fit_tree(
            X = X_train, y = A_train, loss_function = "log_loss",
            regularization = regularization, verbose = FALSE, ...
          )
        }, error = function(e) NULL)
      }

      if (!is.null(tree_e)) {
        all_e_trees[[length(all_e_trees) + 1]] <- tree_e
      }
    }

    # Fit K trees for outcome with CV-selected lambda (2026-05-26)
    fold_indices_control <- lapply(fold_indices, function(idx) idx[idx <= nrow(X_control)])

    for (k in 1:K) {
      test_idx <- fold_indices_control[[k]]
      train_idx <- setdiff(1:nrow(X_control), test_idx)

      X_train <- X_control[train_idx, , drop = FALSE]
      Y_train <- Y_control[train_idx]

      # Select lambda via CV
      outcome_loss <- if (outcome_type == "binary") "log_loss" else "squared_error"
      cv_m0 <- tryCatch({
        optimaltrees::cv_regularization(
          X = X_train, y = Y_train, loss_function = outcome_loss,
          K = 5, refit = TRUE, verbose = FALSE
        )
      }, error = function(e) list(best_lambda = NA_real_, model = NULL))

      # Use CV-selected model or fallback
      tree_m0 <- if (!is.na(cv_m0$best_lambda) && !is.null(cv_m0$model)) {
        cv_m0$model
      } else {
        tryCatch({
          optimaltrees::fit_tree(
            X = X_train, y = Y_train, loss_function = outcome_loss,
            regularization = regularization, verbose = FALSE, ...
          )
        }, error = function(e) NULL)
      }

      if (!is.null(tree_m0)) {
        all_m0_trees[[length(all_m0_trees) + 1]] <- tree_m0
      }
    }
  }

  # Find modal structures
  if (length(all_e_trees) == 0 || length(all_m0_trees) == 0) {
    return(list(success = FALSE, reason = "No trees collected"))
  }

  # Extract structures from model objects
  e_structures <- lapply(all_e_trees, optimaltrees::extract_tree_structure)
  m0_structures <- lapply(all_m0_trees, optimaltrees::extract_tree_structure)

  e_modal <- select_structure_modal(e_structures)
  m0_modal <- select_structure_modal(m0_structures)

  # Convert models to tree structures (nested lists) for averaging
  e_tree_structures <- lapply(all_e_trees, function(model) {
    optimaltrees::get_rashomon_trees(model)[[1]]
  })
  m0_tree_structures <- lapply(all_m0_trees, function(model) {
    optimaltrees::get_rashomon_trees(model)[[1]]
  })

  # Filter to matching trees
  e_filtered <- filter_trees_by_structure(e_tree_structures, e_structures, e_modal$structure)
  m0_filtered <- filter_trees_by_structure(m0_tree_structures, m0_structures, m0_modal$structure)

  list(
    success = TRUE,
    e_trees = e_filtered$matched_trees,
    m0_trees = m0_filtered$matched_trees,
    n_matched_e = e_filtered$n_matched,
    n_matched_m0 = m0_filtered$n_matched,
    n_total_e = length(all_e_trees),
    n_total_m0 = length(all_m0_trees),
    modal_freq_e = e_modal$frequency,
    modal_freq_m0 = m0_modal$frequency,
    structures = list(e = e_modal$structure, m0 = m0_modal$structure)
  )
}

#' Try M-Split with Rashomon Sets (Tier 2)
#'
#' Collect all trees from M×K Rashomon sets, find modal structure, filter and average.
#'
#' @keywords internal
try_msplit_rashomon_sets <- function(X, A, Y, M, K, regularization, outcome_type,
                                     tolerance, verbose, ...) {
  n <- nrow(X)

  Y_control <- Y[A == 0]
  X_control <- X[A == 0, , drop = FALSE]

  all_e_trees <- list()
  all_e_structures <- list()
  all_m0_trees <- list()
  all_m0_structures <- list()

  for (m in 1:M) {
    if (verbose && m %% max(1, M %/% 4) == 0) {
      message(sprintf("M-split Rashomon: Split %d/%d at ε = %.2f", m, M, tolerance))
    }

    fold_indices <- create_folds(n, K, strata = A, seed = NULL)

    # Collect Rashomon trees for propensity
    e_result <- collect_rashomon_trees_at_tolerance(
      X, A, K, fold_indices, tolerance, regularization, "binary", verbose = FALSE, ...
    )
    all_e_trees <- c(all_e_trees, e_result$all_trees)
    all_e_structures <- c(all_e_structures, e_result$all_structures)

    # Collect Rashomon trees for outcome
    fold_indices_control <- lapply(fold_indices, function(idx) idx[idx <= nrow(X_control)])
    m0_result <- collect_rashomon_trees_at_tolerance(
      X_control, Y_control, K, fold_indices_control, tolerance, regularization, outcome_type,
      verbose = FALSE, ...
    )
    all_m0_trees <- c(all_m0_trees, m0_result$all_trees)
    all_m0_structures <- c(all_m0_structures, m0_result$all_structures)
  }

  if (length(all_e_trees) == 0 || length(all_m0_trees) == 0) {
    return(list(success = FALSE, reason = "No Rashomon trees collected"))
  }

  # Find modal structures using pre-extracted structures
  e_modal <- select_structure_modal(all_e_structures)
  m0_modal <- select_structure_modal(all_m0_structures)

  # Filter to matching trees
  e_filtered <- filter_trees_by_structure(all_e_trees, all_e_structures, e_modal$structure)
  m0_filtered <- filter_trees_by_structure(all_m0_trees, all_m0_structures, m0_modal$structure)

  list(
    success = TRUE,
    e_trees = e_filtered$matched_trees,
    m0_trees = m0_filtered$matched_trees,
    n_matched_e = e_filtered$n_matched,
    n_matched_m0 = m0_filtered$n_matched,
    n_total_e = length(all_e_trees),
    n_total_m0 = length(all_m0_trees),
    modal_freq_e = e_modal$frequency,
    modal_freq_m0 = m0_modal$frequency,
    structures = list(e = e_modal$structure, m0 = m0_modal$structure)
  )
}

#' Finalize M-Split Result
#'
#' Average filtered trees, compute ATT, add diagnostics.
#'
#' @keywords internal
finalize_msplit_result <- function(result, X, A, Y, tier, M_used, M_expanded, verbose) {
  n <- nrow(X)

  # Check minimum threshold
  if (result$n_matched_e < 3 || result$n_matched_m0 < 3) {
    return(list(
      theta = NA_real_,
      sigma = NA_real_,
      ci_95 = c(NA_real_, NA_real_),
      e_hat = rep(NA_real_, n),
      m0_hat = rep(NA_real_, n),
      error = sprintf("Insufficient trees: e=%d, m0=%d (need ≥3 each)",
                      result$n_matched_e, result$n_matched_m0),
      tier_used = tier,
      M_used = M_used,
      M_expanded = M_expanded
    ))
  }

  # Issue warning if low agreement but ≥3
  if (result$modal_freq_e < 0.5 || result$modal_freq_m0 < 0.5) {
    low_agreement_msg <- sprintf(
      "Low modal agreement: e=%.1f%%, m0=%.1f%%. Consider increasing M.",
      result$modal_freq_e * 100, result$modal_freq_m0 * 100
    )
    if (verbose) message("WARNING: ", low_agreement_msg)
  }

  # Average trees
  e_averaged <- average_trees(result$e_trees)
  m0_averaged <- average_trees(result$m0_trees)

  # Predict
  e_hat <- predict_from_tree(e_averaged, X)
  m0_hat <- predict_from_tree(m0_averaged, X)

  # Compute ATT
  pi_hat <- mean(A)
  eta <- list(e = e_hat, m0 = m0_hat, m1 = NULL)
  score_at_zero <- psi_att(Y, A, theta = 0, eta, pi_hat)
  sum_a_over_pi <- sum(A / pi_hat)
  theta_hat <- sum(score_at_zero) / sum_a_over_pi

  score_values <- psi_att(Y, A, theta_hat, eta, pi_hat)
  sigma <- att_se(score_values, n)
  ci_95 <- att_ci(theta_hat, sigma, n, level = 0.95)

  list(
    theta = theta_hat,
    sigma = sigma,
    ci_95 = ci_95,
    score_values = score_values,
    e_hat = e_hat,
    m0_hat = m0_hat,
    structures = result$structures,
    averaged_trees = list(e = e_averaged, m0 = m0_averaged),
    n_trees_matched_modal = list(e = result$n_matched_e, m0 = result$n_matched_m0),
    n_trees_total = list(e = result$n_total_e, m0 = result$n_total_m0),
    modal_frequency = list(e = result$modal_freq_e, m0 = result$modal_freq_m0),
    tier_used = tier,
    M_used = M_used,
    M_expanded = M_expanded,
    theory_warning = NULL,  # Tier 2 always uses ε = 0.05
    n = n,
    converged = TRUE
  )
}

#' Estimate ATT with M-Split and Averaged Leaves
#'
#' Estimates the ATT using M-split doubletree to find a modal tree structure
#' across M independent splits, then averages leaf values across all M×K
#' cross-fitted trees to create a single interpretable tree. Uses CV-selected
#' lambda for all tree fits (2026-05-26).
#'
#' @inheritParams estimate_att
#' @param M Number of independent splits. Default 10.
#' @param K Number of folds per split. Default 5.
#' @param regularization Numeric. Fallback tree complexity penalty if CV fails. Default 0.1.
#'
#' @return List with elements:
#'   \item{theta}{ATT point estimate}
#'   \item{sigma}{Standard error}
#'   \item{ci_95}{95\% confidence interval}
#'   \item{e_hat}{Propensity score predictions (from averaged tree)}
#'   \item{m0_hat}{Control outcome predictions (from averaged tree)}
#'   \item{structures}{List with e and m0 modal structures}
#'   \item{averaged_trees}{List with e and m0 averaged trees}
#'   \item{n_trees_matched_modal}{List with number of trees matching modal (e and m0)}
#'   \item{n_trees_total}{List with total trees collected (e and m0)}
#'   \item{modal_frequency}{List with modal frequency (e and m0)}
#'   \item{tier_used}{Character indicating which tier succeeded}
#'   \item{M_used}{Integer showing final M value (may be > input M if expanded)}
#'   \item{M_expanded}{Logical indicating if M was increased}
#'   \item{n_splits}{Number of splits (M from input)}
#'   \item{n_folds}{Number of folds per split (K)}
#'
#' @details
#' This estimator combines the stability of M-split with the interpretability of
#' a single tree, using a three-tier fallback strategy that increases M rather
#' than changing tolerance:
#'
#' \strong{Three-Tier Fallback:}
#' \enumerate{
#'   \item Modal across M×K fold-specific trees (no Rashomon)
#'   \item Modal across M×K Rashomon set trees at ε = 0.05 (theory-compliant)
#'   \item Increase M and retry Tiers 1-2
#' }
#'
#' \strong{Key features:}
#' - Leverages M-split philosophy: more splits provide more diverse trees
#' - Always uses ε = 0.05 for Rashomon (no theory violations)
#' - Dynamically increases M if modal structure doesn't emerge
#' - Minimum 3 trees required for averaging
#'
#' @seealso \code{\link{estimate_att_doubletree_averaged}}, \code{\link{estimate_att_msplit}},
#'   \code{\link{average_trees}}
#'
#' @export
estimate_att_msplit_averaged <- function(X, A, Y, M = 10, K = 5, regularization = 0.1,
                                        outcome_type = c("binary", "continuous"),
                                        rashomon_bound_multiplier = 0.05,
                                        verbose = FALSE, ...) {
  outcome_type <- match.arg(outcome_type)
  n <- nrow(X)
  check_att_data(X, A, Y, outcome_type = outcome_type)
  if (is.matrix(X)) X <- as.data.frame(X)

  M_original <- M

  # Tier 1: Modal across M×K fold-specific trees
  if (verbose) message(sprintf("\n=== Tier 1: Fold-specific with M=%d ===", M))
  tier1_result <- try_msplit_fold_specific(X, A, Y, M, K, regularization, outcome_type, verbose, ...)

  if (tier1_result$success && tier1_result$n_matched_e >= 3 && tier1_result$n_matched_m0 >= 3) {
    if (verbose) message(sprintf("Tier 1 succeeded: e=%d/%d (%.1f%%), m0=%d/%d (%.1f%%)",
                                  tier1_result$n_matched_e, tier1_result$n_total_e,
                                  tier1_result$modal_freq_e * 100,
                                  tier1_result$n_matched_m0, tier1_result$n_total_m0,
                                  tier1_result$modal_freq_m0 * 100))
    return(finalize_msplit_result(tier1_result, X, A, Y,
                                   tier = "tier1_fold_specific",
                                   M_used = M, M_expanded = FALSE, verbose))
  }

  # Tier 2: Modal across M×K Rashomon set trees at ε = 0.05
  if (verbose) message(sprintf("\n=== Tier 2: Rashomon sets at ε=0.05 with M=%d ===", M))
  tier2_result <- try_msplit_rashomon_sets(X, A, Y, M, K, regularization, outcome_type,
                                          tolerance = 0.05, verbose, ...)

  if (tier2_result$success && tier2_result$n_matched_e >= 3 && tier2_result$n_matched_m0 >= 3) {
    if (verbose) message(sprintf("Tier 2 succeeded: e=%d/%d (%.1f%%), m0=%d/%d (%.1f%%)",
                                  tier2_result$n_matched_e, tier2_result$n_total_e,
                                  tier2_result$modal_freq_e * 100,
                                  tier2_result$n_matched_m0, tier2_result$n_total_m0,
                                  tier2_result$modal_freq_m0 * 100))
    return(finalize_msplit_result(tier2_result, X, A, Y,
                                   tier = "tier2_rashomon_0.05",
                                   M_used = M, M_expanded = FALSE, verbose))
  }

  # Tier 3: Increase M and retry
  M_increment <- max(5, M %/% 2)
  M_max <- 50

  while (M < M_max) {
    M_new <- M + M_increment
    if (verbose) message(sprintf("\n=== Tier 3: Increasing M from %d to %d ===", M, M_new))

    # Retry Tier 1 with larger M
    tier1_expanded <- try_msplit_fold_specific(X, A, Y, M_new, K, regularization, outcome_type, verbose, ...)
    if (tier1_expanded$success && tier1_expanded$n_matched_e >= 3 && tier1_expanded$n_matched_m0 >= 3) {
      if (verbose) message(sprintf("Tier 3 (fold-specific) succeeded with M=%d", M_new))
      return(finalize_msplit_result(tier1_expanded, X, A, Y,
                                     tier = "tier1_fold_specific",
                                     M_used = M_new, M_expanded = TRUE, verbose))
    }

    # Retry Tier 2 with larger M
    tier2_expanded <- try_msplit_rashomon_sets(X, A, Y, M_new, K, regularization, outcome_type,
                                              tolerance = 0.05, verbose, ...)
    if (tier2_expanded$success && tier2_expanded$n_matched_e >= 3 && tier2_expanded$n_matched_m0 >= 3) {
      if (verbose) message(sprintf("Tier 3 (Rashomon) succeeded with M=%d", M_new))
      return(finalize_msplit_result(tier2_expanded, X, A, Y,
                                     tier = "tier2_rashomon_0.05",
                                     M_used = M_new, M_expanded = TRUE, verbose))
    }

    M <- M_new
  }

  # All tiers failed
  return(list(
    theta = NA_real_,
    sigma = NA_real_,
    ci_95 = c(NA_real_, NA_real_),
    e_hat = rep(NA_real_, n),
    m0_hat = rep(NA_real_, n),
    error = sprintf("All tiers failed: <3 trees matched modal even with M=%d", M),
    tier_used = "failed",
    M_used = M,
    M_expanded = TRUE,
    n_splits = M_original,
    n_folds = K
  ))
}
