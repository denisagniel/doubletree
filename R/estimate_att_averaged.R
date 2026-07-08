# Averaged-Tree ATT Estimators
#
# Alternative inference methods that average leaf values across cross-fitted trees
# to produce a single interpretable tree while maintaining cross-fit validity.

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
#' @param rashomon_bound_multiplier Numeric or NULL. Rashomon tolerance epsilon_n.
#'   Default: NULL, which uses the theory value log(n)/n (= o(n^{-1/2})) via
#'   \code{optimaltrees::select_epsilon_n(nrow(X))}.
#' @param rashomon_bound_adder Numeric. Additive Rashomon bound. Default: 0
#' @param max_leaves Integer. Maximum number of leaves (optional sieve). Default: NULL
#' @param auto_tune_intersecting Logical. Auto-tune epsilon_n to find intersection?
#'   Default: FALSE. Not valid for inference (post-selection); warns when TRUE.
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
#'   \item Average leaf values across K trees using \code{optimaltrees::average_trees()}
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
  rashomon_bound_multiplier = NULL,
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

  # Resolve Rashomon tolerance: NULL -> theory epsilon_n = log(n)/n (o(n^{-1/2})).
  if (is.null(rashomon_bound_multiplier)) {
    rashomon_bound_multiplier <- optimaltrees::select_epsilon_n(n)
    if (verbose) {
      message("Using theory epsilon_n = log(n)/n = ",
              signif(rashomon_bound_multiplier, 3))
    }
  }
  # Data-adaptive epsilon_n voids the valid-inference guarantee (post-selection).
  if (isTRUE(auto_tune_intersecting)) {
    warning(
      "auto_tune_intersecting = TRUE selects the Rashomon tolerance from the ",
      "data (post-selection) and voids the o(n^{-1/2}) valid-inference ",
      "guarantee; exploratory use only. For inference, keep the fixed theory ",
      "epsilon_n (rashomon_bound_multiplier = NULL).",
      call. = FALSE
    )
  }

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

  # Sample-size-weighted averaging: each fold-refit tree carries an "n_per_leaf"
  # attribute (from refit_structure_on_data) giving the training-observation count
  # per leaf. Weighting by these counts (rather than uniform) gives leaves fitted on
  # more data more influence in the average, matching estimate_att_msplit_averaged.
  # Fall back to uniform weights only if the attribute is missing (defensive).
  leaf_weights <- function(trees) {
    lapply(trees, function(tree) {
      npl <- attr(tree, "n_per_leaf")
      if (is.null(npl)) {
        lv <- optimaltrees::extract_leaf_values(tree)
        setNames(rep(1L, length(lv)), names(lv))
      } else {
        npl
      }
    })
  }
  weights_e <- leaf_weights(e_trees)
  weights_m0 <- leaf_weights(m0_trees)

  e_averaged <- tryCatch({
    optimaltrees::average_trees(e_trees, weights_e)
  }, error = function(e) {
    stop("Failed to average propensity trees: ", e$message, call. = FALSE)
  })

  m0_averaged <- tryCatch({
    optimaltrees::average_trees(m0_trees, weights_m0)
  }, error = function(e) {
    stop("Failed to average outcome trees: ", e$message, call. = FALSE)
  })

  # Predict for ALL observations (no cross-fitting)
  if (verbose) message("\n--- Generating predictions with averaged trees ---")

  # Trees were built on discretized binary features (X_binary).
  # Must apply the same discretization to X before prediction.
  apply_disc <- optimaltrees::apply_discretization
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

  e_hat <- optimaltrees::predict_averaged_tree(e_averaged, X_for_e_pred)
  m0_hat <- optimaltrees::predict_averaged_tree(m0_averaged, X_for_m0_pred)

  # Compute ATT via EIF
  if (verbose) message("\n--- Computing ATT estimate ---")

  # Shared EIF solve (see eif_att_solve in inference.R).
  .att <- eif_att_solve(Y, A, e_hat, m0_hat, n)
  theta_hat <- .att$theta
  score_values <- .att$score_values
  sigma <- .att$sigma
  ci_95 <- .att$ci_95

  if (verbose) {
    message(sprintf("\n=== Results ==="))
    message(sprintf("ATT estimate: %.4f", theta_hat))
    message(sprintf("Standard error: %.4f", sigma))
    message(sprintf("95%% CI: [%.4f, %.4f]", ci_95[1], ci_95[2]))
  }

  # Tolerance multipliers selected by escalation (epsilon_n = c*log(n)/n per nuisance).
  rashomon_c_e  <- if (is.null(nuisance_fits$rashomon_c_e))  NA_real_ else nuisance_fits$rashomon_c_e
  rashomon_c_m0 <- if (is.null(nuisance_fits$rashomon_c_m0)) NA_real_ else nuisance_fits$rashomon_c_m0
  c_vals <- c(rashomon_c_e, rashomon_c_m0)
  epsilon_n_used <- if (all(is.na(c_vals))) rashomon_bound_multiplier else
    max(c_vals, na.rm = TRUE) * (log(n) / n)

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
    epsilon_n = epsilon_n_used,
    rashomon_c_e = rashomon_c_e,
    rashomon_c_m0 = rashomon_c_m0,
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
                                            verbose = FALSE,
                                            outcome_type = c("binary", "continuous")) {
  outcome_type <- match.arg(outcome_type)
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
  # Stage 1: Structure Selection (shared with estimate_att_msplit)
  # ============================================================
  if (verbose) cat("Stage 1: Selecting modal structures...\n")

  modal <- discover_modal_structures(X, A, Y, M = M, K = K, seed_base = seed_base,
                                     outcome_type = outcome_type, verbose = verbose)
  s_star_e <- modal$e
  s_star_m0 <- modal$m0

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

  e_averaged <- optimaltrees::average_trees(trees_e, leaf_counts_e)
  m0_averaged <- optimaltrees::average_trees(trees_m0, leaf_counts_m0)

  # ============================================================
  # Stage 4: Predict and Compute ATT
  # ============================================================
  if (verbose) cat("Stage 4: Computing ATT with averaged tree...\n")

  # Predict for ALL observations using the averaged tree.
  # Trees use discretized binary features; apply_discretization converts X to
  # the same binary feature space before calling predict_from_tree.
  apply_disc <- optimaltrees::apply_discretization
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

  e_hat  <- pmax(pmin(optimaltrees::predict_averaged_tree(e_averaged,  X_for_e_pred),  0.99), 0.01)
  m0_hat <- pmax(pmin(optimaltrees::predict_averaged_tree(m0_averaged, X_for_m0_pred), 0.99), 0.01)

  # Compute ATT via the shared EIF solve (see eif_att_solve in inference.R).
  .att <- eif_att_solve(Y, A, e_hat, m0_hat, n)
  theta_hat <- .att$theta
  score <- .att$score_values
  sigma <- .att$sigma
  ci_95 <- .att$ci_95

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
