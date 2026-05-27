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
#' 2. Stage 2: Refit modal structure M×K times (same as approach 5)
#' 3. Stage 3: **Average leaf values** across all M×K trees → 1 averaged tree
#' 4. Predict all observations using single averaged tree (no cross-fitting)
#'
#' @export
estimate_att_msplit_averaged_new <- function(X, A, Y,
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

    # Fit propensity with adaptive CV
    cv_e <- optimaltrees::cv_regularization_adaptive(
      X = X_train, y = A_train, loss_function = "log_loss",
      K = 5, max_iterations = 10, refit = TRUE, verbose = FALSE
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
      K = 5, max_iterations = 10, refit = TRUE, verbose = FALSE
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

  # Select modal structures
  s_star_e <- select_structure_modal(structures_e)
  s_star_m0 <- select_structure_modal(structures_m0)

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
  tree_idx <- 1

  for (m in seq_len(M)) {
    seed_m <- if (!is.null(seed_base)) seed_base + m else NULL
    folds_m <- create_folds(n, K, strata = A, seed = seed_m)

    for (k in seq_len(K)) {
      test_idx <- which(folds_m == k)
      train_idx <- which(folds_m != k)

      # Refit propensity tree with modal structure
      X_train <- X[train_idx, , drop = FALSE]
      A_train <- A[train_idx]

      tree_e_mk <- optimaltrees::refit_tree_structure(
        structure = s_star_e$structure,
        X_new = X_train,
        y_new = A_train,
        loss_function = "log_loss",
        store_training_data = FALSE,
        discretization_metadata = s_star_e$discretization_metadata
      )

      # Extract as tree structure (for averaging)
      trees_e[[tree_idx]] <- optimaltrees::get_rashomon_trees(tree_e_mk)[[1]]

      # Refit outcome tree with modal structure
      control_idx <- which(A_train == 0)
      Y_train_control <- Y[train_idx][A_train == 0]
      X_train_control <- X_train[control_idx, , drop = FALSE]

      outcome_loss <- if (outcome_type == "binary") "log_loss" else "squared_error"

      tree_m0_mk <- optimaltrees::refit_tree_structure(
        structure = s_star_m0$structure,
        X_new = X_train_control,
        y_new = Y_train_control,
        loss_function = outcome_loss,
        store_training_data = FALSE,
        discretization_metadata = s_star_m0$discretization_metadata
      )

      # Extract as tree structure
      trees_m0[[tree_idx]] <- optimaltrees::get_rashomon_trees(tree_m0_mk)[[1]]

      tree_idx <- tree_idx + 1
    }

    if (verbose && m %% max(1, M %/% 10) == 0) {
      cat(sprintf("  Refitting: %d/%d splits complete\n", m, M))
    }
  }

  # ============================================================
  # Stage 3: Average Leaf Values Across All M×K Trees
  # ============================================================
  if (verbose) cat("Stage 3: Averaging leaf values across M×K trees...\n")

  e_averaged <- average_trees(trees_e)
  m0_averaged <- average_trees(trees_m0)

  # ============================================================
  # Stage 4: Predict and Compute ATT
  # ============================================================
  if (verbose) cat("Stage 4: Computing ATT with averaged tree...\n")

  # Predict for all observations using single averaged tree
  e_hat <- predict_from_tree(e_averaged, X)
  m0_hat <- predict_from_tree(m0_averaged, X)

  # Compute ATT
  pi_hat <- mean(A)
  eta <- list(e = e_hat, m0 = m0_hat, m1 = NULL)
  score <- psi_att(Y, A, theta = 0, eta = eta, pi_hat = pi_hat)
  theta_hat <- sum(score) / sum(A / pi_hat)

  # Standard error
  score_centered <- score - mean(score)
  sigma <- sqrt(mean(score_centered^2) / n)

  # CI
  ci_95 <- theta_hat + c(-1, 1) * qnorm(0.975) * sigma / sqrt(n)

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
