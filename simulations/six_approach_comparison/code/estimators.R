# Estimator Functions for Six-Approach Comparison
# Created: 2026-05-01
#
# Six approaches:
# (i)   Full-sample tree
# (ii)  Standard cross-fit (separate trees)
# (iii) Doubletree (cross-fit intersection)
# (iv)  Doubletree structure + single fit
# (v)   M-split doubletree
# (vi)  M-split structure + single fit

library(optimaltrees)

#' Compute ATT from EIF
#'
#' @param Y Outcome vector
#' @param A Treatment vector (0/1)
#' @param e_hat Propensity score estimates
#' @param m0_hat Outcome predictions for controls
#' @return ATT estimate
compute_att <- function(Y, A, e_hat, m0_hat) {
  n <- length(Y)

  # VALIDATE INPUTS (no silent fixes)
  if (length(e_hat) != n || length(m0_hat) != n) {
    stop("Prediction vectors must have length n=", n,
         ". Got e_hat: ", length(e_hat), ", m0_hat: ", length(m0_hat))
  }

  if (any(is.na(e_hat)) || any(is.na(m0_hat))) {
    stop("Predictions contain NA values. e_hat NAs: ", sum(is.na(e_hat)),
         ", m0_hat NAs: ", sum(is.na(m0_hat)))
  }

  if (any(!is.finite(e_hat)) || any(!is.finite(m0_hat))) {
    stop("Predictions contain non-finite values (Inf/NaN).")
  }

  # Check predictions are in valid range (but allow 0 and 1 - will be clipped)
  if (any(e_hat < 0) || any(e_hat > 1)) {
    stop("Propensity scores outside [0,1]. Got range: [",
         min(e_hat), ", ", max(e_hat), "]")
  }

  if (any(m0_hat < 0) || any(m0_hat > 1)) {
    stop("Outcome probabilities outside [0,1]. Got range: [",
         min(m0_hat), ", ", max(m0_hat), "]")
  }

  # Clip away from boundaries to avoid division by zero in EIF
  # Trees can legitimately predict 0 or 1 (certain predictions)
  n_clipped_e <- sum(e_hat < 0.01 | e_hat > 0.99)
  n_clipped_m0 <- sum(m0_hat < 0.01 | m0_hat > 0.99)

  if (n_clipped_e > 0.1 * n) {
    warning("Clipping ", n_clipped_e, " (", round(100*n_clipped_e/n, 1),
            "%) propensity scores to [0.01, 0.99]")
  }
  if (n_clipped_m0 > 0.1 * n) {
    warning("Clipping ", n_clipped_m0, " (", round(100*n_clipped_m0/n, 1),
            "%) outcome probabilities to [0.01, 0.99]")
  }

  e_hat <- pmax(pmin(e_hat, 0.99), 0.01)
  m0_hat <- pmax(pmin(m0_hat, 0.99), 0.01)

  pi_hat <- mean(A)

  # EIF for ATT
  psi <- (A / pi_hat) * (Y - m0_hat) -
         ((1 - A) * e_hat) / (pi_hat * (1 - e_hat)) * (Y - m0_hat)

  # Final validation of EIF scores
  if (any(!is.finite(psi))) {
    stop("EIF scores contain non-finite values. This should not happen after validation.")
  }

  theta_hat <- mean(psi)
  return(theta_hat)
}

#' Compute standard error from EIF
#'
#' @param Y Outcome vector
#' @param A Treatment vector
#' @param e_hat Propensity estimates
#' @param m0_hat Outcome predictions
#' @param theta_hat ATT estimate
#' @return Standard error
compute_se <- function(Y, A, e_hat, m0_hat, theta_hat) {
  n <- length(Y)

  # VALIDATE INPUTS
  if (length(e_hat) != n || length(m0_hat) != n) {
    stop("Prediction vectors must have length n=", n,
         ". Got e_hat: ", length(e_hat), ", m0_hat: ", length(m0_hat))
  }

  if (any(is.na(e_hat)) || any(is.na(m0_hat))) {
    stop("Predictions contain NA values. e_hat NAs: ", sum(is.na(e_hat)),
         ", m0_hat NAs: ", sum(is.na(m0_hat)))
  }

  if (any(!is.finite(e_hat)) || any(!is.finite(m0_hat))) {
    stop("Predictions contain non-finite values (Inf/NaN).")
  }

  if (any(e_hat < 0) || any(e_hat > 1)) {
    stop("Propensity scores must be in [0,1]. Got range: [",
         min(e_hat), ", ", max(e_hat), "]")
  }

  if (any(m0_hat < 0) || any(m0_hat > 1)) {
    stop("Outcome probabilities must be in [0,1]. Got range: [",
         min(m0_hat), ", ", max(m0_hat), "]")
  }

  # Apply bounds for numerical stability
  n_clipped_e <- sum(e_hat < 0.01 | e_hat > 0.99)
  n_clipped_m0 <- sum(m0_hat < 0.01 | m0_hat > 0.99)

  if (n_clipped_e > 0.1 * n) {
    warning("Clipping ", n_clipped_e, " (", round(100*n_clipped_e/n, 1),
            "%) propensity scores to [0.01, 0.99] in SE calculation")
  }
  if (n_clipped_m0 > 0.1 * n) {
    warning("Clipping ", n_clipped_m0, " (", round(100*n_clipped_m0/n, 1),
            "%) outcome probabilities to [0.01, 0.99] in SE calculation")
  }

  e_hat <- pmax(pmin(e_hat, 0.99), 0.01)
  m0_hat <- pmax(pmin(m0_hat, 0.99), 0.01)

  pi_hat <- mean(A)

  # EIF scores
  psi <- (A / pi_hat) * (Y - m0_hat - theta_hat) -
         ((1 - A) * e_hat) / (pi_hat * (1 - e_hat)) * (Y - m0_hat)

  # Validate EIF scores
  if (any(!is.finite(psi))) {
    stop("EIF scores contain non-finite values in SE calculation.")
  }

  # Variance
  se <- sqrt(mean(psi^2) / n)

  if (!is.finite(se) || se <= 0) {
    stop("Computed SE is invalid: ", se)
  }

  return(se)
}

#' Create K-fold indices
#'
#' @param n Sample size
#' @param K Number of folds
#' @return Vector of fold assignments (1 to K)
create_folds <- function(n, K = 5) {
  folds <- rep(1:K, length.out = n)
  sample(folds)  # Shuffle
}

# ============================================================================
# Approach (i): Full-Sample Tree
# ============================================================================

#' Estimate ATT using full-sample trees
#'
#' Fits trees on all n observations (has overfitting bias)
#'
#' @param X Covariate data.frame
#' @param A Treatment vector
#' @param Y Outcome vector
#' @param regularization Not used (CV selects lambda). Kept for backward compatibility.
#' @return List with theta, se, trees
#' @export
estimate_att_fullsample <- function(X, A, Y, regularization = 0.1) {
  n <- nrow(X)

  # Fit propensity tree on all data with CV-selected lambda
  cv_e <- optimaltrees::cv_regularization(
    X = X,
    y = A,
    loss_function = "log_loss",
    K = 5,
    refit = TRUE,
    verbose = FALSE
  )

  # No fallback - CV must succeed
  if (is.na(cv_e$best_lambda)) {
    stop(
      "CV failed for propensity model in fullsample approach.\n",
      "Possible fixes:\n",
      "  1. Check data quality (enough variation in A?)\n",
      "  2. Try different lambda_grid in cv_regularization()\n",
      "  3. Increase K in cv_regularization() for more stable CV\n",
      "  4. Check for numerical issues (NaN, Inf in X or A)",
      call. = FALSE
    )
  }

  e_tree <- cv_e$model

  e_pred <- predict(e_tree, X, type = "prob")
  if (!is.matrix(e_pred) || ncol(e_pred) != 2) {
    stop("Propensity tree predict() returned unexpected format. Expected 2-column matrix, got: ",
         class(e_pred), " with dims: ", paste(dim(e_pred), collapse="x"))
  }
  e_hat <- e_pred[, 2]  # P(A=1|X)

  # Fit outcome tree on controls (all data) with CV-selected lambda
  control_idx <- A == 0
  cv_m0 <- optimaltrees::cv_regularization(
    X = X[control_idx, , drop = FALSE],
    y = Y[control_idx],
    loss_function = "log_loss",
    K = 5,
    refit = TRUE,
    verbose = FALSE
  )

  # No fallback - CV must succeed
  if (is.na(cv_m0$best_lambda)) {
    stop(
      "CV failed for outcome model in fullsample approach.\n",
      "Possible fixes:\n",
      "  1. Check data quality (enough control units? variation in Y?)\n",
      "  2. Try different lambda_grid in cv_regularization()\n",
      "  3. Increase K in cv_regularization() for more stable CV\n",
      "  4. Check for numerical issues (NaN, Inf in X or Y)",
      call. = FALSE
    )
  }

  m0_tree <- cv_m0$model

  m0_pred <- predict(m0_tree, X, type = "prob")
  if (!is.matrix(m0_pred) || ncol(m0_pred) != 2) {
    stop("Outcome tree predict() returned unexpected format. Expected 2-column matrix, got: ",
         class(m0_pred), " with dims: ", paste(dim(m0_pred), collapse="x"))
  }
  m0_hat <- m0_pred[, 2]  # P(Y=1|A=0,X)

  # Compute ATT
  theta_hat <- compute_att(Y, A, e_hat, m0_hat)

  # SE (use cross-fitted for conservatism)
  # For this approach, we'll compute SE using these same predictions
  # (This is not quite right, but we'll test with bias correction)
  se <- compute_se(Y, A, e_hat, m0_hat, theta_hat)

  list(
    theta = theta_hat,
    se = se,
    e_hat = e_hat,
    m0_hat = m0_hat,
    trees = list(e = e_tree, m0 = m0_tree)
  )
}

# ============================================================================
# Approach (ii): Standard Cross-Fit (Separate Trees)
# ============================================================================

#' Estimate ATT using standard cross-fitting
#'
#' Fits separate tree per fold (no structure intersection)
#'
#' @param X Covariate data.frame
#' @param A Treatment vector
#' @param Y Outcome vector
#' @param K Number of folds
#' @param regularization Not used (CV selects lambda). Kept for backward compatibility.
#' @return List with theta, se
#' @export
estimate_att_crossfit <- function(X, A, Y, K = 5, regularization = 0.1) {
  n <- nrow(X)
  folds <- create_folds(n, K)

  e_hat <- numeric(n)
  m0_hat <- numeric(n)

  for (k in 1:K) {
    train_idx <- folds != k
    test_idx <- folds == k

    # Fit propensity tree on training fold with CV-selected lambda
    cv_e_k <- optimaltrees::cv_regularization(
      X = X[train_idx, , drop = FALSE],
      y = A[train_idx],
      loss_function = "log_loss",
      K = 5,
      refit = TRUE,
      verbose = FALSE
    )

    # No fallback - CV must succeed
    if (is.na(cv_e_k$best_lambda)) {
      stop(
        "CV failed for propensity model in crossfit approach (fold ", k, ").\n",
        "Possible fixes:\n",
        "  1. Check data quality (enough variation in A in this fold?)\n",
        "  2. Try different lambda_grid in cv_regularization()\n",
        "  3. Increase K in cv_regularization() for more stable CV\n",
        "  4. Check for numerical issues (NaN, Inf in X or A)\n",
        "  5. Try different random seed (fold might have unusual split)",
        call. = FALSE
      )
    }

    e_tree_k <- cv_e_k$model

    e_pred_k <- predict(e_tree_k, X[test_idx, , drop = FALSE], type = "prob")
    if (!is.matrix(e_pred_k) || ncol(e_pred_k) != 2) {
      stop("Propensity tree predict() (fold ", k, ") returned unexpected format. Expected 2-column matrix, got: ",
           class(e_pred_k), " with dims: ", paste(dim(e_pred_k), collapse="x"))
    }
    e_hat[test_idx] <- e_pred_k[, 2]  # P(A=1|X)

    # Fit outcome tree on controls in training fold with CV-selected lambda
    control_train_idx <- train_idx & (A == 0)
    cv_m0_k <- optimaltrees::cv_regularization(
      X = X[control_train_idx, , drop = FALSE],
      y = Y[control_train_idx],
      loss_function = "log_loss",
      K = 5,
      refit = TRUE,
      verbose = FALSE
    )

    # No fallback - CV must succeed
    if (is.na(cv_m0_k$best_lambda)) {
      stop(
        "CV failed for outcome model in crossfit approach (fold ", k, ").\n",
        "Possible fixes:\n",
        "  1. Check data quality (enough control units in this fold? variation in Y?)\n",
        "  2. Try different lambda_grid in cv_regularization()\n",
        "  3. Increase K in cv_regularization() for more stable CV\n",
        "  4. Check for numerical issues (NaN, Inf in X or Y)\n",
        "  5. Try different random seed (fold might have unusual split)",
        call. = FALSE
      )
    }

    m0_tree_k <- cv_m0_k$model

    m0_pred_k <- predict(m0_tree_k, X[test_idx, , drop = FALSE], type = "prob")
    if (!is.matrix(m0_pred_k) || ncol(m0_pred_k) != 2) {
      stop("Outcome tree predict() (fold ", k, ") returned unexpected format. Expected 2-column matrix, got: ",
           class(m0_pred_k), " with dims: ", paste(dim(m0_pred_k), collapse="x"))
    }
    m0_hat[test_idx] <- m0_pred_k[, 2]  # P(Y=1|A=0,X)
  }

  # Compute ATT and SE
  theta_hat <- compute_att(Y, A, e_hat, m0_hat)
  se <- compute_se(Y, A, e_hat, m0_hat, theta_hat)

  list(
    theta = theta_hat,
    se = se,
    e_hat = e_hat,
    m0_hat = m0_hat
  )
}

# ============================================================================
# Approach (iii): Doubletree (Cross-Fit Intersection)
# ============================================================================

#' Estimate ATT using doubletree (Rashomon intersection)
#'
#' Uses doubletree package with use_rashomon = TRUE
#'
#' @param X Covariate data.frame
#' @param A Treatment vector
#' @param Y Outcome vector
#' @param K Number of folds
#' @param regularization Not used (CV selects lambda). Kept for backward compatibility.
#' @return List with theta, se, structures
#' @export
estimate_att_doubletree <- function(X, A, Y, K = 5, regularization = 0.1) {
  # Theory-justified Rashomon bound: 2*sqrt(log(n)/n) per manuscript Appendix A.5
  # At n=500: ~0.22; n=1000: ~0.17; n=2000: ~0.12 (decays as required for valid inference)
  eps_n <- 2 * sqrt(log(nrow(X)) / nrow(X))

  # Use doubletree package implementation with Rashomon
  result <- doubletree::estimate_att(
    X = X,
    A = A,
    Y = Y,
    K = K,
    # regularization parameter removed - CV is used by default
    outcome_type = "binary",
    use_rashomon = TRUE,
    rashomon_bound_multiplier = eps_n,
    auto_tune_intersecting = FALSE,  # pure Rashomon: fall back if no intersection at eps_n
    verbose = FALSE
  )

  # Validate result structure
  # NA is acceptable (Rashomon intersection may fail), but Inf/NaN is not
  if (!is.na(result$theta) && !is.finite(result$theta)) {
    stop("doubletree::estimate_att returned non-finite theta: ", result$theta)
  }

  if (!is.na(result$sigma) && (!is.finite(result$sigma) || result$sigma <= 0)) {
    stop("doubletree::estimate_att returned invalid sigma: ", result$sigma)
  }

  # If theta is NA, this indicates Rashomon intersection failed
  # Return with error message for tracking
  if (is.na(result$theta)) {
    return(list(
      theta = NA_real_,
      se = NA_real_,
      e_hat = rep(NA_real_, nrow(X)),
      m0_hat = rep(NA_real_, nrow(X)),
      structures = list(e = NULL, m0 = NULL),
      error = "Rashomon intersection found no common structure"
    ))
  }

  list(
    theta = result$theta,
    se = result$sigma,
    e_hat = if (!is.null(result$nuisance_fits$cf_e)) {
      predict(result$nuisance_fits$cf_e, X, type = "prob")
    } else {
      rep(NA, nrow(X))
    },
    m0_hat = if (!is.null(result$nuisance_fits$cf_m0)) {
      predict(result$nuisance_fits$cf_m0, X, type = "prob")
    } else {
      rep(NA, nrow(X))
    },
    structures = list(
      e = if (!is.null(result$nuisance_fits$cf_e)) result$nuisance_fits$cf_e@n_intersecting else NULL,
      m0 = if (!is.null(result$nuisance_fits$cf_m0)) result$nuisance_fits$cf_m0@n_intersecting else NULL
    ),
    # Also return the full CF objects (for approach iv to access fold_refits)
    cf_e = result$nuisance_fits$cf_e,
    cf_m0 = result$nuisance_fits$cf_m0
  )
}

# ============================================================================
# Approach (iv): Doubletree Structure + Averaged Leaves
# ============================================================================

#' Estimate ATT using doubletree with averaged leaf values
#'
#' Wrapper for doubletree::estimate_att_doubletree_averaged()
#'
#' @param X Covariate data.frame
#' @param A Treatment vector
#' @param Y Outcome vector
#' @param K Number of folds (for structure selection)
#' @param regularization Not used (CV selects lambda). Kept for backward compatibility.
#' @return List with theta, se, structures
#' @export
estimate_att_doubletree_averaged <- function(X, A, Y, K = 5, regularization = 0.1) {
  # Theory-justified Rashomon bound: 2*sqrt(log(n)/n) per manuscript Appendix A.5
  # auto_tune_intersecting=TRUE starts from eps_n and increases only if no intersection found
  eps_n <- 2 * sqrt(log(nrow(X)) / nrow(X))

  result <- doubletree::estimate_att_doubletree_averaged(
    X = X,
    A = A,
    Y = Y,
    K = K,
    regularization = regularization,
    outcome_type = "binary",
    rashomon_bound_multiplier = eps_n,
    auto_tune_intersecting = TRUE,
    verbose = FALSE
  )

  # Convert to simulation format (theta → theta, sigma → se)
  list(
    theta = result$theta,
    se = result$sigma,
    e_hat = result$e_hat,
    m0_hat = result$m0_hat,
    structures = result$structures,
    averaged_trees = if(!is.null(result$averaged_trees)) result$averaged_trees else NULL,
    error = if (!is.null(result$error)) result$error else NULL
  )
}

# ============================================================================
# Approach (v): M-Split Doubletree
# ============================================================================

#' Estimate ATT using M-split doubletree
#'
#' Finds modal structure across M splits, averages predictions.
#'
#' @param X Covariate data.frame
#' @param A Treatment vector
#' @param Y Outcome vector
#' @param M Number of independent splits
#' @param K Number of folds per split
#' @param regularization Not used (CV selects lambda). Kept for backward compatibility.
#' @return List with theta, se, structures, diagnostics
#' @export
estimate_att_msplit <- function(X, A, Y, M = 10, K = 5, regularization = 0.1) {
  result <- doubletree::estimate_att_msplit(
    X = X,
    A = A,
    Y = Y,
    M = M,
    K = K,
    regularization = regularization,
    outcome_type = "binary",
    verbose = FALSE
  )

  list(
    theta = result$theta,
    se = result$sigma,
    e_hat = result$averaged_predictions$e,
    m0_hat = result$averaged_predictions$m0,
    structures = result$structures,
    diagnostics = result$diagnostics
  )
}

# ============================================================================
# Approach (vi): M-Split Structure + M×K Averaged Leaves
# ============================================================================

#' Estimate ATT using M-split with M×K leaf averaging
#'
#' Gets modal structure from M splits, then averages leaf values across
#' all M×K cross-fitted trees (K trees per split, M splits).
#'
#' @param X Covariate data.frame
#' @param A Treatment vector
#' @param Y Outcome vector
#' @param M Number of splits (for structure selection)
#' @param K Number of folds per split
#' @param regularization Not used (CV selects lambda). Kept for backward compatibility.
#' @return List with theta, se, structures, diagnostics
#' @export
estimate_att_msplit_averaged <- function(X, A, Y, M = 10, K = 5, regularization = 0.1) {
  result <- doubletree::estimate_att_msplit_averaged(
    X = X,
    A = A,
    Y = Y,
    M = M,
    K = K,
    outcome_type = "binary",
    seed_base = NULL
  )

  # Convert to simulation format (theta → theta, sigma → se)
  list(
    theta = result$theta,
    se = result$sigma,
    e_hat = result$e_hat,
    m0_hat = result$m0_hat,
    structures = result$structures,
    n_trees_averaged = if(!is.null(result$n_trees_averaged)) result$n_trees_averaged else NULL,
    averaged_trees = if(!is.null(result$averaged_trees)) result$averaged_trees else NULL,
    error = if (!is.null(result$error)) result$error else NULL
  )
}
