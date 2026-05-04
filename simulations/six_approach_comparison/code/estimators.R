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
  psi <- (A / pi_hat) * (Y - m0_hat) +
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
  psi <- (A / pi_hat) * (Y - m0_hat - theta_hat) +
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
#' @param regularization Tree complexity penalty
#' @return List with theta, se, trees
#' @export
estimate_att_fullsample <- function(X, A, Y, regularization = 0.1) {
  n <- nrow(X)

  # Fit propensity tree on all data
  e_tree <- optimaltrees::fit_tree(
    X = X,
    y = A,
    loss_function = "log_loss",
    regularization = regularization,
    verbose = FALSE
  )
  e_pred <- predict(e_tree, X, type = "prob")
  if (!is.matrix(e_pred) || ncol(e_pred) != 2) {
    stop("Propensity tree predict() returned unexpected format. Expected 2-column matrix, got: ",
         class(e_pred), " with dims: ", paste(dim(e_pred), collapse="x"))
  }
  e_hat <- e_pred[, 2]  # P(A=1|X)

  # Fit outcome tree on controls (all data)
  control_idx <- A == 0
  m0_tree <- optimaltrees::fit_tree(
    X = X[control_idx, , drop = FALSE],
    y = Y[control_idx],
    loss_function = "log_loss",
    regularization = regularization,
    verbose = FALSE
  )
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
#' @param regularization Tree complexity penalty
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

    # Fit propensity tree on training fold
    e_tree_k <- optimaltrees::fit_tree(
      X = X[train_idx, , drop = FALSE],
      y = A[train_idx],
      loss_function = "log_loss",
      regularization = regularization,
      verbose = FALSE
    )
    e_pred_k <- predict(e_tree_k, X[test_idx, , drop = FALSE], type = "prob")
    if (!is.matrix(e_pred_k) || ncol(e_pred_k) != 2) {
      stop("Propensity tree predict() (fold ", k, ") returned unexpected format. Expected 2-column matrix, got: ",
           class(e_pred_k), " with dims: ", paste(dim(e_pred_k), collapse="x"))
    }
    e_hat[test_idx] <- e_pred_k[, 2]  # P(A=1|X)

    # Fit outcome tree on controls in training fold
    control_train_idx <- train_idx & (A == 0)
    m0_tree_k <- optimaltrees::fit_tree(
      X = X[control_train_idx, , drop = FALSE],
      y = Y[control_train_idx],
      loss_function = "log_loss",
      regularization = regularization,
      verbose = FALSE
    )
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
#' @param regularization Tree complexity penalty
#' @return List with theta, se, structures
#' @export
estimate_att_doubletree <- function(X, A, Y, K = 5, regularization = 0.1) {
  # Use doubletree package implementation with Rashomon
  result <- doubletree::estimate_att(
    X = X,
    A = A,
    Y = Y,
    K = K,
    regularization = regularization,
    outcome_type = "binary",
    use_rashomon = TRUE,
    rashomon_bound_multiplier = 0.05,
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
      predict(result$nuisance_fits$cf_e, X)
    } else {
      rep(NA, nrow(X))
    },
    m0_hat = if (!is.null(result$nuisance_fits$cf_m0)) {
      predict(result$nuisance_fits$cf_m0, X)
    } else {
      rep(NA, nrow(X))
    },
    structures = list(
      e = if (!is.null(result$nuisance_fits$cf_e)) result$nuisance_fits$cf_e$structure else NULL,
      m0 = if (!is.null(result$nuisance_fits$cf_m0)) result$nuisance_fits$cf_m0$structure else NULL
    )
  )
}

# ============================================================================
# Approach (iv): Doubletree Structure + Single Fit
# ============================================================================

#' Estimate ATT using doubletree structure fitted once
#'
#' Gets structure from cross-fit intersection, but refits on all data
#'
#' @param X Covariate data.frame
#' @param A Treatment vector
#' @param Y Outcome vector
#' @param K Number of folds (for structure selection)
#' @param regularization Tree complexity penalty
#' @return List with theta, se, structures
#' @export
estimate_att_doubletree_singlefit <- function(X, A, Y, K = 5, regularization = 0.1) {
  # Stage 1: Get structures from doubletree
  # Catch errors from doubletree call
  result_doubletree <- tryCatch({
    estimate_att_doubletree(X, A, Y, K, regularization)
  }, error = function(e) {
    return(list(
      theta = NA_real_,
      se = NA_real_,
      e_hat = rep(NA_real_, nrow(X)),
      m0_hat = rep(NA_real_, nrow(X)),
      structures = list(e = NULL, m0 = NULL),
      error = paste("doubletree call failed:", e$message)
    ))
  })

  # Check if doubletree call failed
  if (!is.null(result_doubletree$error)) {
    return(result_doubletree)
  }

  e_structure <- result_doubletree$structures$e
  m0_structure <- result_doubletree$structures$m0

  # Check if structures were successfully found
  if (is.null(e_structure) || is.null(m0_structure)) {
    return(list(
      theta = NA_real_,
      se = NA_real_,
      e_hat = rep(NA_real_, nrow(X)),
      m0_hat = rep(NA_real_, nrow(X)),
      structures = list(e = NULL, m0 = NULL),
      error = "Rashomon intersection failed - no common structure found"
    ))
  }

  # Stage 2: Refit structures on ALL data
  e_refit <- optimaltrees::refit_tree_structure(
    structure = e_structure,
    X = X,
    y = A,
    loss_function = "log_loss"
  )
  e_pred <- predict(e_refit, X, type = "prob")
  if (!is.matrix(e_pred) || ncol(e_pred) != 2) {
    stop("Refitted propensity tree predict() returned unexpected format. Expected 2-column matrix, got: ",
         class(e_pred), " with dims: ", paste(dim(e_pred), collapse="x"))
  }
  e_hat <- e_pred[, 2]  # P(A=1|X)

  control_idx <- A == 0
  m0_refit <- optimaltrees::refit_tree_structure(
    structure = m0_structure,
    X = X[control_idx, , drop = FALSE],
    y = Y[control_idx],
    loss_function = "log_loss"
  )
  m0_pred <- predict(m0_refit, X, type = "prob")
  if (!is.matrix(m0_pred) || ncol(m0_pred) != 2) {
    stop("Refitted outcome tree predict() returned unexpected format. Expected 2-column matrix, got: ",
         class(m0_pred), " with dims: ", paste(dim(m0_pred), collapse="x"))
  }
  m0_hat <- m0_pred[, 2]  # P(Y=1|A=0,X)

  # Compute ATT and SE
  theta_hat <- compute_att(Y, A, e_hat, m0_hat)
  se <- compute_se(Y, A, e_hat, m0_hat, theta_hat)

  list(
    theta = theta_hat,
    se = se,
    e_hat = e_hat,
    m0_hat = m0_hat,
    structures = list(e = e_structure, m0 = m0_structure)
  )
}

# ============================================================================
# Approach (v): M-Split Doubletree
# ============================================================================

#' Estimate ATT using M-split doubletree
#'
#' Finds modal structure across M splits, averages predictions
#'
#' @param X Covariate data.frame
#' @param A Treatment vector
#' @param Y Outcome vector
#' @param M Number of independent splits
#' @param K Number of folds per split
#' @param regularization Tree complexity penalty
#' @return List with theta, se, structures, diagnostics
#' @export
estimate_att_msplit <- function(X, A, Y, M = 10, K = 5, regularization = 0.1) {
  # Use the doubletree package implementation
  # This was already implemented in earlier session
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
# Approach (vi): M-Split Structure + Single Fit
# ============================================================================

#' Estimate ATT using M-split structure fitted once
#'
#' Gets modal structure from M-split, but refits on all data
#'
#' @param X Covariate data.frame
#' @param A Treatment vector
#' @param Y Outcome vector
#' @param M Number of splits (for structure selection)
#' @param K Number of folds per split
#' @param regularization Tree complexity penalty
#' @return List with theta, se, structures, diagnostics
#' @export
estimate_att_msplit_singlefit <- function(X, A, Y, M = 10, K = 5, regularization = 0.1) {
  # Stage 1: Get modal structure from M-split
  result_msplit <- estimate_att_msplit(X, A, Y, M, K, regularization)
  e_structure <- result_msplit$structures$e
  m0_structure <- result_msplit$structures$m0

  # Stage 2: Refit structures on ALL data
  e_refit <- optimaltrees::refit_tree_structure(
    structure = e_structure,
    X = X,
    y = A,
    loss_function = "log_loss"
  )
  e_pred <- predict(e_refit, X, type = "prob")
  if (!is.matrix(e_pred) || ncol(e_pred) != 2) {
    stop("Refitted propensity tree (M-split) predict() returned unexpected format. Expected 2-column matrix, got: ",
         class(e_pred), " with dims: ", paste(dim(e_pred), collapse="x"))
  }
  e_hat <- e_pred[, 2]  # P(A=1|X)

  control_idx <- A == 0
  m0_refit <- optimaltrees::refit_tree_structure(
    structure = m0_structure,
    X = X[control_idx, , drop = FALSE],
    y = Y[control_idx],
    loss_function = "log_loss"
  )
  m0_pred <- predict(m0_refit, X, type = "prob")
  if (!is.matrix(m0_pred) || ncol(m0_pred) != 2) {
    stop("Refitted outcome tree (M-split) predict() returned unexpected format. Expected 2-column matrix, got: ",
         class(m0_pred), " with dims: ", paste(dim(m0_pred), collapse="x"))
  }
  m0_hat <- m0_pred[, 2]  # P(Y=1|A=0,X)

  # Compute ATT and SE
  theta_hat <- compute_att(Y, A, e_hat, m0_hat)
  se <- compute_se(Y, A, e_hat, m0_hat, theta_hat)

  list(
    theta = theta_hat,
    se = se,
    e_hat = e_hat,
    m0_hat = m0_hat,
    structures = list(e = e_structure, m0 = m0_structure),
    diagnostics = result_msplit$diagnostics
  )
}
