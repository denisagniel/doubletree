#' ATT Estimation with Random Forest Nuisances (ranger)
#'
#' Cross-fitted ATT estimator using random forests for both
#' propensity e(X) and outcome m0(X) estimation. Matches the structure
#' of doubletree::estimate_att but uses ranger instead of trees.
#'
#' Paper reference: forest baseline comparison
#'
#' @param X Data.frame of covariates
#' @param A Binary treatment (0/1)
#' @param Y Binary outcome (0/1)
#' @param K Number of folds for cross-fitting (default 5)
#' @param seed Random seed for fold creation
#' @param num.trees Number of trees in forest (default 500)
#' @param mtry Features per split (default sqrt(p) for classification)
#' @param min.node.size Minimum node size (default 5)
#' @param probability Use probability forests (required, always TRUE)
#' @param verbose Print progress messages
#'
#' @return List with components:
#'   - theta: Point estimate of ATT
#'   - sigma: Standard error
#'   - ci: 95% confidence interval (length 2 vector)
#'   - n_treated: Number of treated units
#'   - convergence: Always "converged" (ranger doesn't fail)
#'
#' @details
#' Cross-fitting procedure:
#' 1. Split data into K folds (stratified by treatment)
#' 2. For each fold k:
#'    - Train e(X) on folds != k using ranger (binary classification)
#'    - Predict e(X) on fold k
#'    - Train m0(X) on control units in folds != k (binary classification)
#'    - Predict m0(X) on fold k
#' 3. Compute orthogonal score: ψ = (Y - m0(X)) * A / e(X)
#' 4. θ̂ = E[ψ] / E[A/e(X)]
#' 5. Inference via EIF-based variance estimator
#'
#' @examples
#' source("production/dgps/dgps_smooth.R")
#' d <- generate_dgp_binary_att(400, tau = 0.10, seed = 123)
#' fit <- att_forest(d$X, d$A, d$Y, K = 5, seed = 123)
#' print(fit$theta)
#' print(fit$ci)

att_forest <- function(X, A, Y, K = 5, seed = NULL,
                           num.trees = 500, mtry = NULL, min.node.size = 5,
                           probability = TRUE, verbose = FALSE) {

  # Check required package
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' required for forest-based ATT estimation. Install with: install.packages('ranger')")
  }

  # Input validation
  n <- nrow(X)
  if (length(A) != n || length(Y) != n) {
    stop("X, A, Y must have same number of rows")
  }
  if (!all(A %in% c(0, 1))) {
    stop("A must be binary (0/1)")
  }
  if (!all(Y %in% c(0, 1))) {
    stop("Y must be binary (0/1)")
  }

  n_treated <- sum(A == 1)
  if (n_treated == 0) {
    stop("No treated units")
  }
  if (n_treated == n) {
    stop("No control units")
  }

  # Set mtry to default (sqrt(p) for classification)
  if (is.null(mtry)) {
    mtry <- floor(sqrt(ncol(X)))
  }

  # Create folds (stratified by treatment)
  # Use same logic as doubletree
  if (!is.null(seed)) set.seed(seed)

  # Stratified sampling
  treated_idx <- which(A == 1)
  control_idx <- which(A == 0)

  folds <- integer(n)
  folds[treated_idx] <- sample(rep(1:K, length.out = length(treated_idx)))
  folds[control_idx] <- sample(rep(1:K, length.out = length(control_idx)))

  # Initialize storage for cross-fitted predictions
  e_hat <- numeric(n)    # Propensity scores
  m0_hat <- numeric(n)   # Control outcome predictions

  # Cross-fitting loop
  for (k in 1:K) {
    if (verbose) {
      cat(sprintf("Fold %d/%d...\n", k, K))
    }

    test_idx <- which(folds == k)
    train_idx <- which(folds != k)

    # --- Propensity Score e(X) ---
    # Train on all units in training folds
    train_data_e <- cbind(X[train_idx, , drop = FALSE], A = A[train_idx])

    rf_e <- ranger::ranger(
      A ~ .,
      data = train_data_e,
      num.trees = num.trees,
      mtry = mtry,
      min.node.size = min.node.size,
      probability = TRUE,  # Required for classification
      num.threads = 1,     # Avoid nested parallelism
      verbose = FALSE
    )

    # Predict on test fold (get P(A=1|X))
    pred_e <- predict(rf_e, X[test_idx, , drop = FALSE])
    e_hat[test_idx] <- pred_e$predictions[, 2]  # Column 2 is P(A=1)

    # --- Control Outcome m0(X) ---
    # Train ONLY on control units in training folds
    control_train_idx <- train_idx[A[train_idx] == 0]

    if (length(control_train_idx) < 10) {
      warning(sprintf("Fold %d: Only %d control units in training set",
                      k, length(control_train_idx)))
    }

    train_data_m0 <- cbind(X[control_train_idx, , drop = FALSE],
                           Y = Y[control_train_idx])

    rf_m0 <- ranger::ranger(
      Y ~ .,
      data = train_data_m0,
      num.trees = num.trees,
      mtry = mtry,
      min.node.size = min.node.size,
      probability = TRUE,
      num.threads = 1,
      verbose = FALSE
    )

    # Predict on test fold (get P(Y=1|X, A=0))
    pred_m0 <- predict(rf_m0, X[test_idx, , drop = FALSE])
    m0_hat[test_idx] <- pred_m0$predictions[, 2]  # Column 2 is P(Y=1)
  }

  # Clip propensity scores to avoid extreme weights
  # Use same bounds as doubletree (1e-6, 1-1e-6)
  e_hat <- pmax(pmin(e_hat, 1 - 1e-6), 1e-6)
  m0_hat <- pmax(pmin(m0_hat, 1 - 1e-6), 1e-6)

  # --- Orthogonal Score and Point Estimate ---
  # ATT orthogonal score (Chernozhukov et al. 2018):
  # ψ_i(θ) = (A_i/π)*(Y_i - m0_i - θ) - (1/π)*(e_i*(1 - A_i)/(1 - e_i))*(Y_i - m0_i)
  # where π = P(A=1) = mean(A)

  pi_hat <- mean(A)

  # Solve for θ: sum(ψ_i(θ)) = 0
  # This gives: θ = sum(ψ_i(0)) / sum(A_i/π)
  score_at_zero_term1 <- (A / pi_hat) * (Y - m0_hat)
  score_at_zero_term2 <- (1 / pi_hat) * (e_hat * (1 - A) / (1 - e_hat)) * (Y - m0_hat)
  score_at_zero_term2[!is.finite(score_at_zero_term2)] <- 0  # Safety

  score_at_zero <- score_at_zero_term1 - score_at_zero_term2
  sum_a_over_pi <- sum(A / pi_hat)

  theta <- sum(score_at_zero) / sum_a_over_pi

  # --- Variance Estimation ---
  # Compute score at θ̂ (proper Neyman-orthogonal score)
  score_term1 <- (A / pi_hat) * (Y - m0_hat - theta)
  score_term2 <- (1 / pi_hat) * (e_hat * (1 - A) / (1 - e_hat)) * (Y - m0_hat)
  score_term2[!is.finite(score_term2)] <- 0

  score_values <- score_term1 - score_term2

  # Variance: Var(√n * θ̂) = E[ψ(θ)^2]
  # Standard error: SE(θ̂) = sqrt(Var) / sqrt(n)
  variance <- mean(score_values^2)
  sigma <- sqrt(variance / n)

  # 95% CI (normal approximation)
  ci <- theta + c(-1.96, 1.96) * sigma

  # Return results (same structure as doubletree::estimate_att)
  result <- list(
    theta = theta,
    sigma = sigma,
    ci = ci,
    n_treated = n_treated,
    convergence = "converged",
    method = "forest",
    hyperparams = list(
      num.trees = num.trees,
      mtry = mtry,
      min.node.size = min.node.size
    )
  )

  if (verbose) {
    cat(sprintf("\nForest ATT Results:\n"))
    cat(sprintf("  ATT estimate: %.4f\n", theta))
    cat(sprintf("  Std error:    %.4f\n", sigma))
    cat(sprintf("  95%% CI:       [%.4f, %.4f]\n", ci[1], ci[2]))
  }

  return(result)
}

# Alias for consistency with simulation scripts
method_forest <- att_forest
