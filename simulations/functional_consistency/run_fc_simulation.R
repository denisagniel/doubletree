#' Run Functional Consistency Simulation
#'
#' @param n Sample size
#' @param dgp DGP type: "simple", "complex", or "sparse"
#' @param method Method: "standard_msplit" or "averaged_tree"
#' @param K Number of folds for cross-fitting
#' @param M Number of sample splits
#' @param seed Random seed for replication
#'
#' @return Data frame with one row containing results
run_fc_simulation <- function(n, dgp, method, K, M, seed) {
  set.seed(seed)

  # Generate data according to DGP
  dgp_data <- generate_dgp(n, dgp)
  X <- dgp_data$X
  A <- dgp_data$A
  Y <- dgp_data$Y
  att_true <- dgp_data$att_true

  # Get regularization for this n (tune so n_ℓ ∝ n^β)
  # Use β = 0.4 as reasonable middle ground
  regularization <- 0.1

  # Run estimation
  if (method == "standard_msplit") {
    result <- run_standard_msplit(X, A, Y, M, K, regularization)
  } else if (method == "averaged_tree") {
    result <- run_averaged_tree(X, A, Y, M, K, regularization)
  } else if (method == "pattern_aggregation") {
    result <- run_pattern_aggregation(X, A, Y, M, K, regularization)
  } else {
    stop("Unknown method: ", method)
  }

  # Compute metrics
  bias <- result$att_est - att_true
  coverage <- (att_true >= result$ci_lower && att_true <= result$ci_upper)
  bias_se_ratio <- abs(bias) / result$se
  standardized_bias <- sqrt(n) * bias
  ci_width <- result$ci_upper - result$ci_lower

  # Return results
  data.frame(
    n = n,
    dgp = dgp,
    method = method,
    K = K,
    M = M,
    seed = seed,
    att_true = att_true,
    att_est = result$att_est,
    se = result$se,
    ci_lower = result$ci_lower,
    ci_upper = result$ci_upper,
    bias = bias,
    coverage = coverage,
    bias_se_ratio = bias_se_ratio,
    standardized_bias = standardized_bias,
    ci_width = ci_width,
    max_diff_e = result$max_diff_e,
    max_diff_m0 = result$max_diff_m0,
    n_leaves_e = result$n_leaves_e,
    n_leaves_m0 = result$n_leaves_m0,
    mean_leaf_size_e = result$mean_leaf_size_e,
    mean_leaf_size_m0 = result$mean_leaf_size_m0,
    stringsAsFactors = FALSE
  )
}

#' Generate data according to DGP
generate_dgp <- function(n, dgp) {
  if (dgp == "simple") {
    # Simple: 3 binary covariates, linear effects
    p <- 3
    X <- data.frame(
      x1 = rbinom(n, 1, 0.5),
      x2 = rbinom(n, 1, 0.5),
      x3 = rbinom(n, 1, 0.5)
    )

    # Treatment
    logit_e <- -0.5 + 0.8*X$x1 + 0.6*X$x2 + 0.4*X$x3
    e_true <- plogis(logit_e)
    A <- rbinom(n, 1, e_true)

    # Outcome with HTE
    tau <- 0.3 + 0.4*X$x1 + 0.3*X$x2
    m0 <- 0.2 + 0.3*X$x1 - 0.2*X$x2 + 0.1*X$x3

  } else if (dgp == "complex") {
    # Complex: 4 binary covariates, interactions
    p <- 4
    X <- data.frame(
      x1 = rbinom(n, 1, 0.5),
      x2 = rbinom(n, 1, 0.5),
      x3 = rbinom(n, 1, 0.5),
      x4 = rbinom(n, 1, 0.5)
    )

    # Treatment with interactions
    logit_e <- -0.5 + 0.8*X$x1 + 0.6*X$x2 + 0.4*X$x3 + 0.5*X$x1*X$x2
    e_true <- plogis(logit_e)
    A <- rbinom(n, 1, e_true)

    # Outcome with interactions
    tau <- 0.3 + 0.4*X$x1 + 0.3*X$x2 + 0.2*X$x1*X$x3
    m0 <- 0.2 + 0.3*X$x1 - 0.2*X$x2 + 0.1*X$x3 + 0.15*X$x2*X$x4

  } else if (dgp == "sparse") {
    # Sparse: 5 binary covariates (many patterns, small n_x)
    p <- 5
    X <- data.frame(
      x1 = rbinom(n, 1, 0.5),
      x2 = rbinom(n, 1, 0.5),
      x3 = rbinom(n, 1, 0.5),
      x4 = rbinom(n, 1, 0.5),
      x5 = rbinom(n, 1, 0.5)
    )

    # Treatment
    logit_e <- -0.5 + 0.8*X$x1 + 0.6*X$x2 + 0.4*X$x3
    e_true <- plogis(logit_e)
    A <- rbinom(n, 1, e_true)

    # Outcome (only first 3 covariates matter)
    tau <- 0.3 + 0.4*X$x1 + 0.3*X$x2
    m0 <- 0.2 + 0.3*X$x1 - 0.2*X$x2 + 0.1*X$x3

  } else {
    stop("Unknown DGP: ", dgp)
  }

  # Generate outcome
  m1 <- m0 + tau
  Y <- A * m1 + (1-A) * m0 + rnorm(n, 0, 0.1)

  # True ATT
  att_true <- mean(tau[A==1])

  list(X = X, A = A, Y = Y, att_true = att_true)
}

#' Run standard M-split
run_standard_msplit <- function(X, A, Y, M, K, regularization) {
  n <- nrow(X)

  # Use doubletree::estimate_att_msplit
  msplit_result <- doubletree::estimate_att_msplit(
    X = X, A = A, Y = Y,
    outcome_type = "continuous",
    M = M, K = K,
    regularization = regularization,
    verbose = FALSE
  )

  # Extract results
  list(
    att_est = msplit_result$theta,
    se = msplit_result$sigma,
    ci_lower = msplit_result$ci_95[1],
    ci_upper = msplit_result$ci_95[2],
    max_diff_e = msplit_result$diagnostics$functional_consistency$max_diff_e,
    max_diff_m0 = msplit_result$diagnostics$functional_consistency$max_diff_m0,
    n_leaves_e = msplit_result$diagnostics$n_leaves_e,
    n_leaves_m0 = msplit_result$diagnostics$n_leaves_m0,
    mean_leaf_size_e = n / msplit_result$diagnostics$n_leaves_e,
    mean_leaf_size_m0 = sum(A == 0) / msplit_result$diagnostics$n_leaves_m0
  )
}

#' Run averaged tree approach
run_averaged_tree <- function(X, A, Y, M, K, regularization) {
  n <- nrow(X)

  # First get modal structure using standard M-split
  msplit_result <- doubletree::estimate_att_msplit(
    X = X, A = A, Y = Y,
    outcome_type = "continuous",
    M = M, K = K,
    regularization = regularization,
    verbose = FALSE
  )

  modal_e <- msplit_result$structures$e
  modal_m0 <- msplit_result$structures$m0

  # Collect predictions from all M×K refits on FULL dataset
  e_preds_all <- matrix(NA, nrow = n, ncol = M*K)
  m0_preds_all <- matrix(NA, nrow = n, ncol = M*K)

  idx <- 1
  for (m in 1:M) {
    folds <- sample(rep(1:K, length.out = n))

    for (k in 1:K) {
      train_idx <- which(folds != k)
      X_train <- X[train_idx, ]
      A_train <- A[train_idx]
      Y_train <- Y[train_idx]

      # Refit propensity
      e_tree <- optimaltrees::refit_tree_structure(
        structure = modal_e,
        X_new = X_train,
        y_new = A_train,
        loss_function = "log_loss",
        store_training_data = FALSE
      )

      # Refit outcome
      Y0_train <- Y_train[A_train == 0]
      X0_train <- X_train[A_train == 0, ]

      m0_tree <- optimaltrees::refit_tree_structure(
        structure = modal_m0,
        X_new = X0_train,
        y_new = Y0_train,
        loss_function = "squared_error",
        store_training_data = FALSE
      )

      # Predict on FULL dataset (breaks cross-fitting)
      e_preds_all[, idx] <- predict(e_tree, X, type = "prob")[, 2]
      m0_preds_all[, idx] <- predict(m0_tree, X, type = "response")

      idx <- idx + 1
    }
  }

  # Average predictions
  e_avg <- rowMeans(e_preds_all)
  m0_avg <- rowMeans(m0_preds_all)

  # Compute ATT
  treated_idx <- which(A == 1)
  att_est <- mean(Y[treated_idx]) - mean(m0_avg[treated_idx])

  # EIF-based SE
  n1 <- sum(A)
  eif_components <- (A / (n1/n)) * (Y - m0_avg) -
    ((1 - A) * e_avg / ((1 - e_avg) * (n1/n))) * (Y - m0_avg) +
    m0_avg - att_est

  se_est <- sd(eif_components) / sqrt(n)
  ci_lower <- att_est - 1.96 * se_est
  ci_upper <- att_est + 1.96 * se_est

  # Compute functional consistency
  patterns <- do.call(paste, X)
  max_diff_e <- 0
  max_diff_m0 <- 0

  for (pattern in unique(patterns)) {
    idx_pattern <- which(patterns == pattern)
    if (length(idx_pattern) > 1) {
      max_diff_e <- max(max_diff_e, max(e_avg[idx_pattern]) - min(e_avg[idx_pattern]))
      max_diff_m0 <- max(max_diff_m0, max(m0_avg[idx_pattern]) - min(m0_avg[idx_pattern]))
    }
  }

  list(
    att_est = att_est,
    se = se_est,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    max_diff_e = max_diff_e,
    max_diff_m0 = max_diff_m0,
    n_leaves_e = msplit_result$diagnostics$n_leaves_e,
    n_leaves_m0 = msplit_result$diagnostics$n_leaves_m0,
    mean_leaf_size_e = n / msplit_result$diagnostics$n_leaves_e,
    mean_leaf_size_m0 = sum(A == 0) / msplit_result$diagnostics$n_leaves_m0
  )
}

#' Run pattern aggregation approach
run_pattern_aggregation <- function(X, A, Y, M, K, regularization) {
  n <- nrow(X)

  # Step 1: Get standard M-split cross-fit predictions
  msplit_result <- doubletree::estimate_att_msplit(
    X = X, A = A, Y = Y,
    outcome_type = "continuous",
    M = M, K = K,
    regularization = regularization,
    verbose = FALSE
  )

  # Extract cross-fit predictions (averaged across M splits, but still cross-fit)
  e_crossfit <- msplit_result$averaged_predictions$e
  m0_crossfit <- msplit_result$averaged_predictions$m0

  # Step 2: Identify covariate patterns
  patterns <- do.call(paste, X)
  unique_patterns <- unique(patterns)

  # Step 3: Compute pattern-level averages
  e_pattern <- numeric(n)
  m0_pattern <- numeric(n)

  for (pattern in unique_patterns) {
    idx_pattern <- which(patterns == pattern)
    # Average cross-fit predictions within this pattern
    e_pattern[idx_pattern] <- mean(e_crossfit[idx_pattern])
    m0_pattern[idx_pattern] <- mean(m0_crossfit[idx_pattern])
  }

  # Step 4: Compute ATT using pattern-averaged predictions
  treated_idx <- which(A == 1)
  att_est <- mean(Y[treated_idx]) - mean(m0_pattern[treated_idx])

  # EIF-based SE using pattern predictions
  n1 <- sum(A)
  eif_components <- (A / (n1/n)) * (Y - m0_pattern) -
    ((1 - A) * e_pattern / ((1 - e_pattern) * (n1/n))) * (Y - m0_pattern) +
    m0_pattern - att_est

  se_est <- sd(eif_components) / sqrt(n)
  ci_lower <- att_est - 1.96 * se_est
  ci_upper <- att_est + 1.96 * se_est

  # Step 5: Compute functional consistency (should be perfect by construction)
  max_diff_e <- 0
  max_diff_m0 <- 0

  for (pattern in unique_patterns) {
    idx_pattern <- which(patterns == pattern)
    if (length(idx_pattern) > 1) {
      # By construction, all observations in same pattern have same prediction
      max_diff_e <- max(max_diff_e, max(e_pattern[idx_pattern]) - min(e_pattern[idx_pattern]))
      max_diff_m0 <- max(max_diff_m0, max(m0_pattern[idx_pattern]) - min(m0_pattern[idx_pattern]))
    }
  }

  list(
    att_est = att_est,
    se = se_est,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    max_diff_e = max_diff_e,
    max_diff_m0 = max_diff_m0,
    n_leaves_e = msplit_result$diagnostics$n_leaves_e,
    n_leaves_m0 = msplit_result$diagnostics$n_leaves_m0,
    mean_leaf_size_e = n / msplit_result$diagnostics$n_leaves_e,
    mean_leaf_size_m0 = sum(A == 0) / msplit_result$diagnostics$n_leaves_m0
  )
}
