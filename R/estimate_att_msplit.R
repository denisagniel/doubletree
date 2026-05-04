#' M-Split Doubletree ATT Estimation
#'
#' @description
#' Estimates ATT using M-split algorithm: select modal tree structure across M
#' independent cross-fits, refit that structure on all M splits, average predictions,
#' compute ATT with averaged nuisances.
#'
#' @param X Data.frame or matrix of covariates (binary features for optimaltrees)
#' @param A Integer or numeric vector of treatment (0/1)
#' @param Y Numeric vector of outcome
#' @param M Number of independent sample splits (default 10)
#' @param K Number of cross-validation folds per split (default 5)
#' @param structure_selection "modal" (most frequent), "first", or "lowest_risk" (default "modal")
#' @param seed_base Base seed for reproducibility. Split m uses seed = seed_base + m.
#' @param verbose Logical. Print progress (default TRUE)
#' @param regularization Numeric. Tree complexity penalty (default 0.1)
#' @param outcome_type "binary" or "continuous" (default "binary")
#'
#' @return List with class "msplit_att" containing:
#' \describe{
#'   \item{theta}{Numeric: ATT estimate}
#'   \item{sigma}{Numeric: standard error}
#'   \item{ci_95}{Numeric vector: 95 percent confidence interval}
#'   \item{score_values}{Numeric vector: EIF scores at theta}
#'   \item{structures}{List with e and m0 TreeStructure objects}
#'   \item{predictions_all_splits}{List with e and m0 matrices (n x M)}
#'   \item{averaged_predictions}{List with e and m0 vectors}
#'   \item{diagnostics}{List with structure frequencies, prediction variances, etc.}
#'   \item{M}{Number of splits}
#'   \item{K}{Number of folds}
#'   \item{n}{Sample size}
#'   \item{n_treated}{Number treated}
#'   \item{structure_selection}{Method used}
#' }
#'
#' @details
#' ## Algorithm
#'
#' 1. **Stage 1: Structure Selection**
#'    - Run M independent cross-fits (each with K folds)
#'    - Extract tree structures for e(X) and m0(X) from each split
#'    - Select modal (most frequent) structure for each nuisance
#'
#' 2. **Stage 2: Refit on All M Splits**
#'    - For each split m = 1, ..., M:
#'      - Create K folds
#'      - For each fold k:
#'        - Refit selected e structure on train set
#'        - Refit selected m0 structure on control train set
#'        - Predict on test fold
#'    - Result: n x M matrix of predictions
#'
#' 3. **Stage 3: Average and Compute ATT**
#'    - Average predictions: e_bar = rowMeans(e_matrix)
#'    - Compute EIF score with averaged nuisances
#'    - Estimate theta and standard error
#'
#' ## Diagnostics
#'
#' - Structure frequencies: How often modal structure appears (higher = more stable)
#' - Prediction variance: Var(predictions across M) per observation (lower = more stable)
#' - Functional consistency: max|μ_i - μ_j| for X_i = X_j (should be ~0)
#'
#' @references
#' Theory developed in m-split-theory.tex (4 theorems).
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 400
#' X <- data.frame(x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5))
#' A <- rbinom(n, 1, 0.4)
#' Y <- rbinom(n, 1, 0.3 + 0.2 * A + 0.1 * X$x1)
#'
#' result <- estimate_att_msplit(X, A, Y, M = 10, K = 5, seed_base = 100)
#' print(result)
#' }
#'
#' @export
estimate_att_msplit <- function(X, A, Y,
                                M = 10,
                                K = 5,
                                structure_selection = "modal",
                                seed_base = NULL,
                                verbose = TRUE,
                                regularization = 0.1,
                                outcome_type = "binary") {
  # Validate inputs
  n <- nrow(X)

  if (!is.data.frame(X) && !is.matrix(X)) {
    stop("X must be a data.frame or matrix", call. = FALSE)
  }
  if (length(A) != n) {
    stop("length(A) must equal nrow(X)", call. = FALSE)
  }
  if (length(Y) != n) {
    stop("length(Y) must equal nrow(X)", call. = FALSE)
  }
  if (M < 1) {
    stop("M must be at least 1", call. = FALSE)
  }
  if (K < 2) {
    stop("K must be at least 2", call. = FALSE)
  }

  if (verbose) {
    cat(sprintf("M-split doubletree: M=%d splits, K=%d folds, n=%d\n", M, K, n))
  }

  # ============================================================
  # Stage 1: Structure Selection
  # ============================================================
  if (verbose) cat("Stage 1: Selecting modal structures...\n")

  structures_e <- vector("list", M)
  structures_m0 <- vector("list", M)

  for (m in seq_len(M)) {
    seed_m <- if (!is.null(seed_base)) seed_base + m else NULL

    # Create folds for this split
    folds_m <- create_folds(n, K, strata = A, seed = seed_m)

    # We just need structures, so fit one tree per nuisance (on first fold's training set)
    # Simplified: use fold 1's training set to get a representative structure
    train_idx <- which(folds_m != 1)

    X_train <- X[train_idx, , drop = FALSE]
    A_train <- A[train_idx]
    Y_train <- Y[train_idx]

    # Fit propensity model
    model_e <- optimaltrees::fit_tree(
      X = X_train,
      y = A_train,
      loss_function = "log_loss",
      regularization = regularization,
      store_training_data = TRUE,
      verbose = FALSE
    )
    structures_e[[m]] <- optimaltrees::extract_tree_structure(model_e)

    # Fit outcome model (control outcomes only)
    control_idx <- which(A_train == 0)
    X_control <- X_train[control_idx, , drop = FALSE]
    Y_control <- Y_train[control_idx]

    outcome_loss <- if (outcome_type == "binary") "log_loss" else "squared_error"

    model_m0 <- optimaltrees::fit_tree(
      X = X_control,
      y = Y_control,
      loss_function = outcome_loss,
      regularization = regularization,
      store_training_data = TRUE,
      verbose = FALSE
    )
    structures_m0[[m]] <- optimaltrees::extract_tree_structure(model_m0)

    if (verbose && m %% max(1, M %/% 10) == 0) {
      cat(sprintf("  Extracted structures from %d/%d splits\n", m, M))
    }
  }

  # Select modal structures
  s_star_e <- select_structure_modal(structures_e)
  s_star_m0 <- select_structure_modal(structures_m0)

  if (verbose) {
    cat(sprintf("  Propensity: modal structure appears in %.1f%% of splits\n",
                s_star_e$frequency * 100))
    cat(sprintf("  Outcome: modal structure appears in %.1f%% of splits\n",
                s_star_m0$frequency * 100))
  }

  # ============================================================
  # Stage 2: Refit Fixed Structures on All M Splits
  # ============================================================
  if (verbose) cat("Stage 2: Refitting structures on all M splits...\n")

  predictions_e <- matrix(NA_real_, nrow = n, ncol = M)
  predictions_m0 <- matrix(NA_real_, nrow = n, ncol = M)

  for (m in seq_len(M)) {
    seed_m <- if (!is.null(seed_base)) seed_base + m else NULL
    folds_m <- create_folds(n, K, strata = A, seed = seed_m)

    for (k in seq_len(K)) {
      test_idx <- which(folds_m == k)
      train_idx <- which(folds_m != k)

      # Refit propensity tree
      X_train <- X[train_idx, , drop = FALSE]
      A_train <- A[train_idx]

      tree_e_mk <- optimaltrees::refit_tree_structure(
        s_star_e$structure,
        X_train, A_train,
        loss_function = "log_loss",
        store_training_data = FALSE
      )

      # Refit outcome tree (control outcomes only)
      control_idx <- which(A_train == 0)
      Y_train_control <- Y[train_idx][A_train == 0]
      X_train_control <- X_train[control_idx, , drop = FALSE]

      outcome_loss <- if (outcome_type == "binary") "log_loss" else "squared_error"

      tree_m0_mk <- optimaltrees::refit_tree_structure(
        s_star_m0$structure,
        X_train_control, Y_train_control,
        loss_function = outcome_loss,
        store_training_data = FALSE
      )

      # Predict on test fold
      X_test <- X[test_idx, , drop = FALSE]

      preds_e <- predict(tree_e_mk, X_test, type = "prob")
      predictions_e[test_idx, m] <- preds_e[, 2]  # P(A=1|X)

      if (outcome_loss == "log_loss") {
        preds_m0 <- predict(tree_m0_mk, X_test, type = "prob")
        predictions_m0[test_idx, m] <- preds_m0[, 2]  # P(Y=1|A=0,X)
      } else {
        # Regression: use get_fitted_from_tree directly
        preds_m0 <- optimaltrees::get_fitted_from_tree(tree_m0_mk@trees[[1]], X_test)
        predictions_m0[test_idx, m] <- preds_m0  # E[Y|A=0,X]
      }
    }

    if (verbose && m %% max(1, M %/% 10) == 0) {
      cat(sprintf("  Completed %d/%d splits\n", m, M))
    }
  }

  # Average predictions
  e_bar <- rowMeans(predictions_e, na.rm = FALSE)
  m0_bar <- rowMeans(predictions_m0, na.rm = FALSE)

  # ============================================================
  # Stage 3: Compute ATT with Averaged Nuisances
  # ============================================================
  if (verbose) cat("Stage 3: Computing ATT...\n")

  eta_bar <- list(e = e_bar, m0 = m0_bar, m1 = NULL)
  pi_hat <- mean(A)

  # EIF score
  score <- psi_att(Y, A, theta = 0, eta = eta_bar, pi_hat = pi_hat)
  theta_msplit <- sum(score) / sum(A / pi_hat)

  # Standard error: SE(θ̂) = sqrt(Var[ψ] / n) = sqrt(E[(ψ - E[ψ])²] / n)
  score_centered <- score - mean(score)
  sigma_msplit <- sqrt(mean(score_centered^2) / n)

  # Confidence interval
  ci_95 <- theta_msplit + c(-1, 1) * qnorm(0.975) * sigma_msplit / sqrt(n)

  # ============================================================
  # Diagnostics
  # ============================================================
  diagnostics <- list(
    structure_frequency_e = s_star_e$frequency,
    structure_frequency_m0 = s_star_m0$frequency,
    structure_counts_e = s_star_e$counts,
    structure_counts_m0 = s_star_m0$counts,
    n_leaves_e = s_star_e$structure@n_leaves,
    n_leaves_m0 = s_star_m0$structure@n_leaves,
    max_depth_e = s_star_e$structure@max_depth,
    max_depth_m0 = s_star_m0$structure@max_depth,
    prediction_variance_e = apply(predictions_e, 1, var, na.rm = TRUE),
    prediction_variance_m0 = apply(predictions_m0, 1, var, na.rm = TRUE),
    mean_prediction_variance_e = mean(apply(predictions_e, 1, var, na.rm = TRUE)),
    mean_prediction_variance_m0 = mean(apply(predictions_m0, 1, var, na.rm = TRUE)),
    functional_consistency = compute_functional_consistency(predictions_e, predictions_m0, X)
  )

  # ============================================================
  # Return
  # ============================================================
  result <- list(
    theta = theta_msplit,
    sigma = sigma_msplit,
    ci_95 = ci_95,
    score_values = score,

    # Structures
    structures = list(
      e = s_star_e$structure,
      m0 = s_star_m0$structure
    ),

    # Predictions
    predictions_all_splits = list(
      e = predictions_e,
      m0 = predictions_m0
    ),
    averaged_predictions = list(
      e = e_bar,
      m0 = m0_bar
    ),

    # Diagnostics
    diagnostics = diagnostics,

    # Metadata
    M = M,
    K = K,
    n = n,
    n_treated = sum(A),
    structure_selection = structure_selection,
    outcome_type = outcome_type,
    regularization = regularization
  )

  class(result) <- c("msplit_att", "list")
  result
}

#' Print Method for msplit_att
#'
#' @param x msplit_att object
#' @param ... Additional arguments (ignored)
#' @export
print.msplit_att <- function(x, ...) {
  cat("M-Split Doubletree ATT Estimation\n")
  cat("=================================\n\n")
  cat(sprintf("Estimate:  %.4f\n", x$theta))
  cat(sprintf("Std Error: %.4f\n", x$sigma))
  cat(sprintf("95%% CI:    [%.4f, %.4f]\n", x$ci_95[1], x$ci_95[2]))
  cat(sprintf("\nSample: n=%d, n_treated=%d, M=%d splits, K=%d folds\n",
              x$n, x$n_treated, x$M, x$K))

  cat("\nStructure Selection:\n")
  cat(sprintf("  Propensity: modal frequency = %.1f%%, %d leaves (depth %d)\n",
              x$diagnostics$structure_frequency_e * 100,
              x$diagnostics$n_leaves_e,
              x$diagnostics$max_depth_e))
  cat(sprintf("  Outcome:    modal frequency = %.1f%%, %d leaves (depth %d)\n",
              x$diagnostics$structure_frequency_m0 * 100,
              x$diagnostics$n_leaves_m0,
              x$diagnostics$max_depth_m0))

  cat("\nFunctional Consistency (max|μ̄ᵢ-μ̄ⱼ| for Xᵢ=Xⱼ):\n")
  cat(sprintf("  e:  %.6f\n", x$diagnostics$functional_consistency$max_diff_e))
  cat(sprintf("  m0: %.6f\n", x$diagnostics$functional_consistency$max_diff_m0))

  cat("\nPrediction Variance (mean across observations):\n")
  cat(sprintf("  e:  %.6f\n", x$diagnostics$mean_prediction_variance_e))
  cat(sprintf("  m0: %.6f\n", x$diagnostics$mean_prediction_variance_m0))

  invisible(x)
}
