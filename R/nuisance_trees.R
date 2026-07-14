#' Fit tree with optional CV for regularization selection
#' @param X Data.frame or matrix of features.
#' @param y Vector of outcomes.
#' @param loss_function Character. Loss function for optimaltrees.
#' @param cv_reg Logical. Whether to use CV for regularization selection.
#' @param cv_K Integer. Number of CV folds for lambda selection.
#' @param reg Numeric. Fixed regularization value (used if cv_reg = FALSE).
#' @param discretize_method Character. Discretization method.
#' @param discretize_bins Integer or "adaptive". Number of bins.
#' @param verbose Logical. Print messages.
#' @param fold_id Integer. Current fold ID (for messages).
#' @param model_name Character. Name of model (for messages, e.g. "e", "m0").
#' @param ... Additional arguments passed to optimaltrees functions.
#' @return optimaltrees model object.
#' @noRd
fit_tree_with_cv <- function(X, y, loss_function, cv_reg, cv_K, reg,
                             discretize_method, discretize_bins,
                             verbose, fold_id, model_name, ...) {
  if (cv_reg) {
    cv_result <- optimaltrees::cv_regularization_adaptive(
      X, y,
      loss_function = loss_function,
      discretize_method = discretize_method,
      discretize_bins = discretize_bins,
      K = cv_K,
      max_iterations = 10,
      refit = TRUE,
      verbose = FALSE,
      ...
    )
    if (verbose) {
      message("  Fold ", fold_id, " ", model_name,
              ": selected lambda = ", round(cv_result$best_lambda, 5),
              " (", cv_result$iterations, " iter, ",
              if (cv_result$converged) "converged" else "max iter", ")")
    }
    return(cv_result$model)
  } else {
    return(optimaltrees::fit_tree(
      X, y,
      loss_function = loss_function,
      discretize_method = discretize_method,
      discretize_bins = discretize_bins,
      regularization = reg,
      verbose = verbose,
      ...
    ))
  }
}

#' Fit nuisance trees for one cross-fitting fold
#'
#' Fits propensity e(X) and outcome regression m0(X) on training data
#' (all rows not in fold_id) using optimaltrees. Propensity uses log_loss; outcome
#' trees use log_loss (binary Y) or squared_error (continuous Y). Used internally by estimate_att.
#'
#' @param X Data.frame or matrix of covariates (binary 0/1).
#' @param A Integer vector of treatment (0/1).
#' @param Y Numeric vector of outcome (binary 0/1 or continuous).
#' @param fold_id Integer. Fold index (1..K) to leave out; training = fold_indices != fold_id.
#' @param fold_indices Integer vector of length nrow(X) with fold assignment.
#' @param outcome_type Character. "binary" or "continuous"; determines loss for m0.
#' @param regularization Numeric. Penalty per leaf for optimaltrees. Default 0.1.
#' @param cv_regularization Logical. If TRUE, use cross-validation to select lambda. Default FALSE.
#' @param cv_K Integer. Number of CV folds for lambda selection. Default 5.
#' @param verbose Logical. Passed to optimaltrees. Default FALSE.
#' @param discretize_method Character. "quantiles" for quantile-based discretization.
#' @param discretize_bins Integer or "adaptive". Number of bins for discretization.
#' @param ... Additional arguments passed to optimaltrees::fit_tree.
#' @return List with elements e_model, m0_model (optimaltrees fit objects),
#'   discretization_e, discretization_m0 (breaks for test set), outcome_type.
#' @noRd
fit_nuisances_fold <- function(X, A, Y, fold_id, fold_indices, outcome_type = "binary",
                               regularization = 0.1, cv_regularization = FALSE, cv_K = 5,
                               verbose = FALSE,
                               discretize_method = "quantiles", discretize_bins = "adaptive",
                               ...) {
  train_idx <- which(fold_indices != fold_id)
  X_tr <- X[train_idx, , drop = FALSE]
  A_tr <- A[train_idx]
  Y_tr <- Y[train_idx]
  n0 <- sum(A_tr == 0)
  n1 <- sum(A_tr == 1)
  if (n0 == 0) {
    stop("No control units (A=0) in training set for fold ", fold_id, "; use fewer folds or check data.")
  }
  if (n1 == 0) {
    stop("No treated units (A=1) in training set for fold ", fold_id, "; use fewer folds or check data.")
  }
  loss_outcome <- if (outcome_type == "continuous") "squared_error" else "log_loss"

  # Fit propensity model (optimaltrees handles discretization with threshold encoding)
  e_model <- fit_tree_with_cv(
    X_tr, A_tr,
    loss_function = "log_loss",
    cv_reg = cv_regularization,
    cv_K = cv_K,
    reg = regularization,
    discretize_method = discretize_method,
    discretize_bins = discretize_bins,
    verbose = verbose,
    fold_id = fold_id,
    model_name = "e",
    ...
  )

  # Check for model fitting failure using n_trees instead of tree_json
  # (tree_json may be NULL for models with discretization, which is OK)
  # optimaltrees is now fully S7, so use @ accessor
  n_trees_e <- e_model@n_trees

  if (is.null(n_trees_e) || n_trees_e == 0) {
    stop(
      "TreeFARMS model fitting failed for propensity model in fold ", fold_id, ".\n",
      "This indicates:\n",
      "  - No valid trees found within regularization constraints\n",
      "  - Or data quality issues (extreme values, collinearity)\n",
      "\n",
      "Diagnostic info:\n",
      "  Sample size: ", nrow(X_tr), "\n",
      "  Features: ", ncol(X_tr), "\n",
      "  Regularization: ", regularization, "\n",
      "\n",
      "Suggested actions:\n",
      "  - Increase regularization to allow simpler trees\n",
      "  - Check for data anomalies or perfect separation\n",
      "  - Verify features have variation\n",
      call. = FALSE
    )
  }

  # Fit m0 model on controls (optimaltrees handles discretization with threshold encoding)
  X_tr_controls <- X_tr[A_tr == 0, , drop = FALSE]
  m0_model <- fit_tree_with_cv(
    X_tr_controls, Y_tr[A_tr == 0],
    loss_function = loss_outcome,
    cv_reg = cv_regularization,
    cv_K = cv_K,
    reg = regularization,
    discretize_method = discretize_method,
    discretize_bins = discretize_bins,
    verbose = verbose,
    fold_id = fold_id,
    model_name = "m0",
    ...
  )

  # Check for model fitting failure using n_trees instead of tree_json
  # (tree_json may be NULL for models with discretization, which is OK)
  # optimaltrees is now fully S7, so use @ accessor
  n_trees_m0 <- m0_model@n_trees

  if (is.null(n_trees_m0) || n_trees_m0 == 0) {
    stop(
      "TreeFARMS model fitting failed for control outcome model (m0) in fold ", fold_id, ".\n",
      "This indicates:\n",
      "  - No valid trees found within regularization constraints\n",
      "  - Or data quality issues (extreme values, collinearity)\n",
      "\n",
      "Diagnostic info:\n",
      "  Sample size: ", nrow(X_tr_controls), "\n",
      "  Features: ", ncol(X_tr_controls), "\n",
      "  Regularization: ", regularization, "\n",
      "\n",
      "Suggested actions:\n",
      "  - Increase regularization to allow simpler trees\n",
      "  - Check for data anomalies or perfect separation\n",
      "  - Verify features have variation\n",
      call. = FALSE
    )
  }

  # Issue #31: Nuisance quality diagnostics can now be re-enabled.
  # predict() was fixed in Batch 1 to handle discretization automatically.
  # Diagnostics would compare in-sample predictions to training data.
  # Currently disabled for simplicity - can be added if needed.

  list(
    e_model = e_model,
    m0_model = m0_model,
    outcome_type = outcome_type
    # Note: discretization metadata now stored in model$discretization_metadata (optimaltrees)
  )
}

#' Predict nuisances for given rows using a fold's fitted models
#'
#' Returns e(X), m0(X) for the requested rows. Propensity e from
#' predict(..., type = "prob"); m0 from type = "prob" (binary) or default predict (continuous).
#'
#' @param models List from fit_nuisances_fold (e_model, m0_model, outcome_type).
#' @param X Full covariate matrix (we subset to fold_rows).
#' @param fold_rows Integer vector of row indices to predict for.
#' @return List with elements e, m0 (numeric vectors of length length(fold_rows)).
#' @noRd
predict_nuisances_fold <- function(models, X, fold_rows) {
  if (length(fold_rows) == 0) {
    return(list(e = numeric(0), m0 = numeric(0)))
  }

  X_sub <- X[fold_rows, , drop = FALSE]
  outcome_type <- if (!is.null(models$outcome_type)) models$outcome_type else "binary"

  # Issue #29: Validate outcome_type
  if (!outcome_type %in% c("binary", "continuous")) {
    stop("Invalid outcome_type: ", outcome_type,
         ". Must be 'binary' or 'continuous'.", call. = FALSE)
  }

  # Note: discretization is now handled automatically by predict() using
  # the model's stored discretization metadata. No need to manually discretize.

  # Predict propensity scores (predict() will apply discretization if needed)
  pe <- predict(models$e_model, X_sub, type = "prob")

  # Validate propensity prediction format
  if (!is.matrix(pe) || ncol(pe) != 2) {
    stop("Propensity model predict() returned unexpected format. ",
         "Expected 2-column matrix (class 0/1 probabilities), got: ",
         if (is.matrix(pe)) paste0("matrix with ", ncol(pe), " columns")
         else paste0(class(pe), " (not a matrix)"),
         call. = FALSE)
  }
  e_vec <- pe[, 2L]

  # Predict control outcomes
  if (outcome_type == "continuous") {
    pm0 <- predict(models$m0_model, X_sub)
    if (!is.numeric(pm0) || length(pm0) != nrow(X_sub)) {
      stop("Control outcome model predict() returned unexpected format. ",
           "Expected numeric vector of length ", nrow(X_sub), ", got: ",
           class(pm0), " with length ", length(pm0),
           call. = FALSE)
    }
    m0_vec <- as.numeric(pm0)
  } else {
    pm0 <- predict(models$m0_model, X_sub, type = "prob")
    if (!is.matrix(pm0) || ncol(pm0) != 2) {
      stop("Control outcome model predict() returned unexpected format. ",
           "Expected 2-column matrix (class 0/1 probabilities), got: ",
           if (is.matrix(pm0)) paste0("matrix with ", ncol(pm0), " columns")
           else paste0(class(pm0), " (not a matrix)"),
           call. = FALSE)
    }
    m0_vec <- pm0[, 2L]
  }
  list(e = e_vec, m0 = m0_vec)
}

#' Build fold-specific nuisance predictions for all observations
#'
#' For each observation i, get eta^{(-k(i))}(X_i) using the models fitted for fold k(i).
#' Propensity scores are clamped to [e_min, e_max] at prediction time to ensure
#' numerical stability in downstream score computation.
#'
#' @param nuisance_fits List of length K; each element is list(e_model, m0_model).
#' @param X Data.frame or matrix of covariates.
#' @param fold_indices Integer vector of length n (fold id per row).
#' @param e_min Lower bound for propensity clamping.
#' @param e_max Upper bound for propensity clamping.
#' @return List with elements e, m0 (each numeric vector of length n).
#' @noRd
get_fold_specific_eta <- function(nuisance_fits, X, fold_indices,
                                   e_min = .PROPENSITY_LOWER_BOUND,
                                   e_max = .PROPENSITY_UPPER_BOUND) {
  n <- nrow(X)
  K <- length(nuisance_fits)
  e_out <- numeric(n)
  m0_out <- numeric(n)
  for (k in seq_len(K)) {
    idx_k <- which(fold_indices == k)
    if (length(idx_k) == 0) next
    pred <- predict_nuisances_fold(nuisance_fits[[k]], X, idx_k)
    e_out[idx_k] <- pred$e
    m0_out[idx_k] <- pred$m0
  }

  # Clamp propensities at prediction time for numerical stability
  e_clamped <- pmax(e_min, pmin(e_max, e_out))

  list(e = e_clamped, m0 = m0_out)
}

#' Safe Rashomon fitting with automatic fallback
#' @param X Data.frame or matrix of features.
#' @param y Vector of outcomes.
#' @param K Integer. Number of folds.
#' @param loss_function Character. Loss function for optimaltrees.
#' @param regularization Numeric. Regularization parameter.
#' @param rashomon_bound_multiplier Numeric. Rashomon tolerance.
#' @param rashomon_bound_adder Numeric. Additive Rashomon bound.
#' @param max_leaves Integer or NULL. Maximum leaves constraint.
#' @param fold_indices Integer vector. Fold assignments.
#' @param auto_tune_intersecting Logical. Auto-tune epsilon_n.
#' @param verbose Logical. Print messages.
#' @param nuisance_name Character. Name for error messages ("propensity" or "control outcome").
#' @param ... Additional arguments passed to cross_fitted_rashomon.
#' @return cross_fitted_rashomon object or NULL if fallback needed.
#' @noRd
NULL

#' Escalate the Rashomon tolerance until the cross-fold intersection is non-empty
#'
#' Keeps lambda FIXED (at the CV-selected value) and widens only the Rashomon bound
#' multiplier epsilon_n = c * (log(n)/n), stepping c upward over `c_grid` until the
#' K-fold intersection is non-empty (or the grid is exhausted). Motivation: the margin
#' theorem guarantees tau* is in the intersection with probability -> 1, so an empty
#' intersection is a finite-sample artifact; widening epsilon_n resolves it.
#'
#' \strong{Validity caveat.} For any FIXED c, epsilon_n = c*log(n)/n is still
#' o(n^{-1/2}) asymptotically, so escalation is asymptotically inference-valid. But
#' c is chosen HERE FROM THE DATA (the first c yielding a non-empty intersection),
#' which is a post-selection step, and at realistic n a large c (up to 1000) makes
#' epsilon_n large in finite samples. Escalation therefore TRADES the fixed-epsilon_n
#' finite-sample validity guarantee for a non-empty intersection. This function does
#' NOT cap c and does NOT silently fall back to fold-specific trees; it RETURNS the
#' selected c so downstream analysis can evaluate empirically whether large c degrades
#' coverage (see the escalation coverage sweep in the arbitration simulation). Callers
#' opt into this behavior via escalate_intersection = TRUE; it is not the default.
#'
#' @return List with `fit` (CFRashomon or NULL if empty even at c_max) and `c`
#'   (the multiplier that first yielded a non-empty intersection, or NA if none).
#' @keywords internal
# Default escalation grid for the Rashomon tolerance multiplier c (epsilon_n =
# c*log(n)/n): geometric from the theory value c=1 up to c=1000.
.RASHOMON_C_GRID <- c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1000)

escalate_rashomon_intersection <- function(X, y, K, loss_function, regularization,
                                            fold_indices, max_leaves, verbose,
                                            nuisance_name, max_depth = 4L,
                                            c_grid = .RASHOMON_C_GRID,
                                            ...) {
  n <- nrow(X)
  eps_base <- optimaltrees::select_epsilon_n(n)   # = log(n)/n
  for (c_val in c_grid) {
    fit <- safe_rashomon_fit(
      X, y, K, loss_function, regularization,
      rashomon_bound_multiplier = c_val * eps_base, rashomon_bound_adder = 0,
      max_leaves = max_leaves, fold_indices = fold_indices,
      auto_tune_intersecting = FALSE,   # we escalate c ourselves; lambda stays fixed
      verbose = verbose, nuisance_name = nuisance_name, max_depth = max_depth, ...
    )
    if (!is.null(fit)) {
      if (verbose && c_val > 1) {
        message("  ", nuisance_name, ": intersection non-empty at c = ", c_val,
                " (epsilon_n = ", signif(c_val * eps_base, 4), ")")
      }
      return(list(fit = fit, c = c_val))
    }
  }
  if (verbose) {
    message("  ", nuisance_name, ": intersection EMPTY even at c = ",
            max(c_grid), " (epsilon_n = ", signif(max(c_grid) * eps_base, 4), ")")
  }
  list(fit = NULL, c = NA_real_)
}

safe_rashomon_fit <- function(X, y, K, loss_function, regularization,
                              rashomon_bound_multiplier, rashomon_bound_adder,
                              max_leaves, fold_indices, auto_tune_intersecting,
                              verbose, nuisance_name, ...) {
  tryCatch({
    out <- optimaltrees::cross_fitted_rashomon(
      X, y,
      K = K,
      loss_function = loss_function,
      regularization = regularization,
      rashomon_bound_multiplier = rashomon_bound_multiplier,
      rashomon_bound_adder = rashomon_bound_adder,
      rashomon_ignore_trivial_extensions = FALSE,  # CRITICAL: Keep all trivial extensions for fold intersection
      max_leaves = max_leaves,
      fold_indices = fold_indices,
      auto_tune_intersecting = auto_tune_intersecting,
      verbose = verbose,
      ...
    )
    if (out@n_intersecting > 0) out else NULL
  }, error = function(e) {
    msg <- conditionMessage(e)

    # This function is called inside escalate_rashomon_intersection's c-grid loop,
    # where a NULL return is the EXPECTED signal for "no intersection at this
    # tolerance; try the next c". So most errors here are benign non-convergence and
    # must NOT stop or warn (that would be noisy across the grid).
    #
    # Genuine data/parameter errors, however, will recur at every c and should stop
    # immediately with a clear message. We detect them by message stem. This is a
    # heuristic proxy (optimaltrees does not throw classed conditions); broadened to
    # cover the validation stems its checks actually emit. If optimaltrees later
    # exposes condition classes, switch to inherits()-based dispatch here.
    data_error_pattern <- paste(
      "must be", "must have", "must equal", "must contain", "invalid", "cannot",
      "not allowed", "required", "length of", "no control units", "no treated",
      sep = "|"
    )
    if (grepl(data_error_pattern, msg, ignore.case = TRUE)) {
      stop("Error fitting ", nuisance_name, " model: ", msg,
           "\nCheck your data and parameters.", call. = FALSE)
    }

    # Intersection empty or optimization issues: fallback (expected during escalation)
    if (verbose) {
      message("Rashomon fitting failed for ", nuisance_name, ": ", msg,
              "\nFalling back to fold-specific trees.")
    }
    return(NULL)
  })
}

#' Fit nuisances via Rashomon intersection (cross_fitted_rashomon)
#'
#' Fits e(X), m0(X) using optimaltrees::cross_fitted_rashomon with the same
#' fold assignment as the ATT estimation. When a nuisance has no intersecting trees (n_intersecting == 0),
#' falls back to single-tree-per-fold for that nuisance via fit_nuisances_fold.
#'
#' @param X Data.frame or matrix of covariates (binary 0/1).
#' @param A Integer vector of treatment (0/1).
#' @param Y Numeric vector of outcome (binary 0/1 or continuous).
#' @param fold_indices Integer vector of length nrow(X) with fold id per row (1..K).
#' @param outcome_type Character. "binary" or "continuous"; determines loss for m0.
#' @param regularization Numeric. Passed to optimaltrees. Default 0.1.
#' @param cv_regularization Logical. If TRUE, use CV to select lambda for each nuisance. Default FALSE.
#' @param cv_K Integer. Number of CV folds for lambda selection. Default 5.
#' @param verbose Logical. Passed to optimaltrees. Default FALSE.
#' @param rashomon_bound_multiplier,rashomon_bound_adder,max_leaves,auto_tune_intersecting Passed to cross_fitted_rashomon.
#' @param escalate_intersection Logical. If TRUE (and no explicit rashomon_bound_multiplier),
#'   widen epsilon_n = c*log(n)/n over .RASHOMON_C_GRID until the cross-fold intersection
#'   is non-empty (see escalate_rashomon_intersection). Trades the fixed-epsilon_n validity
#'   guarantee for a non-empty intersection; opt-in, not the default.
#' @param ... Additional arguments passed to optimaltrees::cross_fitted_rashomon.
#' @return List with cf_e, cf_m0 (cf_rashomon or NULL if fallback), fallback_fits (list of K per-fold fits or NULL), outcome_type.
#' @noRd
fit_nuisances_rashomon <- function(X, A, Y, fold_indices, outcome_type = "binary",
                                   regularization = 0.1, cv_regularization = FALSE, cv_K = 5,
                                   verbose = FALSE,
                                   rashomon_bound_multiplier = NULL, rashomon_bound_adder = 0,
                                   max_leaves = NULL, auto_tune_intersecting = FALSE,
                                   escalate_intersection = FALSE,
                                   discretize_method = "quantiles", discretize_bins = "adaptive",
                                   max_depth = 4L,
                                   ...) {
  n <- nrow(X)
  # Did the caller pass an explicit tolerance? Capture BEFORE resolving NULL, since
  # an explicit multiplier pins a single fixed tolerance (no escalation) even when
  # escalate_intersection = TRUE (see c_grid selection below).
  explicit_multiplier <- !is.null(rashomon_bound_multiplier)
  # Defensive resolution: theory epsilon_n = log(n)/n if not supplied by caller.
  if (is.null(rashomon_bound_multiplier)) {
    rashomon_bound_multiplier <- optimaltrees::select_epsilon_n(n)
  }
  K <- max(fold_indices, na.rm = TRUE)
  loss_outcome <- if (outcome_type == "continuous") "squared_error" else "log_loss"

  # If CV requested, select regularization for each nuisance
  if (cv_regularization) {
    if (verbose) message("Selecting regularization via adaptive CV...")

    # Cap lambda at 20x theory value to prevent over-regularization to stumps.
    # Without this cap, cv_regularization_adaptive can select lambda >> log(n)/n,
    # making the stump optimal and yielding degenerate (constant) nuisance predictions.
    # 20x (not 15x) ensures cap > floor at all n up to 2000: cap=20*log(n)/n, floor=log(n)/n.
    # Consistent with the cap used in estimate_att_msplit_averaged.
    max_lambda_cap <- (log(n) / n) * 20

    # Adaptive CV for propensity
    cv_e <- optimaltrees::cv_regularization_adaptive(X, A, loss_function = "log_loss",
                                         K = cv_K, max_iterations = 10,
                                         refit = FALSE, verbose = FALSE,
                                         max_lambda = max_lambda_cap, ...)
    reg_e <- cv_e$best_lambda
    if (verbose) {
      message("  e: selected lambda = ", round(reg_e, 5),
              " (", cv_e$iterations, " iter, ",
              if (cv_e$converged) "converged" else "max iter", ")")
    }

    # Adaptive CV for m0
    idx0 <- which(A == 0)
    X0 <- X[idx0, , drop = FALSE]
    Y0 <- Y[idx0]
    n0 <- length(idx0)
    max_lambda_cap_m0 <- (log(n0) / n0) * 20
    cv_m0 <- optimaltrees::cv_regularization_adaptive(X0, Y0, loss_function = loss_outcome,
                                          K = cv_K, max_iterations = 10,
                                          refit = FALSE, verbose = FALSE,
                                          max_lambda = max_lambda_cap_m0, ...)
    reg_m0 <- cv_m0$best_lambda
    if (verbose) {
      message("  m0: selected lambda = ", round(reg_m0, 5),
              " (", cv_m0$iterations, " iter, ",
              if (cv_m0$converged) "converged" else "max iter", ")")
    }

    # Floor lambda at theory minimum log(n)/n for Rashomon fitting.
    # CV may select a very small lambda (optimal for prediction), but Rashomon
    # intersection needs a minimum lambda to keep Rashomon sets tractable.
    # Theory minimum is the smallest lambda where any regularization is applied.
    # max_depth=4L bounds the Rashomon set size from above; this floor bounds it below.
    # Using log(n)/n (not sqrt(log(n)/n)) avoids floor-cap inversion at n>=1680.
    theory_floor_e  <- log(n)  / n
    theory_floor_m0 <- log(n0) / n0
    reg_e  <- max(reg_e,  theory_floor_e)
    reg_m0 <- max(reg_m0, theory_floor_m0)
    if (verbose) {
      if (reg_e  > cv_e$best_lambda)  message("  e:  lambda floored to theory scale: ", round(reg_e,  5))
      if (reg_m0 > cv_m0$best_lambda) message("  m0: lambda floored to theory scale: ", round(reg_m0, 5))
    }
  } else {
    # Use fixed regularization
    reg_e <- regularization
    reg_m0 <- regularization
    idx0 <- which(A == 0)
    X0 <- X[idx0, , drop = FALSE]
    Y0 <- Y[idx0]
    n0 <- length(idx0)
  }

  # Rashomon tolerance policy. Two distinct behaviors, selected by the caller:
  #
  #  (1) FIXED tolerance (default; escalate_intersection = FALSE). Use the single
  #      tolerance epsilon_n = rashomon_bound_multiplier * log(n)/n -- the theory
  #      value log(n)/n when the multiplier was NULL, or the caller's explicit value.
  #      This is o(n^{-1/2}) for any fixed multiplier and carries the valid-inference
  #      guarantee (manuscript Cor. "Rashomon tolerance without the intersection
  #      trade-off"). If the intersection is empty at this tolerance, we do NOT widen
  #      it; the caller falls back to fold-specific trees (still valid). The c-grid is
  #      the single point {multiplier}, so escalate_rashomon_intersection does one fit.
  #
  #  (2) ESCALATED tolerance (escalate_intersection = TRUE, and NO explicit multiplier).
  #      Keep lambda fixed and widen epsilon_n = c * log(n)/n over .RASHOMON_C_GRID
  #      (c = 1, 2, ..., 1000) until the cross-fold intersection is non-empty, recording
  #      the selected c. An empty theory-tolerance intersection is a finite-sample
  #      artifact (margin theorem), and escalation forces a single shared tree instead
  #      of abandoning it. IMPORTANT: for c > 1 the tolerance is NO LONGER guaranteed
  #      o(n^{-1/2}), so this TRADES the fixed-epsilon_n validity guarantee for a
  #      non-empty intersection. It is an opt-in practical/study device whose finite-
  #      sample coverage must be validated empirically (see the escalation coverage
  #      sweep in simulations/six-approach-arbitration), NOT a silent default.
  #
  # An explicit multiplier always pins a single fixed tolerance (behavior 1) even when
  # escalate_intersection = TRUE, so the documented single-tolerance API is preserved.
  use_escalation <- isTRUE(escalate_intersection) && !explicit_multiplier
  n_full <- nrow(X)
  c_grid_e <- if (use_escalation) .RASHOMON_C_GRID else
    rashomon_bound_multiplier / optimaltrees::select_epsilon_n(n_full)

  esc_e <- escalate_rashomon_intersection(
    X, A, K, "log_loss", reg_e, fold_indices, max_leaves,
    verbose, "propensity", max_depth = max_depth, c_grid = c_grid_e, ...
  )
  cf_e <- esc_e$fit
  rashomon_c_e <- esc_e$c

  idx0 <- which(A == 0)
  n0 <- length(idx0)
  if (n0 == 0) stop("No control units (A=0); cannot fit m0.")

  X0 <- X[idx0, , drop = FALSE]
  Y0 <- Y[idx0]
  fold0 <- fold_indices[idx0]

  c_grid_m0 <- if (use_escalation) .RASHOMON_C_GRID else
    rashomon_bound_multiplier / optimaltrees::select_epsilon_n(n0)

  esc_m0 <- escalate_rashomon_intersection(
    X0, Y0, K, loss_outcome, reg_m0, fold0, max_leaves,
    verbose, "control outcome", max_depth = max_depth, c_grid = c_grid_m0, ...
  )
  cf_m0 <- esc_m0$fit
  rashomon_c_m0 <- esc_m0$c

  # Efficient fallback: only refit failed models
  fallback_fits <- NULL
  need_e_fallback <- is.null(cf_e)
  need_m0_fallback <- is.null(cf_m0)

  if (need_e_fallback || need_m0_fallback) {
    if (verbose) {
      failed_models <- c(if (need_e_fallback) "propensity" else NULL,
                        if (need_m0_fallback) "control outcome" else NULL)
      message("Rashomon intersection failed for: ", paste(failed_models, collapse = ", "),
              ". Using fold-specific trees for all nuisances.")
    }
    fallback_fits <- vector("list", K)
    for (k in seq_len(K)) {
      fallback_fits[[k]] <- fit_nuisances_fold(X, A, Y, fold_id = k, fold_indices = fold_indices,
                                              outcome_type = outcome_type, regularization = regularization,
                                              cv_regularization = cv_regularization, cv_K = cv_K,
                                              verbose = verbose, ...)
    }
  }

  list(
    cf_e = cf_e,
    cf_m0 = cf_m0,
    fallback_fits = fallback_fits,
    outcome_type = outcome_type,
    # Rashomon tolerance multipliers selected by escalation (epsilon_n = c*log(n)/n).
    # NA if the intersection was empty even at c_max (fell back to fold-specific).
    rashomon_c_e = rashomon_c_e,
    rashomon_c_m0 = rashomon_c_m0
  )
}

#' Fold-specific eta from Rashomon fits (and optional fallback)
#'
#' Builds e(X), m0(X) of length n using predict(..., fold_indices) on
#' cf_rashomon objects, or fallback per-fold fits when a nuisance had n_intersecting == 0.
#' Propensity scores are clamped to [e_min, e_max] at prediction time.
#'
#' @param rashomon_list Return value of fit_nuisances_rashomon.
#' @param X Data.frame or matrix of covariates.
#' @param fold_indices Integer vector of length n (fold id per row).
#' @param e_min Lower bound for propensity clamping.
#' @param e_max Upper bound for propensity clamping.
#' @return List with elements e, m0 (each numeric vector of length n).
#' @noRd
get_fold_specific_eta_rashomon <- function(rashomon_list, X, fold_indices,
                                           e_min = .PROPENSITY_LOWER_BOUND,
                                           e_max = .PROPENSITY_UPPER_BOUND) {
  n <- nrow(X)
  e_out <- numeric(n)
  m0_out <- numeric(n)
  outcome_type <- rashomon_list$outcome_type

  if (!is.null(rashomon_list$cf_e) && rashomon_list$cf_e@n_intersecting > 0) {
    pe <- predict(rashomon_list$cf_e, X, fold_indices = fold_indices, type = "prob")
    if (!is.matrix(pe) || ncol(pe) < 2L) {
      stop(
        "predict() for propensity returned unexpected type or shape.\n",
        "  Expected: matrix with 2 columns (class probabilities)\n",
        "  Got: class = ", paste(class(pe), collapse = "/"),
        ", ncol = ", if (is.matrix(pe)) ncol(pe) else "N/A",
        call. = FALSE
      )
    }
    e_out <- pe[, 2L]
  } else if (!is.null(rashomon_list$fallback_fits)) {
    eta_fb <- get_fold_specific_eta(rashomon_list$fallback_fits, X, fold_indices, e_min, e_max)
    e_out <- eta_fb$e
  }

  if (!is.null(rashomon_list$cf_m0) && rashomon_list$cf_m0@n_intersecting > 0) {
    if (outcome_type == "continuous") {
      m0_out <- as.numeric(predict(rashomon_list$cf_m0, X, fold_indices = fold_indices))
    } else {
      pm0 <- predict(rashomon_list$cf_m0, X, fold_indices = fold_indices, type = "prob")
      if (!is.matrix(pm0) || ncol(pm0) < 2L) {
        stop(
          "predict() for outcome returned unexpected type or shape.\n",
          "  Expected: matrix with 2 columns (class probabilities)\n",
          "  Got: class = ", paste(class(pm0), collapse = "/"),
          ", ncol = ", if (is.matrix(pm0)) ncol(pm0) else "N/A",
          call. = FALSE
        )
      }
      m0_out <- pm0[, 2L]
    }
  } else if (!is.null(rashomon_list$fallback_fits)) {
    eta_fb <- get_fold_specific_eta(rashomon_list$fallback_fits, X, fold_indices, e_min, e_max)
    m0_out <- eta_fb$m0
  }

  # Clamp propensities at prediction time for numerical stability
  e_clamped <- pmax(e_min, pmin(e_max, e_out))

  list(e = e_clamped, m0 = m0_out)
}
