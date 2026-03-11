#' Fit nuisance trees for one cross-fitting fold
#'
#' Fits propensity e(X) and outcome regression m0(X) on training data
#' (all rows not in fold_id) using treefarmr. Propensity uses log_loss; outcome
#' trees use log_loss (binary Y) or squared_error (continuous Y). Used internally by dml_att.
#'
#' @param X Data.frame or matrix of covariates (binary 0/1).
#' @param A Integer vector of treatment (0/1).
#' @param Y Numeric vector of outcome (binary 0/1 or continuous).
#' @param fold_id Integer. Fold index (1..K) to leave out; training = fold_indices != fold_id.
#' @param fold_indices Integer vector of length nrow(X) with fold assignment.
#' @param outcome_type Character. "binary" or "continuous"; determines loss for m0.
#' @param regularization Numeric. Penalty per leaf for treefarmr. Default 0.1.
#' @param cv_regularization Logical. If TRUE, use cross-validation to select lambda. Default FALSE.
#' @param cv_K Integer. Number of CV folds for lambda selection. Default 5.
#' @param verbose Logical. Passed to treefarmr. Default FALSE.
#' @param discretize_method Character. "quantiles" for quantile-based discretization.
#' @param discretize_bins Integer or "adaptive". Number of bins for discretization.
#' @param ... Additional arguments passed to optimaltrees::fit_tree.
#' @return List with elements e_model, m0_model (treefarmr fit objects),
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

  # Fit propensity model (treefarmr handles discretization with threshold encoding)
  if (cv_regularization) {
    cv_e <- optimaltrees::cv_regularization(
      X_tr, A_tr,  # Pass continuous features directly
      loss_function = "log_loss",
      discretize_method = discretize_method,
      discretize_bins = discretize_bins,
      K = cv_K, refit = TRUE, verbose = FALSE, ...
    )
    e_model <- cv_e$model
    if (verbose) message("  Fold ", fold_id, " e: selected lambda = ", round(cv_e$best_lambda, 5))
  } else {
    e_model <- optimaltrees::fit_tree(
      X_tr, A_tr,  # Pass continuous features directly
      loss_function = "log_loss",
      discretize_method = discretize_method,
      discretize_bins = discretize_bins,
      regularization = regularization, verbose = verbose, ...
    )
  }

  # Check for model fitting failure (should not happen with model_limit=0)
  if (is.null(e_model$model$tree_json) && is.null(e_model$model$result_data)) {
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

  # Fit m0 model on controls (treefarmr handles discretization with threshold encoding)
  X_tr_controls <- X_tr[A_tr == 0, , drop = FALSE]
  if (cv_regularization) {
    cv_m0 <- optimaltrees::cv_regularization(
      X_tr_controls, Y_tr[A_tr == 0],  # Pass continuous features directly
      loss_function = loss_outcome,
      discretize_method = discretize_method,
      discretize_bins = discretize_bins,
      K = cv_K, refit = TRUE, verbose = FALSE, ...
    )
    m0_model <- cv_m0$model
    if (verbose) message("  Fold ", fold_id, " m0: selected lambda = ", round(cv_m0$best_lambda, 5))
  } else {
    m0_model <- optimaltrees::fit_tree(
      X_tr_controls, Y_tr[A_tr == 0],  # Pass continuous features directly
      loss_function = loss_outcome,
      discretize_method = discretize_method,
      discretize_bins = discretize_bins,
      regularization = regularization, verbose = verbose, ...
    )
  }

  # Check for model fitting failure (should not happen with model_limit=0)
  if (is.null(m0_model$model$tree_json) && is.null(m0_model$model$result_data)) {
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

  # NOTE: Nuisance quality diagnostics temporarily disabled
  # TODO: Re-enable after predict() properly handles continuous→binary discretization
  # The predict function currently expects binary features matching training data,
  # but with automatic discretization we pass continuous features.
  # Need to either:
  # 1. Have predict() auto-discretize using stored metadata, OR
  # 2. Store discretized X_train in model object, OR
  # 3. Implement diagnostics using fold-out predictions instead of in-sample

  list(
    e_model = e_model,
    m0_model = m0_model,
    outcome_type = outcome_type
    # Note: discretization metadata now stored in model$discretization_metadata (treefarmr)
  )
}

#' Apply discretization metadata to new data
#'
#' Helper to discretize test data using training discretization metadata.
#' @param X_new Data.frame of continuous features
#' @param metadata Discretization metadata from treefarmr model
#' @return Data.frame with binary features
#' @noRd
apply_discretization_metadata <- function(X_new, metadata) {
  if (is.null(metadata) || is.null(metadata$features)) {
    return(X_new)
  }

  binary_cols_list <- list()

  for (col_name in names(X_new)) {
    x <- X_new[[col_name]]
    feat_meta <- metadata$features[[col_name]]

    if (is.null(feat_meta)) {
      binary_cols_list[[col_name]] <- data.frame(x)
      names(binary_cols_list[[col_name]]) <- col_name
      next
    }

    if (feat_meta$type == "binary") {
      binary_cols_list[[col_name]] <- data.frame(x)
      names(binary_cols_list[[col_name]]) <- col_name
    } else if (feat_meta$type == "binary_converted") {
      x_binary <- as.numeric(x == max(feat_meta$original_values))
      binary_cols_list[[col_name]] <- data.frame(x_binary)
      names(binary_cols_list[[col_name]]) <- col_name
    } else if (feat_meta$type == "constant") {
      x_binary <- rep(0, length(x))
      binary_cols_list[[col_name]] <- data.frame(x_binary)
      names(binary_cols_list[[col_name]]) <- feat_meta$new_names
    } else if (feat_meta$type == "continuous") {
      thresholds <- feat_meta$thresholds
      new_names <- feat_meta$new_names
      threshold_cols <- lapply(seq_along(thresholds), function(i) {
        as.integer(x <= thresholds[i])
      })
      binary_cols <- as.data.frame(threshold_cols)
      names(binary_cols) <- new_names
      binary_cols_list[[col_name]] <- binary_cols
    }
  }

  do.call(cbind, binary_cols_list)
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

  # Apply discretization to test data using training metadata
  # treefarmr stores discretization metadata in model$discretization
  e_metadata <- models$e_model$discretization
  m0_metadata <- models$m0_model$discretization

  X_sub_e <- if (!is.null(e_metadata)) {
    apply_discretization_metadata(X_sub, e_metadata)
  } else {
    X_sub  # Already binary
  }

  X_sub_m0 <- if (!is.null(m0_metadata)) {
    apply_discretization_metadata(X_sub, m0_metadata)
  } else {
    X_sub  # Already binary
  }

  # Predict on discretized data
  pe <- predict(models$e_model, X_sub_e, type = "prob")
  e_vec <- if (is.matrix(pe)) pe[, 2L] else rep(0.5, nrow(X_sub_e))

  if (outcome_type == "continuous") {
    pm0 <- predict(models$m0_model, X_sub_m0)
    m0_vec <- as.numeric(if (is.matrix(pm0)) pm0[, 1L] else pm0)
  } else {
    pm0 <- predict(models$m0_model, X_sub_m0, type = "prob")
    m0_vec <- if (is.matrix(pm0)) pm0[, 2L] else rep(0.5, nrow(X_sub_m0))
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

#' Fit nuisances via Rashomon intersection (cross_fitted_rashomon)
#'
#' Fits e(X), m0(X) using optimaltrees::cross_fitted_rashomon with the same
#' fold assignment as DML. When a nuisance has no intersecting trees (n_intersecting == 0),
#' falls back to single-tree-per-fold for that nuisance via fit_nuisances_fold.
#'
#' @param X Data.frame or matrix of covariates (binary 0/1).
#' @param A Integer vector of treatment (0/1).
#' @param Y Numeric vector of outcome (binary 0/1 or continuous).
#' @param fold_indices Integer vector of length nrow(X) with fold id per row (1..K).
#' @param outcome_type Character. "binary" or "continuous"; determines loss for m0.
#' @param regularization Numeric. Passed to treefarmr. Default 0.1.
#' @param cv_regularization Logical. If TRUE, use CV to select lambda for each nuisance. Default FALSE.
#' @param cv_K Integer. Number of CV folds for lambda selection. Default 5.
#' @param verbose Logical. Passed to treefarmr. Default FALSE.
#' @param rashomon_bound_multiplier,rashomon_bound_adder,max_leaves,auto_tune_intersecting Passed to cross_fitted_rashomon.
#' @param ... Additional arguments passed to optimaltrees::cross_fitted_rashomon.
#' @return List with cf_e, cf_m0 (cf_rashomon or NULL if fallback), fallback_fits (list of K per-fold fits or NULL), outcome_type.
#' @noRd
fit_nuisances_rashomon <- function(X, A, Y, fold_indices, outcome_type = "binary",
                                   regularization = 0.1, cv_regularization = FALSE, cv_K = 5,
                                   verbose = FALSE,
                                   rashomon_bound_multiplier = 0.05, rashomon_bound_adder = 0,
                                   max_leaves = NULL, auto_tune_intersecting = FALSE, ...) {
  n <- nrow(X)
  K <- max(fold_indices, na.rm = TRUE)
  loss_outcome <- if (outcome_type == "continuous") "squared_error" else "log_loss"

  # If CV requested, select regularization for each nuisance
  if (cv_regularization) {
    if (verbose) message("Selecting regularization via CV...")

    # CV for propensity
    cv_e <- optimaltrees::cv_regularization(X, A, loss_function = "log_loss",
                                         K = cv_K, refit = FALSE, verbose = FALSE, ...)
    reg_e <- cv_e$best_lambda
    if (verbose) message("  e: selected lambda = ", round(reg_e, 5))

    # CV for m0
    idx0 <- which(A == 0)
    X0 <- X[idx0, , drop = FALSE]
    Y0 <- Y[idx0]
    cv_m0 <- optimaltrees::cv_regularization(X0, Y0, loss_function = loss_outcome,
                                          K = cv_K, refit = FALSE, verbose = FALSE, ...)
    reg_m0 <- cv_m0$best_lambda
    if (verbose) message("  m0: selected lambda = ", round(reg_m0, 5))
  } else {
    # Use fixed regularization
    reg_e <- regularization
    reg_m0 <- regularization
  }

  cf_e <- tryCatch({
    out <- optimaltrees::cross_fitted_rashomon(
      X, A, K = K, loss_function = "log_loss", regularization = reg_e,
      rashomon_bound_multiplier = rashomon_bound_multiplier, rashomon_bound_adder = rashomon_bound_adder,
      max_leaves = max_leaves, fold_indices = fold_indices,
      auto_tune_intersecting = auto_tune_intersecting, verbose = verbose, ...
    )
    if (out$n_intersecting > 0) out else NULL
  }, error = function(e) {
    # Distinguish error types for better diagnostics
    msg <- conditionMessage(e)

    # Data/parameter errors: should stop immediately
    if (grepl("must be|invalid|cannot", msg, ignore.case = TRUE)) {
      stop("Error fitting propensity model: ", msg,
           "\nCheck your data and parameters.", call. = FALSE)
    }

    # Intersection empty: acceptable, use fallback
    if (grepl("intersection.*empty", msg, ignore.case = TRUE)) {
      if (verbose) {
        message("No intersecting trees found for propensity model. ",
                "Falling back to fold-specific trees.")
      }
      return(NULL)  # Fallback will handle
    }

    # Optimization/convergence issues: warn and fallback
    if (verbose) {
      warning("Rashomon fitting failed for propensity: ", msg,
              "\nFalling back to fold-specific trees.",
              call. = FALSE, immediate. = TRUE)
    }
    return(NULL)
  })

  idx0 <- which(A == 0)
  n0 <- length(idx0)
  if (n0 == 0) stop("No control units (A=0); cannot fit m0.")

  cf_m0 <- tryCatch({
    X0 <- X[idx0, , drop = FALSE]
    Y0 <- Y[idx0]
    fold0 <- fold_indices[idx0]
    out <- optimaltrees::cross_fitted_rashomon(
      X0, Y0, K = K, loss_function = loss_outcome, regularization = reg_m0,
      rashomon_bound_multiplier = rashomon_bound_multiplier, rashomon_bound_adder = rashomon_bound_adder,
      max_leaves = max_leaves, fold_indices = fold0,
      auto_tune_intersecting = auto_tune_intersecting, verbose = verbose, ...
    )
    if (out$n_intersecting > 0) out else NULL
  }, error = function(e) {
    # Distinguish error types for better diagnostics
    msg <- conditionMessage(e)

    # Data/parameter errors: should stop immediately
    if (grepl("must be|invalid|cannot", msg, ignore.case = TRUE)) {
      stop("Error fitting control outcome model (m0): ", msg,
           "\nCheck your data and parameters.", call. = FALSE)
    }

    # Intersection empty: acceptable, use fallback
    if (grepl("intersection.*empty", msg, ignore.case = TRUE)) {
      if (verbose) {
        message("No intersecting trees found for control outcome model. ",
                "Falling back to fold-specific trees.")
      }
      return(NULL)  # Fallback will handle
    }

    # Optimization/convergence issues: warn and fallback
    if (verbose) {
      warning("Rashomon fitting failed for control outcome model: ", msg,
              "\nFalling back to fold-specific trees.",
              call. = FALSE, immediate. = TRUE)
    }
    return(NULL)
  })

  need_fallback <- is.null(cf_e) || is.null(cf_m0)
  fallback_fits <- NULL
  if (need_fallback) {
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
    outcome_type = outcome_type
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

  if (!is.null(rashomon_list$cf_e) && rashomon_list$cf_e$n_intersecting > 0) {
    pe <- predict(rashomon_list$cf_e, X, fold_indices = fold_indices, type = "prob")
    e_out <- if (is.matrix(pe)) pe[, 2L] else rep(0.5, n)
  } else if (!is.null(rashomon_list$fallback_fits)) {
    eta_fb <- get_fold_specific_eta(rashomon_list$fallback_fits, X, fold_indices, e_min, e_max)
    e_out <- eta_fb$e
  }

  if (!is.null(rashomon_list$cf_m0) && rashomon_list$cf_m0$n_intersecting > 0) {
    if (outcome_type == "continuous") {
      m0_out <- as.numeric(predict(rashomon_list$cf_m0, X, fold_indices = fold_indices))
    } else {
      pm0 <- predict(rashomon_list$cf_m0, X, fold_indices = fold_indices, type = "prob")
      m0_out <- if (is.matrix(pm0)) pm0[, 2L] else rep(0.5, n)
    }
  } else if (!is.null(rashomon_list$fallback_fits)) {
    eta_fb <- get_fold_specific_eta(rashomon_list$fallback_fits, X, fold_indices, e_min, e_max)
    m0_out <- eta_fb$m0
  }

  # Clamp propensities at prediction time for numerical stability
  e_clamped <- pmax(e_min, pmin(e_max, e_out))

  list(e = e_clamped, m0 = m0_out)
}
