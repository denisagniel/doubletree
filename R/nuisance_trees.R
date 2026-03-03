#' Fit nuisance trees for one cross-fitting fold
#'
#' Fits propensity e(X), and outcome regressions m0(X), m1(X) on training data
#' (all rows not in fold_id) using treefarmr. Propensity uses log_loss; outcome
#' trees use log_loss (binary Y) or squared_error (continuous Y). Used internally by dml_att.
#'
#' @param X Data.frame or matrix of covariates (binary 0/1).
#' @param A Integer vector of treatment (0/1).
#' @param Y Numeric vector of outcome (binary 0/1 or continuous).
#' @param fold_id Integer. Fold index (1..K) to leave out; training = fold_indices != fold_id.
#' @param fold_indices Integer vector of length nrow(X) with fold assignment.
#' @param outcome_type Character. "binary" or "continuous"; determines loss for m0, m1.
#' @param regularization Numeric. Penalty per leaf for treefarmr. Default 0.1.
#' @param cv_regularization Logical. If TRUE, use cross-validation to select lambda. Default FALSE.
#' @param cv_K Integer. Number of CV folds for lambda selection. Default 5.
#' @param verbose Logical. Passed to treefarmr. Default FALSE.
#' @param ... Additional arguments passed to treefarmr::fit_tree.
#' @return List with elements e_model, m0_model, m1_model (treefarmr fit objects).
#'   If no A=0 or no A=1 in training set, the corresponding model is NULL (caller should handle).
#' @noRd
fit_nuisances_fold <- function(X, A, Y, fold_id, fold_indices, outcome_type = "binary",
                               regularization = 0.1, cv_regularization = FALSE, cv_K = 5,
                               verbose = FALSE, ...) {
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

  # Fit propensity model
  if (cv_regularization) {
    cv_e <- treefarmr::cv_regularization(X_tr, A_tr, loss_function = "log_loss",
                                         K = cv_K, refit = TRUE, verbose = FALSE, ...)
    e_model <- cv_e$model
    if (verbose) message("  Fold ", fold_id, " e: selected lambda = ", round(cv_e$best_lambda, 5))
  } else {
    e_model <- treefarmr::fit_tree(X_tr, A_tr, loss_function = "log_loss",
                                   regularization = regularization, verbose = verbose, ...)
  }

  # Fit m0 model
  if (cv_regularization) {
    cv_m0 <- treefarmr::cv_regularization(X_tr[A_tr == 0, , drop = FALSE], Y_tr[A_tr == 0],
                                          loss_function = loss_outcome, K = cv_K,
                                          refit = TRUE, verbose = FALSE, ...)
    m0_model <- cv_m0$model
    if (verbose) message("  Fold ", fold_id, " m0: selected lambda = ", round(cv_m0$best_lambda, 5))
  } else {
    m0_model <- treefarmr::fit_tree(X_tr[A_tr == 0, , drop = FALSE], Y_tr[A_tr == 0],
                                    loss_function = loss_outcome, regularization = regularization,
                                    verbose = verbose, ...)
  }

  # Fit m1 model
  if (cv_regularization) {
    cv_m1 <- treefarmr::cv_regularization(X_tr[A_tr == 1, , drop = FALSE], Y_tr[A_tr == 1],
                                          loss_function = loss_outcome, K = cv_K,
                                          refit = TRUE, verbose = FALSE, ...)
    m1_model <- cv_m1$model
    if (verbose) message("  Fold ", fold_id, " m1: selected lambda = ", round(cv_m1$best_lambda, 5))
  } else {
    m1_model <- treefarmr::fit_tree(X_tr[A_tr == 1, , drop = FALSE], Y_tr[A_tr == 1],
                                    loss_function = loss_outcome, regularization = regularization,
                                    verbose = verbose, ...)
  }

  # Nuisance quality diagnostics
  # Check propensity model quality
  e_pred <- predict(e_model, X_tr, type = "prob")
  e_vals <- if (is.matrix(e_pred)) e_pred[, 2L] else rep(0.5, nrow(X_tr))

  # Check for degenerate predictions
  if (sd(e_vals) < 1e-6) {
    warning("Fold ", fold_id, ": Propensity model predictions nearly constant ",
            "(sd < 1e-6). Model may not have converged or data may lack signal. ",
            "Consider different regularization or checking data quality.",
            call. = FALSE)
  }

  # Check for extreme predictions
  prop_extreme <- mean(e_vals < 0.01 | e_vals > 0.99)
  if (prop_extreme > 0.1) {
    warning("Fold ", fold_id, ": ", round(100 * prop_extreme, 1),
            "% of propensity predictions are extreme (<0.01 or >0.99). ",
            "This may indicate separation or poor overlap. ",
            "Consider larger regularization or checking for rare covariate patterns.",
            call. = FALSE)
  }

  # Check for predictions outside [0,1] (shouldn't happen but defensive)
  if (any(e_vals < 0 | e_vals > 1)) {
    stop("Fold ", fold_id, ": Propensity predictions outside [0,1]. ",
         "This indicates a bug in treefarmr or corrupted model output.",
         call. = FALSE)
  }

  # Similar checks for outcome models if continuous
  if (outcome_type == "continuous") {
    m0_pred <- predict(m0_model, X_tr[A_tr == 0, , drop = FALSE])
    m1_pred <- predict(m1_model, X_tr[A_tr == 1, , drop = FALSE])

    # Check for constant predictions
    if (sd(m0_pred) < 1e-6) {
      warning("Fold ", fold_id, ": Control outcome model predictions nearly constant",
              call. = FALSE)
    }
    if (sd(m1_pred) < 1e-6) {
      warning("Fold ", fold_id, ": Treated outcome model predictions nearly constant",
              call. = FALSE)
    }
  }

  list(e_model = e_model, m0_model = m0_model, m1_model = m1_model, outcome_type = outcome_type)
}

#' Predict nuisances for given rows using a fold's fitted models
#'
#' Returns e(X), m0(X), m1(X) for the requested rows. Propensity e from
#' predict(..., type = "prob"); m0, m1 from type = "prob" (binary) or default predict (continuous).
#'
#' @param models List from fit_nuisances_fold (e_model, m0_model, m1_model, outcome_type).
#' @param X Full covariate matrix (we subset to fold_rows).
#' @param fold_rows Integer vector of row indices to predict for.
#' @return List with elements e, m0, m1 (numeric vectors of length length(fold_rows)).
#' @noRd
predict_nuisances_fold <- function(models, X, fold_rows) {
  if (length(fold_rows) == 0) {
    return(list(e = numeric(0), m0 = numeric(0), m1 = numeric(0)))
  }
  X_sub <- X[fold_rows, , drop = FALSE]
  outcome_type <- if (!is.null(models$outcome_type)) models$outcome_type else "binary"
  pe <- predict(models$e_model, X_sub, type = "prob")
  e_vec <- if (is.matrix(pe)) pe[, 2L] else rep(0.5, nrow(X_sub))
  if (outcome_type == "continuous") {
    pm0 <- predict(models$m0_model, X_sub)
    pm1 <- predict(models$m1_model, X_sub)
    m0_vec <- as.numeric(if (is.matrix(pm0)) pm0[, 1L] else pm0)
    m1_vec <- as.numeric(if (is.matrix(pm1)) pm1[, 1L] else pm1)
  } else {
    pm0 <- predict(models$m0_model, X_sub, type = "prob")
    pm1 <- predict(models$m1_model, X_sub, type = "prob")
    m0_vec <- if (is.matrix(pm0)) pm0[, 2L] else rep(0.5, nrow(X_sub))
    m1_vec <- if (is.matrix(pm1)) pm1[, 2L] else rep(0.5, nrow(X_sub))
  }
  list(e = e_vec, m0 = m0_vec, m1 = m1_vec)
}

#' Build fold-specific nuisance predictions for all observations
#'
#' For each observation i, get eta^{(-k(i))}(X_i) using the models fitted for fold k(i).
#'
#' @param nuisance_fits List of length K; each element is list(e_model, m0_model, m1_model).
#' @param X Data.frame or matrix of covariates.
#' @param fold_indices Integer vector of length n (fold id per row).
#' @return List with elements e, m0, m1 (each numeric vector of length n).
#' @noRd
get_fold_specific_eta <- function(nuisance_fits, X, fold_indices) {
  n <- nrow(X)
  K <- length(nuisance_fits)
  e_out <- numeric(n)
  m0_out <- numeric(n)
  m1_out <- numeric(n)
  for (k in seq_len(K)) {
    idx_k <- which(fold_indices == k)
    if (length(idx_k) == 0) next
    pred <- predict_nuisances_fold(nuisance_fits[[k]], X, idx_k)
    e_out[idx_k] <- pred$e
    m0_out[idx_k] <- pred$m0
    m1_out[idx_k] <- pred$m1
  }
  list(e = e_out, m0 = m0_out, m1 = m1_out)
}

#' Fit nuisances via Rashomon intersection (cross_fitted_rashomon)
#'
#' Fits e(X), m0(X), m1(X) using treefarmr::cross_fitted_rashomon with the same
#' fold assignment as DML. When a nuisance has no intersecting trees (n_intersecting == 0),
#' falls back to single-tree-per-fold for that nuisance via fit_nuisances_fold.
#'
#' @param X Data.frame or matrix of covariates (binary 0/1).
#' @param A Integer vector of treatment (0/1).
#' @param Y Numeric vector of outcome (binary 0/1 or continuous).
#' @param fold_indices Integer vector of length nrow(X) with fold id per row (1..K).
#' @param outcome_type Character. "binary" or "continuous"; determines loss for m0, m1.
#' @param regularization Numeric. Passed to treefarmr. Default 0.1.
#' @param cv_regularization Logical. If TRUE, use CV to select lambda for each nuisance. Default FALSE.
#' @param cv_K Integer. Number of CV folds for lambda selection. Default 5.
#' @param verbose Logical. Passed to treefarmr. Default FALSE.
#' @param rashomon_bound_multiplier,rashomon_bound_adder,max_leaves,auto_tune_intersecting Passed to cross_fitted_rashomon.
#' @param ... Additional arguments passed to treefarmr::cross_fitted_rashomon.
#' @return List with cf_e, cf_m0, cf_m1 (cf_rashomon or NULL if fallback), fallback_fits (list of K per-fold fits or NULL), outcome_type.
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
    cv_e <- treefarmr::cv_regularization(X, A, loss_function = "log_loss",
                                         K = cv_K, refit = FALSE, verbose = FALSE, ...)
    reg_e <- cv_e$best_lambda
    if (verbose) message("  e: selected lambda = ", round(reg_e, 5))

    # CV for m0
    idx0 <- which(A == 0)
    X0 <- X[idx0, , drop = FALSE]
    Y0 <- Y[idx0]
    cv_m0 <- treefarmr::cv_regularization(X0, Y0, loss_function = loss_outcome,
                                          K = cv_K, refit = FALSE, verbose = FALSE, ...)
    reg_m0 <- cv_m0$best_lambda
    if (verbose) message("  m0: selected lambda = ", round(reg_m0, 5))

    # CV for m1
    idx1 <- which(A == 1)
    X1 <- X[idx1, , drop = FALSE]
    Y1 <- Y[idx1]
    cv_m1 <- treefarmr::cv_regularization(X1, Y1, loss_function = loss_outcome,
                                          K = cv_K, refit = FALSE, verbose = FALSE, ...)
    reg_m1 <- cv_m1$best_lambda
    if (verbose) message("  m1: selected lambda = ", round(reg_m1, 5))
  } else {
    # Use fixed regularization
    reg_e <- regularization
    reg_m0 <- regularization
    reg_m1 <- regularization
  }

  cf_e <- tryCatch({
    out <- treefarmr::cross_fitted_rashomon(
      X, A, K = K, loss_function = "log_loss", regularization = reg_e,
      rashomon_bound_multiplier = rashomon_bound_multiplier, rashomon_bound_adder = rashomon_bound_adder,
      max_leaves = max_leaves, fold_indices = fold_indices,
      auto_tune_intersecting = auto_tune_intersecting, verbose = verbose, ...
    )
    if (out$n_intersecting > 0) out else NULL
  }, error = function(e) NULL)

  idx0 <- which(A == 0)
  idx1 <- which(A == 1)
  n0 <- length(idx0)
  n1 <- length(idx1)
  if (n0 == 0) stop("No control units (A=0); cannot fit m0.")
  if (n1 == 0) stop("No treated units (A=1); cannot fit m1.")

  cf_m0 <- tryCatch({
    X0 <- X[idx0, , drop = FALSE]
    Y0 <- Y[idx0]
    fold0 <- fold_indices[idx0]
    out <- treefarmr::cross_fitted_rashomon(
      X0, Y0, K = K, loss_function = loss_outcome, regularization = reg_m0,
      rashomon_bound_multiplier = rashomon_bound_multiplier, rashomon_bound_adder = rashomon_bound_adder,
      max_leaves = max_leaves, fold_indices = fold0,
      auto_tune_intersecting = auto_tune_intersecting, verbose = verbose, ...
    )
    if (out$n_intersecting > 0) out else NULL
  }, error = function(e) NULL)

  cf_m1 <- tryCatch({
    X1 <- X[idx1, , drop = FALSE]
    Y1 <- Y[idx1]
    fold1 <- fold_indices[idx1]
    out <- treefarmr::cross_fitted_rashomon(
      X1, Y1, K = K, loss_function = loss_outcome, regularization = reg_m1,
      rashomon_bound_multiplier = rashomon_bound_multiplier, rashomon_bound_adder = rashomon_bound_adder,
      max_leaves = max_leaves, fold_indices = fold1,
      auto_tune_intersecting = auto_tune_intersecting, verbose = verbose, ...
    )
    if (out$n_intersecting > 0) out else NULL
  }, error = function(e) NULL)

  need_fallback <- is.null(cf_e) || is.null(cf_m0) || is.null(cf_m1)
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
    cf_m1 = cf_m1,
    fallback_fits = fallback_fits,
    outcome_type = outcome_type
  )
}

#' Fold-specific eta from Rashomon fits (and optional fallback)
#'
#' Builds e(X), m0(X), m1(X) of length n using predict(..., fold_indices) on
#' cf_rashomon objects, or fallback per-fold fits when a nuisance had n_intersecting == 0.
#'
#' @param rashomon_list Return value of fit_nuisances_rashomon.
#' @param X Data.frame or matrix of covariates.
#' @param fold_indices Integer vector of length n (fold id per row).
#' @return List with elements e, m0, m1 (each numeric vector of length n).
#' @noRd
get_fold_specific_eta_rashomon <- function(rashomon_list, X, fold_indices) {
  n <- nrow(X)
  e_out <- numeric(n)
  m0_out <- numeric(n)
  m1_out <- numeric(n)
  outcome_type <- rashomon_list$outcome_type

  if (!is.null(rashomon_list$cf_e) && rashomon_list$cf_e$n_intersecting > 0) {
    pe <- predict(rashomon_list$cf_e, X, fold_indices = fold_indices, type = "prob")
    e_out <- if (is.matrix(pe)) pe[, 2L] else rep(0.5, n)
  } else if (!is.null(rashomon_list$fallback_fits)) {
    eta_fb <- get_fold_specific_eta(rashomon_list$fallback_fits, X, fold_indices)
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
    eta_fb <- get_fold_specific_eta(rashomon_list$fallback_fits, X, fold_indices)
    m0_out <- eta_fb$m0
  }

  if (!is.null(rashomon_list$cf_m1) && rashomon_list$cf_m1$n_intersecting > 0) {
    if (outcome_type == "continuous") {
      m1_out <- as.numeric(predict(rashomon_list$cf_m1, X, fold_indices = fold_indices))
    } else {
      pm1 <- predict(rashomon_list$cf_m1, X, fold_indices = fold_indices, type = "prob")
      m1_out <- if (is.matrix(pm1)) pm1[, 2L] else rep(0.5, n)
    }
  } else if (!is.null(rashomon_list$fallback_fits)) {
    eta_fb <- get_fold_specific_eta(rashomon_list$fallback_fits, X, fold_indices)
    m1_out <- eta_fb$m1
  }

  list(e = e_out, m0 = m0_out, m1 = m1_out)
}
