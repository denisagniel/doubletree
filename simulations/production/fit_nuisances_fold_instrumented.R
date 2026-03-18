# Instrumented fit_nuisances_fold to find exact hang location

fit_nuisances_fold_instrumented <- function(X, A, Y, fold_id, fold_indices, outcome_type = "binary",
                               regularization = 0.1, cv_regularization = FALSE, cv_K = 5,
                               verbose = FALSE,
                               discretize_method = "quantiles", discretize_bins = "adaptive",
                               ...) {

  cat("[NUISANCE] Starting fold", fold_id, "\n"); flush.console()

  train_idx <- which(fold_indices != fold_id)
  X_tr <- X[train_idx, , drop = FALSE]
  A_tr <- A[train_idx]
  Y_tr <- Y[train_idx]

  cat("[NUISANCE] Train set: n=", length(train_idx), "\n"); flush.console()

  n0 <- sum(A_tr == 0)
  n1 <- sum(A_tr == 1)

  cat("[NUISANCE] n0=", n0, ", n1=", n1, "\n"); flush.console()

  if (n0 == 0) {
    stop("No control units (A=0) in training set for fold ", fold_id, "; use fewer folds or check data.")
  }
  if (n1 == 0) {
    stop("No treated units (A=1) in training set for fold ", fold_id, "; use fewer folds or check data.")
  }

  loss_outcome <- if (outcome_type == "continuous") "squared_error" else "log_loss"
  cat("[NUISANCE] loss_outcome=", loss_outcome, "\n"); flush.console()

  # Fit propensity model
  cat("[NUISANCE] Fitting propensity model (e_model)...\n"); flush.console()

  if (cv_regularization) {
    cat("[NUISANCE]   Using cv_regularization...\n"); flush.console()
    cv_e <- optimaltrees::cv_regularization(
      X_tr, A_tr,
      loss_function = "log_loss",
      discretize_method = discretize_method,
      discretize_bins = discretize_bins,
      K = cv_K, refit = TRUE, verbose = FALSE, ...
    )
    e_model <- cv_e$model
    cat("[NUISANCE]   cv_regularization complete\n"); flush.console()
  } else {
    cat("[NUISANCE]   Using fit_tree with regularization=", regularization, "...\n"); flush.console()
    e_model <- optimaltrees::fit_tree(
      X_tr, A_tr,
      loss_function = "log_loss",
      discretize_method = discretize_method,
      discretize_bins = discretize_bins,
      regularization = regularization, verbose = verbose, ...
    )
    cat("[NUISANCE]   fit_tree complete\n"); flush.console()
  }

  cat("[NUISANCE] Propensity model fitted, n_trees=", e_model$n_trees, "\n"); flush.console()

  # Check for model fitting failure
  if (is.null(e_model$n_trees) || e_model$n_trees == 0) {
    stop("TreeFARMS model fitting failed for propensity model in fold ", fold_id, call. = FALSE)
  }

  # Fit m0 model on controls
  cat("[NUISANCE] Fitting outcome|control model (m0_model)...\n"); flush.console()
  X_tr_controls <- X_tr[A_tr == 0, , drop = FALSE]
  cat("[NUISANCE]   Control subset: n=", nrow(X_tr_controls), "\n"); flush.console()

  if (cv_regularization) {
    cat("[NUISANCE]   Using cv_regularization...\n"); flush.console()
    cv_m0 <- optimaltrees::cv_regularization(
      X_tr_controls, Y_tr[A_tr == 0],
      loss_function = loss_outcome,
      discretize_method = discretize_method,
      discretize_bins = discretize_bins,
      K = cv_K, refit = TRUE, verbose = FALSE, ...
    )
    m0_model <- cv_m0$model
    cat("[NUISANCE]   cv_regularization complete\n"); flush.console()
  } else {
    cat("[NUISANCE]   Using fit_tree with regularization=", regularization, "...\n"); flush.console()
    m0_model <- optimaltrees::fit_tree(
      X_tr_controls, Y_tr[A_tr == 0],
      loss_function = loss_outcome,
      discretize_method = discretize_method,
      discretize_bins = discretize_bins,
      regularization = regularization, verbose = verbose, ...
    )
    cat("[NUISANCE]   fit_tree complete\n"); flush.console()
  }

  cat("[NUISANCE] Outcome|control model fitted, n_trees=", m0_model$n_trees, "\n"); flush.console()

  # Check for model fitting failure
  if (is.null(m0_model$n_trees) || m0_model$n_trees == 0) {
    stop("TreeFARMS model fitting failed for control outcome model (m0) in fold ", fold_id, call. = FALSE)
  }

  cat("[NUISANCE] Fold", fold_id, "complete\n"); flush.console()

  list(
    e_model = e_model,
    m0_model = m0_model,
    outcome_type = outcome_type
  )
}
