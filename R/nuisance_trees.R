#' Fit nuisance trees for one cross-fitting fold
#'
#' Fits propensity e(X), and outcome regressions m0(X), m1(X) on training data
#' (all rows not in fold_id) using treefarmr with log_loss. Used internally by dml_att.
#'
#' @param X Data.frame or matrix of covariates (binary 0/1).
#' @param A Integer vector of treatment (0/1).
#' @param Y Integer vector of outcome (binary 0/1).
#' @param fold_id Integer. Fold index (1..K) to leave out; training = fold_indices != fold_id.
#' @param fold_indices Integer vector of length nrow(X) with fold assignment.
#' @param regularization Numeric. Penalty per leaf for treefarmr. Default 0.1.
#' @param verbose Logical. Passed to treefarmr. Default FALSE.
#' @param ... Additional arguments passed to treefarmr::fit_tree.
#' @return List with elements e_model, m0_model, m1_model (treefarmr fit objects).
#'   If no A=0 or no A=1 in training set, the corresponding model is NULL (caller should handle).
#' @noRd
fit_nuisances_fold <- function(X, A, Y, fold_id, fold_indices, regularization = 0.1, verbose = FALSE, ...) {
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
  e_model <- treefarmr::fit_tree(X_tr, A_tr, loss_function = "log_loss", regularization = regularization, verbose = verbose, ...)
  m0_model <- treefarmr::fit_tree(X_tr[A_tr == 0, , drop = FALSE], Y_tr[A_tr == 0], loss_function = "log_loss", regularization = regularization, verbose = verbose, ...)
  m1_model <- treefarmr::fit_tree(X_tr[A_tr == 1, , drop = FALSE], Y_tr[A_tr == 1], loss_function = "log_loss", regularization = regularization, verbose = verbose, ...)
  list(e_model = e_model, m0_model = m0_model, m1_model = m1_model)
}

#' Predict nuisances for given rows using a fold's fitted models
#'
#' Returns e(X), m0(X), m1(X) for the requested rows using the three models from one fold.
#' Uses P(A=1) and P(Y=1) from treefarmr predict(..., type = "prob") (second column).
#'
#' @param models List from fit_nuisances_fold (e_model, m0_model, m1_model).
#' @param X Full covariate matrix (we subset to fold_rows).
#' @param fold_rows Integer vector of row indices to predict for.
#' @return List with elements e, m0, m1 (numeric vectors of length length(fold_rows)).
#' @noRd
predict_nuisances_fold <- function(models, X, fold_rows) {
  if (length(fold_rows) == 0) {
    return(list(e = numeric(0), m0 = numeric(0), m1 = numeric(0)))
  }
  X_sub <- X[fold_rows, , drop = FALSE]
  pe <- predict(models$e_model, X_sub, type = "prob")
  e_vec <- if (is.matrix(pe)) pe[, 2L] else rep(0.5, nrow(X_sub))
  pm0 <- predict(models$m0_model, X_sub, type = "prob")
  m0_vec <- if (is.matrix(pm0)) pm0[, 2L] else rep(0.5, nrow(X_sub))
  pm1 <- predict(models$m1_model, X_sub, type = "prob")
  m1_vec <- if (is.matrix(pm1)) pm1[, 2L] else rep(0.5, nrow(X_sub))
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
