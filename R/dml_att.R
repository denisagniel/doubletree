#' DML estimator for the Average Treatment effect on the Treated (ATT)
#'
#' Estimates the ATT using double machine learning with optimal decision trees
#' (via treefarmr) for the nuisance functions e(X), m0(X), m1(X). Uses log-loss
#' only; binary outcome Y required. Continuous outcome is not implemented yet.
#'
#' @param X Data.frame or matrix of covariates. Must be binary (0/1) for treefarmr.
#' @param A Integer or numeric vector of treatment (0/1).
#' @param Y Integer or numeric vector of outcome (0/1). Binary only; continuous Y not supported.
#' @param K Number of cross-fitting folds. Default 5.
#' @param regularization Numeric. Tree complexity penalty passed to treefarmr. Default 0.1.
#' @param stratified Logical. If TRUE (default), fold assignment is stratified by A.
#' @param seed Optional. Random seed for fold creation.
#' @param verbose Logical. Passed to treefarmr. Default FALSE.
#' @param ... Additional arguments passed to treefarmr::fit_tree.
#' @return List with elements: theta (point estimate), sigma (estimated SE), ci_95 (Wald 95% CI),
#'   score_values (influence at theta), nuisance_fits (per-fold models), fold_indices, n, K.
#' @references Manuscript equation (2) for the orthogonal score; Chernozhukov et al. for DML.
#' @export
dml_att <- function(X, A, Y, K = 5, regularization = 0.1, stratified = TRUE, seed = NULL, verbose = FALSE, ...) {
  check_dml_att_data(X, A, Y)
  if (is.matrix(X)) X <- as.data.frame(X)
  n <- nrow(X)
  fold_indices <- create_folds(n, K, strata = if (stratified) A else NULL, seed = seed)

  nuisance_fits <- vector("list", K)
  for (k in seq_len(K)) {
    nuisance_fits[[k]] <- fit_nuisances_fold(X, A, Y, fold_id = k, fold_indices = fold_indices, regularization = regularization, verbose = verbose, ...)
  }

  eta <- get_fold_specific_eta(nuisance_fits, X, fold_indices)
  pi_hat <- mean(A)

  # Closed form: psi(theta) = psi(0) - theta*(A/pi), so sum(psi(theta)) = 0 => theta = sum(psi(0)) / sum(A/pi).
  sum_a_over_pi <- sum(A / pi_hat)
  if (sum_a_over_pi < 1e-10) stop("No treated units (sum(A) ~ 0).")
  score_at_zero <- psi_att(Y, A, theta = 0, eta, pi_hat)
  theta <- sum(score_at_zero) / sum_a_over_pi

  score_values <- psi_att(Y, A, theta, eta, pi_hat)
  sigma_sq <- dml_att_variance(score_values, n)
  sigma <- sqrt(sigma_sq)
  ci_95 <- dml_att_ci(theta, sigma, n, level = 0.95)

  list(
    theta = theta,
    sigma = sigma,
    ci_95 = ci_95,
    score_values = score_values,
    nuisance_fits = nuisance_fits,
    fold_indices = fold_indices,
    n = n,
    K = K
  )
}
