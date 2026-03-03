#' DML estimator for the Average Treatment effect on the Treated (ATT)
#'
#' Estimates the ATT using double machine learning with optimal decision trees
#' (via treefarmr) for the nuisance functions e(X), m0(X), m1(X). Binary outcome
#' uses log-loss for all nuisances; continuous outcome uses log-loss for propensity
#' and squared_error for m0, m1 (requires treefarmr to support squared_error).
#'
#' When \code{use_rashomon = TRUE}, nuisances are fit via
#' \code{treefarmr::cross_fitted_rashomon}: one interpretable tree per nuisance
#' (intersection of Rashomon sets across folds) with fold-specific refits for
#' valid DML. The same K and fold assignment are used for Rashomon and the score.
#'
#' @param X Data.frame or matrix of covariates. Must be binary (0/1) for treefarmr.
#' @param A Integer or numeric vector of treatment (0/1).
#' @param Y Numeric vector of outcome. Binary (0/1) when outcome_type is "binary"; any numeric when "continuous".
#' @param K Number of cross-fitting folds. Default 5.
#' @param outcome_type Character. "binary" (default) or "continuous". Continuous requires treefarmr squared_error loss for m0, m1.
#' @param regularization Numeric. Tree complexity penalty passed to treefarmr. Default 0.1. Ignored if \code{cv_regularization = TRUE}.
#' @param cv_regularization Logical. If TRUE, use cross-validation to select regularization parameter for each nuisance function. Default FALSE.
#' @param cv_K Integer. Number of folds for cross-validation of regularization. Default 5. Only used if \code{cv_regularization = TRUE}.
#' @param stratified Logical. If TRUE (default), fold assignment is stratified by A.
#' @param seed Optional. Random seed for fold creation.
#' @param verbose Logical. Passed to treefarmr. Default FALSE.
#' @param use_rashomon Logical. If TRUE, fit nuisances via \code{treefarmr::cross_fitted_rashomon} (one interpretable tree per nuisance via intersection + refit per fold). Default FALSE (single tree per fold).
#' @param rashomon_bound_multiplier Numeric. Rashomon tolerance \eqn{\varepsilon_n} passed to
#'   \code{cross_fitted_rashomon} when \code{use_rashomon = TRUE}. Default: 0.05.
#'   For theoretically justified values, use \code{treefarmr::select_epsilon_n(nrow(X))}.
#'   Recommended: \eqn{\varepsilon_n = c\sqrt{(\log n)/n}} for \eqn{c \in \{1,2,3\}}.
#'   See manuscript Appendix A.5.
#' @param rashomon_bound_adder Numeric. Additive Rashomon bound (not recommended for DML).
#'   Default: 0.
#' @param max_leaves Optional integer. Passed to \code{cross_fitted_rashomon} when \code{use_rashomon = TRUE}. Restricts Rashomon set to trees with at most this many leaves.
#' @param auto_tune_intersecting Logical. Passed to \code{cross_fitted_rashomon} when \code{use_rashomon = TRUE}. If TRUE, tune until at least one intersecting tree is found. Default FALSE.
#' @param ... Additional arguments passed to treefarmr (\code{fit_tree} when \code{use_rashomon = FALSE}, \code{cross_fitted_rashomon} when \code{use_rashomon = TRUE}).
#' @return List with elements: theta (point estimate), sigma (estimated SE), ci_95 (Wald 95% CI),
#'   score_values (influence at theta), nuisance_fits (per-fold models or Rashomon list), fold_indices, n, K.
#' @references Manuscript equation (2) for the orthogonal score; Chernozhukov et al. for DML.
#' @examples
#' \dontrun{
#' # Binary outcome (default)
#' library(treefarmr)  # Required dependency
#' set.seed(42)
#' n <- 300
#' X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5))
#' A <- rbinom(n, 1, plogis(0.5 * X$X1 - 0.2))
#' Y <- rbinom(n, 1, 0.3 + 0.2 * X$X1 + 0.15 * A)
#'
#' # Estimate ATT using DML with tree-based nuisances
#' fit <- dml_att(X, A, Y, K = 5, regularization = 0.1)
#' print(fit$theta)   # Point estimate
#' print(fit$ci_95)   # 95% Wald confidence interval
#'
#' # With Rashomon DML (one interpretable tree per nuisance)
#' # Recommended: Use theory-justified epsilon_n
#' epsilon_n <- treefarmr::select_epsilon_n(n, method = "fixed", c = 2)
#' fit_rashomon <- dml_att(X, A, Y, K = 5, use_rashomon = TRUE,
#'                        rashomon_bound_multiplier = epsilon_n)
#' print(fit_rashomon$theta)
#'
#' # Alternative: Use default (0.05) for quick exploratory analysis
#' fit_rashomon_quick <- dml_att(X, A, Y, K = 5, use_rashomon = TRUE)
#' print(fit_rashomon_quick$theta)
#' }
#' @export
dml_att <- function(X, A, Y, K = 5, outcome_type = c("binary", "continuous"),
                   regularization = 0.1, cv_regularization = FALSE, cv_K = 5,
                   stratified = TRUE, seed = NULL, verbose = FALSE,
                   use_rashomon = FALSE, rashomon_bound_multiplier = 0.05,
                   rashomon_bound_adder = 0, max_leaves = NULL,
                   auto_tune_intersecting = FALSE, ...) {
  outcome_type <- match.arg(outcome_type)
  check_dml_att_data(X, A, Y, outcome_type = outcome_type)
  if (is.matrix(X)) X <- as.data.frame(X)
  n <- nrow(X)

  # Validate parameters
  if (!is.numeric(K) || length(K) != 1 || K < 2) {
    stop("K must be a single integer >= 2, got: ", K, call. = FALSE)
  }

  n_treated <- sum(A == 1)
  n_control <- sum(A == 0)
  if (n_treated < K) {
    stop("Insufficient treated units for K-fold cross-fitting. ",
         "Need at least K=", K, " treated units, got: ", n_treated, call. = FALSE)
  }
  if (n_control < K) {
    stop("Insufficient control units for K-fold cross-fitting. ",
         "Need at least K=", K, " control units, got: ", n_control, call. = FALSE)
  }

  if (!cv_regularization && (!is.numeric(regularization) || length(regularization) != 1 || regularization <= 0)) {
    stop("regularization must be a single positive numeric value, got: ",
         regularization, call. = FALSE)
  }

  if (cv_regularization && (!is.numeric(cv_K) || length(cv_K) != 1 || cv_K < 2)) {
    stop("cv_K must be a single integer >= 2 when cv_regularization = TRUE, got: ",
         cv_K, call. = FALSE)
  }

  if (!is.numeric(rashomon_bound_multiplier) || length(rashomon_bound_multiplier) != 1 || rashomon_bound_multiplier < 0) {
    stop("rashomon_bound_multiplier must be a single non-negative numeric value, got: ",
         rashomon_bound_multiplier, call. = FALSE)
  }

  fold_indices <- create_folds(n, K, strata = if (stratified) A else NULL, seed = seed)

  if (use_rashomon) {
    nuisance_fits <- fit_nuisances_rashomon(X, A, Y, fold_indices, outcome_type = outcome_type,
                                           regularization = regularization,
                                           cv_regularization = cv_regularization, cv_K = cv_K,
                                           verbose = verbose,
                                           rashomon_bound_multiplier = rashomon_bound_multiplier,
                                           rashomon_bound_adder = rashomon_bound_adder,
                                           max_leaves = max_leaves,
                                           auto_tune_intersecting = auto_tune_intersecting, ...)
    eta <- get_fold_specific_eta_rashomon(nuisance_fits, X, fold_indices)
  } else {
    nuisance_fits <- vector("list", K)
    for (k in seq_len(K)) {
      nuisance_fits[[k]] <- fit_nuisances_fold(X, A, Y, fold_id = k, fold_indices = fold_indices,
                                              outcome_type = outcome_type,
                                              regularization = regularization,
                                              cv_regularization = cv_regularization, cv_K = cv_K,
                                              verbose = verbose, ...)
    }
    eta <- get_fold_specific_eta(nuisance_fits, X, fold_indices)
  }

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
