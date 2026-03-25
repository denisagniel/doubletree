#' Estimate the Average Treatment Effect on the Treated (ATT)
#'
#' Estimates the ATT using efficient influence function-based estimation with
#' cross-fitting and optimal decision trees (via optimaltrees) for the nuisance
#' functions e(X) and m0(X). This is a doubly robust, semiparametric estimator.
#' Binary outcomes use log-loss for both nuisances; continuous outcomes use
#' log-loss for propensity and squared_error for m0.
#'
#' When \code{use_rashomon = TRUE}, nuisances are fit via
#' \code{optimaltrees::cross_fitted_rashomon}: one interpretable tree per nuisance
#' (e, m0) via intersection of Rashomon sets across folds with fold-specific refits
#' for valid cross-fitted estimation. The same K and fold assignment are used for
#' Rashomon and the score.
#'
#' @param X Data.frame or matrix of covariates. Must be binary (0/1) for optimaltrees.
#' @param A Integer or numeric vector of treatment (0/1).
#' @param Y Numeric vector of outcome. Binary (0/1) when outcome_type is "binary"; any numeric when "continuous".
#' @param K Number of cross-fitting folds. Default 5.
#' @param outcome_type Character. "binary" (default) or "continuous". Continuous requires optimaltrees squared_error loss for m0, m1.
#' @param regularization Numeric. Tree complexity penalty passed to optimaltrees. Default 0.1. Ignored if \code{cv_regularization = TRUE}.
#' @param cv_regularization Logical. If TRUE, use cross-validation to select
#'   regularization parameter \eqn{\lambda} separately for each nuisance function
#'   (e, m0). If FALSE (default), use fixed \code{regularization} value.
#'
#'   \strong{When to use TRUE:} You don't know the right penalty or want robustness
#'   across varied data structures. Adds computational cost (nested CV) but improves
#'   model selection.
#'
#'   \strong{When to use FALSE:} You have a theory-justified choice (e.g., from
#'   \code{optimaltrees::cv_regularization()} on pilot data) or want speed. Fixed
#'   \code{regularization} is faster and reproducible.
#'
#'   \strong{Theory:} Manuscript recommends \eqn{\lambda \propto (\log n)/n} for
#'   minimax-optimal trees. See \code{optimaltrees::cv_regularization()} for automatic
#'   selection implementing this rate.
#' @param cv_K Integer. Number of folds for cross-validation of regularization. Default 5. Only used if \code{cv_regularization = TRUE}.
#' @param stratified Logical. If TRUE (default), fold assignment is stratified by A.
#' @param seed Optional. Random seed for fold creation.
#' @param verbose Logical. Passed to optimaltrees. Default FALSE.
#' @param use_rashomon Logical. If TRUE, fit nuisances via \code{optimaltrees::cross_fitted_rashomon} (one interpretable tree per nuisance via intersection + refit per fold). Default FALSE (single tree per fold).
#' @param rashomon_bound_multiplier Numeric. Rashomon tolerance \eqn{\varepsilon_n}
#'   controlling the size of the Rashomon set (trees with penalized risk
#'   \eqn{\le (1 + \varepsilon) \cdot \text{best}}). Default: 0.05 (for quick exploration).
#'
#'   \strong{Recommended:} Use theory-justified value via
#'   \code{optimaltrees::select_epsilon_n(nrow(X), method = "fixed", c = 2)}.
#'   This sets \eqn{\varepsilon_n = c\sqrt{(\log n)/n}}, which satisfies the
#'   rate requirement \eqn{o(n^{-1/2})} for the 2-nuisance ATT score (manuscript Appendix A.5).
#'
#'   \strong{Trade-off:} Smaller \eqn{\varepsilon_n} yields trees closer to optimal
#'   but higher risk of empty intersection. Larger \eqn{\varepsilon_n} facilitates
#'   intersection (interpretability gain) but includes more sub-optimal trees.
#'   Typical range: 0.02 (tight) to 0.10 (loose).
#' @param rashomon_bound_adder Numeric. Additive Rashomon bound (not recommended for cross-fitted estimation).
#'   Default: 0.
#' @param max_leaves Optional integer. Passed to \code{cross_fitted_rashomon} when \code{use_rashomon = TRUE}. Restricts Rashomon set to trees with at most this many leaves.
#' @param auto_tune_intersecting Logical. If TRUE, automatically increase
#'   \code{rashomon_bound_multiplier} until at least one tree structure appears in
#'   the intersection across all K folds. Default: FALSE.
#'
#'   \strong{Use with caution:} Automatically selecting \eqn{\varepsilon_n} based on
#'   intersection status is a heuristic that may yield arbitrarily large tolerance.
#'   If intersection is empty for reasonable \eqn{\varepsilon_n \le 0.2}, this
#'   indicates substantial cross-fold heterogeneity. Consider falling back to
#'   fold-specific trees (\code{use_rashomon = FALSE}) instead.
#' @param discretize_method Character. Method for discretizing continuous features.
#'   Default: "quantiles" (theory-recommended, do not override unless you have good reason).
#'   Uses threshold encoding (k bins → k-1 features) for computational efficiency.
#' @param discretize_bins Integer or "adaptive". Number of bins for discretization.
#'   Default: "adaptive" (theory-recommended, do not override unless you have good reason).
#'   Uses b_n = max(2, ceiling(log(n)/3)) as suggested by nonparametric theory
#'   for optimal bias-variance tradeoff. Threshold encoding: k bins → k-1 binary features.
#' @param ... Additional arguments passed to optimaltrees (\code{fit_tree} when \code{use_rashomon = FALSE}, \code{cross_fitted_rashomon} when \code{use_rashomon = TRUE}).
#' @return List with elements: theta (point estimate), sigma (estimated SE), ci_95 (Wald 95% CI),
#'   score_values (influence at theta), nuisance_fits (per-fold models or Rashomon list), fold_indices, n, K,
#'   converged (logical; TRUE if rashomon intersection succeeded or if use_rashomon=FALSE),
#'   epsilon_n (numeric; rashomon_bound_multiplier if use_rashomon=TRUE, NA otherwise).
#' @references Manuscript equation (2) for the orthogonal score.
#' @examples
#' \dontrun{
#' # Decision guide for key parameters:
#'
#' # 1. epsilon_n (rashomon_bound_multiplier):
#' #    - Use theory: epsilon_n <- optimaltrees::select_epsilon_n(n, c = 2)
#' #    - Quick exploratory: rashomon_bound_multiplier = 0.05 (default)
#'
#' # 2. regularization:
#' #    - Use CV if unsure: cv_regularization = TRUE
#' #    - Fixed if known: regularization = 0.1 (or from pilot CV)
#'
#' # 3. Rashomon vs fold-specific:
#' #    - Rashomon (use_rashomon = TRUE): interpretability, single tree/nuisance
#' #    - Fold-specific (FALSE): robustness, no intersection requirement
#'
#' # Recommended workflow for new dataset:
#' library(optimaltrees)  # Required dependency
#' set.seed(42)
#' n <- 300
#' X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5))
#' A <- rbinom(n, 1, plogis(0.5 * X$X1 - 0.2))
#' Y <- rbinom(n, 1, 0.3 + 0.2 * X$X1 + 0.15 * A)
#'
#' # Theory-justified epsilon_n
#' epsilon_n <- optimaltrees::select_epsilon_n(n, method = "fixed", c = 2)
#'
#' # Recommended: Rashomon with CV regularization
#' fit <- estimate_att(
#'   X, A, Y,
#'   K = 5,
#'   use_rashomon = TRUE,
#'   rashomon_bound_multiplier = epsilon_n,
#'   cv_regularization = TRUE  # Auto-select lambda
#' )
#' print(fit$theta)   # Point estimate
#' print(fit$ci_95)   # 95% Wald confidence interval
#'
#' # Alternative: Quick exploratory analysis
#' fit_quick <- estimate_att(X, A, Y, K = 5, regularization = 0.1)
#' print(fit_quick$theta)
#' }
#' @export
estimate_att <- function(X, A, Y, K = 5, outcome_type = c("binary", "continuous"),
                   regularization = 0.1, cv_regularization = FALSE, cv_K = 5,
                   stratified = TRUE, seed = NULL, verbose = FALSE,
                   use_rashomon = FALSE, rashomon_bound_multiplier = 0.05,
                   rashomon_bound_adder = 0, max_leaves = NULL,
                   auto_tune_intersecting = FALSE,
                   discretize_method = "quantiles",
                   discretize_bins = "adaptive",
                   ...) {
  outcome_type <- match.arg(outcome_type)
  check_att_data(X, A, Y, outcome_type = outcome_type)
  if (is.matrix(X)) X <- as.data.frame(X)
  n <- nrow(X)

  # Validate critical parameters only (R's type coercion handles the rest)

  # K must be integer >= 2
  if (!is.numeric(K) || length(K) != 1 || K < 2) {
    stop("K must be a single integer >= 2, got: ", K, call. = FALSE)
  }
  if (K != as.integer(K)) {
    stop("K must be an integer, got: ", K, call. = FALSE)
  }

  # max_leaves must be positive integer if provided
  if (!is.null(max_leaves)) {
    if (!is.numeric(max_leaves) || length(max_leaves) != 1 || max_leaves < 1) {
      stop("max_leaves must be NULL or a single positive integer, got: ", max_leaves, call. = FALSE)
    }
    if (max_leaves != as.integer(max_leaves)) {
      stop("max_leaves must be an integer, got: ", max_leaves, call. = FALSE)
    }
  }

  # Valid discretize_method
  valid_methods <- c("quantiles", "median")
  if (!discretize_method %in% valid_methods) {
    stop("discretize_method must be one of: ", paste(valid_methods, collapse = ", "),
         ", got: ", discretize_method, call. = FALSE)
  }

  # Consolidated validation: check for sufficient treated/control units
  # Need at least 2 units per fold for meaningful cross-fitting
  n_treated <- sum(A == 1)
  n_control <- sum(A == 0)
  min_per_fold <- 2

  if (n_treated < K * min_per_fold) {
    stop("Insufficient treated units for K=", K, " fold cross-fitting. ",
         "Need at least ", K * min_per_fold, " treated units, got: ", n_treated, ". ",
         "Either reduce K or collect more data.",
         call. = FALSE)
  }

  if (n_control < K * min_per_fold) {
    stop("Insufficient control units for K=", K, " fold cross-fitting. ",
         "Need at least ", K * min_per_fold, " control units, got: ", n_control, ". ",
         "Either reduce K or collect more data.",
         call. = FALSE)
  }

  # Validate treatment proportion (pi_hat) is in (0,1)
  # This is redundant with the sample size checks above, but serves as explicit validation
  pi_hat <- mean(A)
  if (pi_hat <= 0 || pi_hat >= 1) {
    stop("Invalid treatment proportion: pi_hat = ", pi_hat,
         ". This should not happen after sample size validation.",
         call. = FALSE)
  }

  # Regularization must be positive if not using CV
  if (!cv_regularization) {
    if (!is.numeric(regularization) || length(regularization) != 1 || regularization <= 0) {
      stop("regularization must be a single positive numeric value, got: ",
           regularization, call. = FALSE)
    }
  }

  # cv_K must be integer >= 2 if using CV regularization
  if (cv_regularization) {
    if (!is.numeric(cv_K) || length(cv_K) != 1 || cv_K < 2) {
      stop("cv_K must be a single integer >= 2 when cv_regularization = TRUE, got: ",
           cv_K, call. = FALSE)
    }
    if (cv_K != as.integer(cv_K)) {
      stop("cv_K must be an integer, got: ", cv_K, call. = FALSE)
    }
  }

  if (!is.numeric(rashomon_bound_multiplier) || length(rashomon_bound_multiplier) != 1 || rashomon_bound_multiplier < 0) {
    stop("rashomon_bound_multiplier must be a single non-negative numeric value, got: ",
         rashomon_bound_multiplier, call. = FALSE)
  }

  # Issue #21: Validate rashomon_bound_adder
  if (!is.numeric(rashomon_bound_adder) || length(rashomon_bound_adder) != 1 || rashomon_bound_adder < 0) {
    stop("rashomon_bound_adder must be a single non-negative numeric value, got: ",
         rashomon_bound_adder, call. = FALSE)
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
                                           auto_tune_intersecting = auto_tune_intersecting,
                                           discretize_method = discretize_method,
                                           discretize_bins = discretize_bins, ...)
    eta <- get_fold_specific_eta_rashomon(nuisance_fits, X, fold_indices)
  } else {
    nuisance_fits <- vector("list", K)
    for (k in seq_len(K)) {
      nuisance_fits[[k]] <- fit_nuisances_fold(X, A, Y, fold_id = k, fold_indices = fold_indices,
                                              outcome_type = outcome_type,
                                              regularization = regularization,
                                              cv_regularization = cv_regularization, cv_K = cv_K,
                                              verbose = verbose,
                                              discretize_method = discretize_method,
                                              discretize_bins = discretize_bins, ...)
    }
    eta <- get_fold_specific_eta(nuisance_fits, X, fold_indices)
  }

  # pi_hat already computed during validation above
  # Closed form: psi(theta) = psi(0) - theta*(A/pi), so sum(psi(theta)) = 0 => theta = sum(psi(0)) / sum(A/pi).
  sum_a_over_pi <- sum(A / pi_hat)
  score_at_zero <- psi_att(Y, A, theta = 0, eta, pi_hat)
  theta <- sum(score_at_zero) / sum_a_over_pi

  score_values <- psi_att(Y, A, theta, eta, pi_hat)
  sigma <- att_se(score_values, n)
  ci_95 <- att_ci(theta, sigma, n, level = 0.95)

  # Add predictions to nuisance_fits for diagnostics
  nuisance_fits$propensity <- eta$e
  nuisance_fits$outcome_control <- eta$m0

  # Add convergence information for rashomon method
  if (use_rashomon) {
    # Rashomon converged if both models have intersecting trees (no fallback)
    converged <- !is.null(nuisance_fits$cf_e) &&
                 !is.null(nuisance_fits$cf_m0) &&
                 nuisance_fits$cf_e@n_intersecting > 0 &&
                 nuisance_fits$cf_m0@n_intersecting > 0
    epsilon_n <- rashomon_bound_multiplier
  } else {
    # Non-rashomon always converges (uses fold-specific trees)
    converged <- TRUE
    epsilon_n <- NA_real_
  }

  list(
    theta = theta,
    sigma = sigma,
    ci_95 = ci_95,
    score_values = score_values,
    nuisance_fits = nuisance_fits,
    fold_indices = fold_indices,
    n = n,
    K = K,
    converged = converged,
    epsilon_n = epsilon_n
  )
}
