#' Cross-Fitted Estimate of the Average Treatment Effect on the Treated (ATT)
#'
#' Estimates the ATT by K-fold cross-fitting: for each fold \eqn{k}, a single
#' optimal decision tree per nuisance (propensity \eqn{e(X)} and control outcome
#' \eqn{m_0(X)}) is fit on \eqn{D^{(-k)}} via optimaltrees and used to predict the
#' held-out fold \eqn{k}. The resulting out-of-sample nuisance predictions are fed
#' to the efficient influence function (EIF) solver, giving a doubly robust,
#' semiparametric estimator with a plain Wald interval. Binary outcomes use
#' log-loss for both nuisances; continuous outcomes use log-loss for the
#' propensity and squared_error for \eqn{m_0}.
#'
#' Both the tree STRUCTURE and the leaf values are learned per fold, so each
#' fold's nuisance fit is orthogonal to the fold it predicts. There is no
#' Rashomon-set structure intersection and no shared tree across folds; for that
#' variant see \code{\link{estimate_att_rashomon}}.
#'
#' @param X Data.frame or matrix of covariates. Must be binary (0/1) for optimaltrees.
#' @param A Integer or numeric vector of treatment (0/1).
#' @param Y Numeric vector of outcome. Binary (0/1) when outcome_type is "binary"; any numeric when "continuous".
#' @param K Number of cross-fitting folds. Default 5.
#' @param outcome_type Character. "binary" (default) or "continuous". Continuous requires optimaltrees squared_error loss for m0, m1.
#' @param regularization Numeric. Tree complexity penalty passed to optimaltrees. Default 0.1.
#'   Only used if \code{cv_regularization = FALSE}. For most applications, use
#'   \code{cv_regularization = TRUE} (default) for data-adaptive selection.
#' @param cv_regularization Logical. If TRUE (default), use cross-validation to select
#'   regularization parameter \eqn{\lambda} separately for each nuisance function
#'   (e, m0) using a theory-driven grid centered on \eqn{(\log n)/n}. If FALSE, use
#'   fixed \code{regularization} value.
#'
#'   \strong{When to use TRUE (recommended):} You don't know the right penalty or want
#'   robustness across varied data structures. Uses theory-driven grid:
#'   \eqn{(\log n)/n \times [0.25, 0.5, 1, 2, 4]}. Adds computational cost (nested CV)
#'   but improves model selection and inference quality.
#'
#'   \strong{When to use FALSE:} You have a theory-justified fixed value (e.g., from
#'   \code{optimaltrees::cv_regularization()} on pilot data) or need maximum speed.
#'   Set \code{cv_regularization = FALSE} only when you have strong theoretical
#'   justification for a specific \eqn{\lambda} value.
#'
#'   \strong{Theory:} Manuscript recommends \eqn{\lambda \propto (\log n)/n} for
#'   minimax-optimal trees. CV automatically implements this recommendation.
#' @param cv_K Integer. Number of folds for cross-validation of regularization. Default 5. Only used if \code{cv_regularization = TRUE}.
#' @param stratified Logical. If TRUE (default), fold assignment is stratified by A.
#' @param seed Optional. Random seed for fold creation.
#' @param verbose Logical. Passed to optimaltrees. Default FALSE.
#' @param max_depth Integer. Maximum GOSDT search depth for the nuisance fits
#'   (\code{0L} = unlimited). Default \code{4L}, matching
#'   \code{\link{estimate_att_rashomon}} so the two paths are symmetric. This bounds
#'   the branch-and-bound \emph{search space} GOSDT explores, not the leaf count of the
#'   returned tree -- those are orthogonal: a low-leaf-count tree can still come from an
#'   unbounded-depth search (leaf count and depth are not tightly coupled; a tree
#'   with \code{L} leaves can be as deep as \code{L - 1}). With continuous
#'   covariates, discretization yields many threshold-indicator features, so an
#'   unbounded search must still explore/rule out combinatorially many deep
#'   candidate splits before any leaf-count regularization gets a chance to reject
#'   them -- a tractability guard, not a substitute for leaf-budget control (see
#'   \code{bisect_lambda_to_budget()} in optimaltrees for the latter). Set
#'   \code{max_depth = 0L} for unlimited depth (not recommended with continuous covariates).
#' @param discretize_method Character. Method for discretizing continuous features.
#'   Default: "quantiles" (theory-recommended, do not override unless you have good reason).
#'   Uses threshold encoding (k bins → k-1 features) for computational efficiency.
#' @param discretize_bins Integer or "adaptive". Number of bins for discretization.
#'   Default: "adaptive" (theory-recommended, do not override unless you have good reason).
#'   Uses b_n = max(2, ceiling(log(n)/3)) as suggested by nonparametric theory
#'   for optimal bias-variance tradeoff. Threshold encoding: k bins → k-1 binary features.
#' @param ... Additional arguments passed to optimaltrees via \code{fit_tree}.
#'
#' @details
#' \strong{Regularization Selection:} By default (\code{cv_regularization = TRUE}),
#' the regularization parameter \eqn{\lambda} is selected via 5-fold cross-validation
#' on each training fold, using a theory-driven grid: \eqn{(\log n / n) \times [0.25, 0.5, 1, 2, 4]}.
#' This implements the manuscript's recommendation that \eqn{\lambda \propto (\log n)/n} for
#' minimax-optimal trees. Fixed regularization (\code{cv_regularization = FALSE}) should
#' only be used when you have strong theoretical justification for a specific value.
#'
#' @return List with elements: theta (point estimate); sigma (Wald standard error);
#'   ci_95 (Wald 95\% confidence interval); ci_95_wald (identical to \code{ci_95} on
#'   this path, retained so callers can compare display and Wald intervals uniformly
#'   across estimators); score_values (influence at theta); nuisance_fits (per-fold
#'   models, plus \code{propensity} and \code{outcome_control} out-of-sample
#'   predictions); fold_indices; n; K.
#' @references Manuscript equation (2) for the orthogonal score.
#' @examples
#' \dontrun{
#' # Decision guide for key parameters:
#'
#' # regularization:
#' #    - Default (recommended): cv_regularization = TRUE (data-adaptive)
#' #    - Fixed only when theory-justified: cv_regularization = FALSE, regularization = 0.05
#'
#' # Recommended workflow for new dataset:
#' library(optimaltrees)  # Required dependency
#' set.seed(42)
#' n <- 300
#' X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5))
#' A <- rbinom(n, 1, plogis(0.5 * X$X1 - 0.2))
#' Y <- rbinom(n, 1, 0.3 + 0.2 * X$X1 + 0.15 * A)
#'
#' # Default: CV-selected lambda (recommended)
#' fit1 <- estimate_att_crossfit(X, A, Y, K = 5)
#' print(fit1$theta)   # Point estimate
#' print(fit1$ci_95)   # 95\% Wald confidence interval
#'
#' # Alternative: Fixed lambda (when theory-justified)
#' fit2 <- estimate_att_crossfit(
#'   X, A, Y,
#'   K = 5,
#'   cv_regularization = FALSE,
#'   regularization = 0.05
#' )
#' print(fit2$theta)
#' }
#' @export
estimate_att_crossfit <- function(X, A, Y, K = 5, outcome_type = c("binary", "continuous"),
                   regularization = 0.1, cv_regularization = TRUE, cv_K = 5,
                   stratified = TRUE, seed = NULL, verbose = FALSE,
                   max_depth = 4L,
                   discretize_method = "quantiles",
                   discretize_bins = "adaptive",
                   ...) {
  outcome_type <- match.arg(outcome_type)
  check_att_data(X, A, Y, outcome_type = outcome_type)
  if (is.matrix(X)) X <- as.data.frame(X)
  n <- nrow(X)

  # Bound the GOSDT *search space* (0L = unlimited). This is orthogonal to
  # leaf-count/interpretability budgeting (see bisect_lambda_to_budget() in
  # optimaltrees): leaf count and depth are not tightly coupled, so a low-leaf-count
  # tree can still come from an unbounded-depth search, and continuous-covariate
  # discretization means an unbounded search must explore combinatorially many deep
  # candidate splits before the leaf-count penalty ever gets to reject them. This path
  # (fit_nuisances_fold) was previously uncapped, so continuous-covariate DGPs blew up
  # and its nuisance trees were deeper than the Rashomon twin's -- muddying comparisons.
  # Default 4L matches estimate_att_rashomon (see MEMORY estimate-att-depth-cap-asymmetry).
  # Callers can pass max_depth = 0L to restore unlimited depth.
  if (!is.numeric(max_depth) || length(max_depth) != 1 || is.na(max_depth) || max_depth < 0) {
    stop("max_depth must be a single non-negative integer (0 = unlimited), got: ",
         max_depth, call. = FALSE)
  }
  max_depth <- as.integer(max_depth)

  # Validate critical parameters only (R's type coercion handles the rest)

  # K must be integer >= 2
  if (!is.numeric(K) || length(K) != 1 || K < 2) {
    stop("K must be a single integer >= 2, got: ", K, call. = FALSE)
  }
  if (K != as.integer(K)) {
    stop("K must be an integer, got: ", K, call. = FALSE)
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

  fold_indices <- create_folds(n, K, strata = if (stratified) A else NULL, seed = seed)

  nuisance_fits <- vector("list", K)
  for (k in seq_len(K)) {
    # max_depth flows via ... -> fit_tree_with_cv -> cv_regularization_adaptive/fit_tree,
    # bounding GOSDT depth on this cross-fit path (was previously uncapped).
    nuisance_fits[[k]] <- fit_nuisances_fold(X, A, Y, fold_id = k, fold_indices = fold_indices,
                                            outcome_type = outcome_type,
                                            regularization = regularization,
                                            cv_regularization = cv_regularization, cv_K = cv_K,
                                            verbose = verbose,
                                            max_depth = max_depth,
                                            discretize_method = discretize_method,
                                            discretize_bins = discretize_bins, ...)
  }
  eta <- get_fold_specific_eta(nuisance_fits, X, fold_indices)

  # Shared EIF solve (closed form theta = sum(psi(0)) / sum(A/pi); see eif_att_solve).
  .att <- eif_att_solve(Y, A, eta$e, eta$m0, n)
  theta <- .att$theta
  score_values <- .att$score_values
  sigma <- .att$sigma
  ci_95 <- .att$ci_95              # plain Wald interval (used as-is on this path)
  ci_95_wald <- .att$ci_95

  # No honesty correction is needed here: both the structure AND the leaves are fit
  # out-of-fold, so each fold's nuisance is orthogonal to the fold it predicts and the
  # plain Wald interval covers. (This estimator IS the "fully fold-specific twin" that
  # estimate_att_rashomon() constructs to bias-correct its shared-structure estimate.)

  # Add predictions to nuisance_fits for diagnostics
  nuisance_fits$propensity <- eta$e
  nuisance_fits$outcome_control <- eta$m0

  list(
    theta = theta,
    sigma = sigma,                 # Wald SE
    ci_95 = ci_95,                 # Wald CI
    ci_95_wald = ci_95_wald,       # same interval, named for cross-estimator comparability
    score_values = score_values,
    nuisance_fits = nuisance_fits,
    fold_indices = fold_indices,
    n = n,
    K = K
  )
}
