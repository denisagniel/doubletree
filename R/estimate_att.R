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
#' @param use_rashomon Logical. If TRUE, fit nuisances via \code{optimaltrees::cross_fitted_rashomon} (one interpretable tree per nuisance via intersection + refit per fold). Default FALSE (single tree per fold).
#' @param rashomon_bound_multiplier Numeric or NULL. Rashomon tolerance
#'   \eqn{\varepsilon_n} controlling the size of the Rashomon set (trees with
#'   penalized risk \eqn{\le (1 + \varepsilon_n) \cdot \text{best}}).
#'   \strong{Default: NULL}, which uses the theory-justified fixed value
#'   \eqn{\varepsilon_n = \log(n)/n} via \code{optimaltrees::select_epsilon_n(nrow(X))}.
#'
#'   \strong{Theory:} A fixed, deterministic \eqn{\varepsilon_n \asymp \log(n)/n}
#'   is \eqn{o(n^{-1/2})} and, under the structural-margin condition, yields both a
#'   non-empty cross-fold intersection and valid inference with zero intersection
#'   overhead (manuscript Corollary "Rashomon tolerance without the intersection
#'   trade-off"). Do not enlarge \eqn{\varepsilon_n} to force an intersection; if
#'   the intersection is empty at the theory value, the method falls back to
#'   fold-specific trees (still valid, one tree per fold instead of one shared tree).
#'   An explicit numeric is honored but only fixed, deterministic choices of order
#'   \eqn{\log(n)/n} carry the validity guarantee.
#' @param rashomon_bound_adder Numeric. Additive Rashomon bound (not recommended for cross-fitted estimation).
#'   Default: 0.
#' @param max_leaves Optional integer. Passed to \code{cross_fitted_rashomon} when \code{use_rashomon = TRUE}. Restricts Rashomon set to trees with at most this many leaves.
#' @param auto_tune_intersecting Logical. If TRUE, increase
#'   \code{rashomon_bound_multiplier} until a structure appears in the intersection
#'   across all K folds. Default: FALSE.
#'
#'   \strong{Not valid for inference.} Selecting \eqn{\varepsilon_n} from the data
#'   (post-selection) is not covered by the valid-inference theory and can inflate
#'   \eqn{\varepsilon_n} beyond \eqn{o(n^{-1/2})}, voiding the CLT; a \code{warning}
#'   is emitted when TRUE. Use it only for exploration. For inference, keep the
#'   fixed theory \eqn{\varepsilon_n} (leave \code{rashomon_bound_multiplier = NULL})
#'   and, if the intersection is empty, fall back to fold-specific trees
#'   (\code{use_rashomon = FALSE}).
#' @param escalate_intersection Logical. If TRUE (and \code{rashomon_bound_multiplier}
#'   is left NULL), widen the Rashomon tolerance \eqn{\varepsilon_n = c \cdot \log(n)/n}
#'   over an escalating grid of \eqn{c} (1, 2, ..., 1000) until the cross-fold
#'   intersection is non-empty, rather than pinning the single theory value. Default FALSE.
#'   Distinct from \code{auto_tune_intersecting}: escalation keeps \eqn{\lambda} fixed and
#'   steps the tolerance multiplier at the doubletree level, recording the selected \eqn{c}.
#'
#'   \strong{Trades the fixed-\eqn{\varepsilon_n} validity guarantee for a non-empty
#'   intersection.} The selected \eqn{c} is data-dependent (post-selection) and can be
#'   large at realistic \eqn{n}, so finite-sample coverage is not guaranteed and must be
#'   validated empirically (see the escalation coverage sweep in the arbitration
#'   simulation). Opt-in practical/study device; the returned \code{rashomon_c_e} /
#'   \code{rashomon_c_m0} record the selected multiplier per nuisance. An explicit
#'   \code{rashomon_bound_multiplier} always pins a single fixed tolerance, overriding it.
#' @param max_depth Integer. Maximum GOSDT tree depth for the nuisance fits
#'   (\code{0L} = unlimited). Default \code{4L}, applied to BOTH the Rashomon and the
#'   plain cross-fit paths so they are symmetric. Bounding depth prevents the
#'   continuous-covariate blow-up (unbounded GOSDT on many thresholds) and keeps the
#'   plain cross-fit nuisances from being deeper than their Rashomon twin. Set
#'   \code{max_depth = 0L} for unlimited depth (not recommended with continuous covariates).
#' @param discretize_method Character. Method for discretizing continuous features.
#'   Default: "quantiles" (theory-recommended, do not override unless you have good reason).
#'   Uses threshold encoding (k bins → k-1 features) for computational efficiency.
#' @param discretize_bins Integer or "adaptive". Number of bins for discretization.
#'   Default: "adaptive" (theory-recommended, do not override unless you have good reason).
#'   Uses b_n = max(2, ceiling(log(n)/3)) as suggested by nonparametric theory
#'   for optimal bias-variance tradeoff. Threshold encoding: k bins → k-1 binary features.
#' @param ... Additional arguments passed to optimaltrees (\code{fit_tree} when \code{use_rashomon = FALSE}, \code{cross_fitted_rashomon} when \code{use_rashomon = TRUE}).
#'
#' @details
#' \strong{Regularization Selection:} By default (\code{cv_regularization = TRUE}),
#' the regularization parameter \eqn{\lambda} is selected via 5-fold cross-validation
#' on each training fold, using a theory-driven grid: \eqn{(\log n / n) \times [0.25, 0.5, 1, 2, 4]}.
#' This implements the manuscript's recommendation that \eqn{\lambda \propto (\log n)/n} for
#' minimax-optimal trees. Fixed regularization (\code{cv_regularization = FALSE}) should
#' only be used when you have strong theoretical justification for a specific value.
#'
#' @return List with elements: theta (point estimate); sigma (estimated SE; the
#'   fully-fold-specific twin's SE when \code{use_rashomon = TRUE}, the plain Wald SE
#'   otherwise); ci_95 (when \code{use_rashomon = TRUE}, the HONEST bias-aware 95\% CI
#'   built from the fully-fold-specific twin -- the shared intersection structure is not
#'   orthogonal to each fold, so its Wald SE undercovers; otherwise the plain Wald 95\% CI);
#'   ci_95_wald (the plain Wald interval of the display estimate, always); theta_crossfit,
#'   sigma_crossfit, delta, se_delta (fully-fold-specific twin, bias diagnostic
#'   \eqn{\hat\theta - \hat\theta_{twin}}, and its SE; NA when \code{use_rashomon = FALSE});
#'   score_values (influence at theta); nuisance_fits (per-fold models or Rashomon list);
#'   fold_indices; n; K; converged (logical; TRUE if rashomon intersection succeeded or if
#'   use_rashomon=FALSE); epsilon_n (numeric; rashomon_bound_multiplier if use_rashomon=TRUE,
#'   NA otherwise).
#' @references Manuscript equation (2) for the orthogonal score.
#' @examples
#' \dontrun{
#' # Decision guide for key parameters:
#'
#' # 1. epsilon_n (rashomon_bound_multiplier):
#' #    - Default NULL uses the theory value optimaltrees::select_epsilon_n(nrow(X)) = log(n)/n
#' #    - Override with an explicit fixed numeric only if you have a reason
#'
#' # 2. regularization:
#' #    - Default (recommended): cv_regularization = TRUE (data-adaptive)
#' #    - Fixed only when theory-justified: cv_regularization = FALSE, regularization = 0.05
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
#' # Default: theory epsilon_n = log(n)/n (rashomon_bound_multiplier = NULL) +
#' # CV-selected lambda (both recommended)
#' fit1 <- estimate_att(
#'   X, A, Y,
#'   K = 5,
#'   use_rashomon = TRUE
#'   # rashomon_bound_multiplier = NULL -> optimaltrees::select_epsilon_n(nrow(X))
#'   # cv_regularization = TRUE is the default
#' )
#' print(fit1$theta)   # Point estimate
#' print(fit1$ci_95)   # 95\% Wald confidence interval
#'
#' # Alternative: Fixed lambda (when theory-justified)
#' fit2 <- estimate_att(
#'   X, A, Y,
#'   K = 5,
#'   cv_regularization = FALSE,
#'   regularization = 0.05
#' )
#' print(fit2$theta)
#' }
#' @export
estimate_att <- function(X, A, Y, K = 5, outcome_type = c("binary", "continuous"),
                   regularization = 0.1, cv_regularization = TRUE, cv_K = 5,
                   stratified = TRUE, seed = NULL, verbose = FALSE,
                   use_rashomon = FALSE, rashomon_bound_multiplier = NULL,
                   rashomon_bound_adder = 0, max_leaves = NULL,
                   auto_tune_intersecting = FALSE,
                   escalate_intersection = FALSE,
                   max_depth = 4L,
                   discretize_method = "quantiles",
                   discretize_bins = "adaptive",
                   ...) {
  outcome_type <- match.arg(outcome_type)
  check_att_data(X, A, Y, outcome_type = outcome_type)
  if (is.matrix(X)) X <- as.data.frame(X)
  n <- nrow(X)

  # Bound GOSDT search depth for BOTH nuisance paths (0L = unlimited). The Rashomon
  # path already defaulted to 4L; the plain cross-fit path (fit_nuisances_fold) was
  # previously uncapped, so continuous-covariate DGPs blew up and its nuisance trees
  # were deeper than the Rashomon twin's -- muddying comparisons. Default 4L makes the
  # two paths symmetric (see MEMORY estimate-att-depth-cap-asymmetry). Callers can pass
  # max_depth = 0L to restore unlimited depth.
  if (!is.numeric(max_depth) || length(max_depth) != 1 || is.na(max_depth) || max_depth < 0) {
    stop("max_depth must be a single non-negative integer (0 = unlimited), got: ",
         max_depth, call. = FALSE)
  }
  max_depth <- as.integer(max_depth)

  # Resolve the Rashomon tolerance. NULL -> theory value epsilon_n = log(n)/n
  # (= o(n^{-1/2})), the fixed, deterministic rate that yields valid inference
  # under the structural-margin condition (manuscript Cor. margin-resolution).
  # EXCEPTION: when escalate_intersection = TRUE and the caller gave no explicit
  # multiplier, leave it NULL so fit_nuisances_rashomon runs the c-grid escalation
  # (widening epsilon_n until the intersection is non-empty) instead of pinning the
  # single theory tolerance. An explicit multiplier is always honored as fixed.
  escalating <- isTRUE(escalate_intersection) && is.null(rashomon_bound_multiplier)
  if (is.null(rashomon_bound_multiplier) && !escalating) {
    rashomon_bound_multiplier <- optimaltrees::select_epsilon_n(n)
    if (verbose && use_rashomon) {
      message("Using theory epsilon_n = log(n)/n = ",
              signif(rashomon_bound_multiplier, 3))
    }
  }

  # Data-adaptive epsilon_n (auto_tune_intersecting) is a post-selection device
  # not covered by the valid-inference theory; warn when combined with Rashomon.
  if (use_rashomon && isTRUE(auto_tune_intersecting)) {
    warning(
      "auto_tune_intersecting = TRUE selects the Rashomon tolerance from the ",
      "data (post-selection) and voids the o(n^{-1/2}) valid-inference ",
      "guarantee; use it for exploration only. For inference, keep the fixed ",
      "theory epsilon_n (rashomon_bound_multiplier = NULL) and fall back to ",
      "use_rashomon = FALSE if the intersection is empty.",
      call. = FALSE
    )
  }

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

  # When escalating, rashomon_bound_multiplier stays NULL (the sentinel that triggers
  # c-grid escalation downstream); validate the numeric form only otherwise.
  if (!escalating &&
      (!is.numeric(rashomon_bound_multiplier) || length(rashomon_bound_multiplier) != 1 || rashomon_bound_multiplier < 0)) {
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
                                           escalate_intersection = escalate_intersection,
                                           max_depth = max_depth,
                                           discretize_method = discretize_method,
                                           discretize_bins = discretize_bins, ...)
    eta <- get_fold_specific_eta_rashomon(nuisance_fits, X, fold_indices)
  } else {
    nuisance_fits <- vector("list", K)
    for (k in seq_len(K)) {
      # max_depth flows via ... -> fit_tree_with_cv -> cv_regularization_adaptive/fit_tree,
      # bounding GOSDT depth on the plain cross-fit path (was previously uncapped).
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
  }

  # Shared EIF solve (closed form theta = sum(psi(0)) / sum(A/pi); see eif_att_solve).
  .att <- eif_att_solve(Y, A, eta$e, eta$m0, n)
  theta <- .att$theta
  score_values <- .att$score_values
  sigma <- .att$sigma
  ci_95 <- .att$ci_95              # plain Wald interval (used as-is on the non-Rashomon path)
  ci_95_wald <- .att$ci_95

  # HONEST (bias-aware) CI for the shared-Rashomon path. The shared intersection
  # structure "saw" every fold, so its per-fold Wald SE underestimates the estimate's
  # true spread (Phase-A 2026-07-15: se/emp_sd falls with n -> undercoverage). We pair
  # the shared estimate with the FULLY fold-specific twin (per-fold structure AND
  # leaves), whose delta = theta_shared - theta_twin captures that selection variance,
  # and report the AK honest interval widened by a conservative bound B = |delta| +
  # z*se_delta with se_delta = 0 (the tightest interval consistent with the coverage
  # guarantee; see honest_ci). This restores coverage >= 0.97 empirically. On the
  # fully-fold-specific path (use_rashomon = FALSE) the twin IS the estimator, so no
  # honesty correction is needed and the plain Wald ci_95 stands.
  theta_crossfit <- NA_real_; sigma_crossfit <- NA_real_
  delta <- NA_real_; se_delta <- NA_real_
  if (use_rashomon) {
    eta_twin <- get_fully_foldspecific_twin(
      X, A, Y, fold_indices, outcome_type = outcome_type,
      regularization = regularization, cv_regularization = cv_regularization,
      cv_K = cv_K, verbose = verbose, max_depth = max_depth,
      discretize_method = discretize_method, discretize_bins = discretize_bins, ...)
    .att_twin <- eif_att_solve(Y, A, eta_twin$e, eta_twin$m0, n)
    theta_crossfit <- .att_twin$theta
    sigma_crossfit <- .att_twin$sigma
    delta <- theta - theta_crossfit
    # se_delta = 0: raw |delta| bound (widening; a real se_delta only inflates B).
    se_delta <- 0
    hon <- honest_ci(theta, sigma_crossfit, delta, se_delta = se_delta, level = 0.95)
    sigma <- sigma_crossfit        # reported SE is the valid twin's SE
    ci_95 <- hon$ci                # reported CI is the honest bias-aware interval

    if (verbose) {
      message(sprintf("Shared-Rashomon ATT (display): %.4f", theta))
      message(sprintf("Fully-fold-specific twin ATT:  %.4f  (SE %.4f)",
                      theta_crossfit, sigma_crossfit))
      message(sprintf("Honest 95%% CI (bias-aware):    [%.4f, %.4f]  (cv %.2f, B %.4f)",
                      ci_95[1], ci_95[2], hon$cv, hon$B))
    }
  }

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
    # Tolerance multipliers selected by escalation (epsilon_n = c*log(n)/n per nuisance).
    rashomon_c_e  <- if (is.null(nuisance_fits$rashomon_c_e))  NA_real_ else nuisance_fits$rashomon_c_e
    rashomon_c_m0 <- if (is.null(nuisance_fits$rashomon_c_m0)) NA_real_ else nuisance_fits$rashomon_c_m0
    # Report the actual escalated tolerance (max over nuisances), not the input NULL.
    c_vals <- c(rashomon_c_e, rashomon_c_m0)
    epsilon_n <- if (all(is.na(c_vals))) NA_real_ else max(c_vals, na.rm = TRUE) * (log(n) / n)
  } else {
    # Non-rashomon always converges (uses fold-specific trees)
    converged <- TRUE
    epsilon_n <- NA_real_
    rashomon_c_e <- NA_real_
    rashomon_c_m0 <- NA_real_
  }

  list(
    theta = theta,
    sigma = sigma,                 # honest twin SE (Rashomon) or Wald SE (fold-specific)
    ci_95 = ci_95,                 # honest bias-aware CI (Rashomon) or Wald CI (fold-specific)
    ci_95_wald = ci_95_wald,       # plain Wald interval of the display estimate (always)
    theta_crossfit = theta_crossfit,  # fully-fold-specific twin (NA if use_rashomon=FALSE)
    sigma_crossfit = sigma_crossfit,
    delta = delta,                 # theta_shared - theta_twin (bias diagnostic)
    se_delta = se_delta,
    score_values = score_values,
    nuisance_fits = nuisance_fits,
    fold_indices = fold_indices,
    n = n,
    K = K,
    converged = converged,
    epsilon_n = epsilon_n,
    rashomon_c_e = rashomon_c_e,
    rashomon_c_m0 = rashomon_c_m0
  )
}
