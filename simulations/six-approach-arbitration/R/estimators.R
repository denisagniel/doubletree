# =============================================================================
# estimators.R -- the seven ATT estimators under comparison
# =============================================================================
# Six original approaches + Alternative A (single_tree). Dispatched on
# config$method. Each returns a NAMED LIST with the shared contract
#   estimate, std_error, ci_lower, ci_upper, converged
# plus (where available) goal-(ii) fields:
#   theta_crossfit, se_crossfit, delta, delta_over_se, intersection_nonempty
# so run_one.R can record single-tree AND cross-fit-twin coverage and the
# fidelity diagnostic. Missing fields are returned as NA.
#
# CORRECTED vs the old bespoke harness (six_approach_comparison/code/estimators.R):
# Rashomon approaches use the fixed theory tolerance
# rashomon_bound_multiplier = NULL -> optimaltrees::select_epsilon_n(n) = log(n)/n
# (= o(n^{-1/2})), and auto_tune_intersecting = FALSE. The old
# eps_n = 2*sqrt(log n/n) + auto_tune=TRUE config is invalidated (see the
# structural-margin resolution, manuscript Cor.).
#
# ESCALATION: with escalate_intersection = FALSE (default), an empty theory-tolerance
# intersection falls back to fold-specific trees; the tolerance is NOT widened and
# rashomon_c_* == 1. Escalation (widening eps_n = c*log(n)/n until the intersection is
# non-empty) is OPT-IN via escalate_intersection = TRUE -- threaded from config$escalate
# so the escalation coverage sweep can measure whether large c degrades coverage. It
# trades the fixed-eps_n validity guarantee for a non-empty intersection.
# =============================================================================

SIM_K        <- 5L    # cross-fitting folds
SIM_M        <- 10L   # modal splits (msplit approaches)
SIM_K_MSPLIT <- 5L    # folds per split (msplit approaches)

# Contract skeleton: fill what an approach provides, NA otherwise.
# rashomon_c_e / rashomon_c_m0: the Rashomon tolerance multiplier selected by
# escalation for each nuisance (epsilon_n = c * log(n)/n). c=1 is the theory value;
# larger c means the theory-tolerance intersection was empty and had to be widened.
# Recorded so /data-analysis can evaluate whether large c degrades coverage.
# ci_lower/ci_upper default to the symmetric Wald interval estimate +/- z*std_error,
# but callers may pass an explicit (possibly asymmetric) interval -- e.g. the honest
# bias-aware CI of the averaged estimators, which is NOT estimate +/- z*std_error.
.result <- function(estimate, std_error, converged,
                    ci_lower = NULL, ci_upper = NULL,
                    theta_crossfit = NA_real_, se_crossfit = NA_real_,
                    delta = NA_real_, delta_over_se = NA_real_,
                    intersection_nonempty = NA,
                    rashomon_c_e = NA_real_, rashomon_c_m0 = NA_real_) {
  z <- stats::qnorm(0.975)
  if (is.null(ci_lower)) ci_lower <- estimate - z * std_error
  if (is.null(ci_upper)) ci_upper <- estimate + z * std_error
  list(
    estimate  = estimate,
    std_error = std_error,
    ci_lower  = ci_lower,
    ci_upper  = ci_upper,
    converged = converged,
    theta_crossfit = theta_crossfit,
    se_crossfit    = se_crossfit,
    delta          = delta,
    delta_over_se  = delta_over_se,
    intersection_nonempty = intersection_nonempty,
    rashomon_c_e   = rashomon_c_e,
    rashomon_c_m0  = rashomon_c_m0
  )
}

# Dispatch. One branch per method in GRID$method. `escalate` (default FALSE) is the
# Rashomon-tolerance escalation flag from config$escalate; only the intersection-based
# methods consume it (others ignore it).
estimate <- function(data, config) {
  X <- data$X; A <- data$A; Y <- data$Y
  escalate <- isTRUE(config$escalate)
  switch(
    config$method,
    full            = .est_full(X, A, Y),
    crossfit        = .est_crossfit(X, A, Y),
    doubletree      = .est_doubletree(X, A, Y, escalate = escalate),
    dt_averaged     = .est_dt_averaged(X, A, Y, escalate = escalate),
    msplit          = .est_msplit(X, A, Y),
    msplit_averaged = .est_msplit_averaged(X, A, Y),
    single_tree     = .est_single_tree(X, A, Y, escalate = escalate),
    stop(sprintf("Unknown method: '%s'", config$method))  # no silent fallback
  )
}

# --- 1. Full-sample: one tree per nuisance, in-sample predictions (baseline) --
# Structure chosen by CV on all n (no honesty) -> structural in-sample bias.
.est_full <- function(X, A, Y) {
  n <- nrow(X)
  cv_e <- optimaltrees::cv_regularization_adaptive(
    X = X, y = A, loss_function = "log_loss", K = SIM_K,
    max_lambda = 20 * log(n) / n, refit = TRUE, verbose = FALSE,
    max_depth = 4L)   # match Rashomon path; bound GOSDT search on continuous X
  if (is.na(cv_e$best_lambda)) stop("full: CV failed for propensity.", call. = FALSE)
  e_hat <- predict(cv_e$model, X, type = "prob")[, 2]

  idx0 <- which(A == 0); n0 <- length(idx0)
  cv_m0 <- optimaltrees::cv_regularization_adaptive(
    X = X[idx0, , drop = FALSE], y = Y[idx0], loss_function = "log_loss",
    K = SIM_K, max_lambda = 20 * log(n0) / n0, refit = TRUE, verbose = FALSE,
    max_depth = 4L)   # match Rashomon path; bound GOSDT search on continuous X
  if (is.na(cv_m0$best_lambda)) stop("full: CV failed for outcome.", call. = FALSE)
  m0_hat <- predict(cv_m0$model, X, type = "prob")[, 2]

  a <- .att(Y, A, e_hat, m0_hat)
  .result(a$theta, a$sigma, converged = TRUE)
}

# --- 2. Standard cross-fit: K separate trees, out-of-sample (valid, no 1 tree) -
# max_depth = 4L matches .est_full and the Rashomon path. As of 2026-07-14 this is the
# estimate_att DEFAULT (the package now caps both nuisance paths at 4L), so passing it
# is redundant -- kept explicit to pin the comparison against future default changes.
# Bounding depth avoids the continuous-covariate blow-up (unbounded GOSDT: n=2000/
# continuous >24 GB, ~1300 s/unit vs ~150 s capped) and keeps the plain cross-fit trees
# from being deeper than the Rashomon twin.
.est_crossfit <- function(X, A, Y) {
  r <- doubletree::estimate_att(X, A, Y, K = SIM_K, use_rashomon = FALSE,
                                max_depth = 4L, verbose = FALSE)
  .result(r$theta, r$sigma, converged = isTRUE(r$converged))
}

# --- 3. Doubletree: Rashomon-intersection structure, cross-fit leaves (twin) --
# As of Phase B (2026-07-15) estimate_att(use_rashomon=TRUE) reports an HONEST bias-aware
# CI built from the FULLY fold-specific twin (the shared intersection structure is not
# orthogonal to each fold, so its Wald SE undercovers). We pass that honest interval as
# ci_lower/ci_upper and record the twin + delta for the 3-way coverage comparison.
.est_doubletree <- function(X, A, Y, escalate = FALSE) {
  r <- doubletree::estimate_att(
    X, A, Y, K = SIM_K, use_rashomon = TRUE,
    rashomon_bound_multiplier = NULL, auto_tune_intersecting = FALSE,
    escalate_intersection = escalate,
    verbose = FALSE)
  .result(r$theta, r$sigma, converged = isTRUE(r$converged),
          ci_lower = r$ci_95[1], ci_upper = r$ci_95[2],
          theta_crossfit = r$theta_crossfit, se_crossfit = r$sigma_crossfit,
          delta = r$delta,
          intersection_nonempty = isTRUE(r$converged),
          rashomon_c_e = r$rashomon_c_e, rashomon_c_m0 = r$rashomon_c_m0)
}

# --- 4. Doubletree-averaged: intersection structure, averaged leaves, 1 tree --
# Point estimate = averaged (biased) tree; reported CI = honest bias-aware CI built
# from the cross-fit twin (see estimate_att_doubletree_averaged). We record BOTH so
# data-analysis can compare naive averaged-tree, twin, and honest-CI coverage.
.est_dt_averaged <- function(X, A, Y, escalate = FALSE) {
  r <- tryCatch(
    doubletree::estimate_att_doubletree_averaged(
      X, A, Y, K = SIM_K, outcome_type = "binary",
      rashomon_bound_multiplier = NULL, auto_tune_intersecting = FALSE,
      escalate_intersection = escalate,
      verbose = FALSE),
    error = function(e) NULL)
  if (is.null(r)) return(.result(NA_real_, NA_real_, converged = FALSE,
                                 intersection_nonempty = FALSE))
  # estimate/std_error/ci are the DISPLAY point + honest bias-aware CI; twin fields
  # carry the valid cross-fit estimator + fidelity diagnostic.
  .result(r$theta, r$sigma, converged = isTRUE(r$converged),
          ci_lower = r$ci_95[1], ci_upper = r$ci_95[2],
          theta_crossfit = r$theta_crossfit, se_crossfit = r$sigma_crossfit,
          delta = r$delta, delta_over_se = r$delta_over_se,
          intersection_nonempty = TRUE,
          rashomon_c_e = r$rashomon_c_e, rashomon_c_m0 = r$rashomon_c_m0)
}

# --- 5. M-split: modal structure, cross-fit predictions (twin of 6) ----------
# As of Phase B (2026-07-15) estimate_att_msplit reports an HONEST bias-aware CI from the
# fully-fold-specific twin (the modal structure is selected across all splits, same non-
# orthogonality as the Rashomon path). r$sigma is now the twin SE; r$sigma_wald is the old
# modal-averaged Wald SE (kept for reference). Pass the honest interval + twin/delta.
.est_msplit <- function(X, A, Y) {
  r <- doubletree::estimate_att_msplit(X, A, Y, M = SIM_M, K = SIM_K_MSPLIT,
                                       verbose = FALSE)
  .result(r$theta, r$sigma, converged = isTRUE(r$converged),
          ci_lower = r$ci_95[1], ci_upper = r$ci_95[2],
          theta_crossfit = r$theta_crossfit, se_crossfit = r$sigma_crossfit,
          delta = r$delta)
}

# --- 6. M-split-averaged: modal structure, averaged leaves, 1 tree -----------
# Like dt_averaged: point estimate = averaged (biased) tree; reported CI = honest
# bias-aware CI from the cross-fit twin (se_delta from the M splits). No Rashomon
# intersection, so escalate does not apply here.
.est_msplit_averaged <- function(X, A, Y) {
  r <- tryCatch(
    doubletree::estimate_att_msplit_averaged(X, A, Y, M = SIM_M, K = SIM_K_MSPLIT,
                                             outcome_type = "binary", verbose = FALSE),
    error = function(e) NULL)
  if (is.null(r)) return(.result(NA_real_, NA_real_, converged = FALSE))
  .result(r$theta, r$sigma, converged = isTRUE(r$converged),
          ci_lower = r$ci_95[1], ci_upper = r$ci_95[2],
          theta_crossfit = r$theta_crossfit, se_crossfit = r$sigma_crossfit,
          delta = r$delta, delta_over_se = r$delta_over_se)
}

# --- 7. Alternative A: honest single tree (intersection struct, all-n leaves) -
# Point estimate = the interpretable single tree (goal i); reported CI = HONEST bias-aware
# interval from the fully-fold-specific twin (Phase B 2026-07-15; the single tree's Wald SE
# undercovers). Twin + delta recorded (goal ii).
.est_single_tree <- function(X, A, Y, escalate = FALSE) {
  r <- tryCatch(
    doubletree::estimate_att_single_tree(
      X, A, Y, K = SIM_K, outcome_type = "binary",
      rashomon_bound_multiplier = NULL, escalate_intersection = escalate,
      inference = "single", verbose = FALSE),
    error = function(e) NULL)
  if (is.null(r)) {
    # Empty intersection (margin fails) -> no single tree; record as non-converged.
    return(.result(NA_real_, NA_real_, converged = FALSE,
                   intersection_nonempty = FALSE))
  }
  # r$theta = single-tree point estimate; r$sigma = twin SE; r$ci_95 = honest interval.
  .result(r$theta, r$sigma, converged = isTRUE(r$converged),
          ci_lower = r$ci_95[1], ci_upper = r$ci_95[2],
          theta_crossfit = r$theta_crossfit, se_crossfit = r$sigma_crossfit,
          delta = r$delta, delta_over_se = r$delta_over_se,
          intersection_nonempty = TRUE,
          rashomon_c_e = r$rashomon_c_e, rashomon_c_m0 = r$rashomon_c_m0)
}

# --- shared EIF ATT (for the full-sample baseline) ---------------------------
.att <- function(Y, A, e_hat, m0_hat) {
  n <- length(Y)
  e_hat <- pmin(pmax(e_hat, 0.01), 0.99)
  pi_hat <- mean(A)
  eta <- list(e = e_hat, m0 = m0_hat, m1 = NULL)
  s0 <- doubletree::psi_att(Y, A, theta = 0, eta = eta, pi_hat = pi_hat)
  theta <- sum(s0) / sum(A / pi_hat)
  sv <- doubletree::psi_att(Y, A, theta = theta, eta = eta, pi_hat = pi_hat)
  list(theta = theta, sigma = doubletree::att_se(sv, n))
}

# NOTE: true_value() is defined in dgp.R (true ATT = 0.15 for all DGPs).
