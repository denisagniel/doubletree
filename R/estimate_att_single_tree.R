#' Estimate ATT with a single honest tree per nuisance (Alternative A)
#'
#' @description
#' Produces \strong{one interpretable tree per nuisance} (\eqn{e(X)} and
#' \eqn{m_0(X)}) while retaining valid inference, exploiting the structural-margin
#' result (manuscript Theorem "Exact structure recovery"; Corollary "Rashomon
#' tolerance without the intersection trade-off").
#'
#' The estimator selects a tree \emph{structure} honestly via the cross-fold
#' Rashomon intersection (the structure that is near-optimal in every fold), then
#' refits that structure's leaf values \strong{once on all \eqn{n}} observations.
#' Under the margin condition the intersection structure \eqn{\tau^\ast} is
#' recovered exactly on every fold with probability tending to one, so it is
#' asymptotically deterministic; the only remaining in-sample bias is the
#' leaf-value term \eqn{O(L_n / n) = o(n^{-1/2})} (negligible under the DML
#' sparsity condition \eqn{\bar\alpha > s/2}). Hence a single all-\eqn{n} tree per
#' nuisance yields \eqn{\sqrt{n}(\hat\theta - \theta_0) \rightsquigarrow N(0,\sigma^2)}
#' with the ordinary EIF variance --- this is \strong{goal (i)}: a single-tree
#' estimator that is itself valid.
#'
#' The function also returns the \strong{cross-fit twin} (the same structure with
#' fold-specific leaves and out-of-sample predictions, i.e. the Approach-3
#' "doubletree" estimator) and the \strong{fidelity diagnostic}
#' \eqn{\hat\delta = \hat\theta_{\text{single}} - \hat\theta_{\text{cf}}}. This is
#' \strong{goal (ii)}: \eqn{\hat\delta} measures whether the displayed single tree
#' faithfully represents the valid estimator. Note the bias-corrected point
#' estimate \eqn{\hat\theta_{\text{single}} - \hat\delta \equiv \hat\theta_{\text{cf}}};
#' the diagnostic's value is in flagging \eqn{|\hat\delta| \gtrsim \mathrm{SE}}
#' (the single tree is misleading; prefer the cross-fit inference / fold-specific
#' fallback) versus \eqn{|\hat\delta| \ll \mathrm{SE}} (faithful summary).
#'
#' @inheritParams estimate_att_doubletree_averaged
#' @param inference Character. Which point estimate + CI to report in
#'   \code{theta}/\code{ci_95}: \code{"single"} (goal i; the all-\eqn{n} single
#'   tree) or \code{"crossfit"} (goal ii; the valid twin, with the single tree
#'   reported only for display). Default \code{"single"}. Both estimates are always
#'   returned in the result regardless of this choice.
#'
#' @return A list with:
#'   \item{theta, sigma, ci_95}{Point estimate, SE, and 95\% Wald CI for the chosen
#'     \code{inference} target.}
#'   \item{theta_single, sigma_single, ci_95_single}{Single all-\eqn{n} tree
#'     estimator (goal i).}
#'   \item{theta_crossfit, sigma_crossfit, ci_95_crossfit}{Cross-fit twin
#'     (Approach 3) estimator (goal ii).}
#'   \item{delta}{\eqn{\hat\theta_{\text{single}} - \hat\theta_{\text{cf}}}, the
#'     fidelity diagnostic.}
#'   \item{delta_over_se}{\eqn{\hat\delta / \mathrm{SE}_{\text{cf}}}; large absolute
#'     values indicate the single tree misrepresents the valid estimator.}
#'   \item{tree_e, tree_m0}{The single interpretable tree per nuisance (nested
#'     lists), leaf values fit on all \eqn{n} (controls only for \code{tree_m0}).}
#'   \item{converged}{Logical; TRUE if the Rashomon intersection was non-empty for
#'     both nuisances (so a genuine single shared structure exists).}
#'   \item{epsilon_n, n, K}{Tolerance used, sample size, folds.}
#'
#' @details
#' When the intersection is empty for a nuisance (margin fails: near-tied
#' structures, weak overlap), there is no single shared structure and the function
#' errors with guidance to fall back to \code{estimate_att(use_rashomon = FALSE)}
#' (fold-specific trees; valid but not single-tree). This mirrors the honest
#' fallback in the theory: do not force a single tree when the margin does not hold.
#'
#' @seealso \code{\link{estimate_att}} (fold-specific / cross-fit),
#'   \code{\link{estimate_att_doubletree_averaged}} (single averaged tree, Approach 4).
#' @export
estimate_att_single_tree <- function(
  X, A, Y,
  K = 5,
  outcome_type = c("binary", "continuous"),
  cv_regularization = TRUE,
  cv_K = 5,
  regularization = 0.1,
  stratified = TRUE,
  seed = NULL,
  verbose = FALSE,
  rashomon_bound_multiplier = NULL,
  rashomon_bound_adder = 0,
  max_leaves = NULL,
  escalate_intersection = FALSE,
  discretize_method = "quantiles",
  discretize_bins = "adaptive",
  inference = c("single", "crossfit"),
  ...
) {
  outcome_type <- match.arg(outcome_type)
  inference <- match.arg(inference)
  check_att_data(X, A, Y, outcome_type = outcome_type)
  if (is.matrix(X)) X <- as.data.frame(X)
  n <- nrow(X)

  # Resolve tolerance to the theory value log(n)/n unless supplied. When escalating
  # (and no explicit multiplier), keep NULL so fit_nuisances_rashomon runs the c-grid
  # escalation instead of a single fixed tolerance.
  escalating <- isTRUE(escalate_intersection) && is.null(rashomon_bound_multiplier)
  if (is.null(rashomon_bound_multiplier) && !escalating) {
    rashomon_bound_multiplier <- optimaltrees::select_epsilon_n(n)
    if (verbose) {
      message("Using theory epsilon_n = log(n)/n = ",
              signif(rashomon_bound_multiplier, 3))
    }
  }

  fold_indices <- create_folds(n, K, strata = if (stratified) A else NULL, seed = seed)

  # One Rashomon fit yields BOTH the shared structure (for the single tree) and
  # the fold-specific cross-fit refits (for the valid twin).
  if (verbose) message("\n--- Fitting nuisances with Rashomon intersection ---")
  nuisance_fits <- fit_nuisances_rashomon(
    X = X, A = A, Y = Y, fold_indices = fold_indices,
    outcome_type = outcome_type, regularization = regularization,
    cv_regularization = cv_regularization, cv_K = cv_K, verbose = verbose,
    rashomon_bound_multiplier = rashomon_bound_multiplier,
    rashomon_bound_adder = rashomon_bound_adder,
    max_leaves = max_leaves, auto_tune_intersecting = FALSE,
    escalate_intersection = escalate_intersection,
    discretize_method = discretize_method, discretize_bins = discretize_bins,
    ...
  )

  cf_e <- nuisance_fits$cf_e
  cf_m0 <- nuisance_fits$cf_m0

  # A single shared structure requires a non-empty intersection for each nuisance.
  empty_e <- is.null(cf_e) || cf_e@n_intersecting == 0
  empty_m0 <- is.null(cf_m0) || cf_m0@n_intersecting == 0
  if (empty_e || empty_m0) {
    which_empty <- paste(c(if (empty_e) "propensity", if (empty_m0) "outcome"),
                         collapse = " and ")
    eps_label <- if (is.null(rashomon_bound_multiplier)) "the escalated tolerance grid"
                 else signif(rashomon_bound_multiplier, 3)
    stop(
      "Rashomon intersection empty for ", which_empty, " at epsilon_n = ",
      eps_label, ": no single shared structure exists ",
      "(the structural margin does not hold on this data).\n",
      "  Fall back to estimate_att(use_rashomon = FALSE) for a valid (but ",
      "fold-specific, not single-tree) estimate.",
      call. = FALSE
    )
  }

  # ---- Goal (i): single tree per nuisance, leaves refit on ALL n ----
  # The intersecting STRUCTURE is shared across folds; take fold 1's refit as the
  # structure template, then re-estimate its leaves on all n (controls for m0).
  if (verbose) message("\n--- Refitting shared structure on all n (single tree) ---")

  struct_e <- extract_k_trees_from_rashomon(cf_e)[[1]]
  struct_m0 <- extract_k_trees_from_rashomon(cf_m0)[[1]]

  apply_disc <- optimaltrees::apply_discretization

  # Propensity: refit on all n. Discretize X with the propensity's metadata.
  Xb_e <- if (!is.null(cf_e@disc_metadata)) apply_disc(X, cf_e@disc_metadata) else X
  tree_e <- optimaltrees::refit_structure_on_data(struct_e, Xb_e, A,
                                                  allow_partial_leaves = TRUE)

  # Control outcome: refit on control units only (m0 = E[Y | A=0, X]).
  idx0 <- which(A == 0)
  if (length(idx0) == 0) stop("No control units (A=0); cannot fit m0.", call. = FALSE)
  Xb_m0_all <- if (!is.null(cf_m0@disc_metadata)) apply_disc(X, cf_m0@disc_metadata) else X
  tree_m0 <- optimaltrees::refit_structure_on_data(
    struct_m0, Xb_m0_all[idx0, , drop = FALSE], Y[idx0],
    allow_partial_leaves = TRUE
  )

  # Predict all n with the single trees. predict_from_tree returns the leaf
  # prediction directly (P(Y=1|X) for binary, leaf mean for continuous).
  e_single <- optimaltrees::predict_averaged_tree(tree_e, Xb_e)
  e_single <- pmin(pmax(e_single, .PROPENSITY_LOWER_BOUND), .PROPENSITY_UPPER_BOUND)
  m0_single <- optimaltrees::predict_averaged_tree(tree_m0, Xb_m0_all)

  att_single <- .att_from_eta(Y, A, e_single, m0_single, n)

  # ---- Goal (ii): cross-fit twin (Approach 3) from the SAME fit ----
  if (verbose) message("\n--- Cross-fit twin (out-of-sample predictions) ---")
  eta_cf <- get_fold_specific_eta_rashomon(nuisance_fits, X, fold_indices)
  att_cf <- .att_from_eta(Y, A, eta_cf$e, eta_cf$m0, n)

  # Fidelity diagnostic.
  delta <- att_single$theta - att_cf$theta
  delta_over_se <- if (att_cf$sigma > 0) delta / att_cf$sigma else NA_real_

  if (verbose) {
    message(sprintf("\n=== Results ==="))
    message(sprintf("Single-tree ATT : %.4f  (SE %.4f)", att_single$theta, att_single$sigma))
    message(sprintf("Cross-fit ATT   : %.4f  (SE %.4f)", att_cf$theta, att_cf$sigma))
    message(sprintf("delta = single - cf : %.4f  (delta/SE_cf = %.2f)", delta, delta_over_se))
    if (is.finite(delta_over_se) && abs(delta_over_se) > 1) {
      message("  NOTE: |delta| > SE_cf -- the single tree may misrepresent the ",
              "valid estimator; prefer the cross-fit inference.")
    }
  }

  # Tolerance multipliers selected by escalation (epsilon_n = c*log(n)/n per nuisance).
  rashomon_c_e  <- if (is.null(nuisance_fits$rashomon_c_e))  NA_real_ else nuisance_fits$rashomon_c_e
  rashomon_c_m0 <- if (is.null(nuisance_fits$rashomon_c_m0)) NA_real_ else nuisance_fits$rashomon_c_m0
  c_vals <- c(rashomon_c_e, rashomon_c_m0)
  epsilon_n_used <- if (all(is.na(c_vals))) rashomon_bound_multiplier else
    max(c_vals, na.rm = TRUE) * (log(n) / n)

  chosen <- if (inference == "single") att_single else att_cf
  list(
    theta = chosen$theta, sigma = chosen$sigma, ci_95 = chosen$ci_95,
    theta_single = att_single$theta, sigma_single = att_single$sigma,
    ci_95_single = att_single$ci_95,
    theta_crossfit = att_cf$theta, sigma_crossfit = att_cf$sigma,
    ci_95_crossfit = att_cf$ci_95,
    delta = delta, delta_over_se = delta_over_se,
    tree_e = tree_e, tree_m0 = tree_m0,
    converged = TRUE,
    epsilon_n = epsilon_n_used,
    rashomon_c_e = rashomon_c_e,
    rashomon_c_m0 = rashomon_c_m0,
    inference = inference,
    n = n, K = K
  )
}

#' ATT point estimate, SE, and 95\% CI from plugged-in nuisances
#'
#' Thin wrapper kept for the single-tree call sites; delegates to the shared
#' \code{\link{eif_att_solve}} (inference.R). Returns only theta/sigma/ci_95.
#' @noRd
.att_from_eta <- function(Y, A, e_hat, m0_hat, n) {
  res <- eif_att_solve(Y, A, e_hat, m0_hat, n)
  list(theta = res$theta, sigma = res$sigma, ci_95 = res$ci_95)
}
