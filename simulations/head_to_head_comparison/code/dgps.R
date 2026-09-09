# ============================================================
# dgps.R
# Study: head_to_head_comparison  (doubletree, S5)
# Spec:  quality_reports/specs/2026-09-08_head-to-head-comparison.md  (see §3, §4)
#
# WHAT IS NEW HERE AND WHAT IS NOT
#
# Exactly ONE new DGP: DGP-A, the shared-interaction gate-swapped construction
# of spec §3. Everything else this study runs is reused UNMODIFIED by sourcing
# another study's file (see code/common.R):
#
#   DGP-A2 = F1        partition_recovery_clt/code/dgps.R :: dgp_spec_grid_sparse()
#   DGP-B  = eps-dial  honest_inference_sparsity_failure/code/dgps_delta_dial.R
#                        :: dgp_spec_delta_dial()
#   DGP-C  = ST3       partition_recovery_clt/code/dgps.R :: dgp_spec_weak_overlap()
#
# Those three files are READ ONLY for this study. The spec is explicit about
# this ("Reuse, don't reimplement" in §3 twice, and "not by editing
# partition_recovery_clt/code/dgps.R"), and so is the delegation brief.
#
# ------------------------------------------------------------
# WHY DGP-A HAD TO BE WRITTEN AT ALL
# ------------------------------------------------------------
#
# The first version of this study planned to reuse F1 verbatim for claim C1
# (doubletree pays no interpretability-for-efficiency cost under sparsity, while
# a linear working model WITHOUT the right basis is genuinely misspecified). Spec
# §3's "Correction" records what happened when that plan was checked instead of
# assumed: F1's population bilinear remainder under main-effects-only
# misspecification is EXACTLY ZERO (1.7e-18). So is the rejected "complex" DGP's
# (-5e-20). Both would have reported "GLM main-effects: unbiased" for C1 -- a
# null result caused by the DGP, discoverable only after a full implementation
# and a full run.
#
# The mechanism is general, not a calibration slip. A population OLS
# main-effects fit that includes the shared confounder X1 has a residual with
# E[residual | X1] = 0 exactly (OLS orthogonality). If the "other" interacting
# covariate of one nuisance is conditionally independent of the "other"
# interacting covariate of the OTHER nuisance given X1 -- true of F1
# (e_0 on (X1, X2), mu_0 on (X1, X3), X2 ⊥ X3 | X1) -- then the bilinear
# remainder's cross term integrates to exactly zero REGARDLESS OF EFFECT SIZE.
# No amount of turning up e_gap or mu_gap fixes it.
#
# DGP-A defeats that argument by construction: e_0 and mu_0 depend on the SAME
# covariate pair (X1, X2), with the GATE SWAPPED between them --
#
#   e_0's dependence on X2 activates only when X1 = 1
#   mu_0's dependence on X1 activates only when X2 = 1
#
# -- so there is no covariate conditionally independent of the other nuisance's
# "other" covariate to average the cross term away. Both residual functions land
# on the SAME (X1, X2) interaction contrast, and their product does not cancel.
#
#   e_0(X1, X2) = e_lo                    if X1 = 0
#                 e_mid                   if X1 = 1, X2 = 0
#                 e_mid + e_gap           if X1 = 1, X2 = 1
#
#   mu_0(X1, X2) = mu_lo                  if X2 = 0
#                  mu_mid                 if X1 = 0, X2 = 1
#                  mu_mid + mu_gap        if X1 = 1, X2 = 1
#
#   mu_1 = mu_0 + tau ;  X3, X4, X5 pure noise ;  p = 5
#
# Constants are F1's own (e_lo = 0.20, e_mid = 0.45, e_gap = 0.30,
# mu_lo = 0.10, mu_mid = 0.50, mu_gap = 0.40, tau = 0.08) so the two DGPs differ
# in DEPENDENCE STRUCTURE and nothing else -- the magnitudes are held fixed on
# purpose, because "GLM is biased here and not there" must not be readable as an
# effect-size difference.
#
# THE NUMBER THIS FILE EXISTS TO PRODUCE, AND THE ASSERTION THAT GUARDS IT.
# Spec §3 verified by hand and numerically that DGP-A's population bilinear
# remainder under main-effects-only projections is 0.01445 -- NOT zero -- which
# at pi = P(A = 1) = 0.4 implies an asymptotic bias of 0.01445/0.4 = 0.0361 in
# theta-hat, a 45% relative bias against theta_0 = tau = 0.08.
# `verify_main_effects_remainder()` recomputes that from this file's own formulas
# and ERRORS if it comes out near zero. That assertion is the whole reason DGP-A
# exists: a silent ~0 here means the gate-swap was implemented wrong and the
# study would once again report a null result for C1 for a reason that has
# nothing to do with any estimator.
#
# EXACT TREE-SUFFICIENCY. e_0 takes 3 distinct values over the (X1, X2) grid
# (a 3-leaf tree splitting X1 then X2 within X1 = 1 represents it exactly); mu_0
# likewise with the split order reversed. leaf_budget = 4 (the full 2x2 grid)
# removes any dependence on which split order the global optimiser happens to
# prefer, matching F1's own convention. ass:sparsity therefore HOLDS exactly,
# which is what makes R1 a test of C1 rather than of C2.
#
# BOUNDED OUTCOME. Y is binary, so prop:selection-rate's and thm:main's
# boundedness preconditions hold. Asserted in `new_dgp_spec()` (S1's, reused).
#
# theta_0 IS EXACT. tau is additive on the probability scale and never clipped
# (mu_mid + mu_gap + tau = 0.98 < 1), so theta_0 = tau; it is nevertheless
# computed from the closed-form cell probabilities by `new_dgp_spec()` and the
# identity is checked below rather than assumed.
# ============================================================

# ---- 1 population specification: DGP-A -----------------------------------

#' Population spec for the shared-interaction, gate-swapped DGP (DGP-A)
#'
#' Depends on \code{make_cell_grid()} and \code{new_dgp_spec()} from S1's
#' \code{partition_recovery_clt/code/} (sourced by this study's
#' \code{common.R}); it deliberately does not reimplement either.
#'
#' @param p Number of binary covariates, \code{>= 5}. \code{X1} and \code{X2}
#'   drive BOTH nuisances -- that shared pair is the construction (see this
#'   file's header) -- and \code{X3..Xp} are pure noise. At least three noise
#'   coordinates are kept so \code{p = 5} matches F1, DGP-B and this project's
#'   other DGPs: a \code{p} difference between the C1 DGP and the C3 DGP would
#'   otherwise confound the comparison with a dimension change.
#' @param e_lo,e_mid,e_gap Propensity leaf values. \code{e_0 = e_lo} when
#'   \code{X1 = 0} (so \code{X2} is INERT there); \code{e_mid} when
#'   \code{X1 = 1, X2 = 0}; \code{e_mid + e_gap} when \code{X1 = 1, X2 = 1}.
#'   \code{X1} is the gate.
#' @param mu_lo,mu_mid,mu_gap Control-outcome leaf values, with the gate SWAPPED:
#'   \code{mu_0 = mu_lo} when \code{X2 = 0} (so \code{X1} is INERT there);
#'   \code{mu_mid} when \code{X1 = 0, X2 = 1}; \code{mu_mid + mu_gap} when
#'   \code{X1 = 1, X2 = 1}. \code{X2} is the gate.
#' @param tau Additive treatment effect on the probability scale.
#' @param leaf_budget Lbar. 4 by default; see the header note on split order.
#' @param id,label Identifiers carried into the results.
#' @return A DGP spec list; pass to \code{build_dgp()} (S1's).
dgp_spec_shared_interaction <- function(p = 5L,
                                        e_lo = 0.20, e_mid = 0.45, e_gap = 0.30,
                                        mu_lo = 0.10, mu_mid = 0.50, mu_gap = 0.40,
                                        tau = 0.08,
                                        leaf_budget = 4L,
                                        id = "A", label = "DGP-A shared-interaction, gate-swapped") {
  stopifnot(p >= 5L, p == as.integer(p))
  p <- as.integer(p)
  e_hi <- e_mid + e_gap
  mu_hi <- mu_mid + mu_gap
  e_vals <- c(e_lo, e_mid, e_hi)
  mu_vals <- c(mu_lo, mu_mid, mu_hi)
  if (any(e_vals <= 0) || any(e_vals >= 1)) {
    stop("Propensity leaf values must lie strictly inside (0, 1); got ",
         paste(signif(e_vals, 4), collapse = ", "), ".", call. = FALSE)
  }
  if (any(mu_vals <= 0) || any(mu_vals >= 1)) {
    stop("Control-outcome leaf values must lie strictly inside (0, 1); got ",
         paste(signif(mu_vals, 4), collapse = ", "), ".", call. = FALSE)
  }
  if (length(unique(e_vals)) != 3L || length(unique(mu_vals)) != 3L) {
    stop("The three leaf values of each nuisance must be distinct, otherwise ",
         "the intended 3-leaf gate-swapped structure collapses.", call. = FALSE)
  }
  if (mu_hi + tau > 1) {
    stop("mu_mid + mu_gap + tau = ", signif(mu_hi + tau, 4), " > 1, so mu_1 ",
         "would need clipping and tau(x) would not be constant. Reduce tau or ",
         "mu_gap.", call. = FALSE)
  }

  cells <- make_cell_grid(p)
  # X_k iid Bernoulli(0.5), so every cell has probability 2^-p.
  cell_prob <- rep(1 / nrow(cells), nrow(cells))

  e0 <- shared_interaction_e0(cells$X1, cells$X2, e_lo, e_mid, e_hi)
  mu0 <- shared_interaction_mu0(cells$X1, cells$X2, mu_lo, mu_mid, mu_hi)
  mu1 <- mu0 + tau

  spec <- new_dgp_spec(
    id = id, label = label, p = p, leaf_budget = leaf_budget,
    cells = cells, cell_prob = cell_prob,
    e0 = e0, mu0 = mu0, mu1 = mu1, tau = tau,
    draw = function(n) draw_shared_interaction(n, p, e_lo, e_mid, e_hi,
                                               mu_lo, mu_mid, mu_hi, tau),
    # BOTH nuisances declare BOTH X1 and X2. That is the machine-readable
    # statement of the construction: not merely "confounded" (X1 shared, as in
    # F1) but "shares the whole interacting PAIR", which is what defeats the
    # exact-cancellation argument of spec §3's correction. S1's
    # check_variable_roles(), run by build_dgp(), verifies on the exact cell grid
    # that each nuisance depends on exactly these and nothing else -- so a
    # mis-keyed gate cannot pass as a self-consistent DGP with a quietly wrong
    # remainder.
    e_vars = c("X1", "X2"),
    mu_vars = c("X1", "X2"),
    params = list(e_lo = e_lo, e_mid = e_mid, e_gap = e_gap, e_hi = e_hi,
                  mu_lo = mu_lo, mu_mid = mu_mid, mu_gap = mu_gap,
                  mu_hi = mu_hi, tau = tau, p = p,
                  gate_e = "X1", gate_mu = "X2")
  )

  # theta_0 = tau identically (tau additive, mu_1 unclipped). new_dgp_spec()
  # computes theta_0 from the cell probabilities; the identity is CHECKED here
  # rather than assumed, because a clipped mu_1 would break it silently.
  if (abs(spec$theta0 - tau) > 1e-12) {
    stop("theta_0 = ", signif(spec$theta0, 12), " != tau = ", tau,
         "; tau(x) is not constant, so mu_1 must have been clipped.",
         call. = FALSE)
  }
  spec
}

#' DGP-A's propensity, gated on X1
#'
#' Factored out so the population spec and the sample draw cannot drift apart:
#' both call this. (S1's \code{check_spec_vs_draw()} verifies the agreement per
#' cell anyway; sharing the function makes the check redundant rather than
#' load-bearing, which is the safer ordering.)
shared_interaction_e0 <- function(x1, x2, e_lo, e_mid, e_hi) {
  ifelse(x1 == 0L, e_lo, ifelse(x2 == 0L, e_mid, e_hi))
}

#' DGP-A's control outcome, gated on X2 (the SWAP)
shared_interaction_mu0 <- function(x1, x2, mu_lo, mu_mid, mu_hi) {
  ifelse(x2 == 0L, mu_lo, ifelse(x1 == 0L, mu_mid, mu_hi))
}

#' One sample from DGP-A
#'
#' Draw order (X columns, then A, then Y0, Y1) is fixed so a given seed
#' reproduces a given dataset exactly. Returns the DETERMINISTIC nuisance
#' vectors, which \code{att_oracle()} consumes directly and
#' \code{check_spec_vs_draw()} compares against the population spec at 1e-12.
#' \code{Y0} is returned as well so the confounding diagnostic can compute
#' \code{cor(A, Y0)} from data rather than from a comment (S1's lesson).
draw_shared_interaction <- function(n, p, e_lo, e_mid, e_hi,
                                    mu_lo, mu_mid, mu_hi, tau) {
  X <- as.data.frame(
    lapply(seq_len(p), function(k) as.integer(stats::rbinom(n, 1L, 0.5)))
  )
  names(X) <- paste0("X", seq_len(p))
  e <- shared_interaction_e0(X$X1, X$X2, e_lo, e_mid, e_hi)
  A <- as.integer(stats::rbinom(n, 1L, e))
  m0 <- shared_interaction_mu0(X$X1, X$X2, mu_lo, mu_mid, mu_hi)
  m1 <- m0 + tau
  Y0 <- as.integer(stats::rbinom(n, 1L, m0))
  Y1 <- as.integer(stats::rbinom(n, 1L, m1))
  list(X = X, A = A, Y = A * Y1 + (1L - A) * Y0,
       e_true = e, mu0_true = m0, Y0 = Y0, Y1 = Y1)
}

# ---- 2 the main-effects projection and the bilinear remainder -------------

#' Population weighted main-effects-only projection of a nuisance
#'
#' The population OLS fit of \code{gamma0} on \code{(1, X1, ..., Xp)} under cell
#' weights \code{w} -- i.e. what an analyst's main-effects working model
#' converges to, with no interaction terms available. Computed exactly over the
#' \eqn{2^p} cells; no data, no sampling.
#'
#' This is a LINEAR-PROBABILITY projection, not a logistic one, which is what
#' spec §3's verified 0.01445 is computed under and therefore what
#' \code{verify_main_effects_remainder()} must reproduce. The GLM arms
#' (\code{att_linear()}) fit a LOGIT main-effects model, whose population limit
#' differs slightly; the remainder computed here is a DGP DIAGNOSTIC establishing
#' that the misspecification is non-degenerate, not a prediction of any
#' particular arm's realised bias. That distinction is stated because conflating
#' them would invite reading a small numerical discrepancy as an implementation
#' error.
#'
#' @param gamma0 True nuisance value per cell, in canonical cell order.
#' @param cells The cell grid (\code{spec$cells}).
#' @param w Cell weights. \code{spec$cell_prob} for the propensity (fit on all
#'   rows); \code{spec$nu} for the control outcome if the CONTROL-weighted
#'   projection is wanted (what a controls-only fit actually targets).
#' @return Numeric vector of fitted values per cell.
main_effects_projection <- function(gamma0, cells, w) {
  stopifnot(length(gamma0) == nrow(cells), length(w) == nrow(cells))
  stopifnot(all(is.finite(gamma0)), all(is.finite(w)), all(w >= 0), sum(w) > 0)
  d <- cbind(cells, .gamma0 = gamma0)
  fit <- stats::lm(.gamma0 ~ ., data = d, weights = w)
  as.numeric(stats::fitted(fit))
}

#' DGP-A's population bilinear remainder under main-effects-only misspecification
#'
#' Computes prop:bilinear's remainder
#'   \eqn{E[(1 - e_0)\{odds(e_0) - odds(\hat e)\}(\mu_0 - \hat m_0)]}
#' at \eqn{\hat e = e_{lin}}, \eqn{\hat m_0 = \mu_{lin}}, the population
#' main-effects projections. Delegates the arithmetic to S2's
#' \code{population_bias()} (reused, not reimplemented -- spec §3), which also
#' returns the implied \code{theta} limit and cor:width's error norms.
#'
#' Reports the remainder under BOTH weightings of the mu projection:
#'   \code{P}  -- the uniform/population weighting, which is what spec §3's
#'               verified 0.01445 is computed under and what
#'               \code{verify_main_effects_remainder()} asserts on;
#'   \code{nu} -- the control weighting \eqn{\nu_x = P_x\{1 - e_0(x)\}}, which is
#'               what a controls-only outcome model actually targets. Reported as
#'               a sensitivity so the assertion is not mistaken for the only
#'               defensible reading.
#'
#' @param spec A \code{dgp_spec_shared_interaction()} result (built or unbuilt).
#' @param clip Propensity clip forwarded to \code{population_bias()}.
#' @return A list with the two projections and, per weighting, the
#'   \code{population_bias()} result.
main_effects_bilinear_remainder <- function(spec, clip = c(0.01, 0.99)) {
  e_lin <- main_effects_projection(spec$e0, spec$cells, spec$cell_prob)
  mu_lin_P <- main_effects_projection(spec$mu0, spec$cells, spec$cell_prob)
  mu_lin_nu <- main_effects_projection(spec$mu0, spec$cells, spec$nu)

  # A projection can leave (0, 1) in principle; if it does, odds(e_lin) is not
  # defined and the clip below is doing real work rather than being inert. Say so
  # rather than let it pass.
  out_of_range <- range(e_lin)
  list(
    e_lin = e_lin,
    mu_lin_P = mu_lin_P,
    mu_lin_nu = mu_lin_nu,
    e_lin_range = out_of_range,
    e_lin_needs_clip = out_of_range[[1]] <= clip[[1]] || out_of_range[[2]] >= clip[[2]],
    P = population_bias(spec, e_lin, mu_lin_P, clip = clip),
    nu = population_bias(spec, e_lin, mu_lin_nu, clip = clip)
  )
}

#' MANDATORY GUARD: DGP-A's remainder must be the verified NONZERO value
#'
#' Spec §3 verified, by hand and numerically, that DGP-A's population bilinear
#' remainder under main-effects-only projections is \strong{0.01445} -- implying
#' an asymptotic bias of \eqn{0.01445/\pi = 0.0361} at \eqn{\pi = 0.4}, a 45\%
#' relative bias against \eqn{\theta_0 = \tau = 0.08}.
#'
#' THIS CHECK IS THE POINT OF THE DGP. Two DGPs previously proposed for claim C1
#' -- F1 reused verbatim, and \code{single_tree_coverage}'s "complex" DGP -- were
#' each found to have a remainder of EXACTLY ZERO (1.7e-18 and -5e-20), for a
#' structural reason no amount of effect-size tuning repairs (this file's
#' header). Either would have reported "GLM main-effects: unbiased" for C1, a
#' null result caused entirely by the DGP. If the value computed here comes back
#' near zero, the gate-swap is implemented wrong and NOTHING downstream is
#' informative; the study must stop here.
#'
#' @param spec A \code{dgp_spec_shared_interaction()} result.
#' @param expected Verified target from spec §3.
#' @param tol Absolute tolerance against \code{expected}. \code{1e-4} is much
#'   tighter than any plausible transcription error and much looser than
#'   floating-point noise on a 32-cell exact sum.
#' @param floor Minimum magnitude below which the remainder is treated as
#'   "cancelled" and the DGP is rejected outright, with the cancellation
#'   diagnosis in the message rather than a bare tolerance failure.
#' @return Invisibly, a one-row data.frame of the verified quantities.
verify_main_effects_remainder <- function(spec, expected = 0.01445,
                                          tol = 1e-4, floor = 1e-3) {
  rem <- main_effects_bilinear_remainder(spec)
  b <- rem$P$bilinear_form
  bias <- rem$P$bias

  if (abs(b) < floor) {
    stop("DGP-A's population bilinear remainder is ", signif(b, 4),
         ", i.e. CANCELLED to within ", floor, ".\n",
         "This is the exact failure mode spec §3's correction documents: F1 and ",
         "the 'complex' DGP both gave EXACTLY zero here, so a main-effects GLM ",
         "would be reported as unbiased and claim C1 would get a null result ",
         "caused by the DGP, not by any estimator.\n",
         "The gate-swap is not doing what it must: check that e_0 gates on X1 ",
         "(X2 inert when X1 = 0) while mu_0 gates on X2 (X1 inert when X2 = 0), ",
         "and that BOTH nuisances use the SAME pair (X1, X2). Do not proceed.",
         call. = FALSE)
  }
  if (abs(b - expected) > tol) {
    stop("DGP-A's population bilinear remainder is ", signif(b, 6),
         ", but spec §3 verified ", expected, " (tol ", tol, ").\n",
         "The construction or the constants have drifted from the spec. Either ",
         "fix them, or re-derive and update the spec's verified value -- do not ",
         "loosen the tolerance.", call. = FALSE)
  }
  # prop:bilinear's identity: pi * bias == bilinear_form. An algebraic
  # cross-check on population_bias() itself, free to compute.
  if (abs(rem$P$bias - rem$P$bias_from_bilinear) > 1e-10) {
    stop("population_bias(): bias = ", signif(rem$P$bias, 10),
         " but bilinear_form/pi = ", signif(rem$P$bias_from_bilinear, 10),
         "; prop:bilinear's identity fails, so one of the two is miscomputed.",
         call. = FALSE)
  }
  if (rem$e_lin_needs_clip) {
    stop("The main-effects propensity projection leaves the clip range ",
         "(range [", signif(rem$e_lin_range[[1]], 4), ", ",
         signif(rem$e_lin_range[[2]], 4), "]), so the reported remainder is a ",
         "CLIPPED quantity and not comparable to spec §3's value.", call. = FALSE)
  }

  invisible(data.frame(
    dgp = spec$id,
    theta0 = spec$theta0,
    pi_pop = rem$P$pi_pop,
    bilinear_remainder_P = b,
    implied_bias_P = bias,
    relative_bias_P = bias / spec$theta0,
    bilinear_remainder_nu = rem$nu$bilinear_form,
    implied_bias_nu = rem$nu$bias,
    relative_bias_nu = rem$nu$bias / spec$theta0,
    D_w = rem$P$D_w, D_mu = rem$P$D_mu,
    expected = expected, tol = tol,
    stringsAsFactors = FALSE
  ))
}

#' Empirical confounding diagnostics for DGP-A, from an actual draw
#'
#' S1's lesson applied literally: a documented confounding fix that never reached
#' the running code made that study's central metric pairing vacuous, and it was
#' caught only by computing \code{cor(A, Y0)} and
#' \code{E[Y0 | A = 1] - E[Y0 | A = 0]} from data rather than reading comments.
#'
#' Also reports the four conditional means spec §3 verified by hand
#' (\code{E[e_0 | X1] = 0.20 / 0.60}, \code{E[mu_0 | X1] = 0.30 / 0.50}), since
#' those are the statement that X1 moves BOTH nuisances.
#'
#' @param spec A \code{dgp_spec_shared_interaction()} result.
#' @param n Draw size for the empirical half.
#' @param seed Seed for the draw.
#' @return A one-row data.frame.
confounding_check_shared_interaction <- function(spec, n = 200000L,
                                                 seed = 20260908L) {
  set.seed(seed)
  d <- spec$draw(n)
  P <- spec$cell_prob
  e0 <- spec$e0
  mu0 <- spec$mu0
  x1 <- spec$cells$X1
  wmean <- function(v, w) sum(w * v) / sum(w)
  data.frame(
    dgp = spec$id, n = n,
    emp_cor_A_Y0 = stats::cor(d$A, d$Y0),
    emp_EY0_treated = mean(d$Y0[d$A == 1L]),
    emp_EY0_control = mean(d$Y0[d$A == 0L]),
    emp_gap = mean(d$Y0[d$A == 1L]) - mean(d$Y0[d$A == 0L]),
    pop_Ee0_x1_0 = wmean(e0[x1 == 0L], P[x1 == 0L]),
    pop_Ee0_x1_1 = wmean(e0[x1 == 1L], P[x1 == 1L]),
    pop_Emu0_x1_0 = wmean(mu0[x1 == 0L], P[x1 == 0L]),
    pop_Emu0_x1_1 = wmean(mu0[x1 == 1L], P[x1 == 1L]),
    pop_Emu0_treated = sum(P * e0 * mu0) / sum(P * e0),
    pop_Emu0_control = sum(P * (1 - e0) * mu0) / sum(P * (1 - e0)),
    pi_pop = sum(P * e0),
    stringsAsFactors = FALSE
  )
}
