# ============================================================================
# DEPRECATED: This file is deprecated as of doubletree 0.0.0.9000
#
# REASON: the manuscript no longer scopes doubletree's validity claims to the
# off-grid, continuum-threshold-recovery setting (formerly "Instantiation 2").
# The companion two-stage tree-fitting algorithm this estimator depends on
# (optimaltrees::fit_twostage()) has an unresolved correctness issue in its
# off-grid search-radius scaling and an unaudited core theorem, and is not
# validated at production leaf budgets.
#
# REPLACEMENT: estimate_att() for grid-exact sparsity, or
# estimate_att_crossfit() when structural sparsity does not plausibly hold.
#
# STATUS: Kept for reference only. Not exported, not tested, not documented
# in the manuscript. Do not use in new code.
# ============================================================================

#' Default minimum leaf mass for the two-stage estimator (internal)
#'
#' @description
#' Resolves \code{m_n} for \code{\link{estimate_att_twostage}} when the
#' caller leaves it \code{NULL}. Unlike \code{\link{estimate_att}}'s
#' \code{m_n = 1L} default -- licensed only by Assumption
#' \code{ass:construct}(d) ("$m_n/n\\to0$ and $\\bar L m_n\\lesssim n$"),
#' both satisfied trivially at \eqn{m_n = 1} -- \code{theory.tex}
#' Proposition \code{prop:boundary-recovery} carries the \strong{full}
#' rate conditions of eq.~\code{eq:stage-a-rates} as a hypothesis,
#' including (R4) \eqn{m_n\to\infty}, \eqn{m_n\lambda_n\to\infty} and (R5)
#' \eqn{m_n/\sqrt n\to\infty}. At \eqn{m_n=1} neither holds and
#' Theorem \code{thm:main-P} (via \code{prop:P-twostage}) does not apply.
#'
#' \eqn{m_n=\lceil\sqrt{n\log\log n}\rceil} satisfies (R5) exactly
#' (\eqn{m_n/\sqrt n=\sqrt{\log\log n}\to\infty}) and is compatible with
#' (R4)'s \eqn{m_n r_n/n\to0} given (R1)'s \eqn{r_n=o(\sqrt{n/\log r_n})}
#' (Remark \code{rem:stage-a-rate-compatibility} works this out for the
#' analogous \eqn{\sqrt n\ll m_n\ll n/r_n} regime). It is Oracle-vetted, not
#' independently proof-audited beyond that cross-check.
#'
#' @param n Sample size the floor is computed for. For the outcome
#'   (control-only) fit, pass the control count, not the full \code{n}.
#' @return Integer minimum leaf mass.
#' @noRd
.twostage_default_m_n <- function(n) {
  # log(log(n)) is undefined/negative for tiny n; the asymptotic regime this
  # formula targets is meaningless below that range anyway, and the
  # leaf-budget feasibility check downstream will reject unworkably small n
  # regardless. Floor log(log(n)) at a small positive constant rather than
  # letting NaN propagate silently.
  loglogn <- log(log(pmax(n, 20)))
  as.integer(ceiling(sqrt(n * pmax(loglogn, 0.1))))
}

#' Default penalty for the two-stage estimator, per nuisance (internal)
#'
#' @description
#' Resolves \code{lambda_n} for one nuisance fit of
#' \code{\link{estimate_att_twostage}} when the caller leaves it
#' \code{NULL}. Two things make this different from
#' \code{\link{estimate_att}}'s single shared \code{lambda_n =
#' log(n)/n} default:
#'
#' \strong{The rate condition is reversed and ladder-relative.} Condition
#' (R3) of eq.~\code{eq:stage-a-rates} requires
#' \eqn{\max(1/r_n,\sqrt{\log(r_n)/n})\ll\lambda_n\ll1} -- \eqn{\lambda_n}
#' must \emph{dominate} the grid-approximation error, the reverse of the
#' fixed-grid intuition (Remark \code{rem:penalty-regime-tension}).
#' \code{fit_twostage()} applies one \code{lambda_n} across its whole
#' \code{m_ladder} of resolutions, not one per rung, so the default must
#' dominate at every rung. (R3)'s lower bound is decreasing in \eqn{r_n},
#' so dominating at the \emph{coarsest} rung, \eqn{r_0=\min(\code{m\_ladder})},
#' dominates at every finer one.
#'
#' \strong{The bound is dimensionless; the solver's objective is not.}
#' \code{optimaltrees::fit_tree()}'s penalised risk is mean-normalised, so
#' \eqn{\lambda_n} competes with \eqn{\mathrm{Var}(y)} directly (a
#' propensity fit on \eqn{y=A\in\{0,1\}} and an outcome fit on continuous
#' \eqn{Y} live on entirely different scales). Multiplying (R3)'s
#' dimensionless bound by \code{stats::var(y)} converts it onto the scale
#' the solver actually minimises.
#'
#' @param y Response the fit is for (\code{A} for the propensity tree,
#'   \code{Y0} for the controls-only outcome tree).
#' @param n Sample size for this fit (full \code{n} for propensity,
#'   control count for the outcome fit).
#' @param m_ladder The same resolution ladder passed to
#'   \code{optimaltrees::fit_twostage()}.
#' @param kappa Positive multiplier. Default \code{1}; do not raise casually
#'   -- \code{fit_twostage()}'s own default \code{lambda_n = 0.1} already
#'   over-penalises a typical propensity fit into a single-leaf tree
#'   (verified: a split reducing \eqn{\mathrm{Var}(A)} by e.g. 0.0225 cannot
#'   clear a 0.1-per-leaf cost), so \code{kappa} should shrink this bound
#'   toward, not away from, its lower edge.
#' @return Numeric penalty, on the response's own variance scale.
#' @noRd
.twostage_default_lambda_n <- function(y, n, m_ladder, kappa = 1) {
  r0 <- min(as.integer(m_ladder))
  kappa * stats::var(y) * max(1 / r0, sqrt(log(r0) / n))
}

#' Warn when a covariate's scale makes the Stage-B search interval
#' meaningless (internal)
#'
#' @description
#' \code{optimaltrees::fit_twostage()} only surfaces its own
#' scale-mismatch warning when \code{verbose = TRUE}
#' (\code{fit_twostage.R} gates it on \code{isTRUE(verbose)}), yet the
#' Stage-B search radius \eqn{\rho_n=M_n/r_n} is computed in \strong{raw}
#' covariate units regardless of \code{verbose}: an off-unit-scale
#' coordinate silently makes Stage B either a no-op (radius too small to
#' move any threshold) or effectively unrestricted (radius wider than the
#' whole support). Because this affects whether Definition
#' \code{def:stageb-interval}'s hypotheses hold at all, this check runs
#' unconditionally, independent of \code{verbose}.
#'
#' @param X Covariate data.frame.
#' @return Invisible \code{NULL}; warns (does not error) on out-of-range
#'   columns.
#' @noRd
.twostage_check_covariate_scale <- function(X) {
  ranges <- vapply(X, function(col) {
    if (is.numeric(col)) diff(range(col)) else NA_real_
  }, numeric(1))
  bad <- names(X)[!is.na(ranges) & (ranges > 10 | ranges < 0.1)]
  if (length(bad) > 0) {
    warning(
      "Column(s) ", paste(bad, collapse = ", "), " have range far from ",
      "unit scale. The Stage-B search interval rho_n = M_n/r_n (Definition ",
      "def:stageb-interval) is computed in RAW covariate units, so an ",
      "off-scale coordinate makes Stage B either a no-op (radius too small ",
      "to move any threshold) or effectively unrestricted. Rescale the ",
      "affected column(s) toward [0, 1], or pass an explicit M_n sized for ",
      "this column's own scale.",
      call. = FALSE, immediate. = TRUE
    )
  }
  invisible(NULL)
}

#' Full-Sample ATT with Two-Stage Off-Grid Trees (Instantiation 2, continuum sparsity)
#'
#' @description
#' Estimates the ATT exactly as specified in \code{inst/paper/theory.tex}
#' \strong{Instantiation 2} ("continuum sparsity and off-grid threshold
#' refinement"): both nuisance trees are fit by \code{optimaltrees}'s
#' two-stage procedure -- an exact grid search (Stage~A) that recovers each
#' tree's \emph{topology}, followed by an exact off-grid scan (Stage~B) that
#' refines each threshold within a shrinking search interval -- rather than
#' the single-stage, pre-specified-cutpoint search
#' \code{\link{estimate_att}} uses. As in \code{estimate_att()}, both trees
#' and all leaf values are computed from all \eqn{n} observations: there is
#' \strong{no sample splitting, no cross-fitting}. The resulting in-sample
#' nuisance predictions are fed to the unchanged efficient influence
#' function (EIF) solver, giving the plain Wald interval whose asymptotic
#' validity is Theorem \code{thm:main-P} / Corollary \code{cor:variance-P},
#' via Proposition \code{prop:P-twostage} (the two-stage guarantees imply
#' Condition~(P)).
#'
#' Use this estimator instead of \code{\link{estimate_att}} when the true
#' nuisances are exactly representable by \emph{some} tree of the declared
#' leaf budget, but at threshold values that need not lie on any
#' pre-specified grid (continuum sparsity, \code{ass:sparsity-jump}) --
#' i.e., when grid-exact sparsity is implausible but tree-shaped structure
#' at continuum thresholds is not. Requires at least one continuous
#' covariate; for purely binary \code{X} there is nothing for Stage~B to
#' refine and \code{estimate_att()} is the right tool.
#'
#' @details
#' \strong{Both nuisances use squared-error loss (\code{ass:construct}(a)),
#' unlike \code{estimate_att()}.} \code{ass:construct}(a) requires squared
#' error for both nuisances so that excess risk equals squared \eqn{L_2}
#' distance to the truth, with loss-to-norm constant \eqn{1} and no
#' positivity condition -- this is a genuine requirement of the theory this
#' function implements, not merely a limitation of
#' \code{optimaltrees::fit_twostage()} (which currently supports
#' \code{loss_function = "squared_error"} only, disclosed as v1,
#' regression-only). \code{estimate_att()}'s \code{log_loss} default for
#' the propensity tree is licensed there by a separate structure-search
#' argument (leaf-value invariance under any proper scoring rule, for a
#' \emph{fixed} partition) that this function does not need or use. There
#' is accordingly no \code{outcome_type}/\code{propensity_loss} argument
#' here: \code{Y} may be binary (0/1) or continuous, always fit by squared
#' error.
#'
#' \strong{Rate conditions (\code{eq:stage-a-rates}).} Unlike
#' \code{estimate_att()}'s \code{m_n = 1L} default -- licensed only by
#' \code{ass:construct}(d) -- \code{prop:boundary-recovery} (and hence
#' \code{prop:P-twostage} / \code{thm:main-P}) carries the \strong{full}
#' joint rate conditions (R1)-(R5) as a hypothesis, including \eqn{m_n\to
#' \infty} and \eqn{m_n/\sqrt n\to\infty}. \code{m_n = NULL} (default)
#' resolves to \eqn{\lceil\sqrt{n\log\log n}\rceil} per nuisance (control
#' count for the outcome fit), which satisfies these; see
#' \code{.twostage_default_m_n}. Likewise \code{lambda_n = NULL} (default)
#' resolves \strong{separately for each nuisance} to a penalty on that
#' response's own variance scale satisfying (R3) at the coarsest resolution
#' in \code{m_ladder}; see \code{.twostage_default_lambda_n}. If you supply
#' \code{lambda_n} yourself, the \emph{same} value is used for both
#' nuisances -- only sensible if you have already accounted for the scale
#' difference between \eqn{\mathrm{Var}(A)} and \eqn{\mathrm{Var}(Y_0)}
#' yourself; when in doubt, leave it \code{NULL}.
#'
#' \strong{Shared vs. per-nuisance tuning.} \code{leaf_budget},
#' \code{depth_budget}, \code{m_ladder}, and \code{M_n} are shared across
#' both nuisance fits -- this is the theory-faithful choice, not merely API
#' parity with \code{estimate_att()}: \code{prop:twostage-superconsistency}'s
#' conclusion and the event \eqn{E_n=E_{e,n}\cap E_{\mu,n}} are stated
#' \emph{jointly} over both nuisances at one \eqn{(r_n,M_n)} sequence, and
#' \code{ass:budget} gives a single fixed \eqn{\bar L} bounding both. Only
#' \code{lambda_n} (and, by default, \code{m_n}) genuinely differ by
#' nuisance, because they are compared against each response's own
#' variance/sample size.
#'
#' \strong{Covariate scale.} The Stage-B search radius
#' \eqn{\rho_n=M_n/r_n} (\code{def:stageb-interval}) is computed in raw
#' covariate units. This function warns (independent of \code{verbose})
#' when any continuous column's range is far from unit scale, since
#' \code{optimaltrees::fit_twostage()} only surfaces its own version of
#' this warning when \code{verbose = TRUE}.
#'
#' \strong{No dedicated group/control-count mechanism is used.} The
#' outcome tree is fit on the control subset directly
#' (\code{X[A == 0, ]}, \code{Y[A == 0]}), exactly as in
#' \code{estimate_att()}, so a plain \code{m_n} floor on that fit is
#' already a control-count floor as \code{ass:solver-twostage}(g)
#' requires. \code{optimaltrees::fit_twostage()}'s \code{group}/
#' \code{group_value}/\code{m_n_group} mechanism (for enforcing a subgroup
#' floor while fitting on the \emph{full} sample) is not used and is
#' blocked from \code{...} below, along with \code{discretize_bins} (owned
#' by \code{m_ladder}), \code{loss_function}, \code{max_depth}, and
#' \code{min_leaf_n} (owned by \code{m_n}), to prevent an inadvertent
#' override that would silently change what is actually being estimated.
#'
#' \strong{No automatic fallback.} Exactly as \code{estimate_att()}
#' documents for itself: this function reports diagnostics
#' (\code{certified_e}, \code{stop_reason_e}, etc.) but never switches
#' estimator on the basis of them.
#'
#' \strong{Wall time.} \code{time_budget} and \code{fit_time_limit} are
#' each applied \emph{per} \code{optimaltrees::fit_twostage()} call, and
#' this function makes two calls (propensity, outcome) -- total wall time
#' can be up to twice either budget.
#'
#' @param X Data.frame or matrix of covariates. Must include at least one
#'   continuous (non-binary) column; for purely binary \code{X} use
#'   \code{\link{estimate_att}} instead.
#' @param A Integer or numeric vector of treatment (0/1).
#' @param Y Numeric vector of outcome. May be binary (0/1) or continuous;
#'   always fit by squared error (see Details).
#' @param leaf_budget Integer, \strong{required, no default}, shared by both
#'   nuisances. The leaf budget \eqn{\bar L} of \code{ass:budget}, exactly
#'   as in \code{\link{estimate_att}}.
#' @param depth_budget \code{NULL} (default), \code{"full"}, or a positive
#'   integer, shared by both nuisances. The depth budget \eqn{\bar D} of
#'   \code{ass:sparsity-jump}(S), which requires
#'   \eqn{\lceil\log_2\bar L\rceil\le\bar D}. \code{NULL} resolves to the
#'   minimal compatible \eqn{\bar D=\lceil\log_2\bar L\rceil}, which
#'   satisfies (S) but is not "full compliance" in
#'   \code{optimaltrees::fit_twostage()}'s own internal sense (its returned
#'   \code{certified_full_class_*} will be \code{FALSE} by design at this
#'   default -- this does not mean the fit failed, see
#'   \code{certified_full_class_e}/\code{certified_full_class_m0} below).
#'   Forwarded unchanged to \code{optimaltrees::fit_twostage()}.
#' @param m_ladder Numeric vector of at least two ascending grid
#'   resolutions, shared by both nuisances. The theory's \eqn{r_n}
#'   (\code{eq:stage-a-rates}), realised here as an escalating ladder:
#'   \code{optimaltrees::fit_twostage()} climbs it until two adjacent
#'   resolutions agree on topology (or the ladder is exhausted). Default
#'   \code{c(16, 32, 64, 128)}, \code{optimaltrees}'s own default.
#'   \strong{Ladder-appropriateness for a given \eqn{n} is the caller's
#'   responsibility}: this function warns when
#'   \eqn{m_n\cdot\max(\code{m\_ladder})/n} is not small, since that is
#'   condition (R4)'s \eqn{m_nr_n/n\to0} visibly failing at the finest
#'   rung.
#' @param lambda_n \code{NULL} (default, resolved separately per nuisance;
#'   see Details) or a single positive numeric applied to both nuisances.
#' @param M_n \code{NULL} (default) or a positive numeric, shared by both
#'   nuisances. The Stage-B interval-scale constant of
#'   \code{def:stageb-interval}; \code{NULL} resolves to
#'   \code{optimaltrees::fit_twostage()}'s own default
#'   (\eqn{\sqrt{\text{resolution}}} per rung), which is a documented
#'   package default, not theory-derived (the theory only requires
#'   \eqn{M_n\to\infty}, \eqn{M_n/r_n\to0}).
#' @param m_n \code{NULL} (default, resolved separately per nuisance from
#'   the full/control sample size; see Details) or a single positive
#'   integer applied to both nuisances.
#' @param time_budget Numeric, seconds, per \code{fit_twostage()} call
#'   (see "Wall time" in Details). Default \code{3600}.
#' @param fit_time_limit Numeric, seconds, per individual tree fit inside
#'   the C++ solver, per \code{fit_twostage()} call. Default \code{600}.
#' @param verbose Logical. Forwarded to \code{optimaltrees}. Default
#'   \code{FALSE}. The covariate-scale check in Details runs regardless of
#'   this setting.
#' @param ... Additional arguments forwarded to
#'   \code{optimaltrees::fit_twostage()} for \emph{both} the \eqn{e} and
#'   \eqn{\mu} fits. \code{discretize_bins}, \code{loss_function},
#'   \code{max_depth}, \code{min_leaf_n}, \code{group}, \code{group_value},
#'   and \code{m_n_group} are reserved (owned by the arguments above) and
#'   will raise an error if supplied here.
#'
#' @return A list with elements:
#'   \item{theta, sigma, ci_95, ci_95_wald, score_values, n}{As in
#'     \code{\link{estimate_att}}: point estimate, plug-in Wald standard
#'     error (\code{eq:vhat}, mildly downward biased -- see
#'     "finite-sample bias order" below), 95\% Wald interval (twice,
#'     for cross-estimator naming parity), per-observation influence
#'     values, and sample size.}
#'   \item{nuisance_fits}{List with \code{e_fit}, \code{m0_fit} (the full
#'     \code{optimaltrees_twostage_fit} objects returned by
#'     \code{optimaltrees::fit_twostage()}) and the in-sample
#'     \code{propensity} (clamped) and \code{outcome_control} predictions.}
#'   \item{leaf_budget, depth_budget, m_ladder, M_n}{Shared tuning inputs,
#'     echoed back.}
#'   \item{lambda_n_e, lambda_n_m0, m_n_e, m_n_m0}{The resolved per-nuisance
#'     tuning values actually used.}
#'   \item{certified_e, certified_m0}{Logical. \code{TRUE} iff that fit's
#'     winning rung is proven exact over the declared problem \emph{and}
#'     topology-stable across two rungs. \code{FALSE} with a non-\code{NULL}
#'     model is the \strong{common} case (it requires two completed
#'     rungs agreeing) -- report it, never auto-switch, exactly as
#'     \code{estimate_att()}'s diagnostics are documented.}
#'   \item{certified_full_class_e, certified_full_class_m0}{Logical.
#'     \code{certified_* && depth_sufficient}; \code{FALSE} by design at
#'     the default \code{depth_budget = NULL} (see that argument).}
#'   \item{stop_reason_e, stop_reason_m0}{Character. Why each fit's ladder
#'     stopped: \code{"topology_stable"} (success) or one of
#'     \code{"stage_a_truncated"}, \code{"stage_a_deadline"},
#'     \code{"stage_b_binary_split"}, \code{"budget_exhausted"},
#'     \code{"ladder_exhausted"}.}
#'   \item{m_used_e, m_used_m0}{The grid resolution (rung of \code{m_ladder})
#'     the winning model actually came from.}
#'   \item{n_leaves_e, n_leaves_m0}{Realised leaf counts of the refined
#'     trees (read from \code{fit_twostage()}'s \code{n_leaves_collapsed}
#'     field -- a legacy name; under this pipeline it is simply the
#'     refined tree's own leaf count, since collapsing is off).}
#'
#' \strong{Finite-sample bias order.} As in \code{estimate_att()}, the
#' in-sample plug-in \code{sigma} is mildly downward biased. Remark after
#' Corollary \code{cor:variance-twostage} gives a \emph{sharper} order here
#' than the generic \eqn{(2\bar L+1)/n}: \eqn{(L_{0,e}+L_{0,\mu}+1)/n} at
#' the \emph{realised} leaf counts. Callers wanting the finite-sample
#' rescaling can use \code{n_leaves_e}/\code{n_leaves_m0} directly:
#' \code{sigma * sqrt(n / (n - (n_leaves_e + n_leaves_m0 + 1)))}.
#'
#' \strong{Reading a failed or over-penalised fit.} If either fit's
#' \code{$model} is \code{NULL} (ladder exhausted without a usable model),
#' this function raises a classed condition \code{doubletree_twostage_no_model}
#' with the failing fit's \code{stop_reason} and full \code{ladder_topologies}
#' log attached, so callers can program against it (e.g. widen
#' \code{m_ladder}, raise \code{time_budget}/\code{fit_time_limit}, or --
#' if \code{stop_reason} is \code{"stage_b_binary_split"} -- switch to
#' \code{estimate_att()}, since that specific refusal means the design is
#' effectively binary). Separately, \eqn{\lambda_n\bar L<\Delta_j/2}'s upper
#' bound on the penalty squeeze (\code{rem:penalty-regime-tension}) is not
#' checkable from a single fit; its observable symptom is
#' \code{budget_slack > 0} together with \code{lambda_binding = FALSE} in
#' the corresponding \code{nuisance_fits$*_fit$ladder_topologies}.
#'
#' @references
#' \code{inst/paper/theory.tex}: Section \code{sec:twostage} (two-stage
#' estimator), Definition \code{cond:P} / \code{def:stageb-interval},
#' Assumption \code{ass:sparsity-jump} (continuum sparsity and jump
#' regularity), Assumption \code{ass:solver-twostage} (two-stage solver
#' fidelity), eq.~\code{eq:stage-a-rates} (rate conditions (R1)-(R5)),
#' Proposition \code{prop:boundary-recovery} (Stage-A topology recovery),
#' Proposition \code{prop:search-interval} (search-interval validity),
#' Proposition \code{prop:twostage-superconsistency} (simultaneous
#' threshold superconsistency), Proposition \code{prop:P-twostage} (the
#' two-stage guarantees imply Condition~(P)), Theorem \code{thm:main-P}
#' (CLT under Condition~(P)), Corollary \code{cor:variance-P} / eq.
#' \code{eq:vhat} (variance estimation).
#'
#' @seealso \code{\link{estimate_att}} for the Instantiation-1 (grid-exact)
#'   flagship estimator; \code{\link{estimate_att_crossfit}} for the
#'   cross-fitting fallback that requires neither form of sparsity.
#'
#' @examples
#' \dontrun{
#' library(optimaltrees)
#' set.seed(42)
#' n <- 600
#' # Continuum sparsity: e and m0 each depend on ONE continuous covariate
#' # through a single off-grid threshold -- not aligned to any particular
#' # pre-specified cutpoint doubletree's own fixed grid would use.
#' X <- data.frame(X1 = runif(n), X2 = runif(n))
#' A <- rbinom(n, 1, ifelse(X$X1 > 0.37, 0.65, 0.35))
#' Y <- rbinom(n, 1, 0.25 + 0.30 * (X$X2 > 0.62) + 0.15 * A)
#'
#' fit <- estimate_att_twostage(X, A, Y, leaf_budget = 2L)
#' fit$theta
#' fit$ci_95
#'
#' # Sparsity-proxy diagnostics (report, never auto-switch):
#' c(fit$certified_e, fit$certified_m0)
#' c(fit$stop_reason_e, fit$stop_reason_m0)
#' }
#' @noRd
estimate_att_twostage <- function(X, A, Y, leaf_budget, depth_budget = NULL,
                                   m_ladder = c(16, 32, 64, 128),
                                   lambda_n = NULL, M_n = NULL, m_n = NULL,
                                   time_budget = 3600, fit_time_limit = 600,
                                   verbose = FALSE, ...) {
  check_att_data(X, A, Y, outcome_type = "continuous")
  if (is.matrix(X)) X <- as.data.frame(X)
  n <- nrow(X)

  # -- Reserved names: owned by named arguments above, not user-overridable
  # via ... without silently changing what is being estimated. --------------
  reserved <- c("discretize_bins", "loss_function", "max_depth", "min_leaf_n",
                "group", "group_value", "m_n_group")
  dots <- list(...)
  clash <- intersect(names(dots), reserved)
  if (length(clash) > 0) {
    stop("estimate_att_twostage() reserves ", paste(clash, collapse = ", "),
         " (owned by m_ladder/m_n/leaf_budget or not used by this ",
         "estimator's design -- see ?estimate_att_twostage Details). ",
         "Do not pass them via ...", call. = FALSE)
  }

  # -- leaf_budget: Lbar, shared, required, no default (as in estimate_att()). -
  if (missing(leaf_budget) || is.null(leaf_budget)) {
    stop("`leaf_budget` is required and has no default. It is the leaf ",
         "budget Lbar of theory.tex Assumption ass:budget, shared by both ",
         "nuisances: a FIXED, analyst-chosen structural parameter, so no ",
         "package-wide default is meaningful. Supply e.g. leaf_budget = 2L.",
         call. = FALSE)
  }
  if (!is.numeric(leaf_budget) || length(leaf_budget) != 1 ||
      is.na(leaf_budget) || leaf_budget < 1 ||
      leaf_budget != as.integer(leaf_budget)) {
    stop("leaf_budget must be a single positive integer, got: ", leaf_budget,
         call. = FALSE)
  }
  leaf_budget <- as.integer(leaf_budget)

  # -- depth_budget: Dbar, shared. NULL/"full"/positive integer, validated
  # lightly here; fit_twostage() re-validates and is authoritative. ----------
  if (!is.null(depth_budget) &&
      !(identical(depth_budget, "full") ||
        (is.numeric(depth_budget) && length(depth_budget) == 1 &&
         !is.na(depth_budget) && depth_budget >= 1 &&
         depth_budget == as.integer(depth_budget)))) {
    stop("depth_budget must be NULL, \"full\", or a single positive ",
         "integer, got: ", depth_budget, call. = FALSE)
  }

  # -- m_ladder: shared. At least two ascending resolutions. -----------------
  if (!is.numeric(m_ladder) || length(m_ladder) < 2 || any(is.na(m_ladder)) ||
      any(m_ladder < 1)) {
    stop("m_ladder must be a numeric vector of at least two positive grid ",
         "resolutions (theory.tex's r_n ladder), got: ",
         paste(m_ladder, collapse = ", "), call. = FALSE)
  }

  # -- M_n: shared. NULL, or a single positive numeric. -----------------------
  if (!is.null(M_n) && (!is.numeric(M_n) || length(M_n) != 1 ||
                        is.na(M_n) || M_n <= 0)) {
    stop("M_n must be NULL or a single positive number, got: ", M_n,
         call. = FALSE)
  }

  # -- m_n: NULL (per-nuisance default) or a single positive integer applied
  # to both. ------------------------------------------------------------------
  if (!is.null(m_n) && (!is.numeric(m_n) || length(m_n) != 1 ||
                        is.na(m_n) || m_n < 1 ||
                        m_n != as.integer(m_n))) {
    stop("m_n must be NULL or a single positive integer, got: ", m_n,
         call. = FALSE)
  }

  # -- lambda_n: NULL (per-nuisance default) or a single positive number
  # applied to both. -----------------------------------------------------------
  if (!is.null(lambda_n) && (!is.numeric(lambda_n) || length(lambda_n) != 1 ||
                             !is.finite(lambda_n) || lambda_n <= 0)) {
    stop("lambda_n must be NULL or a single positive finite number, got: ",
         lambda_n, call. = FALSE)
  }

  # -- Requires at least one continuous covariate (see roxygen Description). -
  is_binary_col <- vapply(X, function(col) {
    (is.numeric(col) || is.logical(col)) && all(col %in% c(0, 1))
  }, logical(1))
  if (all(is_binary_col)) {
    stop("estimate_att_twostage() requires at least one continuous ",
         "(non-binary) covariate -- with purely binary X there is nothing ",
         "for Stage B's off-grid refinement to do. Use estimate_att(), ",
         "which is exactly the right tool for grid-exact sparsity on ",
         "binary covariates.", call. = FALSE)
  }
  .twostage_check_covariate_scale(X)

  # -- Sufficient treated/control mass, mirroring estimate_att(). ------------
  n_treated <- sum(A == 1)
  n_control <- sum(A == 0)
  if (n_treated < 1) {
    stop("No treated units (A = 1); the ATT is not defined.", call. = FALSE)
  }
  if (n_control < 1) {
    stop("No control units (A = 0); the control-outcome tree cannot be fit.",
         call. = FALSE)
  }

  # == Per-nuisance rate-condition defaults (theory.tex eq:stage-a-rates). ===
  # Propensity fit uses the full sample; outcome fit uses controls only, so
  # its n-dependent defaults are computed from n_control, not n.
  m_n_e  <- if (is.null(m_n)) .twostage_default_m_n(n) else as.integer(m_n)
  m_n_m0 <- if (is.null(m_n)) .twostage_default_m_n(n_control) else as.integer(m_n)
  lambda_n_e  <- if (is.null(lambda_n)) {
    .twostage_default_lambda_n(A, n, m_ladder)
  } else lambda_n
  lambda_n_m0 <- if (is.null(lambda_n)) {
    .twostage_default_lambda_n(Y[A == 0], n_control, m_ladder)
  } else lambda_n

  if (n_control < leaf_budget * m_n_m0) {
    stop("Insufficient control units for leaf_budget = ", leaf_budget,
         " at the resolved outcome-fit m_n = ", m_n_m0, " (theory.tex's ",
         "rate condition (R5), m_n/sqrt(n_control) -> infinity, not the ",
         "trivial m_n = 1 floor estimate_att() uses). The control-outcome ",
         "tree is fit on controls only, so this requires at least ",
         "leaf_budget * m_n = ", leaf_budget * m_n_m0, " control units, ",
         "got: ", n_control, ". Reduce leaf_budget, or collect more data; ",
         "m_n cannot safely be lowered without leaving the rate regime ",
         "Theorem thm:main-P requires.", call. = FALSE)
  }
  if (n < leaf_budget * m_n_e) {
    stop("Insufficient observations for leaf_budget = ", leaf_budget,
         " at the resolved propensity-fit m_n = ", m_n_e, ": requires at ",
         "least ", leaf_budget * m_n_e, " observations, got: ", n, ".",
         call. = FALSE)
  }
  r_max <- max(m_ladder)
  if (m_n_e * r_max / n > 0.5 || m_n_m0 * r_max / n_control > 0.5) {
    warning(
      "m_n * max(m_ladder) / n is not small (propensity: ",
      signif(m_n_e * r_max / n, 3), "; outcome: ",
      signif(m_n_m0 * r_max / n_control, 3), "). This is condition (R4)'s ",
      "m_n*r_n/n -> 0 visibly failing at the finest rung of m_ladder for ",
      "this sample size. Consider a coarser (smaller-max) m_ladder, or a ",
      "larger n, before trusting the finest rung's fit.",
      call. = FALSE, immediate. = TRUE
    )
  }

  pi_hat <- mean(A)
  if (pi_hat <= 0 || pi_hat >= 1) {
    stop("Invalid treatment proportion: pi_hat = ", pi_hat,
         ". This should not happen after sample size validation.",
         call. = FALSE)
  }

  # == Two-stage nuisance fits: BOTH on all n observations (no splitting). ===
  common_args <- list(leaf_budget = leaf_budget, depth_budget = depth_budget,
                      m_ladder = m_ladder, M_n = M_n,
                      time_budget = time_budget,
                      fit_time_limit = fit_time_limit, verbose = verbose)

  e_fit <- do.call(optimaltrees::fit_twostage, c(
    list(X = X, y = A, lambda_n = lambda_n_e, m_n = m_n_e), common_args, dots
  ))
  .twostage_stop_if_no_model(e_fit, which = "propensity")

  X0 <- X[A == 0, , drop = FALSE]
  Y0 <- Y[A == 0]
  m0_fit <- do.call(optimaltrees::fit_twostage, c(
    list(X = X0, y = Y0, lambda_n = lambda_n_m0, m_n = m_n_m0), common_args, dots
  ))
  .twostage_stop_if_no_model(m0_fit, which = "outcome")

  # == In-sample predictions for ALL n rows, exactly as estimate_att(). ======
  e_hat_raw <- predict(e_fit$model, X)
  e_hat <- pmax(.PROPENSITY_LOWER_BOUND, pmin(.PROPENSITY_UPPER_BOUND, e_hat_raw))
  m0_hat <- predict(m0_fit$model, X)

  # == Shared EIF solve (cor:variance-P / eq:vhat), reused UNCHANGED. ========
  .att <- eif_att_solve(Y, A, e_hat, m0_hat, n)

  nuisance_fits <- list(
    e_fit = e_fit,
    m0_fit = m0_fit,
    propensity = e_hat,
    outcome_control = m0_hat
  )

  list(
    theta = .att$theta,
    sigma = .att$sigma,
    ci_95 = .att$ci_95,
    ci_95_wald = .att$ci_95,
    score_values = .att$score_values,
    nuisance_fits = nuisance_fits,
    n = n,
    leaf_budget = leaf_budget,
    depth_budget = depth_budget,
    m_ladder = m_ladder,
    M_n = M_n,
    lambda_n_e = lambda_n_e,
    lambda_n_m0 = lambda_n_m0,
    m_n_e = m_n_e,
    m_n_m0 = m_n_m0,
    certified_e = e_fit$certified,
    certified_m0 = m0_fit$certified,
    certified_full_class_e = e_fit$certified_full_class,
    certified_full_class_m0 = m0_fit$certified_full_class,
    stop_reason_e = e_fit$stop_reason,
    stop_reason_m0 = m0_fit$stop_reason,
    m_used_e = e_fit$m_used,
    m_used_m0 = m0_fit$m_used,
    n_leaves_e = e_fit$n_leaves_collapsed,
    n_leaves_m0 = m0_fit$n_leaves_collapsed
  )
}

#' Raise a classed condition when a two-stage fit produced no usable model (internal)
#'
#' @description
#' \code{optimaltrees::fit_twostage()}'s own refusal conditions
#' (\code{optimaltrees_refine_infeasible}, \code{optimaltrees_stage_b_binary_split})
#' are caught internally by its own ladder loop and never propagate to the
#' caller -- only a \code{NULL} \code{$model} (ladder exhausted without a
#' winning rung) ever surfaces here. Raised as a classed condition, with the
#' full \code{ladder_topologies} log attached, so callers can program
#' against \code{stop_reason} rather than parsing a message string.
#'
#' @param fit The \code{optimaltrees_twostage_fit} object to check.
#' @param which Character, \code{"propensity"} or \code{"outcome"}, for the
#'   error message.
#' @return Invisible \code{NULL} if \code{fit$model} is not \code{NULL}.
#' @noRd
.twostage_stop_if_no_model <- function(fit, which) {
  if (!is.null(fit$model)) return(invisible(NULL))
  advice <- switch(fit$stop_reason,
    "budget_exhausted" = "raise time_budget/fit_time_limit",
    "stage_a_deadline" = "raise time_budget/fit_time_limit",
    "stage_b_binary_split" = paste(
      "Stage A repeatedly split on a binary passthrough column -- this ",
      "design is effectively binary; use estimate_att() instead"),
    "ladder_exhausted" = "widen or coarsen m_ladder",
    "no advice available for this stop_reason"
  )
  cond <- structure(
    class = c("doubletree_twostage_no_model", "error", "condition"),
    list(
      message = paste0(
        "estimate_att_twostage()'s ", which, " fit produced no usable ",
        "model (fit_twostage()$model is NULL). stop_reason = \"",
        fit$stop_reason, "\". Suggested next step: ", advice, ". Full ",
        "per-rung log is attached as the condition's `ladder_topologies` ",
        "field."
      ),
      call = sys.call(-1),
      ladder_topologies = fit$ladder_topologies,
      stop_reason = fit$stop_reason
    )
  )
  stop(cond)
}
