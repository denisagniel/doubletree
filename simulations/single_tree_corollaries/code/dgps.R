# ============================================================
# dgps.R
# Study: single_tree_corollaries  (doubletree)
# Spec:  quality_reports/specs/2026-09-09_single-outcome-tree-corollaries.md  (§3, §4)
#
# WHAT THIS FILE IS
#
# Four population specifications (Regimes A, B, C, E) plus the population
# variance machinery the study's gates are read against. Regime D has NO DGP of
# its own -- spec §3 folds `rem:single-tree-se` into A-C as an extra metric,
# because the naive/correct SE contrast needs exactly the leaf-wise quantities
# already computed for the primary comparison.
#
# EVERY CONSTANT BELOW IS THE SPEC'S, NOT A CALIBRATION CHOICE MADE HERE.
# Spec §4's gates are PRE-COMPUTED numbers (Regime B's V_tau/V = 0.539,
# Regime C's 1.979, Regime C's naive/correct 0.311, the manuscript's own
# 54.2x two-point-leaf ratio). `verify_population_gates()` at the bottom of this
# file recomputes them from these formulas and ERRORS on any drift. A loosened
# constant here silently invalidates a pre-registered gate, which is the one
# failure mode this file exists to make impossible.
#
# ------------------------------------------------------------
# THE TWO POPULATION VARIANCES, AND WHY BOTH ARE NEEDED
# ------------------------------------------------------------
#
# manuscript.tex eq:score, ATT:
#   psi(O; theta, eta) = pi^{-1}[ A{Y - mu(X) - theta} - (1-A) w(X){Y - mu(X)} ]
#   w(x) = e(x) / (1 - e(x)),   pi = P(A = 1)
#
# The treated-arm and control-arm terms are orthogonal (A(1-A) = 0), so
#
#   Var(psi) = pi^{-2} { Tr + C },
#     Tr = E[ e_0(X) { (tau(X) - theta_0)^2 + sigma_1^2 } ]      (treated arm)
#     C  = E[ (1 - e_0(X)) w(X)^2 sigma_0^2(X) ]                 (control arm)
#
# The UNTIED (efficient, two-tree) estimator uses w = e_0/(1 - e_0), giving
#   C_untied = E[ f(e_0) sigma_0^2 ],   f(e) = e^2 / (1 - e).
# The TIED (single-outcome-tree) estimator, at a FIXED tau in S_mu, uses the
# leaf-averaged propensity ebar(l) = E[e_0 | l], giving
#   C_tied = E[ (1 - e_0) w_tie(l)^2 sigma_0^2 ],  w_tie(l) = ebar(l)/(1 - ebar(l)).
#
# Spec §3 reports V = Tr + C_untied and V_tau = Tr + C_tied WITHOUT the pi^{-2}
# factor (it cancels in the ratio the gates are stated on). This file returns
# BOTH: `*_num` for the spec's own numbers and gate comparisons, and `V_full` /
# `V_tau_full` = `*_num` / pi^2 for the empirical-variance-to-analytic ratios of
# spec §5, where the factor does NOT cancel (Var(theta-hat) ~ V_full / n).
#
# Jensen gives C_tied <= C_untied whenever sigma_0^2 is leaf-constant
# (`cor:single-tree-coarsening`, super-efficiency); Regime C is the spec's
# falsifiable counterexample where anti-correlated heteroskedasticity REVERSES it.
#
# ------------------------------------------------------------
# WHERE THE CATE LIVES, AND WHY IT IS NOT A DEVIATION
# ------------------------------------------------------------
#
# Spec §3's preamble gives the treated arm "an additive CATE plus its own noise
# sigma_1^2", and each regime's text pins e_0 and mu_0 to X1 (and X2) with the
# remaining covariate declared PURE NOISE. Those two statements are about
# different objects and are both honoured here: the CATE is carried by X3, so
#
#   * mu_0 = E[Y | A = 0, X] does not depend on X3, and the m0 tree is fit on
#     CONTROLS ONLY, so X3 is invisible to it -- tau in S_mu is untouched;
#   * e_0 does not depend on X3, so tau in S_e (Regime A) is untouched;
#   * partition recovery is therefore not contaminated by the CATE at all.
#
# Putting the CATE on X1 or X2 instead would make theta_0 depend on the same
# coordinate as the propensity and change Tr's factorisation, for no gain.
#
# Tr is then E[e_0]{Var(tau) + sigma_1^2} = pi{Var(tau) + sigma_1^2} exactly,
# because X3 is independent of (X1, X2). Spec §3 pins Tr NUMERICALLY in Regimes
# B (Tr = 0.5 * C_untied ~ 0.81) and C (Tr ~ 0.15) but leaves the split between
# CATE spread and sigma_1^2 free. This file therefore fixes a mild CATE spread
# and SOLVES sigma_1^2 = Tr_target/pi - Var(tau) in closed form, so Tr hits the
# spec's target EXACTLY rather than approximately. sigma_1^2 is a derived
# quantity here, never a hand-typed decimal.
# ============================================================

# ---- 0 population primitives ---------------------------------------------

#' Untied control-arm integrand weight, f(e) = e^2 / (1 - e)
#'
#' Named because it is the function the manuscript's own worked two-point-leaf
#' example (Gate 1) is stated in terms of, and because the reversal in Regime C
#' is a statement about how f interacts with sigma_0^2 within a leaf.
f_untied <- function(e) e^2 / (1 - e)

#' All 2^p binary cells, canonical order (X1 fastest-varying is NOT used --
#' `expand.grid` order is kept so cell tables read the same way everywhere).
#'
#' @param p Number of binary covariates.
#' @return A data.frame with columns X1..Xp of 0/1 integers.
st_cell_grid <- function(p) {
  stopifnot(length(p) == 1L, p >= 1L, p == as.integer(p))
  g <- expand.grid(rep(list(0:1), as.integer(p)), KEEP.OUT.ATTRS = FALSE)
  names(g) <- paste0("X", seq_len(p))
  g[] <- lapply(g, as.integer)
  g
}

# ---- 1 the DGP spec constructor ------------------------------------------

#' Build one population specification for this study
#'
#' Holds the exact cell-level population (probabilities, both nuisances, both
#' conditional variances, the CATE) plus a sample-draw closure, and derives every
#' population quantity the gates and metrics need. Nothing is estimated here.
#'
#' @param id,label Identifiers carried into the results.
#' @param cells Cell grid (\code{st_cell_grid()} output).
#' @param cell_prob Cell probabilities; must sum to 1.
#' @param e0,mu0,sigma0_sq,tau Per-cell true propensity, control outcome,
#'   control-outcome variance, and CATE.
#' @param sigma1_sq Treated-arm noise variance (scalar; constant by design --
#'   spec §3 gives the treated arm "its own noise sigma_1^2", singular).
#' @param leaf_maps Named list of integer vectors, one per candidate partition
#'   tau, each of length \code{nrow(cells)} giving the leaf index of each cell.
#'   Regimes A-C have exactly one (\code{"true"}); Regime E has the two ties.
#' @param leaf_fun Function mapping a design data.frame to the integer leaf index
#'   under \code{leaf_maps$true}; \code{NULL} for Regime E, which has no single
#'   true partition and therefore no oracle arm (spec §3, Regime E).
#' @param leaf_budget Lbar. Spec §2: set EXACTLY to |leaves(tau)|, never larger,
#'   because a larger budget lets the tree over-split and silently exit S_mu.
#' @param draw \code{function(n)} returning a list with X, A, Y and the
#'   deterministic truths (\code{e_true}, \code{mu0_true}, \code{leaf_true}).
#' @param oracle_arm Logical; whether an oracle-tau arm is meaningful.
#' @return A DGP spec list, with population quantities attached.
new_st_dgp <- function(id, label, cells, cell_prob, e0, mu0, sigma0_sq, tau,
                       sigma1_sq, leaf_maps, leaf_fun, leaf_budget, draw,
                       oracle_arm = TRUE) {
  K <- nrow(cells)
  stopifnot(
    length(cell_prob) == K, length(e0) == K, length(mu0) == K,
    length(sigma0_sq) == K, length(tau) == K,
    length(sigma1_sq) == 1L, sigma1_sq > 0,
    all(cell_prob >= 0), abs(sum(cell_prob) - 1) < 1e-12,
    all(e0 > 0), all(e0 < 1), all(sigma0_sq > 0),
    length(leaf_maps) >= 1L, !is.null(names(leaf_maps))
  )
  for (nm in names(leaf_maps)) {
    if (length(leaf_maps[[nm]]) != K) {
      stop("leaf_maps$", nm, " has length ", length(leaf_maps[[nm]]),
           ", expected one leaf index per cell (", K, ").", call. = FALSE)
    }
  }

  pi_pop <- sum(cell_prob * e0)
  # theta_0 = E[tau(X) | A = 1] = E[e_0 tau] / E[e_0]. Computed, never assumed
  # equal to mean(tau) -- that identity only holds when tau is independent of e_0,
  # which is true here by the X3 construction but is exactly the kind of thing
  # that should be a checked consequence rather than a premise (spec §4 Gate 6).
  theta0 <- sum(cell_prob * e0 * tau) / pi_pop

  # Treated-arm term of Var(psi) * pi^2. Identical for the tied and untied
  # estimators (both use the same mu = mu_0 on tau in S_mu), which is why spec §3
  # can pin it as a single number "Tr" per regime.
  Tr <- sum(cell_prob * e0 * ((tau - theta0)^2 + sigma1_sq))
  C_untied <- sum(cell_prob * f_untied(e0) * sigma0_sq)

  tied <- lapply(leaf_maps, function(lf) st_tied_terms(cell_prob, e0, mu0,
                                                       sigma0_sq, lf))
  names(tied) <- names(leaf_maps)

  spec <- list(
    id = id, label = label, p = ncol(cells), leaf_budget = as.integer(leaf_budget),
    cells = cells, cell_prob = cell_prob,
    e0 = e0, mu0 = mu0, sigma0_sq = sigma0_sq, tau = tau, sigma1_sq = sigma1_sq,
    leaf_maps = leaf_maps, leaf_fun = leaf_fun, draw = draw,
    oracle_arm = oracle_arm,
    pi_pop = pi_pop, theta0 = theta0,
    Tr = Tr, C_untied = C_untied,
    # Untied (efficient) variance. `_num` = the spec §3 quantity; `_full` = the
    # thing an empirical variance is actually compared against.
    V_num = Tr + C_untied,
    V_full = (Tr + C_untied) / pi_pop^2,
    tied = tied,
    n_leaves = vapply(leaf_maps, function(lf) length(unique(lf)), integer(1))
  )

  # leaf_budget MUST equal |leaves(tau)| exactly (spec §2). Asserted, because a
  # budget mismatch is the single most likely way this study reports "the tree
  # missed tau" when in fact the harness allowed it to over-split.
  if (length(unique(spec$n_leaves)) != 1L) {
    stop("DGP '", id, "': candidate partitions have differing leaf counts (",
         paste(spec$n_leaves, collapse = ", "), "); leaf_budget cannot be set ",
         "to |leaves(tau)| unambiguously.", call. = FALSE)
  }
  if (spec$leaf_budget != spec$n_leaves[[1L]]) {
    stop("DGP '", id, "': leaf_budget = ", spec$leaf_budget, " but |leaves(tau)| = ",
         spec$n_leaves[[1L]], ". Spec §2 requires them EQUAL -- a larger budget ",
         "lets the tree over-split and silently exit S_mu.", call. = FALSE)
  }
  # tau in S_mu is what makes theta_tau = theta_0 (Gate 6). Verified per
  # candidate partition on the exact cell grid rather than asserted in a comment.
  for (nm in names(leaf_maps)) {
    if (tied[[nm]]$mu_projection_error > 1e-12) {
      stop("DGP '", id, "', partition '", nm, "': mu_0 is NOT leaf-constant ",
           "(weighted projection error ", signif(tied[[nm]]$mu_projection_error, 4),
           "). tau is not in S_mu, so theta_tau != theta_0 and every variance ",
           "claim in this study is off-target for this regime.", call. = FALSE)
    }
  }
  spec
}

#' Control-arm variance term under a TIED partition, plus its S_mu check
#'
#' \code{w_tie(l) = ebar(l)/(1 - ebar(l))} with \code{ebar(l) = E[e_0 | l]} under
#' the CONTROL-free cell weighting P (the propensity leaf value a tied estimator
#' uses is the leaf's treated fraction, whose population limit is E[e_0 | l]).
#'
#' \code{mu_projection_error} is the P-weighted variance of mu_0 WITHIN leaves;
#' zero iff mu_0 is leaf-constant, i.e. iff tau is in S_mu.
#'
#' @return A list of the per-leaf quantities plus the scalar \code{C_tied}.
st_tied_terms <- function(cell_prob, e0, mu0, sigma0_sq, leaf) {
  leaf <- as.integer(leaf)
  P_leaf <- as.numeric(tapply(cell_prob, leaf, sum))
  ebar <- as.numeric(tapply(cell_prob * e0, leaf, sum)) / P_leaf
  w_tie <- ebar / (1 - ebar)
  idx <- match(leaf, sort(unique(leaf)))
  C_tied <- sum(cell_prob * (1 - e0) * w_tie[idx]^2 * sigma0_sq)

  mubar <- as.numeric(tapply(cell_prob * mu0, leaf, sum)) / P_leaf
  mu_projection_error <- sum(cell_prob * (mu0 - mubar[idx])^2)

  list(
    leaf = leaf, P_leaf = P_leaf, ebar = ebar, w_tie = w_tie,
    C_tied = C_tied, mu_projection_error = mu_projection_error
  )
}

#' Every population quantity the gates and the §5 metrics are read against
#'
#' One row per (dgp, partition). \code{ratio_Vtau_V} is spec §4's gate statistic;
#' \code{naive_over_correct} is Regime D's (spec §4 Gate 5), and
#' \code{implied_coverage_naive} turns it into the coverage number
#' \code{rem:single-tree-se} predicts a naive-SE user would realise.
#'
#' @param spec A \code{new_st_dgp()} result.
#' @return A data.frame, one row per candidate partition.
st_population_table <- function(spec) {
  rows <- lapply(names(spec$tied), function(nm) {
    tt <- spec$tied[[nm]]
    Vt_num <- spec$Tr + tt$C_tied
    naive_ratio <- spec$Tr / Vt_num
    data.frame(
      dgp = spec$id, partition = nm,
      pi_pop = spec$pi_pop, theta0 = spec$theta0,
      Tr = spec$Tr, C_untied = spec$C_untied, C_tied = tt$C_tied,
      V_num = spec$V_num, V_tau_num = Vt_num,
      ratio_Vtau_V = Vt_num / spec$V_num,
      V_full = spec$V_full, V_tau_full = Vt_num / spec$pi_pop^2,
      # rem:single-tree-se: the naive SE keeps ONLY the treated-arm term, so the
      # variance it reports is Tr/pi^2 against a truth of (Tr + C_tied)/pi^2.
      naive_over_correct = naive_ratio,
      se_ratio_naive = sqrt(naive_ratio),
      implied_coverage_naive =
        2 * stats::pnorm(stats::qnorm(0.975) * sqrt(naive_ratio)) - 1,
      n_leaves = length(tt$P_leaf),
      mu_projection_error = tt$mu_projection_error,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# ---- 2 the shared sample draw --------------------------------------------

#' One sample from a cell-level population
#'
#' Draw order (X columns in index order, then A, then the two potential outcomes)
#' is fixed so a given seed reproduces a given dataset byte-for-byte. Returns the
#' DETERMINISTIC truths alongside the data: \code{e_true} and \code{mu0_true} for
#' diagnostics, and \code{leaf_true} because the oracle arm of spec §2 is defined
#' as leaf-wise means on the KNOWN leaves and must not re-derive them from a fit.
#'
#' @param n Sample size.
#' @param cells,cell_prob,e0,mu0,sigma0_sq,tau,sigma1_sq Population, as in
#'   \code{new_st_dgp()}.
#' @param leaf Integer leaf index per cell, or \code{NULL} (Regime E).
#' @param px Marginal P(X_k = 1) per covariate, used only to draw X; the cell
#'   probabilities must be the induced product measure (asserted by the callers).
st_draw <- function(n, cells, cell_prob, e0, mu0, sigma0_sq, tau, sigma1_sq,
                    leaf, px) {
  p <- ncol(cells)
  stopifnot(length(px) == p)
  X <- as.data.frame(lapply(seq_len(p), function(k) {
    as.integer(stats::rbinom(n, 1L, px[[k]]))
  }))
  names(X) <- names(cells)

  # Row -> cell index by binary encoding, so the sample and the population share
  # one indexing convention and cannot drift.
  key_cells <- do.call(paste, c(cells, sep = "-"))
  key_rows <- do.call(paste, c(X, sep = "-"))
  ci <- match(key_rows, key_cells)
  if (anyNA(ci)) stop("A drawn row does not match any cell.", call. = FALSE)

  e_i <- e0[ci]
  A <- as.integer(stats::rbinom(n, 1L, e_i))
  mu0_i <- mu0[ci]
  Y0 <- mu0_i + stats::rnorm(n, 0, sqrt(sigma0_sq[ci]))
  Y1 <- mu0_i + tau[ci] + stats::rnorm(n, 0, sqrt(sigma1_sq))

  list(
    X = X, A = A, Y = A * Y1 + (1L - A) * Y0,
    e_true = e_i, mu0_true = mu0_i, tau_true = tau[ci], Y0 = Y0, Y1 = Y1,
    leaf_true = if (is.null(leaf)) NULL else as.integer(leaf[ci])
  )
}

#' sigma_1^2 that makes Tr hit a target EXACTLY
#'
#' Tr = pi{Var(tau) + sigma_1^2} because the CATE rides on the independent
#' coordinate X3 (see this file's header). Solving rather than typing a decimal
#' is what keeps Regime B's Tr = 0.5 * C_untied and Regime C's Tr = 0.15 exact,
#' and hence keeps spec §4's gate values reproducible to machine precision.
st_sigma1_sq_for_Tr <- function(Tr_target, pi_pop, var_tau) {
  out <- Tr_target / pi_pop - var_tau
  if (!is.finite(out) || out <= 0) {
    stop("No positive sigma_1^2 achieves Tr = ", Tr_target, " at pi = ", pi_pop,
         " with Var(tau) = ", var_tau, ".", call. = FALSE)
  }
  out
}

# ---- 3 Regime A: saturated leaves (cor:single-tree-saturated) -------------

#' Regime A -- tau in S_e AND S_mu, so V_tau = V is a CONSTRUCTION IDENTITY
#'
#' Spec §3, Regime A, verbatim: p = 3, X1 the true driver, X2 and X3 pure noise
#' for BOTH nuisances; e_0 in {0.3, 0.7} and mu_0 in {0, 1} across the two leaves
#' of tau = {X1 = 0}, {X1 = 1}; sigma_0^2 = 1 constant; leaf_budget = 2.
#'
#' Because e_0 is leaf-CONSTANT, ebar(l) = e_0(l) and therefore
#' C_tied == C_untied identically -- V_tau/V = 1 to machine precision (Gate 2).
#' That is a check on this DGP, not on any estimator.
#'
#' The CATE (X3, mild) and sigma_1^2 = 1 are free here: spec §3 pins Tr only in
#' Regimes B and C, and Tr cancels exactly from Regime A's ratio anyway.
dgp_regime_A <- function(e_lo = 0.3, e_hi = 0.7, mu_lo = 0, mu_hi = 1,
                         sigma0_sq = 1, tau_lo = 0, tau_hi = 0.5,
                         sigma1_sq = 1, id = "A") {
  cells <- st_cell_grid(3L)
  px <- c(0.5, 0.5, 0.5)
  cell_prob <- rep(1 / nrow(cells), nrow(cells))

  e0 <- ifelse(cells$X1 == 1L, e_hi, e_lo)
  mu0 <- ifelse(cells$X1 == 1L, mu_hi, mu_lo)
  s0 <- rep(sigma0_sq, nrow(cells))
  tau_x <- ifelse(cells$X3 == 1L, tau_hi, tau_lo)
  leaf <- as.integer(cells$X1) + 1L

  new_st_dgp(
    id = id,
    label = "Regime A: saturated leaves (tau in S_e AND S_mu; V_tau/V = 1 exactly)",
    cells = cells, cell_prob = cell_prob,
    e0 = e0, mu0 = mu0, sigma0_sq = s0, tau = tau_x, sigma1_sq = sigma1_sq,
    leaf_maps = list(true = leaf),
    leaf_fun = function(X) as.integer(X$X1) + 1L,
    leaf_budget = 2L,
    draw = function(n) st_draw(n, cells, cell_prob, e0, mu0, s0, tau_x,
                               sigma1_sq, leaf, px)
  )
}

# ---- 4 Regime B: coarsening + homoskedasticity ----------------------------

#' Regime B -- tau in S_mu but NOT S_e, homoskedastic: super-efficiency
#'
#' Spec §3, Regime B, verbatim: same p = 3 layout; mu_0 depends on X1 only (same
#' 2 leaves as A, so tau in S_mu); e_0 varies WITHIN each X1-leaf via X2, taking
#' 0.2 / 0.8 equally likely and independently of X1, so tau is NOT in S_e;
#' sigma_0^2 = 1 everywhere; leaf_budget = 2.
#'
#' Spec §3's verified numbers, reproduced by this construction:
#'   C_untied = E[f(e_0) sigma_0^2]                        = 1.625
#'   C_tied   = E[(1 - e_0) sigma_0^2] (ebar/(1 - ebar))^2 = 0.500
#'   control-arm ratio                                     = 0.308
#'   Tr = 0.5 * C_untied                                   = 0.8125
#'   V = 2.4375, V_tau = 1.3125, V_tau/V                    = 0.539   (Gate 3)
#'
#' ebar = 0.5 in BOTH leaves (e_0 is independent of X1), so w_tie = 1 and
#' C_tied collapses to E[1 - e_0] = 0.5 -- the whole coarsening effect is the
#' Jensen gap between E[f(e_0)] and f(E[e_0]).
dgp_regime_B <- function(e_lo = 0.2, e_hi = 0.8, mu_lo = 0, mu_hi = 1,
                         sigma0_sq = 1, tau_lo = 0, tau_hi = 0.5,
                         Tr_multiple = 0.5, id = "B") {
  cells <- st_cell_grid(3L)
  px <- c(0.5, 0.5, 0.5)
  cell_prob <- rep(1 / nrow(cells), nrow(cells))

  e0 <- ifelse(cells$X2 == 1L, e_hi, e_lo)      # within-leaf variation via X2
  mu0 <- ifelse(cells$X1 == 1L, mu_hi, mu_lo)   # leaf-constant: tau in S_mu
  s0 <- rep(sigma0_sq, nrow(cells))
  tau_x <- ifelse(cells$X3 == 1L, tau_hi, tau_lo)
  leaf <- as.integer(cells$X1) + 1L

  pi_pop <- sum(cell_prob * e0)
  C_untied <- sum(cell_prob * f_untied(e0) * s0)
  # Spec §3: "Tr set via mild CATE heterogeneity at Tr ~ 0.81 (= 0.5 x C_untied)".
  var_tau <- sum(cell_prob * (tau_x - sum(cell_prob * tau_x))^2)
  sigma1_sq <- st_sigma1_sq_for_Tr(Tr_multiple * C_untied, pi_pop, var_tau)

  new_st_dgp(
    id = id,
    label = "Regime B: coarsening + homoskedasticity (super-efficiency, V_tau/V <= 0.7)",
    cells = cells, cell_prob = cell_prob,
    e0 = e0, mu0 = mu0, sigma0_sq = s0, tau = tau_x, sigma1_sq = sigma1_sq,
    leaf_maps = list(true = leaf),
    leaf_fun = function(X) as.integer(X$X1) + 1L,
    leaf_budget = 2L,
    draw = function(n) st_draw(n, cells, cell_prob, e0, mu0, s0, tau_x,
                               sigma1_sq, leaf, px)
  )
}

# ---- 5 Regime C: heteroskedastic reversal (the falsifiable sub-claim) -----

#' Regime C -- anti-correlated heteroskedasticity REVERSES the ordering
#'
#' Spec §3, Regime C, verbatim, including the mass skew that buys the margin:
#'   P(X1 = 1) = 0.85  (the REVERSAL leaf), P(X1 = 0) = 0.15 (the boring leaf)
#'   boring leaf  : e_0 = 0.5 constant, sigma_0^2 = 1     -- pure additive drag
#'   reversal leaf: X2 drives e_0 in {0.1, 0.8} equally likely, PAIRED with
#'                  sigma_0^2 in {1, 0.01} respectively -- high e_0 with LOW
#'                  variance, the anti-correlation the reversal needs.
#'
#' sigma_0^2 = 0.01, NOT 0. Spec §3 is explicit (Oracle's flag): an exact zero
#' degenerates the plug-in sigma-hat^2 and can perturb tree fitting.
#'
#' Spec §3's verified numbers, reproduced by this construction:
#'   C_untied = 0.15*0.5 + 0.85*0.02156 = 0.09332
#'   C_tied   = 0.15*0.5 + 0.85*0.30191 = 0.33162
#'   within the reversal leaf alone, C_tied/C_untied = 14.0x
#'   Tr ~ 0.15  ->  V = 0.243, V_tau = 0.482, V_tau/V ~ 1.98   (Gate 4)
#'
#' FRAGILITY, STATED HERE BECAUSE IT IS A PROPERTY OF THESE CONSTANTS (spec §3):
#' the 1.98 is sensitive to Tr -- at Tr = 0.3 it falls to ~1.68, and a naive
#' 50/50-mass two-leaf design (tried and rejected during spec drafting) cleared
#' ~1.5 only at Tr = 0 and fell BELOW the gate at Tr = 0.1. The 0.85 skew toward
#' the reversal leaf is what buys the margin. `Tr_target` is exposed as an
#' argument precisely so that sensitivity can be re-measured rather than argued.
dgp_regime_C <- function(p_reversal = 0.85, e_boring = 0.5, e_lo = 0.1, e_hi = 0.8,
                         sigma0_sq_boring = 1, sigma0_sq_lo = 1, sigma0_sq_hi = 0.01,
                         mu_lo = 0, mu_hi = 1, tau_lo = 0, tau_hi = 0.2,
                         Tr_target = 0.15, id = "C") {
  if (sigma0_sq_hi <= 0) {
    stop("sigma0_sq_hi must be strictly positive; spec §3 uses 0.01 and is ",
         "explicit that an exact zero degenerates the plug-in sigma-hat^2.",
         call. = FALSE)
  }
  cells <- st_cell_grid(3L)
  px <- c(p_reversal, 0.5, 0.5)
  cell_prob <- ifelse(cells$X1 == 1L, p_reversal, 1 - p_reversal) * 0.25

  in_rev <- cells$X1 == 1L
  e0 <- ifelse(!in_rev, e_boring, ifelse(cells$X2 == 1L, e_hi, e_lo))
  s0 <- ifelse(!in_rev, sigma0_sq_boring,
               ifelse(cells$X2 == 1L, sigma0_sq_hi, sigma0_sq_lo))
  mu0 <- ifelse(in_rev, mu_hi, mu_lo)   # leaf-constant on X1: tau in S_mu
  tau_x <- ifelse(cells$X3 == 1L, tau_hi, tau_lo)
  leaf <- as.integer(cells$X1) + 1L

  pi_pop <- sum(cell_prob * e0)
  var_tau <- sum(cell_prob * (tau_x - sum(cell_prob * tau_x))^2)
  sigma1_sq <- st_sigma1_sq_for_Tr(Tr_target, pi_pop, var_tau)

  new_st_dgp(
    id = id,
    label = "Regime C: heteroskedastic reversal (falsifiable sub-claim, V_tau/V >= 1.5)",
    cells = cells, cell_prob = cell_prob,
    e0 = e0, mu0 = mu0, sigma0_sq = s0, tau = tau_x, sigma1_sq = sigma1_sq,
    leaf_maps = list(true = leaf),
    leaf_fun = function(X) as.integer(X$X1) + 1L,
    leaf_budget = 2L,
    draw = function(n) st_draw(n, cells, cell_prob, e0, mu0, s0, tau_x,
                               sigma1_sq, leaf, px)
  )
}

# ---- 6 Regime E: non-singleton S_mu (the random-index caveat) -------------

#' Regime E -- TWO distinct 3-leaf partitions both exactly represent mu_0
#'
#' Spec §3, Regime E, verbatim. mu_0 has an OR structure over (X1, X2):
#' mu_0 = a at (0,0) and b otherwise; X3 is noise for both nuisances (and carries
#' the CATE); leaf_budget = 3.
#'
#' The tie is GENUINE, not a duplicate-column artefact:
#'   Tree 1 (split X1, then X2 within X1 = 0): {X1=1} = {10, 11}, {00}, {01}
#'   Tree 2 (split X2, then X1 within X2 = 0): {X2=1} = {01, 11}, {00}, {10}
#' Cell 10 pools with 11 under Tree 1 but is its own leaf under Tree 2, and
#' symmetrically for 01 -- DIFFERENT partitions, both with exactly zero
#' population approximation error for mu_0 (asserted in \code{new_st_dgp()}).
#'
#' e_0 is deliberately ASYMMETRIC across the boundary cells so the tie is
#' consequential: e_0(00) = e_0(11) = 0.5, e_0(01) = 0.2, e_0(10) = 0.7. Tree 1's
#' big leaf averages e_0 over {0.7, 0.5}; Tree 2's over {0.2, 0.5} -- different
#' averages, hence different w_tie, hence psi_tau1 != psi_tau2 and no single
#' fixed variance describes the tied estimator (lem:single-tree-linear's caveat).
#'
#' FITTED ARM ONLY, per spec §3: the point is which tau the FITTING realises, so
#' there is no oracle arm and \code{leaf_fun} is NULL.
dgp_regime_E <- function(mu_a = 0, mu_b = 1, sigma0_sq = 1,
                         e_00 = 0.5, e_01 = 0.2, e_10 = 0.7, e_11 = 0.5,
                         tau_lo = 0, tau_hi = 0.5, sigma1_sq = 1, id = "E") {
  cells <- st_cell_grid(3L)
  px <- c(0.5, 0.5, 0.5)
  cell_prob <- rep(1 / nrow(cells), nrow(cells))

  key <- paste0(cells$X1, cells$X2)
  e0 <- c("00" = e_00, "01" = e_01, "10" = e_10, "11" = e_11)[key]
  mu0 <- ifelse(key == "00", mu_a, mu_b)
  s0 <- rep(sigma0_sq, nrow(cells))
  tau_x <- ifelse(cells$X3 == 1L, tau_hi, tau_lo)

  # The two ties, as CELL groupings. Encoded by (X1, X2) key so the topology
  # names in the results mean exactly what this comment says they mean.
  tie1 <- match(ifelse(key %in% c("10", "11"), "big", key), c("big", "00", "01"))
  tie2 <- match(ifelse(key %in% c("01", "11"), "big", key), c("big", "00", "10"))

  new_st_dgp(
    id = id,
    label = "Regime E: non-singleton S_mu, two genuine 3-leaf ties (fitted arm only)",
    cells = cells, cell_prob = cell_prob,
    e0 = unname(e0), mu0 = mu0, sigma0_sq = s0, tau = tau_x,
    sigma1_sq = sigma1_sq,
    leaf_maps = list(tree1 = as.integer(tie1), tree2 = as.integer(tie2)),
    leaf_fun = NULL,
    leaf_budget = 3L,
    draw = function(n) st_draw(n, cells, cell_prob, unname(e0), mu0, s0, tau_x,
                               sigma1_sq, NULL, px),
    oracle_arm = FALSE
  )
}

# ---- 7 the registry ------------------------------------------------------

ST_DGP_REGISTRY <- list(
  A = dgp_regime_A,
  B = dgp_regime_B,
  C = dgp_regime_C,
  E = dgp_regime_E
)

#' Build one regime's DGP by id
#'
#' Regime D is deliberately absent: spec §3 folds \code{rem:single-tree-se} into
#' Regimes A-C as an additional metric rather than giving it a fourth DGP.
#' Asking for it is a design misunderstanding, so it errors with that sentence
#' rather than with "unknown id".
make_st_dgp <- function(id) {
  if (identical(id, "D")) {
    stop("Regime D has no DGP. Spec §3 folds rem:single-tree-se into Regimes ",
         "A-C as the naive-vs-corrected SE metric, because that contrast needs ",
         "exactly the leaf-wise quantities already computed for the primary ",
         "comparison. Read `naive_over_correct` / `cov_naive` in A, B, C.",
         call. = FALSE)
  }
  if (!id %in% names(ST_DGP_REGISTRY)) {
    stop("Unknown regime '", id, "'; expected one of ",
         paste(names(ST_DGP_REGISTRY), collapse = ", "), ".", call. = FALSE)
  }
  ST_DGP_REGISTRY[[id]]()
}

# ---- 8 the pre-code gates, recomputed ------------------------------------

#' MANDATORY GATE: spec §4's pre-computed numbers, recomputed from these formulas
#'
#' Spec §4's gates were run in R BEFORE this file existed and are pre-registered.
#' They are therefore assertions on the DGP constants, not observations about a
#' run: if a constant drifts, the gate a later result is judged against silently
#' changes. Every value below has a spec-stated target and a tolerance tight
#' enough that only a real change trips it.
#'
#' Gate 1 additionally reproduces the MANUSCRIPT's own worked two-point-leaf
#' example (e_0 in {0.1, 0.8}, sigma_0^2 in {1, 0}) -- the one place this study
#' can check its variance machinery against a number the paper already publishes.
#' It uses sigma_0^2 = 0 exactly, because that is what the manuscript uses; the
#' 0.01 substitution is a DGP concern (tree fitting), not an arithmetic one.
#'
#' @param tol_ratio Absolute tolerance on the spec's ratio values. Spec §3/§4
#'   quote them to THREE decimals (0.539, 1.979, 0.311), so half-ulp rounding is
#'   up to 5e-4 and a 1e-3 tolerance is the tightest value that admits the
#'   spec's own rounding without admitting a real drift. (Regime B computes
#'   0.5384615..., which the spec rounds to 0.539 -- a 5.4e-4 gap that is
#'   rounding, not disagreement.)
#' @return Invisibly, a data.frame of gate rows with pass/fail per gate.
verify_population_gates <- function(tol_ratio = 1e-3) {
  rows <- list()
  add <- function(gate, check, value, target, tol, pass) {
    rows[[length(rows) + 1L]] <<- data.frame(
      gate = gate, check = check, value = value, target = target,
      tol = tol, pass = pass, stringsAsFactors = FALSE
    )
  }

  # Gate 1 -- the manuscript's own worked example.
  e2 <- c(0.1, 0.8); s2 <- c(1, 0); P2 <- c(0.5, 0.5)
  Cu2 <- sum(P2 * f_untied(e2) * s2)
  eb2 <- sum(P2 * e2); w2 <- eb2 / (1 - eb2)
  Ct2 <- sum(P2 * (1 - e2) * s2) * w2^2
  add(1, "manuscript C_untied", Cu2, 0.005556, 1e-5, abs(Cu2 - 0.005556) < 1e-5)
  add(1, "manuscript C_tied", Ct2, 0.301240, 1e-5, abs(Ct2 - 0.301240) < 1e-5)
  add(1, "manuscript ratio (x)", Ct2 / Cu2, 54.2, 0.1, abs(Ct2 / Cu2 - 54.2) < 0.1)

  pop <- lapply(c("A", "B", "C", "E"), function(id) st_population_table(make_st_dgp(id)))
  pop <- do.call(rbind, pop)

  rA <- pop$ratio_Vtau_V[pop$dgp == "A"]
  add(2, "Regime A |V_tau/V - 1|", abs(rA - 1), 0, 1e-12, abs(rA - 1) <= 1e-12)

  rB <- pop$ratio_Vtau_V[pop$dgp == "B"]
  add(3, "Regime B V_tau/V (<= 0.7)", rB, 0.539, tol_ratio,
      abs(rB - 0.539) < tol_ratio && rB <= 0.7)

  rC <- pop$ratio_Vtau_V[pop$dgp == "C"]
  add(4, "Regime C V_tau/V (>= 1.5)", rC, 1.979, tol_ratio,
      abs(rC - 1.979) < tol_ratio && rC >= 1.5)

  nC <- pop$naive_over_correct[pop$dgp == "C"]
  add(5, "Regime C naive/correct (<= 0.6)", nC, 0.311, tol_ratio,
      abs(nC - 0.311) < tol_ratio && nC <= 0.6)
  covC <- pop$implied_coverage_naive[pop$dgp == "C"]
  add(5, "Regime C implied naive coverage", covC, 0.726, 1e-3,
      abs(covC - 0.726) < 1e-3)

  # Gate 6 -- theta_tau = theta_0. mu_tau = mu_0 exactly in every regime (already
  # asserted per partition inside new_st_dgp()); this records the largest
  # projection error across all regimes as the single reportable number.
  worst <- max(pop$mu_projection_error)
  add(6, "max_regime |mu_tau - mu_0|^2 (bias guard)", worst, 0, 1e-12,
      worst <= 1e-12)

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  bad <- out[!out$pass, , drop = FALSE]
  if (nrow(bad)) {
    stop("Population gate(s) FAILED -- spec §4's pre-registered numbers no ",
         "longer follow from this file's constants:\n",
         paste(sprintf("  gate %s: %s = %.8g, target %.8g (tol %.3g)",
                       bad$gate, bad$check, bad$value, bad$target, bad$tol),
               collapse = "\n"),
         "\nDo NOT loosen a tolerance or recalibrate a constant to make this ",
         "pass. Spec §4's gates were computed before any simulation code ",
         "existed; a change here invalidates them.", call. = FALSE)
  }
  invisible(out)
}
