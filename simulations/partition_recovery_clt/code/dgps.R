# ============================================================
# dgps.R
# Study: partition_recovery_clt  (doubletree, S1)
# Spec:  quality_reports/specs/2026-09-01_partition-recovery-clt.md  (see its §3, §7)
#
# WHAT IS NEW HERE AND WHAT IS NOT
#
# Spec §7 is explicit: "Do not build a new DGP family from scratch. Adapt
# dgps_smooth.R / single_tree_coverage/code/dgps.R's existing simple/moderate
# generators, adding the explicit margin-control parameter ST2 needs (not
# currently exposed) and reusing dgps_stress.R's generate_dgp_weak_overlap()
# unmodified for ST3."
#
# So:
#
#   F1 / F1-weak (`dgp_spec_grid_sparse`) is the existing simple/moderate
#     structure -- p iid Bernoulli(0.5) binary covariates, a propensity on
#     (X1, X2), a control outcome on (X1, X3) so the two nuisances SHARE X1,
#     spare noise coordinates, and an additive-on-the-probability-scale
#     treatment effect (mu_1 = mu_0 + tau) -- with exactly three changes, all
#     required by this study:
#
#       (a) LEAF VALUES ARE SET DIRECTLY instead of through a logit link.
#           `generate_dgp_binary_att()` uses e = plogis(0.5 X1 - 0.3 X2), i.e.
#           an ADDITIVE-in-logit propensity, which on binary covariates takes 4
#           distinct values and so needs the full 4-leaf (X1, X2) grid. Two
#           consequences make it unusable as written for this study: S_e is then
#           a SINGLETON at Lbar = 4 (no redundant refinement is affordable, so
#           the non-singleton-S_j behaviour the recovery indicator exists to
#           handle is never exercised), and the implied margin is tiny -- at
#           coefficient 0.3 the leaf values span only [0.38, 0.55], which is
#           precisely why propensity_loss_choice's fits on that family collapse
#           to 1.0-1.5 leaves. A HIERARCHICAL truth (e depends on X2 only when
#           X1 = 1) is representable in 3 leaves, so S_e is genuinely
#           non-singleton, and setting the values directly is what makes the
#           margin a dial rather than an accident of the link function.
#
#       (b) `e_gap` / `mu_gap` are the MARGIN-CONTROL PARAMETERS spec §3's ST2
#           needs and no existing generator exposes. `*_gap` is the distance
#           between the two leaf values that a tree MUST separate to be
#           sufficient; shrinking it makes one specific non-sufficient
#           competitor a near-tie, which is ST2's stated construction. Delta_j is
#           NOT asserted from this formula -- it is computed exactly by
#           enumerate_sufficient_class.R over the whole class. That distinction
#           earns its keep: the two nuisances' tightest competitors do not even
#           have the same SHAPE. For e (uniform weights) it is a merge inside
#           {X1 = 1, Xk = c} for a pure-noise coordinate k, at mass 1/4. For mu
#           (control weights nu_x = p_x{1 - e_0(x)}) it is instead a merge inside
#           {X1 = 1, X2 = 1} -- the region where the propensity is HIGHEST and
#           control mass is therefore most depleted -- at normalised nu mass
#           0.0625/0.6, which is 1.6x CHEAPER than the same merge on a noise
#           coordinate. Both are re-derived by hand in code/review_checks.R R2
#           and agree with the enumeration to 1e-12; a single two-point formula
#           applied to both nuisances would have gotten mu wrong.
#
#       (c) THE TWO NUISANCES SHARE X1, so the DGP is genuinely CONFOUNDED. This
#           is the fix for the defect diagnosed in
#           quality_reports/2026-09-01_partition-recovery-clt-design-audit.md §2
#           and pinned as item 1 of the spec's own §8 addendum. An earlier state
#           of this file put e_0 on (X1, X2) and mu_0 on the DISJOINT (X3, X4),
#           which makes A independent of Y(0) outright: a difference in means is
#           already unbiased and NO amount of partition-recovery failure in
#           either nuisance can move theta_hat. That is fatal to this study's
#           central metric pairing (spec §4) -- recovery is reported ALONGSIDE
#           coverage precisely to separate "Condition (P) held" from "(P) failed
#           but the bias happened to be small," and the second branch is
#           unreachable by construction on an unconfounded DGP. The 18-cell
#           sweep run on the disjoint version demonstrated exactly that: at
#           ST2 n = 200, rec_both = 0.001 (recovery essentially never happened)
#           yet coverage was 0.942, indistinguishable from nominal.
#
#           Gating mu_0 on X1 -- the coordinate that also gates e_0 -- restores
#           confounding: X1 = 1 raises both the propensity and the control
#           outcome, so the naive difference in means is biased upward and
#           SIMULTANEOUS failure of both trees (the product-bias regime double
#           robustness does NOT protect against) is visible in bias and
#           coverage.
#
#           This fix has an intrinsic cost, not a parameter accident:
#           confounding depletes control mass exactly where the propensity is
#           high, and mu_0's fine split lives inside {X1 = 1}, so nu there is
#           thin and Delta_mu shrinks relative to the disjoint version. The
#           lever that buys it back is e's SPREAD, not mu's -- narrowing e keeps
#           control mass in the high-propensity region. code/calibrate_margin.R
#           computes the whole trade-off surface and is the reason the pinned
#           leaf values in common.R are reproducible rather than mystery
#           constants.
#
#   ST3 (`dgp_spec_weak_overlap`) reuses `generate_dgp_weak_overlap()` from
#     simulations/dgps/dgps_stress.R UNMODIFIED for the sample draw. Only its
#     POPULATION counterpart is written out here, because the enumeration needs
#     closed-form cell probabilities and the shipped function returns a sample.
#     Those two descriptions must agree, so `check_spec_vs_draw()` verifies the
#     population values against a large draw from the shipped function rather
#     than trusting the transcription.
#
# BOUNDED OUTCOME (spec §7). Y is binary in every DGP here, so
# prop:selection-rate's boundedness precondition holds. Asserted in
# `new_dgp_spec()` rather than left as a comment.
#
# theta_0 IS EXACT, NOT A SAMPLE QUANTITY. The existing generators return
# `true_att = mean(p1 - p0)` over the REALISED treated units, which is a random
# variable; coverage against a random target is a different (and weaker)
# statement than coverage for theta_0. Here theta_0 is computed in closed form,
#   theta_0 = sum_x p_x e_0(x) {mu_1(x) - mu_0(x)} / sum_x p_x e_0(x),
# from the population cell probabilities.
# ============================================================

# ---- 1 population specification: F1 / F1-weak ----------------------------

#' Population spec for the grid-exact, hierarchically sparse DGP (F1, F1-weak)
#'
#' @param p Number of binary covariates, >= 5. \code{X1, X2} drive the
#'   propensity, \code{X1, X3} the control outcome (so \code{X1} is SHARED and
#'   the DGP is confounded -- see this file's header note (c)), and
#'   \code{X4..Xp} are pure noise. At least one noise coordinate is required:
#'   without it the solver cannot produce a redundant-cutpoint refinement, and
#'   the study would never exercise the case that distinguishes the
#'   currency-correct recovery indicator from a naive same-partition check.
#' @param e_lo,e_mid,e_gap Propensity leaf values.
#'   \code{e_0 = e_lo} when \code{X1 = 0}; \code{e_mid} when
#'   \code{X1 = 1, X2 = 0}; \code{e_mid + e_gap} when \code{X1 = 1, X2 = 1}.
#'   \code{e_gap} is the margin dial: the X2 split inside the \code{X1 = 1}
#'   branch is the split a sufficient partition must make, and \code{e_gap} is
#'   how much it buys.
#' @param mu_lo,mu_mid,mu_gap Control-outcome leaf values, same shape on
#'   \code{(X1, X3)}: \code{mu_lo} when \code{X1 = 0}, \code{mu_mid} when
#'   \code{X1 = 1, X3 = 0}, \code{mu_mid + mu_gap} when \code{X1 = 1, X3 = 1}.
#'   Sharing the \code{X1} gate with \code{e_0} is what makes \code{A} and
#'   \code{Y(0)} dependent; \code{mu_gap} is mu's own margin dial.
#' @param tau Additive treatment effect on the probability scale.
#' @param leaf_budget Lbar.
#' @param id,label Identifiers carried into the results.
#' @return A DGP spec list (see \code{build_dgp()} for the fields consumers use).
dgp_spec_grid_sparse <- function(p = 5L,
                                 e_lo = 0.20, e_mid = 0.45, e_gap = 0.30,
                                 mu_lo = 0.10, mu_mid = 0.50, mu_gap = 0.40,
                                 tau = 0.08,
                                 leaf_budget = 4L,
                                 id = "f1", label = id) {
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
         "the intended 3-leaf structure collapses and S_j changes shape.",
         call. = FALSE)
  }
  # mu_1 = mu_0 + tau must not need clipping: if it did, tau(x) would stop being
  # constant and theta_0 would no longer equal tau. theta_0 is computed exactly
  # below regardless, but the equality is a useful invariant to keep.
  if (mu_hi + tau > 1) {
    stop("mu_mid + mu_gap + tau = ", signif(mu_hi + tau, 4), " > 1, so mu_1 ",
         "would need clipping and tau(x) would not be constant. Reduce tau or ",
         "mu_gap.", call. = FALSE)
  }

  cells <- make_cell_grid(p)
  # X_k iid Bernoulli(0.5), so every cell has probability 2^-p.
  cell_prob <- rep(1 / nrow(cells), nrow(cells))

  e0 <- ifelse(cells$X1 == 0L, e_lo, ifelse(cells$X2 == 0L, e_mid, e_hi))
  # mu_0 is gated on X1, the SAME coordinate that gates e_0. That shared gate is
  # the confounding (header note (c)); its fine split is X3, mu's own coordinate.
  mu0 <- ifelse(cells$X1 == 0L, mu_lo, ifelse(cells$X3 == 0L, mu_mid, mu_hi))
  mu1 <- mu0 + tau

  new_dgp_spec(
    id = id, label = label, p = p, leaf_budget = leaf_budget,
    cells = cells, cell_prob = cell_prob,
    e0 = e0, mu0 = mu0, mu1 = mu1, tau = tau,
    draw = function(n) draw_grid_sparse(n, p, e_lo, e_mid, e_hi,
                                        mu_lo, mu_mid, mu_hi, tau),
    # Declared coordinate roles, checked by check_variable_roles() on every
    # build_dgp() call. Declaring them makes the intended structure
    # machine-checkable instead of a prose claim in a header comment: an earlier
    # state of this file keyed mu_0 on (X3, X4) -- DISJOINT from e_0's (X1, X2) --
    # on BOTH the population and the sample side, so the two agreed with each
    # other, every fit succeeded, every fitted partition landed inside T_Lbar,
    # and the only symptom was that theta_hat could not respond to nuisance
    # misfit at all, because there was no confounding to misfit. Nothing except
    # an intent-level assertion catches that. X1 appearing in BOTH lists is the
    # machine-readable statement that this DGP is confounded.
    e_vars = c("X1", "X2"),
    mu_vars = c("X1", "X3"),
    params = list(e_lo = e_lo, e_mid = e_mid, e_gap = e_gap, e_hi = e_hi,
                  mu_lo = mu_lo, mu_mid = mu_mid, mu_gap = mu_gap,
                  mu_hi = mu_hi, tau = tau, p = p)
  )
}

#' One sample from the F1 / F1-weak DGP
#'
#' Draw order (X columns, then A, then Y0, Y1) is fixed so a given seed
#' reproduces a given dataset exactly.
draw_grid_sparse <- function(n, p, e_lo, e_mid, e_hi, mu_lo, mu_mid, mu_hi, tau) {
  X <- as.data.frame(
    lapply(seq_len(p), function(k) as.integer(stats::rbinom(n, 1L, 0.5)))
  )
  names(X) <- paste0("X", seq_len(p))
  e <- ifelse(X$X1 == 0L, e_lo, ifelse(X$X2 == 0L, e_mid, e_hi))
  A <- as.integer(stats::rbinom(n, 1L, e))
  # Shared X1 gate: must match dgp_spec_grid_sparse()'s population mu0 exactly.
  # check_spec_vs_draw() verifies that per cell rather than trusting this line.
  m0 <- ifelse(X$X1 == 0L, mu_lo, ifelse(X$X3 == 0L, mu_mid, mu_hi))
  m1 <- m0 + tau
  Y0 <- as.integer(stats::rbinom(n, 1L, m0))
  Y1 <- as.integer(stats::rbinom(n, 1L, m1))
  list(X = X, A = A, Y = A * Y1 + (1L - A) * Y0, e_true = e, mu0_true = m0)
}

# ---- 2 population specification: ST3, weak overlap -----------------------

#' Population spec matching \code{generate_dgp_weak_overlap()}
#'
#' The sample draw delegates to the shipped generator in
#' \code{simulations/dgps/dgps_stress.R} with no modification (spec §7). Only
#' the population side is restated here, because the enumeration needs
#' closed-form cell probabilities. The restatement is verified against the
#' shipped function by \code{check_spec_vs_draw()}, which the pilot runs.
#'
#' Structure, from \code{dgps_stress.R}: \code{p = 4} iid Bernoulli(0.5);
#' \code{e = plogis(-2.5 + 4 X1 + 2.5 X2)};
#' \code{p0 = plogis(-0.3 + 0.5 X1 + 0.4 X2)}; \code{p1 = pmin(p0 + tau, 1)}.
#' Both nuisances are additive-in-logit on \code{(X1, X2)} and therefore take 4
#' distinct values, so at \code{Lbar = 4} each sufficient class is the SINGLETON
#' \{full (X1, X2) grid\} -- deliberately the opposite corner from F1's
#' hierarchical, non-singleton case.
dgp_spec_weak_overlap <- function(tau = 0.10, leaf_budget = 4L,
                                  id = "st3_weak_overlap",
                                  label = "ST3 weak overlap") {
  p <- 4L
  cells <- make_cell_grid(p)
  cell_prob <- rep(1 / nrow(cells), nrow(cells))
  e0 <- stats::plogis(-2.5 + 4.0 * cells$X1 + 2.5 * cells$X2)
  mu0 <- stats::plogis(-0.3 + 0.5 * cells$X1 + 0.4 * cells$X2)
  mu1 <- pmin(mu0 + tau, 1)

  new_dgp_spec(
    id = id, label = label, p = p, leaf_budget = leaf_budget,
    cells = cells, cell_prob = cell_prob,
    e0 = e0, mu0 = mu0, mu1 = mu1, tau = tau,
    draw = function(n) {
      # generate_dgp_weak_overlap() is sourced from dgps_stress.R by common.R and
      # called with seed = NULL so the caller's set.seed() governs, exactly as
      # propensity_loss_choice/code/common.R does it.
      d <- generate_dgp_weak_overlap(n = n, tau = tau, seed = NULL)
      list(X = d$X, A = d$A, Y = d$Y, e_true = d$true_e, mu0_true = d$true_m0)
    },
    # Both nuisances are additive-in-logit on (X1, X2) in the shipped generator,
    # so BOTH declare the same two coordinates and there is no noise coordinate.
    e_vars = c("X1", "X2"),
    mu_vars = c("X1", "X2"),
    params = list(tau = tau, p = p,
                  source = "simulations/dgps/dgps_stress.R::generate_dgp_weak_overlap")
  )
}

# ---- 2b guards: the population spec must be what it claims to be ---------

#' Assert each nuisance depends only on the coordinates it declares
#'
#' Checks, on the exact cell grid, that \code{e_0} is constant across every
#' coordinate outside \code{spec$e_vars} and \code{mu_0} across every coordinate
#' outside \code{spec$mu_vars}, and that each DECLARED coordinate is really used.
#' This is an INTENT-level check: a nuisance keyed on the wrong covariate is
#' still a perfectly self-consistent DGP, produces no error anywhere downstream,
#' and shows up only as a margin that is quietly the wrong size. See the note at
#' \code{dgp_spec_grid_sparse()}'s \code{e_vars} for the actual defect this
#' exists to catch.
#'
#' @return Invisibly \code{TRUE}.
check_variable_roles <- function(spec) {
  cells <- spec$cells
  # Range of gamma0 within groups defined by holding every coordinate BUT `v`
  # fixed. Zero range everywhere <=> gamma0 does not depend on `v`.
  range_flipping <- function(gamma0, v) {
    others <- setdiff(names(cells), v)
    key <- do.call(paste, c(lapply(others, function(o) cells[[o]]), sep = "-"))
    tapply(gamma0, key, function(z) diff(range(z)))
  }
  check_one <- function(gamma0, claimed, nm) {
    for (v in setdiff(names(cells), claimed)) {
      if (any(range_flipping(gamma0, v) > 1e-12)) {
        stop(nm, " depends on ", v, ", which is not in its declared variables (",
             paste(claimed, collapse = ", "), "). Either the DGP formula or the ",
             "declared role is wrong; both change the margin.", call. = FALSE)
      }
    }
    for (v in claimed) {
      if (all(range_flipping(gamma0, v) <= 1e-12)) {
        stop(nm, " does NOT depend on its declared variable ", v,
             "; the declared structure is stale.", call. = FALSE)
      }
    }
    invisible(TRUE)
  }
  check_one(spec$e0, spec$e_vars, "e_0")
  check_one(spec$mu0, spec$mu_vars, "mu_0")
  invisible(TRUE)
}

#' Assert the population spec agrees with the sample draw
#'
#' Draws once at large \code{n} and compares, per covariate cell, the
#' generator's own deterministic \code{e_true} / \code{mu0_true} vectors against
#' the spec's closed-form values. Compares the DETERMINISTIC nuisance vectors,
#' not noisy empirical rates, so the tolerance can be tight and a transcription
#' error cannot hide behind sampling noise.
#'
#' Load-bearing for ST3 specifically, where the sample comes from the shipped
#' \code{generate_dgp_weak_overlap()} and the population side is a separate
#' transcription -- but run for every DGP, because a population/sample
#' disagreement would silently invalidate every margin and every recovery
#' verdict, whatever the DGP.
#'
#' @return Invisibly, a data.frame of per-cell comparisons.
check_spec_vs_draw <- function(spec, n = 20000L, tol = 1e-12) {
  set.seed(20260901L)
  d <- spec$draw(n)
  if (is.null(d$e_true) || is.null(d$mu0_true)) {
    stop("DGP '", spec$id, "' draw() returns no e_true/mu0_true, so the ",
         "population spec cannot be checked against it.", call. = FALSE)
  }
  key_of <- function(df) {
    do.call(paste, c(lapply(names(spec$cells), function(v) df[[v]]), sep = "-"))
  }
  samp_key <- key_of(d$X)
  cell_key <- key_of(spec$cells)
  if (!setequal(unique(samp_key), cell_key)) {
    stop("Sampled covariate patterns do not exhaust the ", nrow(spec$cells),
         " cells of the grid at n = ", n, " for DGP '", spec$id, "'.",
         call. = FALSE)
  }
  idx <- match(cell_key, samp_key)
  cmp <- data.frame(
    dgp = spec$id, cell = cell_key,
    e0_spec = spec$e0, e0_drawn = d$e_true[idx],
    mu0_spec = spec$mu0, mu0_drawn = d$mu0_true[idx],
    stringsAsFactors = FALSE
  )
  bad <- abs(cmp$e0_spec - cmp$e0_drawn) > tol |
    abs(cmp$mu0_spec - cmp$mu0_drawn) > tol
  if (any(bad)) {
    print(cmp[bad, ])
    stop("Population spec disagrees with the sample draw for DGP '", spec$id,
         "' on ", sum(bad), " cell(s).", call. = FALSE)
  }
  invisible(cmp)
}

# ---- 3 shared constructor: exact population quantities -------------------

#' Assemble a DGP spec and its exact population quantities
#'
#' Not called directly by study code; \code{dgp_spec_*()} call it.
new_dgp_spec <- function(id, label, p, leaf_budget, cells, cell_prob,
                         e0, mu0, mu1, tau, draw, e_vars, mu_vars, params) {
  n_cells <- nrow(cells)
  stopifnot(length(e0) == n_cells, length(mu0) == n_cells,
            length(mu1) == n_cells, length(cell_prob) == n_cells)
  stopifnot(all(e_vars %in% names(cells)), all(mu_vars %in% names(cells)))
  if (any(e0 <= 0 | e0 >= 1)) {
    stop("e_0 must be strictly inside (0, 1) for positivity (ass:causal).",
         call. = FALSE)
  }
  if (any(mu0 < 0 | mu0 > 1) || any(mu1 < 0 | mu1 > 1)) {
    stop("mu_0 and mu_1 must lie in [0, 1]: Y is binary, hence BOUNDED, which ",
         "is prop:selection-rate's own precondition (spec §7).", call. = FALSE)
  }
  if (abs(sum(cell_prob) - 1) > 1e-12) {
    stop("cell_prob must sum to 1.", call. = FALSE)
  }

  # theta_0 EXACTLY: the treated-weighted average treatment effect.
  p_treated <- sum(cell_prob * e0)
  theta0 <- sum(cell_prob * e0 * (mu1 - mu0)) / p_treated

  list(
    id = id, label = label, p = p, leaf_budget = as.integer(leaf_budget),
    var_names = names(cells),
    cells = cells, cell_prob = cell_prob,
    e0 = e0, mu0 = mu0, mu1 = mu1,
    tau_x = mu1 - mu0,
    theta0 = theta0,
    p_treated = p_treated,
    # Control-weighted measure nu_x = p_x {1 - e_0(x)}: the weighting under
    # which Pi^nu_tau and R^{(mu)} are defined, and the one the mu tree
    # actually sees, since estimate_att() fits it on X[A == 0, ].
    nu = cell_prob * (1 - e0),
    e0_range = range(e0),
    tau = tau,
    draw = draw,
    e_vars = e_vars,
    mu_vars = mu_vars,
    params = params
  )
}

#' Attach the exact sufficient classes and margins to a DGP spec
#'
#' Runs the (moderately expensive, entirely deterministic) enumeration ONCE per
#' DGP. Every replication then reuses the result, so no replication pays for
#' re-enumerating 356 partitions.
#'
#' \code{check_variable_roles()} runs first, unconditionally on every call: it is
#' cheap, and a DGP whose nuisances key on the wrong covariates yields margins of
#' the wrong size with no other symptom, so it must not be an opt-in check that
#' only the pilot performs.
#'
#' @param spec A \code{dgp_spec_*()} result.
#' @return \code{spec} with \code{enumeration}, \code{class_e}, \code{class_mu}
#'   added.
build_dgp <- function(spec) {
  check_variable_roles(spec)
  en <- enumerate_tree_partitions(spec$p, spec$leaf_budget, spec$var_names)
  spec$enumeration <- en
  spec$class_e <- sufficient_class(spec$e0, spec$cell_prob, en)
  spec$class_mu <- sufficient_class(spec$mu0, spec$nu, en)
  spec
}

#' One-line-per-nuisance summary of a built DGP's population structure
#'
#' Everything a table or a reader needs to interpret that DGP's cells, with no
#' simulation involved. \code{lambda_n} comparisons are added later, per \code{n},
#' by the harness: whether \code{lambda_n < Delta_j} decides whether the
#' PENALISED population optimum is even in \code{S_j}, which is a different
#' failure mechanism from sampling noise and must not be conflated with it.
dgp_population_summary <- function(spec) {
  data.frame(
    dgp = spec$id,
    nuisance = c("e", "mu"),
    theta0 = spec$theta0,
    p = spec$p,
    leaf_budget = spec$leaf_budget,
    card_T = spec$enumeration$n_partitions,
    n_sufficient = c(spec$class_e$n_sufficient, spec$class_mu$n_sufficient),
    min_sufficient_leaves = c(spec$class_e$min_sufficient_leaves,
                              spec$class_mu$min_sufficient_leaves),
    n_distinct_values = c(spec$class_e$n_distinct_values,
                          spec$class_mu$n_distinct_values),
    delta_sq = c(spec$class_e$delta_sq, spec$class_mu$delta_sq),
    delta_kl = c(spec$class_e$delta_kl, spec$class_mu$delta_kl),
    argmin_key_kl = c(spec$class_e$argmin_key_kl, spec$class_mu$argmin_key_kl),
    e0_min = spec$e0_range[[1]],
    e0_max = spec$e0_range[[2]],
    stringsAsFactors = FALSE
  )
}
