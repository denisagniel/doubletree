# cor:single-tree -- the single-outcome-tree identity.
#
# Spec: quality_reports/specs/2026-09-09_single-outcome-tree-corollaries.md §1.
# Study: simulations/single_tree_corollaries/ (this file is that study's §1 action
# item, deliberately a UNIT TEST and not a Monte Carlo regime).
#
# THE CLAIM. `cor:single-tree` says doubletree's two-tree ATT estimator, with the
# propensity tree FORCED EQUAL to the outcome tree, collapses to the one-line
# plug-in
#
#     theta_tie = mean( Y[A == 1] - mu_hat[A == 1] ).
#
# It is an ALGEBRAIC identity, true pointwise for any leaf-wise-constant e_hat,
# including at finite n -- not an asymptotic statement. Hence a numerical-agreement
# test at ~1e-10 rather than a simulation regime.
#
# WHY IT IS TRUE, so a failure can be localised rather than merely observed.
# eq:score's control-arm term is  -n_1^{-1} sum_l w_hat(l) sum_{A = 0, l} {Y - mu_hat(l)}.
# When mu_hat(l) is the within-leaf CONTROL mean, the inner sum is exactly zero in
# every leaf, so the whole term vanishes for ANY w_hat that is constant within
# those leaves -- leaving the treated-arm mean above. The identity therefore has
# exactly TWO preconditions, each asserted separately below so a failure names its
# own cause:
#   (P1) e_hat is leaf-constant on mu_hat's OWN partition;
#   (P2) mu_hat(l) is the leaf's CONTROL MEAN.
#
# Note what is NOT a precondition, because it is easy to over-assume and it is the
# reason spec §1 says "for any leaf-wise-constant e_hat": the VALUES of e_hat are
# irrelevant. Setting e_hat(l) = n_1(l)/n(l) is what makes the estimator "the
# outcome tree used as the propensity tree" (the object cor:single-tree is about),
# but the identity holds for any leaf-constant e_hat on the same partition. That is
# demonstrated below, and it is also why the clip of ass:construct(c) cannot break
# the identity -- clipping preserves leaf-constancy. The clip is nevertheless
# asserted inert, because a binding clip WOULD mean e_hat is no longer the forced
# propensity tree and the estimator under test would be a different object.
#
# A MEASURED PACKAGE LIMIT, TESTED RATHER THAN WORKED AROUND. optimaltrees stores
# a leaf's fitted value to six significant decimal digits, so
# estimate_att()'s own `outcome_control` vector satisfies (P3) only to ~1e-6 and
# the identity gap it produces is ~1e-7 -- not 1e-10, for a reason unrelated to
# cor:single-tree. The identity is therefore tested at 1e-10 on the fitted
# PARTITION with leaf values recomputed as exact control means (spec §1's own
# description of mu_hat: "fit on control-only data at a fixed leaf_budget,
# leaf-constant"), and the serialised-prediction gap is tested SEPARATELY at a
# documented 1e-5 so a regression in that rounding surfaces here instead of
# silently loosening the identity.

skip_if_not_installed("optimaltrees")

# ---- helpers (local to this file; the study's harness is not on the test path) --

# Leaf index of every row under a fitted estimate_att() outcome tree. Uses the
# package's own .tree_leaf_paths() -- the identity test must exercise the SAME
# leaf-extraction path the study does, so that a bug there fails here.
tie_fitted_leaf <- function(fit, X) {
  model <- fit$nuisance_fits$m0_model
  meta <- model@discretization_metadata
  Xb <- if (!is.null(meta)) optimaltrees::apply_discretization(X, meta) else X
  paths <- doubletree:::.tree_leaf_paths(model@trees[[1L]], Xb)
  match(paths, sort(unique(paths)))
}

# Everything the identity needs, computed once per dataset.
tie_pieces <- function(X, A, Y, leaf_budget) {
  fit <- doubletree::estimate_att(
    X = X, A = A, Y = Y, leaf_budget = leaf_budget,
    outcome_type = "continuous", verbose = FALSE
  )
  leaf <- tie_fitted_leaf(fit, X)
  lv <- sort(unique(leaf))
  idx <- match(leaf, lv)

  n_l <- as.integer(table(leaf)[as.character(lv)])
  n1_l <- as.integer(tapply(A == 1L, leaf, sum)[as.character(lv)])
  n0_l <- as.integer(tapply(A == 0L, leaf, sum)[as.character(lv)])

  # (P3) leaf values: exact within-leaf control means, on the FITTED partition.
  mu_l <- as.numeric(tapply(Y[A == 0L], leaf[A == 0L], mean)[as.character(lv)])
  mu_hat <- mu_l[idx]
  # (P2) e_hat(l): the leaf's treated fraction = the propensity tree forced equal
  # to the outcome tree.
  e_l <- n1_l / n_l
  e_hat <- e_l[idx]

  list(fit = fit, leaf = leaf, lv = lv, n_l = n_l, n1_l = n1_l, n0_l = n0_l,
       mu_l = mu_l, mu_hat = mu_hat, e_l = e_l, e_hat = e_hat)
}

# The DGPs. Two are the study's own Regimes A and B constants (spec §3) so the
# identity is exercised on the structures the simulation actually uses; two are ad
# hoc, one of them with a 4-leaf interaction structure and a deliberately
# asymmetric propensity, so a bug that happens to cancel in a 2-leaf split cannot
# pass unnoticed. All are continuous-outcome and multi-leaf: a degenerate
# single-leaf tree would satisfy the identity trivially and is excluded by
# `expect_gt(length(lv), 1L)` in every case.
tie_dgps <- list(
  regime_A = list(
    seed = 20260909L, n = 900L, leaf_budget = 2L,
    draw = function(n) {
      X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5),
                      X3 = rbinom(n, 1, 0.5))
      e <- ifelse(X$X1 == 1L, 0.7, 0.3)              # spec §3 Regime A
      A <- rbinom(n, 1, e)
      mu0 <- ifelse(X$X1 == 1L, 1, 0)
      tau <- ifelse(X$X3 == 1L, 0.5, 0)
      list(X = X, A = A, Y = mu0 + A * tau + rnorm(n))
    }
  ),
  regime_B = list(
    seed = 20260910L, n = 900L, leaf_budget = 2L,
    draw = function(n) {
      X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5),
                      X3 = rbinom(n, 1, 0.5))
      e <- ifelse(X$X2 == 1L, 0.8, 0.2)              # spec §3 Regime B: e ind. of X1
      A <- rbinom(n, 1, e)
      mu0 <- ifelse(X$X1 == 1L, 1, 0)
      tau <- ifelse(X$X3 == 1L, 0.5, 0)
      list(X = X, A = A, Y = mu0 + A * tau + 1.25 * rnorm(n))
    }
  ),
  interaction_4leaf = list(
    seed = 20260911L, n = 1200L, leaf_budget = 4L,
    draw = function(n) {
      X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5),
                      X3 = rbinom(n, 1, 0.5))
      # Asymmetric propensity across all four (X1, X2) cells, so w_hat(l) differs
      # in every leaf and no two leaves' control-arm terms can cancel each other.
      e <- c(0.25, 0.45, 0.6, 0.8)[1 + X$X1 + 2 * X$X2]
      A <- rbinom(n, 1, e)
      mu0 <- -1 + 2 * X$X1 - 1.5 * X$X2 + 3 * X$X1 * X$X2
      list(X = X, A = A, Y = mu0 + A * 0.4 + rnorm(n, 0, 0.75))
    }
  ),
  unbalanced_leaves = list(
    seed = 20260912L, n = 1500L, leaf_budget = 2L,
    draw = function(n) {
      # Very unequal leaf masses and unequal within-leaf variances: the identity
      # is pointwise, so it must not depend on balance.
      X <- data.frame(X1 = rbinom(n, 1, 0.85), X2 = rbinom(n, 1, 0.5),
                      X3 = rbinom(n, 1, 0.5))
      e <- ifelse(X$X1 == 1L, 0.35, 0.65)
      A <- rbinom(n, 1, e)
      mu0 <- ifelse(X$X1 == 1L, 2, -2)
      list(X = X, A = A, Y = mu0 + A * 0.3 +
             rnorm(n, 0, ifelse(X$X1 == 1L, 0.5, 1.5)))
    }
  )
)

# ---- the identity ------------------------------------------------------------

test_that("cor:single-tree: the tied plug-in equals the forced two-tree estimator", {
  for (nm in names(tie_dgps)) {
    cfg <- tie_dgps[[nm]]
    set.seed(cfg$seed)
    d <- cfg$draw(cfg$n)
    p <- tie_pieces(d$X, d$A, d$Y, cfg$leaf_budget)

    # Not a degenerate single-leaf tree: the identity would hold trivially.
    expect_gt(length(p$lv), 1L)
    expect_true(all(p$n0_l >= 2L), label = paste0(nm, ": every leaf has controls"))

    # The clip must be INERT: a binding clip would mean e_hat is no longer the
    # forced propensity tree, so the object under test would not be theta_tie.
    expect_true(all(p$e_l > 0.01 & p$e_l < 0.99),
                label = paste0(nm, ": propensity clip is inert"))

    # LHS: the one-line tied plug-in.
    theta_tie <- mean(d$Y[d$A == 1L] - p$mu_hat[d$A == 1L])
    # RHS: doubletree's OWN full two-tree solver, propensity forced = outcome tree.
    forced <- doubletree:::eif_att_solve(d$Y, d$A, p$e_hat, p$mu_hat, length(d$Y))

    expect_equal(theta_tie, forced$theta, tolerance = 1e-10,
                 label = paste0(nm, ": theta_tie"),
                 expected.label = paste0(nm, ": forced two-tree theta"))
  }
})

test_that("cor:single-tree: the mechanism is the control-arm term vanishing leaf-by-leaf", {
  # The identity above could in principle hold by a compensating error. This
  # checks the REASON: the within-leaf control residuals sum to zero, which is
  # what kills eq:score's control-arm term.
  cfg <- tie_dgps$interaction_4leaf
  set.seed(cfg$seed)
  d <- cfg$draw(cfg$n)
  p <- tie_pieces(d$X, d$A, d$Y, cfg$leaf_budget)

  resid0 <- (d$Y - p$mu_hat)[d$A == 0L]
  per_leaf <- tapply(resid0, p$leaf[d$A == 0L], sum)
  expect_lt(max(abs(per_leaf)), 1e-9)

  # And the assembled control-arm contribution to theta is zero to the same order.
  w_hat <- p$e_hat / (1 - p$e_hat)
  ctrl_contrib <- sum(w_hat[d$A == 0L] * resid0) / sum(d$A)
  expect_lt(abs(ctrl_contrib), 1e-9)
})

test_that("the identity holds for ANY leaf-constant e_hat, and breaks off mu_hat's partition", {
  cfg <- tie_dgps$interaction_4leaf
  set.seed(cfg$seed)
  d <- cfg$draw(cfg$n)
  p <- tie_pieces(d$X, d$A, d$Y, cfg$leaf_budget)
  n <- length(d$Y)
  theta_tie <- mean(d$Y[d$A == 1L] - p$mu_hat[d$A == 1L])

  # (P1) HOLDS with deliberately WRONG values: e_hat leaf-constant on mu_hat's own
  # partition but nowhere near the treated fraction, including a non-trivial
  # permutation of the true leaf values. Spec §1's "for any leaf-wise-constant
  # e_hat" predicts the identity is UNAFFECTED. If it were affected, the estimator
  # would secretly depend on the propensity values and the corollary would be
  # false as stated.
  for (e_l_alt in list(rep(0.5, length(p$lv)),
                       rev(p$e_l),
                       pmin(pmax(p$e_l + 0.15, 0.02), 0.98),
                       seq(0.1, 0.9, length.out = length(p$lv)))) {
    e_alt <- e_l_alt[match(p$leaf, p$lv)]
    theta_alt <- doubletree:::eif_att_solve(d$Y, d$A, e_alt, p$mu_hat, n)$theta
    expect_equal(theta_tie, theta_alt, tolerance = 1e-10)
  }

  # (P1) VIOLATED: e_hat leaf-constant on a DIFFERENT partition (a split of X3,
  # which the outcome tree does not use). The control-arm term no longer telescopes
  # to zero, so the identity must FAIL -- otherwise the tests above are vacuous.
  leaf_alt <- as.integer(d$X$X3) + 1L
  e_alt_l <- as.numeric(tapply(d$A, leaf_alt, mean))
  theta_off <- doubletree:::eif_att_solve(d$Y, d$A, e_alt_l[leaf_alt], p$mu_hat, n)$theta
  expect_gt(abs(theta_tie - theta_off), 1e-6)

  # (P2) VIOLATED: the right partition, but leaf values that are not the control
  # means (the ALL-unit leaf means instead -- a plausible implementation slip).
  mu_all_l <- as.numeric(tapply(d$Y, p$leaf, mean)[as.character(p$lv)])
  mu_all <- mu_all_l[match(p$leaf, p$lv)]
  theta_mu_wrong <- doubletree:::eif_att_solve(d$Y, d$A, p$e_hat, mu_all, n)$theta
  expect_gt(abs(mean(d$Y[d$A == 1L] - mu_all[d$A == 1L]) - theta_mu_wrong), 1e-6)
})

test_that("estimate_att()'s serialised leaf values match the exact control means to 1e-5", {
  # optimaltrees rounds a leaf's fitted value to six significant decimal digits.
  # That is the reason the identity is tested on recomputed leaf means rather than
  # on `outcome_control` directly (see this file's header), so the rounding
  # magnitude is pinned here: a regression that made it materially worse would
  # otherwise show up only as an unexplained loosening elsewhere.
  cfg <- tie_dgps$regime_A
  set.seed(cfg$seed)
  d <- cfg$draw(cfg$n)
  p <- tie_pieces(d$X, d$A, d$Y, cfg$leaf_budget)

  serialised <- p$fit$nuisance_fits$outcome_control
  expect_equal(serialised, p$mu_hat, tolerance = 1e-5)
  # And the induced error in theta_tie is of the same small order -- i.e. the
  # substitution changes no result this study reports.
  theta_serialised <- mean(d$Y[d$A == 1L] - serialised[d$A == 1L])
  theta_exact <- mean(d$Y[d$A == 1L] - p$mu_hat[d$A == 1L])
  expect_lt(abs(theta_serialised - theta_exact), 1e-5)
})
