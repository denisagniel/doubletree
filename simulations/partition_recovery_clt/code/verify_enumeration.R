# ============================================================
# verify_enumeration.R
# Study: partition_recovery_clt  (doubletree, S1)
#
# HARD PREREQUISITE for the replication harness (spec §7: "confirm by hand on at
# least one small, fully-inspectable DGP that it correctly classifies a known
# sufficient partition and a known non-sufficient one before trusting it inside
# 1000 replications").
#
# Run:  Rscript simulations/partition_recovery_clt/code/verify_enumeration.R
#       (from the doubletree package root)
#
# WHAT IS VERIFIED, AND AGAINST WHAT
#
#   A. enumerate_tree_partitions() vs an INDEPENDENT brute-force enumerator.
#      The independent route enumerates every SET PARTITION of the 2^p cells
#      with at most Lbar blocks (restricted-growth strings) and filters with a
#      from-scratch tree-realizability predicate that shares no code with the
#      constructive enumerator. Two different algorithms must produce the same
#      set of partitions, for p in {2,3} and Lbar in {2,3,4}. This is what
#      catches an off-by-one in the budget split or a missed topology -- a
#      hand-checked COUNT for one configuration would not.
#
#   B. Sufficiency classification against a fully hand-worked p = 2 DGP. All 8
#      members of T_4 are written out below with the answer derived by hand, and
#      the code must agree on all 8. Critically this includes the two
#      REDUNDANT-REFINEMENT cases (3-leaf refinements of the 2-leaf sufficient
#      partition) which a naive same-partition or same-leaf-count check would
#      score as failures -- test B3 asserts that misclassification explicitly,
#      so the reason for not using the naive check is locked down by a test
#      rather than left as a comment.
#
#   C. Delta_j against closed-form arithmetic done by hand, in both currencies
#      (squared error and KL), under uniform AND non-uniform cell weights. The
#      non-uniform case matters because the mu-nuisance risk is CONTROL-weighted
#      (nu_x = p_x{1 - e_0(x)}), never uniform.
#
#   D. The fitted-model path end to end: fit real trees with optimaltrees and
#      confirm partition_recovered() returns TRUE for a tree that recovers the
#      X1 structure and FALSE for one that collapses to a single leaf.
#
#   E. selection_bound() / implied_c1() round-trip, including both edges at
#      which the reported bound goes uninformative.
#
#   F. The DGP guards: each nuisance depends on exactly the coordinates it
#      declares, the two nuisances SHARE a coordinate so the DGP is genuinely
#      confounded (F1b -- the regression test for the disjoint-support defect
#      that shipped and invalidated 12 cells of a completed sweep), the
#      population spec agrees with the sample draw, theta_0 is exact, and the
#      margin dial actually moves the margin.
# ============================================================

suppressPackageStartupMessages({
  library(testthat)
  library(cli)
})

PKG_ROOT <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(PKG_ROOT, "DESCRIPTION"))) {
  cli::cli_abort("Run from the doubletree package root (no DESCRIPTION here).")
}
CODE <- file.path(PKG_ROOT, "simulations", "partition_recovery_clt", "code")
source(file.path(CODE, "enumerate_sufficient_class.R"))

if (!isNamespaceLoaded("optimaltrees")) {
  .ot_src <- file.path(dirname(PKG_ROOT), "optimaltrees")
  if (dir.exists(.ot_src)) {
    pkgload::load_all(.ot_src, quiet = TRUE, export_all = FALSE)
  } else {
    library(optimaltrees)
  }
}

# ============================================================
# A. Independent brute-force enumerator (shares no code with the constructive one)
# ============================================================

#' Every set partition of 1:m with at most k blocks, as label vectors
#'
#' Restricted-growth strings: label[1] = 1 and label[i] <= 1 + max(label[1:i-1]).
#' Every set partition has exactly one such representation, so this enumerates
#' each partition exactly once with no deduplication step.
all_set_partitions <- function(m, k) {
  out <- list()
  rec <- function(prefix, used) {
    i <- length(prefix) + 1L
    if (i > m) {
      out[[length(out) + 1L]] <<- prefix
      return(invisible())
    }
    for (lab in seq_len(min(used + 1L, k))) {
      rec(c(prefix, lab), max(used, lab))
    }
    invisible()
  }
  rec(integer(0), 0L)
  out
}

#' Is a partition realizable by a decision tree with at most `budget` leaves?
#'
#' Independent recursive predicate. A partition of a region is realizable iff
#' either the whole region is one block, or some free coordinate splits the
#' region without CUTTING any block, and both sides are realizable within a
#' budget split. Written from the definition of a decision tree, deliberately
#' not reusing enumerate_tree_partitions()'s recursion.
tree_realizable <- function(labels, cell_mat, region, free_coords, budget) {
  blocks <- unique(labels[region])
  if (length(blocks) > budget) return(FALSE)
  if (length(blocks) == 1L) return(budget >= 1L)
  for (k in free_coords) {
    lo <- region[cell_mat[region, k] == 0L]
    hi <- region[cell_mat[region, k] == 1L]
    if (length(lo) == 0L || length(hi) == 0L) next
    # A block straddling the split cannot be produced by splitting here.
    if (length(intersect(unique(labels[lo]), unique(labels[hi]))) > 0L) next
    b_lo <- length(unique(labels[lo]))
    b_hi <- length(unique(labels[hi]))
    if (b_lo + b_hi > budget) next
    rem <- free_coords[free_coords != k]
    if (tree_realizable(labels, cell_mat, lo, rem, b_lo) &&
        tree_realizable(labels, cell_mat, hi, rem, b_hi)) {
      return(TRUE)
    }
  }
  FALSE
}

brute_force_class <- function(p, leaf_budget) {
  cells <- make_cell_grid(p)
  cm <- as.matrix(cells)
  m <- nrow(cm)
  cand <- all_set_partitions(m, leaf_budget)
  keep <- vapply(cand, function(lab) {
    tree_realizable(lab, cm, seq_len(m), seq_len(p), leaf_budget)
  }, logical(1))
  sort(vapply(cand[keep], function(lab) paste(canonical_labels(lab), collapse = "."),
              character(1)))
}

test_that("A1. constructive enumeration matches independent brute force", {
  for (p in c(2L, 3L)) {
    for (Lbar in c(2L, 3L, 4L)) {
      got <- sort(enumerate_tree_partitions(p, Lbar)$keys)
      want <- brute_force_class(p, Lbar)
      expect_identical(got, want, info = paste0("p = ", p, ", Lbar = ", Lbar))
    }
  }
})

test_that("A2. |T_Lbar| for p = 2, Lbar = 4 is 8, as hand-enumerated below", {
  # 1 (single leaf) + 2 (one split, X1 or X2) + 4 (root split then one child
  # split: 2 root coords x 2 choices of which child) + 1 (all four singletons,
  # reachable two ways and therefore deduplicated to one partition) = 8.
  expect_equal(enumerate_tree_partitions(2L, 4L)$n_partitions, 8L)
})

test_that("A3. the study's own configuration enumerates cleanly", {
  en5 <- enumerate_tree_partitions(5L, 4L)
  expect_gt(en5$n_partitions, 0L)
  expect_true(all(en5$n_leaves <= 4L))
  expect_false(anyDuplicated(en5$keys) > 0L)
  # Every enumerated partition must be tree-realizable by the independent
  # predicate too. The reverse inclusion at p = 5 is too big to brute force, but
  # this direction is cheap and catches a spurious partition (which would inflate
  # |T_Lbar| in prop:selection-rate's bound).
  cm <- as.matrix(en5$cells)
  ok <- apply(en5$labels, 1L, function(lab) {
    tree_realizable(lab, cm, seq_len(nrow(cm)), seq_len(5L), 4L)
  })
  expect_true(all(ok))
})

# ============================================================
# B. Hand-worked p = 2 DGP
# ============================================================
#
# Cell grid (expand.grid order, X1 varying fastest):
#   cell 1 = (X1=0, X2=0)
#   cell 2 = (X1=1, X2=0)
#   cell 3 = (X1=0, X2=1)
#   cell 4 = (X1=1, X2=1)
#
# Truth: gamma_0 depends on X1 ONLY.  gamma_0 = 0.2 if X1 = 0, 0.8 if X1 = 1,
# i.e. by cell: (0.2, 0.8, 0.2, 0.8).
#
# The 8 members of T_4 and the hand-derived answer for each:
#
#  key        leaves                        gamma_0 per leaf     sufficient?
#  1.1.1.1    {all}                         {.2,.8,.2,.8}        NO
#  1.2.1.2    {X1=0} {X1=1}                 {.2,.2} {.8,.8}      YES  <- reference
#  1.1.2.2    {X2=0} {X2=1}                 {.2,.8} {.2,.8}      NO
#  1.2.3.2    {c1} {X1=1} {c3}              {.2} {.8,.8} {.2}    YES  <- REDUNDANT REFINEMENT
#  1.2.1.3    {X1=0} {c2} {c4}              {.2,.2} {.8} {.8}    YES  <- REDUNDANT REFINEMENT
#  1.2.3.3    {c1} {c2} {X2=1}              {.2} {.8} {.2,.8}    NO
#  1.1.2.3    {X2=0} {c3} {c4}              {.2,.8} {.2} {.8}    NO
#  1.2.3.4    all singletons                {.2}{.8}{.2}{.8}     YES  <- finest refinement
#
# So |S_j| = 4 of 8, and the minimum sufficient leaf count is 2.

HAND_P2 <- list(
  gamma0 = c(0.2, 0.8, 0.2, 0.8),
  expected = c(
    "1.1.1.1" = FALSE,
    "1.2.1.2" = TRUE,
    "1.1.2.2" = FALSE,
    "1.2.3.2" = TRUE,
    "1.2.1.3" = TRUE,
    "1.2.3.3" = FALSE,
    "1.1.2.3" = FALSE,
    "1.2.3.4" = TRUE
  )
)

test_that("B1. sufficiency matches the hand-derived answer for all 8 partitions", {
  en <- enumerate_tree_partitions(2L, 4L)
  expect_setequal(en$keys, names(HAND_P2$expected))
  for (i in seq_len(en$n_partitions)) {
    key <- en$keys[[i]]
    got <- is_sufficient_partition(en$labels[i, ], HAND_P2$gamma0)
    expect_identical(got, unname(HAND_P2$expected[[key]]),
                     info = paste0("partition ", key))
  }
})

test_that("B2. sufficient_class() reports |S_j| = 4 and min 2 sufficient leaves", {
  en <- enumerate_tree_partitions(2L, 4L)
  cl <- sufficient_class(HAND_P2$gamma0, w = rep(0.25, 4L), enumeration = en)
  expect_equal(cl$n_partitions, 8L)
  expect_equal(cl$n_sufficient, 4L)
  expect_equal(cl$min_sufficient_leaves, 2L)
  # S_j is NOT a singleton: the property that forces the currency-correct
  # indicator and forbids a same-partition check.
  expect_gt(cl$n_sufficient, 1L)
  # Excess risk is exactly 0 on S_j and strictly positive off it -- the risk and
  # structural criteria must not disagree.
  expect_true(all(cl$excess_sq[cl$sufficient] < 1e-14))
  expect_true(all(cl$excess_sq[!cl$sufficient] > 0))
  expect_true(all(cl$excess_kl[cl$sufficient] < 1e-14))
  expect_true(all(cl$excess_kl[!cl$sufficient] > 0))
})

test_that("B3. the NAIVE checks the spec forbids really do misclassify", {
  # Locking in the reason enumerate_sufficient_class.R does not implement them.
  en <- enumerate_tree_partitions(2L, 4L)
  reference_key <- "1.2.1.2"     # the 2-leaf sufficient partition
  refinement_key <- "1.2.3.2"    # sufficient, 3 leaves, different shape
  i_ref <- match(reference_key, en$keys)
  i_refine <- match(refinement_key, en$keys)

  # Ground truth from the definition: both are in S_j.
  expect_true(is_sufficient_partition(en$labels[i_ref, ], HAND_P2$gamma0))
  expect_true(is_sufficient_partition(en$labels[i_refine, ], HAND_P2$gamma0))

  # Naive check 1, "same partition as the reference": WRONG on the refinement.
  expect_false(identical(refinement_key, reference_key))

  # Naive check 2, "same leaf count as the reference": also WRONG.
  n_ref <- length(unique(en$labels[i_ref, ]))
  n_refine <- length(unique(en$labels[i_refine, ]))
  expect_equal(n_ref, 2L)
  expect_equal(n_refine, 3L)
  expect_false(n_ref == n_refine)
})

# ============================================================
# C. Delta_j by hand, both currencies, uniform and non-uniform weights
# ============================================================
#
# UNIFORM w = 1/4 per cell.
#
# The two tightest competitors outside S_j are "1.2.3.3" and "1.1.2.3": each has
# exactly ONE leaf that merges a 0.2 cell with a 0.8 cell, the other leaves being
# pure. For a leaf merging two cells of weights w_a, w_b and values g_a, g_b, the
# excess squared risk contributed is (w_a w_b)/(w_a + w_b) * (g_a - g_b)^2
#   = (0.25 * 0.25)/0.5 * 0.6^2 = 0.125 * 0.36 = 0.045
# The other two non-sufficient partitions ("1.1.1.1", "1.1.2.2") merge such a
# pair in BOTH leaves, giving 0.09. So Delta_sq = 0.045 exactly.
#
# KL currency: a leaf merging equal masses of 0.2 and 0.8 has leaf mean 0.5, and
# contributes, weighted by the leaf's own mass 0.5,
#   0.5 * [0.5 KL(0.2||0.5) + 0.5 KL(0.8||0.5)] = 0.5 * [log 2 - H(0.2)]
# since KL(p||0.5) = log 2 - H(p), giving Delta_kl = (log 2 - H(0.2))/2.

binary_entropy <- function(p) -(p * log(p) + (1 - p) * log(1 - p))

test_that("C1. Delta_sq and Delta_kl match closed-form hand arithmetic (uniform w)", {
  en <- enumerate_tree_partitions(2L, 4L)
  cl <- sufficient_class(HAND_P2$gamma0, w = rep(0.25, 4L), enumeration = en)

  expect_equal(cl$delta_sq, 0.045, tolerance = 1e-12)
  expect_equal(cl$delta_kl, (log(2) - binary_entropy(0.2)) / 2, tolerance = 1e-12)

  # The tightest competitor is one of the two single-mixed-leaf partitions.
  expect_true(cl$argmin_key_sq %in% c("1.2.3.3", "1.1.2.3"))
  expect_true(cl$argmin_key_kl %in% c("1.2.3.3", "1.1.2.3"))

  # And the coarsest competitor really is twice as bad, as hand-computed.
  i_root <- match("1.1.1.1", en$keys)
  expect_equal(cl$excess_sq[[i_root]], 0.09, tolerance = 1e-12)
  expect_equal(cl$excess_kl[[i_root]], log(2) - binary_entropy(0.2),
               tolerance = 1e-12)
})

test_that("C2. non-uniform (control-style) weights: hand two-point formula", {
  # w = (0.1, 0.2, 0.3, 0.4) over cells 1..4. This is the shape the mu-nuisance
  # weighting takes: nu_x = p_x {1 - e_0(x)} is never uniform when e_0 varies.
  w <- c(0.1, 0.2, 0.3, 0.4)
  en <- enumerate_tree_partitions(2L, 4L)
  cl <- sufficient_class(HAND_P2$gamma0, w = w, enumeration = en)

  # "1.2.3.3" merges cells 3 (w = .3, g = .2) and 4 (w = .4, g = .8):
  #   excess_sq = (.3 * .4)/(.3 + .4) * .6^2 = (0.12/0.7) * 0.36
  i_a <- match("1.2.3.3", en$keys)
  expect_equal(cl$excess_sq[[i_a]], (0.3 * 0.4 / 0.7) * 0.36, tolerance = 1e-12)

  # "1.1.2.3" merges cells 1 (w = .1, g = .2) and 2 (w = .2, g = .8):
  #   excess_sq = (.1 * .2)/.3 * 0.36 = 0.024
  i_b <- match("1.1.2.3", en$keys)
  expect_equal(cl$excess_sq[[i_b]], (0.1 * 0.2 / 0.3) * 0.36, tolerance = 1e-12)

  # Delta_sq is the smaller of the two, i.e. 0.024, NOT the uniform-weight 0.045.
  expect_equal(cl$delta_sq, 0.024, tolerance = 1e-12)
  expect_identical(cl$argmin_key_sq, "1.1.2.3")

  # Weights are normalised internally, so scaling them must not move Delta.
  cl_scaled <- sufficient_class(HAND_P2$gamma0, w = 7 * w, enumeration = en)
  expect_equal(cl_scaled$delta_sq, cl$delta_sq, tolerance = 1e-14)
  expect_equal(cl_scaled$delta_kl, cl$delta_kl, tolerance = 1e-14)

  # Sufficiency itself must NOT depend on the weights (def:sufficient holds for
  # any strictly positive nu).
  cl_u <- sufficient_class(HAND_P2$gamma0, w = rep(0.25, 4L), enumeration = en)
  expect_identical(cl$sufficient, cl_u$sufficient)
})

test_that("C3. degenerate DGPs are refused loudly, not scored silently", {
  en <- enumerate_tree_partitions(2L, 4L)
  # gamma_0 constant -> every partition sufficient -> Delta undefined.
  expect_error(sufficient_class(rep(0.3, 4L), rep(0.25, 4L), en), "EVERY partition")
  # At Lbar = 2 with a truth varying in both coordinates, S_j is empty.
  en2 <- enumerate_tree_partitions(2L, 2L)
  expect_error(sufficient_class(c(0.1, 0.4, 0.6, 0.9), rep(0.25, 4L), en2),
               "S_j is EMPTY")
})

# ============================================================
# D. Fitted-model path, end to end against the real solver
# ============================================================

test_that("D. partition_recovered() scores real fitted trees correctly", {
  set.seed(20260901L)
  n <- 4000L
  X <- data.frame(X1 = rbinom(n, 1L, 0.5), X2 = rbinom(n, 1L, 0.5))
  e0 <- ifelse(X$X1 == 1L, 0.8, 0.2)
  A <- rbinom(n, 1L, e0)

  en <- enumerate_tree_partitions(2L, 4L)
  cl <- sufficient_class(HAND_P2$gamma0, w = rep(0.25, 4L), enumeration = en)

  # (i) Penalty large enough to kill every split. The population gain from the
  # X1 split is log(2) - H(0.2) = 0.193 on the mean-log-loss scale, so a penalty
  # of 0.5 per leaf forces a single leaf -- a KNOWN non-sufficient partition.
  fit_flat <- optimaltrees::fit_tree(
    X, A, loss_function = "log_loss", regularization = 0.5,
    worker_limit = 1L, verbose = FALSE
  )
  rec_flat <- partition_recovered(fit_flat, cl)
  expect_equal(rec_flat$n_leaves_grid, 1L)
  expect_false(rec_flat$recovered)
  expect_equal(rec_flat$excess_sq, 0.09, tolerance = 1e-6)
  expect_true(rec_flat$in_enumerated_class)

  # (ii) Penalty well below the X1 gain but above the (population-zero) X2 gain:
  # the fit must recover the X1 structure.
  fit_ok <- optimaltrees::fit_tree(
    X, A, loss_function = "log_loss", regularization = 0.05,
    worker_limit = 1L, verbose = FALSE
  )
  rec_ok <- partition_recovered(fit_ok, cl)
  expect_true(rec_ok$recovered)
  expect_equal(rec_ok$excess_sq, 0, tolerance = 1e-14)
  expect_true(rec_ok$in_enumerated_class)

  # (iii) Vanishing penalty: the solver is free to spend leaves on the
  # population-irrelevant X2. Whatever it does, the verdict must still be
  # "recovered", because a refinement of a sufficient partition is sufficient.
  # The fitted-path counterpart of test B3.
  fit_fine <- optimaltrees::fit_tree(
    X, A, loss_function = "log_loss", regularization = 1e-9,
    worker_limit = 1L, verbose = FALSE
  )
  rec_fine <- partition_recovered(fit_fine, cl)
  expect_true(rec_fine$recovered)
  expect_equal(rec_fine$excess_sq, 0, tolerance = 1e-14)

  # A non-binary X must be refused rather than silently mis-assigned.
  Xc <- data.frame(X1 = runif(n), X2 = rbinom(n, 1L, 0.5))
  fit_cont <- optimaltrees::fit_tree(
    Xc, A, loss_function = "log_loss", regularization = 0.05,
    worker_limit = 1L, verbose = FALSE
  )
  expect_error(partition_recovered(fit_cont, cl), "all_binary")
})

# ============================================================
# E. The bound helpers
# ============================================================

test_that("E. selection_bound() and implied_c1() are mutual inverses", {
  card_T <- 8L
  delta <- 0.045
  n <- c(200L, 1000L, 5000L)
  # c1 must put c1 * n * Delta^2 inside (log |T_Lbar|, ~700) across the whole n
  # grid for the round-trip to be well posed: below log|T| the bound saturates at
  # 1, above ~700 exp() underflows to exactly 0. c1 = 20 gives exponents
  # 8.1 / 40.5 / 202.5. Both edges are pinned separately below, because they are
  # the two ways the reported bound column can go uninformative in the real sweep.
  c1 <- 20
  b <- selection_bound(n, delta, card_T, c1 = c1)
  expect_true(all(b > 0 & b < 1))
  expect_equal(implied_c1(b, n, delta, card_T), rep(c1, length(n)),
               tolerance = 1e-8)
  expect_true(all(diff(b) < 0))
  # Edge 1: a tiny margin makes the bound vacuous; capped at 1, not |T_Lbar|.
  expect_equal(selection_bound(1L, 1e-9, card_T, c1 = 1), 1)
  # Edge 2: a large n * Delta^2 underflows the bound to exactly 0. Still a valid
  # bound, but implied_c1() cannot be inverted from it.
  b_tiny <- selection_bound(100000L, 0.5, card_T, c1 = 1)
  expect_equal(b_tiny, 0)
  expect_true(is.na(implied_c1(b_tiny, 100000L, 0.5, card_T)))
  # A zero observed failure rate carries no information about c1.
  expect_true(is.na(implied_c1(0, 1000L, delta, card_T)))
})

# ============================================================
# F. DGP guards
# ============================================================

source(file.path(CODE, "dgps.R"))
source(file.path(PKG_ROOT, "simulations", "dgps", "dgps_stress.R"))

test_that("F1. check_variable_roles() accepts the real DGPs", {
  s <- dgp_spec_grid_sparse(p = 5L, id = "f1")
  expect_true(check_variable_roles(s))
  expect_identical(s$e_vars, c("X1", "X2"))
  expect_identical(s$mu_vars, c("X1", "X3"))
  expect_true(check_variable_roles(dgp_spec_weak_overlap()))
})

test_that("F1b. every DGP is CONFOUNDED: shared support, and A not indep of Y(0)", {
  # THE REGRESSION TEST for the defect diagnosed in the design audit §2 and
  # pinned as item 1 of the spec's §8 addendum. F1 originally put e_0 on
  # (X1, X2) and mu_0 on the DISJOINT (X3, X4), which makes A independent of
  # Y(0): a difference in means is already unbiased and NO amount of
  # partition-recovery failure can move theta_hat. Nothing errored -- an 18-cell
  # x 1000-rep sweep ran to completion on it and reported nominal coverage at
  # ST2 n = 200 while recovery had essentially never happened (rec_both = 0.001).
  # Only an intent-level assertion catches that, so here it is, in both
  # currencies: structural (shared covariate) and behavioural (E[Y(0)|A] moves).
  ey0_given_a <- function(s, a) {
    wa <- if (a == 1L) s$e0 else 1 - s$e0
    sum(s$cell_prob * wa * s$mu0) / sum(s$cell_prob * wa)
  }
  for (s in list(dgp_spec_grid_sparse(p = 5L, id = "f1"),
                 dgp_spec_grid_sparse(p = 5L, e_gap = 0.20, mu_gap = 0.20,
                                      id = "f1_weak"),
                 dgp_spec_weak_overlap())) {
    expect_gt(length(intersect(s$e_vars, s$mu_vars)), 0L)
    expect_gt(abs(ey0_given_a(s, 1L) - ey0_given_a(s, 0L)), 0.01)
  }

  # And the defect itself, reconstructed, must fail BOTH criteria -- so this test
  # would have caught the version that shipped.
  disjoint <- dgp_spec_grid_sparse(p = 5L, id = "old_f1")
  disjoint$mu0 <- ifelse(disjoint$cells$X3 == 0L, 0.10,
                         ifelse(disjoint$cells$X4 == 0L, 0.50, 0.90))
  disjoint$mu_vars <- c("X3", "X4")
  expect_length(intersect(disjoint$e_vars, disjoint$mu_vars), 0L)
  expect_equal(ey0_given_a(disjoint, 1L), ey0_given_a(disjoint, 0L),
               tolerance = 1e-14)
})

test_that("F2. check_variable_roles() catches a nuisance keyed on the wrong covariate", {
  # An earlier state of dgps.R keyed mu_0 on (X1, X3) while DECLARING (X3, X4),
  # consistently on both the population and the sample side. Nothing failed: the
  # two agreed, every fit succeeded, every fitted partition landed inside
  # T_Lbar, and the only symptom was Delta_mu coming out the wrong size.
  # Reproduced here in the current declaration's terms -- mu_0 is now legitimately
  # gated on X1, so the off-declaration coordinate to probe with is X2 -- so the
  # guard cannot regress.
  s <- dgp_spec_grid_sparse(p = 5L, id = "f1")
  broken <- s
  broken$mu0 <- ifelse(s$cells$X2 == 0L, 0.10,
                       ifelse(s$cells$X3 == 0L, 0.50, 0.90))
  expect_error(check_variable_roles(broken), "mu_0 depends on X2")

  # And the stale-declaration direction: a mu_0 that ignores a variable it claims.
  stale <- s
  stale$mu0 <- ifelse(s$cells$X3 == 0L, 0.10, 0.90)   # never uses the X1 gate
  expect_error(check_variable_roles(stale),
               "does NOT depend on its declared variable X1")
})

test_that("F3. check_spec_vs_draw() catches a population/sample mismatch", {
  s <- dgp_spec_grid_sparse(p = 5L, id = "f1")
  expect_s3_class(check_spec_vs_draw(s, n = 20000L), "data.frame")
  # Perturb the population side only; the draw is untouched.
  bad <- s
  bad$mu0 <- s$mu0 + 0.05
  expect_error(check_spec_vs_draw(bad, n = 20000L), "disagrees with the sample draw")
})

test_that("F4. theta_0 is exact, Y is bounded, and clipping is refused", {
  s <- dgp_spec_grid_sparse(p = 5L, tau = 0.08, id = "f1")
  expect_equal(s$theta0, 0.08, tolerance = 1e-14)
  expect_true(all(abs(s$tau_x - 0.08) < 1e-14))
  # Y bounded: spec §7's precondition for prop:selection-rate.
  expect_true(all(s$mu0 >= 0 & s$mu0 <= 1))
  expect_true(all(s$mu1 >= 0 & s$mu1 <= 1))
  # A tau that would need clipping is refused rather than silently moving theta_0
  # away from tau.
  expect_error(dgp_spec_grid_sparse(p = 5L, mu_gap = 0.4, tau = 0.3),
               "would need clipping")
  # ST3's theta_0 is likewise exact.
  expect_equal(dgp_spec_weak_overlap(tau = 0.10)$theta0, 0.10, tolerance = 1e-14)
})

test_that("F5. the margin dial really moves the margin, monotonically", {
  # ST2's whole construction rests on this: shrinking *_gap must shrink Delta_j.
  en <- enumerate_tree_partitions(5L, 4L)
  deltas <- vapply(c(0.40, 0.24, 0.14, 0.08), function(g) {
    s <- dgp_spec_grid_sparse(p = 5L, e_gap = g, mu_gap = g, id = "g")
    sufficient_class(s$e0, s$cell_prob, en)$delta_kl
  }, numeric(1))
  expect_true(all(diff(deltas) < 0))
  expect_true(all(deltas > 0))
})

cli::cli_alert_success("verify_enumeration.R: all checks passed.")

# Reported outside test_that() because testthat captures messages emitted inside
# it. These are the class sizes and margins the study's tables cite.
local({
  for (p in c(4L, 5L, 6L)) {
    en <- enumerate_tree_partitions(p, 4L)
    cli::cli_inform("  |T_4| at p = {p}: {en$n_partitions} partitions")
  }
  en5 <- enumerate_tree_partitions(5L, 4L)
  # The two shapes the study actually runs (F1 and ST2), stated as leaf values
  # rather than as a single "gap" dial, because F1 and ST2 differ only in the gap
  # while e_lo/e_mid (the shared X1 gate) are held fixed. code/calibrate_margin.R
  # is the authority on these numbers and on the lambda_n crossings; this is a
  # smoke print so a DGP change cannot pass verification silently.
  shapes <- list(
    f1      = list(e_gap = 0.30, mu_gap = 0.40),
    f1_weak = list(e_gap = 0.20, mu_gap = 0.20)
  )
  for (nm in names(shapes)) {
    g <- shapes[[nm]]
    s <- dgp_spec_grid_sparse(p = 5L, e_gap = g$e_gap, mu_gap = g$mu_gap, id = nm)
    ce <- sufficient_class(s$e0, s$cell_prob, en5)
    cm <- sufficient_class(s$mu0, s$nu, en5)
    ey0 <- vapply(c(1L, 0L), function(a) {
      wa <- if (a == 1L) s$e0 else 1 - s$e0
      sum(s$cell_prob * wa * s$mu0) / sum(s$cell_prob * wa)
    }, numeric(1))
    cli::cli_inform(
      "  {nm}: Delta_e(kl) = {signif(ce$delta_kl, 5)}, Delta_mu(kl) = {signif(cm$delta_kl, 5)}, |S_e| = {ce$n_sufficient}, |S_mu| = {cm$n_sufficient}, E[Y0|A=1] = {signif(ey0[[1]], 4)} vs E[Y0|A=0] = {signif(ey0[[2]], 4)}"
    )
  }
})
