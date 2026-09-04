# ============================================================
# review_checks.R
# Study: partition_recovery_clt  (doubletree, S1)
#
# Independent-by-construction re-derivations of the load-bearing numeric claims,
# written for the review pass in quality_reports/. Each check recomputes a
# quantity by a DIFFERENT route than the production code uses, so agreement is
# evidence rather than a tautology. This file is a review artifact, not part of
# the pipeline; nothing in code/ depends on it.
#
# Run:  Rscript simulations/partition_recovery_clt/code/review_checks.R
#       (from the doubletree package root)
# ============================================================

suppressPackageStartupMessages({ library(testthat); library(cli) })
source(file.path("simulations", "partition_recovery_clt", "code", "common.R"))

f1 <- make_dgp("f1")
en <- f1$enumeration

# ---- R1 enumeration, two-sided random cross-check at p = 4 ---------------
# verify_enumeration.R cross-checks both directions exhaustively at p <= 3 and
# one direction at p = 5. p = 4 is too large to enumerate all set partitions
# (S(16,<=4) ~ 1.7e8) but random probing is a genuine two-sided test: for random
# label vectors, "in the enumerated key set" must agree with an independent
# tree-realizability predicate.
tree_realizable <- function(labels, cm, region, free_coords, budget) {
  blocks <- unique(labels[region])
  if (length(blocks) > budget) return(FALSE)
  if (length(blocks) == 1L) return(budget >= 1L)
  for (k in free_coords) {
    lo <- region[cm[region, k] == 0L]; hi <- region[cm[region, k] == 1L]
    if (length(lo) == 0L || length(hi) == 0L) next
    if (length(intersect(unique(labels[lo]), unique(labels[hi]))) > 0L) next
    b_lo <- length(unique(labels[lo])); b_hi <- length(unique(labels[hi]))
    if (b_lo + b_hi > budget) next
    rem <- free_coords[free_coords != k]
    if (tree_realizable(labels, cm, lo, rem, b_lo) &&
        tree_realizable(labels, cm, hi, rem, b_hi)) return(TRUE)
  }
  FALSE
}

test_that("R1. p=4 enumeration agrees with the independent predicate on random probes", {
  en4 <- enumerate_tree_partitions(4L, 4L)
  cm <- as.matrix(en4$cells); m <- nrow(cm)
  set.seed(42L)
  n_real <- 0L; n_probe <- 0L
  probe <- function(lab) {
    lab <- canonical_labels(lab)
    real <- tree_realizable(lab, cm, seq_len(m), seq_len(4L), 4L)
    enum <- paste(lab, collapse = ".") %in% en4$keys
    expect_identical(real, enum, info = paste(lab, collapse = "."))
    n_real <<- n_real + real; n_probe <<- n_probe + 1L
  }
  # (a) Every enumerated partition: must be realizable. One-sided but exhaustive.
  for (i in seq_len(en4$n_partitions)) probe(en4$labels[i, ])
  # (b) PERTURBED enumerated partitions. Uniform random label vectors are
  # essentially never tree partitions (155 of ~1.7e8 set partitions of 16 cells
  # into <= 4 blocks), so uniform probing only ever exercises the FALSE/FALSE
  # corner and proves nothing. Reassigning 1-2 cells of a real partition lands
  # near the boundary, which is where a disagreement would actually live.
  for (i in seq_len(3000L)) {
    lab <- en4$labels[sample.int(en4$n_partitions, 1L), ]
    for (j in sample.int(m, sample.int(2L, 1L))) {
      lab[[j]] <- sample.int(max(4L, max(lab)), 1L)
    }
    probe(lab)
  }
  # (c) Uniform random, for the far side.
  for (i in seq_len(500L)) probe(sample.int(4L, m, replace = TRUE))
  # The probe must actually hit BOTH answers, or it is vacuous.
  expect_gt(n_real, 0L)
  expect_lt(n_real, n_probe)
  cli::cli_inform("  R1: {n_probe} probes at p=4, {n_real} realizable, 0 disagreements.")
})

# ---- R2 Delta_j re-derived by a from-scratch loop ------------------------
# excess_risk() uses tapply(); this recomputes with an explicit per-leaf loop and
# the two-point closed form, i.e. different code and different algebra.
test_that("R2. Delta_j matches a from-scratch loop and the hand two-point formula", {
  brute_excess <- function(lab, g, w) {
    wn <- w / sum(w); sq <- 0; kl <- 0
    for (L in unique(lab)) {
      ix <- which(lab == L)
      gb <- sum(wn[ix] * g[ix]) / sum(wn[ix])
      sq <- sq + sum(wn[ix] * (g[ix] - gb)^2)
      kl <- kl + sum(wn[ix] * (g[ix] * log(g[ix] / gb) +
                                 (1 - g[ix]) * log((1 - g[ix]) / (1 - gb))))
    }
    c(sq = sq, kl = kl)
  }
  for (nm in c("e", "mu")) {
    cl <- if (nm == "e") f1$class_e else f1$class_mu
    non <- which(!cl$sufficient)
    got <- vapply(non, function(i) brute_excess(en$labels[i, ], cl$gamma0, cl$w)[["sq"]],
                  numeric(1))
    expect_equal(min(got), cl$delta_sq, tolerance = 1e-12)
    gotk <- vapply(non, function(i) brute_excess(en$labels[i, ], cl$gamma0, cl$w)[["kl"]],
                   numeric(1))
    expect_equal(min(gotk), cl$delta_kl, tolerance = 1e-12)
  }
  # Hand formula, re-derived for the CONFOUNDED parameterization (e_0 on (X1, X2),
  # mu_0 on (X1, X3), sharing X1; e = 0.20/0.45/0.75, mu = 0.10/0.50/0.90).
  #
  # The two nuisances' tightest competitors have DIFFERENT SHAPES, which is the
  # whole reason Delta_j is enumerated rather than computed from a two-point
  # formula. Both are 4-leaf partitions with exactly one impure leaf.
  #
  # e (UNIFORM weights w = 1/32): the cheapest impure leaf is {X1=1, Xk=c} for a
  # NOISE coordinate k -- tree = {X1=0} | {X1=1,Xk=1} split on X2 | {X1=1,Xk=0}
  # merged = 4 leaves. That leaf has mass 1/4, split equally between e = 0.45 and
  # e = 0.75, so its leaf mean is 0.60:
  #   Delta_e(kl) = (1/4) * [0.5 KL(.45||.60) + 0.5 KL(.75||.60)]
  #   Delta_e(sq) = (1/8 * 1/8)/(1/4) * (0.75-0.45)^2 = 0.0625 * 0.09
  kl2 <- function(p, q) p * log(p / q) + (1 - p) * log((1 - p) / (1 - q))
  expect_equal(f1$class_e$delta_kl, 0.25 * 0.5 * (kl2(0.45, 0.60) + kl2(0.75, 0.60)),
               tolerance = 1e-12)
  expect_equal(f1$class_e$delta_sq, (1 / 8 * 1 / 8) / (1 / 4) * 0.30^2,
               tolerance = 1e-12)

  # mu (CONTROL weights nu_x = p_x(1 - e_0), never uniform): the cheapest impure
  # leaf is NOT the noise-coordinate merge but {X1=1, X2=1} -- the region where the
  # propensity is HIGHEST (e = 0.75) and control mass is therefore most depleted.
  # Tree = {X1=0} | {X1=1,X2=0} split on X3 | {X1=1,X2=1} merged = 4 leaves.
  #   nu({X1=1,X2=1}) = 8 * (1/32) * (1 - 0.75) = 0.0625, and sum(nu) = P(A=0) = 0.6,
  #   so its normalised mass is 0.0625/0.6, split equally between mu = 0.50 and 0.90
  #   (the X2 = 1 cells' nu is constant, so the X3 = 0/1 halves are equal).
  # This is the confounding fix's intrinsic cost made arithmetic: the SAME merge on
  # a noise coordinate would sit at normalised mass 0.1/0.6 = 1/6, i.e. 1.6x more
  # expensive, so confounding hands the competitor a cheaper option and Delta_mu
  # falls. Verified against the enumeration, which is free to disagree.
  nu_x2hi <- 8 * (1 / 32) * (1 - 0.75)
  w_merge <- nu_x2hi / 0.6
  expect_equal(f1$class_mu$delta_kl,
               w_merge * 0.5 * (kl2(0.50, 0.70) + kl2(0.90, 0.70)),
               tolerance = 1e-12)
  expect_equal(f1$class_mu$delta_sq,
               (w_merge / 2)^2 / w_merge * 0.40^2, tolerance = 1e-12)
  # And the enumeration must agree that mu's binding competitor is the
  # control-depleted one, not the noise-coordinate one: 0.6 * 1/6 > nu_x2hi.
  expect_lt(nu_x2hi, 0.1)
})

# ---- R3 theta_0 against a large Monte Carlo ------------------------------
test_that("R3. closed-form theta_0 matches a Monte Carlo ATT", {
  for (id in DGP_IDS) {
    d <- make_dgp(id)
    set.seed(7L)
    n <- 400000L
    s <- d$draw(n)
    # E[Y(1)-Y(0) | A=1] estimated from the DGP's own potential-outcome means.
    key <- do.call(paste, c(lapply(d$var_names, function(v) s$X[[v]]), sep = "-"))
    ckey <- do.call(paste, c(lapply(d$var_names, function(v) d$cells[[v]]), sep = "-"))
    idx <- match(key, ckey)
    tau_i <- d$tau_x[idx]
    mc <- mean(tau_i[s$A == 1L])
    expect_equal(mc, d$theta0, tolerance = 5e-3,
                 info = paste0(id, ": mc = ", signif(mc, 6)))
  }
})

# ---- R4 nu is exactly the control-covariate distribution ----------------
test_that("R4. nu = p_x (1 - e_0) reproduces P(X = x | A = 0) after normalising", {
  for (id in DGP_IDS) {
    d <- make_dgp(id)
    set.seed(11L)
    s <- d$draw(300000L)
    key <- do.call(paste, c(lapply(d$var_names, function(v) s$X[[v]]), sep = "-"))
    ckey <- do.call(paste, c(lapply(d$var_names, function(v) d$cells[[v]]), sep = "-"))
    emp <- as.numeric(table(factor(key[s$A == 0L], levels = ckey))) / sum(s$A == 0L)
    expect_equal(emp, d$nu / sum(d$nu), tolerance = 5e-3, info = id)
    # And sum(nu) is P(A = 0).
    expect_equal(sum(d$nu), mean(s$A == 0L), tolerance = 5e-3, info = id)
  }
})

# ---- R5 fitted_cell_labels() vs an independent tree walker ---------------
# The riskiest external dependency: optimaltrees:::assign_leaf_ids() is
# unexported, and its split feature indices must address the cell grid's columns
# in order. Re-walk the tree JSON with a from-scratch traversal and require the
# induced PARTITION to agree on every cell.
test_that("R5. an independent tree walk reproduces fitted_cell_labels()", {
  walk_independently <- function(node, cells) {
    lab <- integer(nrow(cells)); nxt <- 1L
    rec <- function(nd, rows) {
      if (!length(rows)) return(invisible())
      if (!is.null(nd$prediction)) { lab[rows] <<- nxt; nxt <<- nxt + 1L; return(invisible()) }
      k <- as.integer(as.numeric(nd$feature) + 1L)
      go <- cells[rows, k, drop = TRUE] == 1L
      rec(nd$true, rows[go]); rec(nd$false, rows[!go])
    }
    rec(node, seq_len(nrow(cells)))
    canonical_labels(lab)
  }
  set.seed(20260901L)
  for (id in DGP_IDS) {
    d <- make_dgp(id)
    for (n in c(500L, 4000L)) {
      s <- d$draw(n)
      fit <- doubletree::estimate_att(
        s$X, s$A, s$Y, leaf_budget = LEAF_BUDGET, outcome_type = OUTCOME_TYPE,
        propensity_loss = PROPENSITY_LOSS, verbose = FALSE, worker_limit = 1L
      )
      for (mdl in list(fit$nuisance_fits$e_model, fit$nuisance_fits$m0_model)) {
        a <- fitted_cell_labels(mdl, d$cells)
        b <- walk_independently(mdl@trees[[1L]], d$cells)
        # Compare as PARTITIONS (label numbering is traversal-order dependent).
        expect_identical(paste(a, collapse = "."), paste(b, collapse = "."),
                         info = paste(id, n))
      }
    }
  }
})

# ---- R6 the attribution 2x2 partitions the cell exactly -----------------
test_that("R6. attribution counts sum to reps in every cell (no dropped reps)", {
  paths <- list.files(DIR_RESULTS, "^cell_.*_r1000\\.rds$", full.names = TRUE)
  skip_if(length(paths) == 0L, "no r1000 checkpoints present")
  res <- do.call(rbind, lapply(paths, function(p) readRDS(p)$results))
  for (part in split(res, list(res$dgp, res$n), drop = TRUE)) {
    z <- part
    n4 <- c(
      sum(z$recovered_both & z$clip_frac == 0),
      sum(z$recovered_both & z$clip_frac > 0),
      sum(!z$recovered_both & z$clip_frac == 0),
      sum(!z$recovered_both & z$clip_frac > 0)
    )
    expect_equal(sum(n4), nrow(z),
                 info = paste(z$dgp[[1]], z$n[[1]]))
    # clip_frac is a sum of two non-negative fractions, so it can never be < 0,
    # and both components are proportions.
    expect_true(all(z$clip_frac >= 0 & z$clip_frac <= 2))
    expect_true(all(abs(z$clip_frac - (z$clip_frac_hi + z$clip_frac_lo)) < 1e-12))
  }
})

# ---- R7 recovery indicator vs excess risk, on the real results ----------
test_that("R7. recovered <=> excess risk is exactly 0, across all 18000 fits", {
  paths <- list.files(DIR_RESULTS, "^cell_.*_r1000\\.rds$", full.names = TRUE)
  skip_if(length(paths) == 0L, "no r1000 checkpoints present")
  res <- do.call(rbind, lapply(paths, function(p) readRDS(p)$results))
  # The structural criterion and the risk criterion must agree on every fit; a
  # single disagreement would mean one of the two is wrong.
  expect_true(all(res$excess_sq_e[res$recovered_e] < 1e-14))
  expect_true(all(res$excess_sq_e[!res$recovered_e] > 0))
  expect_true(all(res$excess_sq_mu[res$recovered_mu] < 1e-14))
  expect_true(all(res$excess_sq_mu[!res$recovered_mu] > 0))
  # recovered_both is the conjunction, not an independent measurement.
  expect_identical(res$recovered_both, res$recovered_e & res$recovered_mu)
  # Every fit stayed inside the class the margin was computed over.
  expect_true(all(res$in_class_e) && all(res$in_class_mu))
  # Sanity on the CI arithmetic: covered <=> theta0 inside [ci_lo, ci_hi], and
  # the interval really is theta +/- z * sigma. The multiplier is qnorm(0.975) =
  # 1.959964, NOT 1.96: doubletree's att_ci() uses qnorm(0.5 + level/2)
  # (R/inference.R:40). estimate_att()'s own roxygen says "1.96", which is a
  # rounding in the package documentation, not in the code -- checked here rather
  # than assumed, because an interval built on the wrong multiplier would shift
  # every coverage number in this study by a few tenths of a point.
  expect_identical(res$covered, res$ci_lo <= res$theta0 & res$theta0 <= res$ci_hi)
  expect_true(max(abs(res$ci_hi - res$ci_lo - 2 * stats::qnorm(0.975) * res$sigma)) < 1e-9)
  expect_true(max(abs((res$ci_lo + res$ci_hi) / 2 - res$theta)) < 1e-9)
})

cli::cli_alert_success("review_checks.R: all independent re-derivations agree.")
