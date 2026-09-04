# ============================================================
# calibrate_margin.R
# Study: partition_recovery_clt  (doubletree, S1)
# Spec:  quality_reports/specs/2026-09-01_partition-recovery-clt.md  §8 addendum
#        ("New required deliverable")
#
# DETERMINISTIC. Closed-form population enumeration only: no replications, no
# model fits, no random numbers anywhere. Running this file twice gives
# byte-identical output.
#
# WHY IT EXISTS. The leaf values pinned in common.R's make_dgp() are the single
# most consequential design choice in this study -- they set Delta_j, and
# Delta_j versus lambda_n = log(n)/n decides, per cell, whether a selection
# failure is STRUCTURAL (the penalised population optimum is outside S_j, so
# failure needs no sampling noise at all) or genuine noise. The spec's §8
# addendum and the design audit
# (quality_reports/2026-09-01_partition-recovery-clt-design-audit.md §3) both
# quote specific margin numbers. Without this script those numbers are mystery
# constants that rot silently the moment the DGP or the n grid changes. With it
# they are reproducible in ~20 seconds.
#
# WHAT IT PRINTS
#
#   1  CONFOUNDING CHECK. E[Y(0) | A = 1] vs E[Y(0) | A = 0] in closed form, and
#      the resulting naive difference-in-means bias, for every DGP. This is the
#      deterministic test of the audit §2 / addendum item 1 defect: on the old
#      disjoint-covariate F1 these two were EQUAL, the naive bias was 0, and no
#      partition-recovery failure could move theta_hat -- which made this study's
#      central metric pairing (spec §4) vacuous in its second branch.
#
#   2  ADOPTED PARAMETERS. Every (Delta_sq, Delta_kl) the study actually runs on,
#      read from common.R's make_dgp() rather than re-typed here, each with its
#      lambda_n-crossing n -- both the exact real-valued crossing and the first
#      cell of N_GRID that clears it.
#
#   3  THE F1 SHAPE SURFACE. Delta_e / Delta_mu over candidate propensity shapes
#      at fixed mu, showing the trade-off the confounding fix creates: sharing
#      the X1 gate depletes control mass exactly where the propensity is high,
#      and mu's fine split lives inside {X1 = 1}, so a WIDE e shrinks Delta_mu.
#      The lever is e's spread, not mu's.
#
#   4  THE ST2 GAP SURFACE. Delta_j over candidate near-tie gaps on the adopted
#      F1 shape, with each gap's crossing n, so the ST2 choice is visible as a
#      choice between "penalty-dominated at every n in the grid" (uninformative:
#      empirical failure pinned near 1 everywhere, nothing for
#      prop:selection-rate's decay shape to show) and "clean everywhere"
#      (uninformative the other way: no visible near-tie effect at all).
#
# Run:  Rscript --vanilla simulations/partition_recovery_clt/code/calibrate_margin.R
#       (from the doubletree package root)
#
# Gates (all defaulted; a bare invocation prints everything):
#   PRC_CAL_SURFACES  "0" to print only sections 1-2 (adopted params)   [1]
# ============================================================

source(file.path("simulations", "partition_recovery_clt", "code", "common.R"))

SHOW_SURFACES <- !identical(Sys.getenv("PRC_CAL_SURFACES", "1"), "0")
cli::cli_inform("gates: PRC_CAL_SURFACES={SHOW_SURFACES}")

# ---- 1 helpers -----------------------------------------------------------

#' Exact n at which lambda_n = log(n)/n falls below `delta`
#'
#' log(n)/n is strictly decreasing for n >= 3, so the crossing is unique there.
#' Returned as a real number (not rounded to a grid point) so two margins that
#' both "clear at n >= 1000" can still be ordered against each other.
#'
#' @param delta A margin. Must be positive.
#' @param upper Search ceiling; returns \code{Inf} if the crossing is beyond it.
#' @return The real-valued crossing n, or \code{Inf}.
lambda_crossing_n <- function(delta, upper = 1e9) {
  stopifnot(length(delta) == 1L, is.finite(delta), delta > 0)
  f <- function(n) log(n) / n - delta
  if (f(upper) > 0) return(Inf)
  if (f(3) < 0) return(3)
  stats::uniroot(f, interval = c(3, upper), tol = 1e-6)$root
}

#' First sample size in the study's own n grid whose lambda_n clears `delta`
#'
#' @return An integer from \code{N_GRID}, or \code{NA_integer_} if none does.
first_clean_n <- function(delta, n_grid = N_GRID) {
  ok <- log(n_grid) / n_grid < delta
  if (!any(ok)) return(NA_integer_)
  as.integer(n_grid[[which(ok)[[1]]]])
}

#' Closed-form confounding diagnostic for a DGP spec
#'
#' Everything from the population cell probabilities; no draw. \code{naive_bias}
#' is the bias of the difference in means as an estimator of theta_0, which is
#' EXACTLY ZERO when A is independent of Y(0) -- the signature of the disjoint
#' covariate defect this study had to fix.
confounding_summary <- function(spec) {
  p <- spec$cell_prob
  e <- spec$e0
  data.frame(
    dgp = spec$id,
    e_vars = paste(spec$e_vars, collapse = "+"),
    mu_vars = paste(spec$mu_vars, collapse = "+"),
    shared_vars = paste(intersect(spec$e_vars, spec$mu_vars), collapse = "+"),
    theta0 = spec$theta0,
    # E[Y(0) | A = a] = sum_x P(x) P(A = a | x) mu_0(x) / P(A = a)
    EY0_treated = sum(p * e * spec$mu0) / sum(p * e),
    EY0_control = sum(p * (1 - e) * spec$mu0) / sum(p * (1 - e)),
    # E[Y | A = 1] - E[Y | A = 0] - theta_0
    naive_bias = (sum(p * e * spec$mu1) / sum(p * e) -
                    sum(p * (1 - e) * spec$mu0) / sum(p * (1 - e))) - spec$theta0,
    e0_min = min(e), e0_max = max(e),
    stringsAsFactors = FALSE
  )
}

#' One row per (dgp, nuisance): margins in both currencies plus crossing n
margin_table <- function(spec) {
  cls <- list(e = spec$class_e, mu = spec$class_mu)
  do.call(rbind, lapply(names(cls), function(j) {
    cl <- cls[[j]]
    data.frame(
      dgp = spec$id, nuisance = j,
      card_T = spec$enumeration$n_partitions,
      n_suff = cl$n_sufficient,
      min_suff_leaves = cl$min_sufficient_leaves,
      n_values = cl$n_distinct_values,
      delta_sq = cl$delta_sq,
      delta_kl = cl$delta_kl,
      crossing_n_exact = lambda_crossing_n(cl$delta_kl),
      clean_from_n = first_clean_n(cl$delta_kl),
      tightest_competitor = cl$argmin_key_kl,
      stringsAsFactors = FALSE
    )
  }))
}

# lambda_n over the study's own grid, printed once so every crossing claim below
# can be checked by eye against it.
lambda_grid <- data.frame(n = N_GRID, lambda_n = log(N_GRID) / N_GRID)

# ---- 2 report: adopted parameters ---------------------------------------

cli::cli_h1("lambda_n = log(n)/n over the study's n grid")
print(lambda_grid, digits = 5, row.names = FALSE)

specs <- lapply(DGP_IDS, make_dgp)
names(specs) <- DGP_IDS

cli::cli_h1("1. Confounding check (closed form; naive_bias == 0 means A indep Y(0))")
conf <- do.call(rbind, lapply(specs, confounding_summary))
rownames(conf) <- NULL
print(conf, digits = 5, row.names = FALSE)
for (id in DGP_IDS) {
  shared <- intersect(specs[[id]]$e_vars, specs[[id]]$mu_vars)
  if (length(shared) == 0L) {
    cli::cli_abort(c(
      "DGP {.val {id}} has DISJOINT nuisance supports, so A is independent of Y(0).",
      i = "This is the defect recorded in the design audit §2 and the spec's §8 addendum item 1: no partition-recovery failure can move theta_hat, which makes spec §4's metric pairing vacuous. Give the two nuisances a shared covariate."
    ))
  }
}
cli::cli_alert_success("every DGP has at least one shared nuisance covariate")

cli::cli_h1("2. Adopted parameters: margins and lambda_n crossings")
adopted <- do.call(rbind, lapply(specs, margin_table))
rownames(adopted) <- NULL
print(adopted, digits = 5, row.names = FALSE)

cli::cli_inform(c(
  "i" = "Read `clean_from_n` as: from this n onward the PENALISED population optimum lies inside S_j, so a selection failure there is sampling noise. Below it, failure is structural.",
  "i" = "`n_suff` > 1 is what the currency-correct recovery indicator exists for; `n_suff` == 1 (ST3) makes lem:uniform's within-S_j variation check vacuous, not failed."
))

#' One line per DGP naming the numeric parameters actually in force
#'
#' Zips names to values from the SAME filtered list, so a non-numeric entry
#' (st3's `source` string) cannot shift the labels off by one.
print_adopted_params <- function(spec) {
  num <- spec$params[vapply(spec$params, is.numeric, logical(1))]
  cli::cli_inform(paste0(
    spec$id, ": ",
    paste(names(num), signif(unlist(num), 4), sep = " = ", collapse = ", ")
  ))
}
for (id in DGP_IDS) print_adopted_params(specs[[id]])

if (!SHOW_SURFACES) {
  cli::cli_alert_success("calibrate_margin.R done (surfaces suppressed).")
} else {

# ---- 3 surfaces: why these parameters and not others --------------------

# One enumeration, reused across every candidate. The enumeration depends only on
# (p, Lbar), never on the leaf values, so recomputing it per candidate would
# multiply the cost of this script by ~40 for no change in the answer.
EN5 <- enumerate_tree_partitions(5L, 4L)

#' Margins for one candidate (e-shape, mu-shape) without re-enumerating
#'
#' Named (not an inline lambda) so a single candidate can be re-run in a live
#' session: candidate_margins(0.20, 0.45, 0.30, 0.10, 0.50, 0.40).
candidate_margins <- function(e_lo, e_mid, e_gap, mu_lo, mu_mid, mu_gap) {
  s <- dgp_spec_grid_sparse(
    p = 5L, e_lo = e_lo, e_mid = e_mid, e_gap = e_gap,
    mu_lo = mu_lo, mu_mid = mu_mid, mu_gap = mu_gap,
    tau = 0.08, leaf_budget = LEAF_BUDGET, id = "candidate"
  )
  check_variable_roles(s)
  ce <- sufficient_class(s$e0, s$cell_prob, EN5)
  cm <- sufficient_class(s$mu0, s$nu, EN5)
  data.frame(
    e = sprintf("%.2f/%.2f/%.2f", e_lo, e_mid, e_mid + e_gap),
    mu = sprintf("%.2f/%.2f/%.2f", mu_lo, mu_mid, mu_mid + mu_gap),
    e_gap = e_gap, mu_gap = mu_gap,
    delta_sq_e = ce$delta_sq, delta_kl_e = ce$delta_kl,
    delta_sq_mu = cm$delta_sq, delta_kl_mu = cm$delta_kl,
    clean_from_e = first_clean_n(ce$delta_kl),
    clean_from_mu = first_clean_n(cm$delta_kl),
    cross_e = lambda_crossing_n(ce$delta_kl),
    cross_mu = lambda_crossing_n(cm$delta_kl),
    n_suff_e = ce$n_sufficient, n_suff_mu = cm$n_sufficient,
    naive_bias = confounding_summary(s)$naive_bias,
    stringsAsFactors = FALSE
  )
}

cli::cli_h1("3. F1 shape surface: e's spread is the lever on Delta_mu")
# mu held at the widest spread available subject to mu_hi + tau <= 1; e narrowed
# progressively. Narrower e keeps control mass inside {X1 = 1}, where mu's fine
# split lives, so Delta_mu RISES as e narrows even though mu is untouched.
f1_shapes <- list(
  c(0.10, 0.50, 0.40),   # widest e: the pre-fix pinned shape
  c(0.15, 0.475, 0.35),
  c(0.20, 0.45, 0.30),   # the audit's §3(a) recommendation
  c(0.25, 0.425, 0.25),
  c(0.30, 0.40, 0.20)
)
f1_surface <- do.call(rbind, lapply(f1_shapes, function(v) {
  candidate_margins(v[[1]], v[[2]], v[[3]], 0.10, 0.50, 0.40)
}))
print(f1_surface[, c("e", "mu", "delta_kl_e", "delta_kl_mu", "clean_from_e",
                     "clean_from_mu", "cross_e", "cross_mu", "n_suff_e",
                     "n_suff_mu", "naive_bias")],
      digits = 4, row.names = FALSE)

cli::cli_h1("4. ST2 gap surface on the ADOPTED F1 shape")
# Both gaps move together, as spec §3's ST2 construction states ("structurally
# identical to F1 except one candidate partition outside S_j is deliberately
# built to have population risk within a small, explicit distance"). The adopted
# F1 e-shape is read back from make_dgp() so this surface cannot drift out of
# sync with the F1 cell the study actually runs.
f1_prm <- specs[["f1"]]$params
st2_gaps <- c(0.34, 0.30, 0.26, 0.22, 0.18, 0.14, 0.10)
st2_surface <- do.call(rbind, lapply(st2_gaps, function(g) {
  # e_mid must stay put while e_gap shrinks, so the X1 gate difference is held
  # fixed and only the X2 split -- the split the margin is measured on -- narrows.
  candidate_margins(f1_prm$e_lo, f1_prm$e_mid, g,
                    f1_prm$mu_lo, f1_prm$mu_mid, g)
}))
print(st2_surface[, c("e_gap", "delta_kl_e", "delta_kl_mu", "clean_from_e",
                      "clean_from_mu", "cross_e", "cross_mu")],
      digits = 4, row.names = FALSE)

cli::cli_inform(c(
  "i" = paste(
    "An informative ST2 needs the failure rate to CROSS inside the study's own n",
    "grid: a gap whose {.field clean_from_*} is NA for both nuisances is",
    "penalty-dominated everywhere (failure pinned near 1, nothing for the",
    "decay-shape comparison to show); a gap clean from the smallest n has no",
    "visible near-tie effect at all."
  ),
  "i" = "The adopted ST2 gap is the one in common.R's make_dgp(); its row above is the justification."
))

cli::cli_alert_success("calibrate_margin.R done.")
}
