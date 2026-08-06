# ============================================================
# harness.R
# Study: 2026-08-04_single-tree-coverage  (doubletree)
# Purpose: Shared simulation harness for the coverage + structure-recovery study.
#   Defines the per-rep seed scheme, the local-feasibility guard, and run_one_rep()
#   (the single-replication driver that wraps estimate_att_single_tree and records
#   every endpoint). Sourced by run_pilot.R and run_coverage.R AFTER dgps.R and
#   oracle_theta_star.R (which supply THETA0, N_GRID, DGP_NAMES).
# Not run directly.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# Base study seed. Set once at the top of each runner for reproducibility; every
# replication then draws its own deterministic seed via make_seed() (independent,
# reproducible, parallel-safe -- the standard simulation pattern).
BASE_SEED       <- 20260804L
MC_SEED_MODULUS <- 2147483647L   # 2^31 - 1 (valid 32-bit set.seed range)

# --- Local feasibility guard --------------------------------------------------
# Cells that OOM locally under the CURRENT optimaltrees Rashomon enumeration
# (TreeFARMS holds the entire Rashomon set in memory). Empirically (guarded probe,
# 2026-08-04, 16 GB host): the `continuous` DGP explodes the propensity Rashomon
# set (continuous covariates -> many quantile-threshold binary features -> ~1e5
# trees, >10 GB) at EVERY n, killing the probe's parent process. The binary DGPs
# (simple/moderate/complex) have a FIXED feature count (binary covariates create no
# discretization bins). That bounds the Rashomon TREE COUNT but NOT total memory:
# TreeFARMS cost is ~ |Rashomon set| x n, because it retains per-sample structures
# for every tree in the set. So a SINGLE call stays ~0.3 GB / 10-27 s to n = 10000,
# but that native RSS is not reclaimable by gc() and ACCUMULATES across calls within
# one process -- and each rep makes K = 5 folds x 2 nuisances = 10 calls.
#
# CORRECTION (2026-08-05, supersedes the earlier "binary DGPs stay <0.3 GB to
# n = 10000" reading of this note): `complex` at n = 10000 with 5 reps per process
# (~50 calls) was SIGKILLed repeatedly on 2026-08-04 and the OOM killer escalated to
# the driver process, taking the 16 GB host down. Per-CALL cost was never the
# binding constraint; per-PROCESS accumulation is. n = 10000 is therefore run one rep
# per subprocess under a driver-side RSS budget (see pilot_complex_driver.R).
# This matches the sister study's INFEASIBLE_CELLS (six-approach-arbitration,
# sacct 2026-08-03: 29-32 GB OOM for single_tree x continuous, all n).
#
# These cells require a high-memory cluster (see setup-cluster-simulations). They
# are SKIPPED locally and recorded as status = "skipped_infeasible_local" (distinct
# from converged = FALSE, which is an empty intersection on a cell that DID run).
LOCAL_INFEASIBLE_DGPS <- c("continuous")

#' Is a (dgp, n) cell feasible to run on the local (16 GB) host?
#'
#' @param dgp DGP name.
#' @param n Sample size (reserved for future finer rules; currently DGP-level).
#' @return Logical scalar.
is_feasible_local <- function(dgp, n) {
  !(dgp %in% LOCAL_INFEASIBLE_DGPS)
}

#' Deterministic per-(dgp, n, rep) seed, independent of cell ordering
#'
#' Distinct (dgp, n, rep) triples get distinct seeds as long as rep_id <= max_reps.
#' Arithmetic is done in doubles then reduced mod 2^31-1 (base-R integer arithmetic
#' overflows to NA past 2^31), giving a valid set.seed() argument.
#'
#' @param dgp,n,rep_id Cell + replication identity.
#' @param dgp_names,n_grid Reference orderings (default study grid).
#' @param max_reps Replication headroom per cell (seed spacing).
#' @param base Base study seed.
#' @return Integer seed in [0, 2^31-1).
make_seed <- function(dgp, n, rep_id, dgp_names = DGP_NAMES, n_grid = N_GRID,
                      max_reps = 1e6, base = BASE_SEED) {
  di <- match(dgp, dgp_names) - 1L
  ni <- match(n, n_grid) - 1L
  offset <- ((as.double(di) * length(n_grid) + ni) * max_reps) + rep_id
  as.integer((base + offset) %% MC_SEED_MODULUS)
}

# --- Structural / projection diagnostics -------------------------------------
# Added 2026-08-05. Rationale: the study previously recorded structure recovery
# only as `converged` (non-empty Rashomon intersection), and read a "recovery
# plateaus at 0.5-0.6 on complex" signal off it. That is the WRONG event to
# measure. Validity does not require recovering the true partition tau*; it
# requires the selected partition to be FINE ENOUGH to represent the truth. If
# tau_hat REFINES tau*, the true nuisance is constant within every selected leaf,
# so the population cell value equals the truth for ANY weighting, the projection
# error is exactly 0, and the working-model estimand equals theta_0 regardless of
# WHICH refinement was chosen. Only a COARSENING (a leaf that merges cells with
# genuinely different truth) breaks validity. Exact recovery and over-refinement
# are both benign, and the old metric cannot tell them apart from a coarsening.
#
# These functions therefore measure the theory's actual quantity: the L2 error of
# the population projection of the TRUE nuisance onto the selected partition.
# Because the DGPs return e_true and mu0_true exactly, this is computable without
# any asymptotic approximation.

#' Approximation error of the population projection onto a selected partition
#'
#' The cell values are the within-leaf means of the TRUE nuisance (not of the
#' data), i.e. an exact evaluation of the population projection this partition
#' induces. The residual is then evaluated over ALL observations, because the
#' relevant norm in the doubly-robust bias bound is L2(P) over the full covariate
#' distribution -- not over the subsample defining the cell means.
#'
#' @param truth Numeric vector of TRUE nuisance values, length n.
#' @param leaf Character vector of leaf identities, length n.
#' @param fit_on Logical vector selecting the observations that DEFINE the cell
#'   means. Propensity: all n (e_hat averages A over every unit in the leaf).
#'   Outcome: controls only (m0_hat averages Y over A == 0), which is what
#'   realizes the (1 - e_0)-weighted projection of mu_0. Default all TRUE.
#' @return List: delta (L2 projection error), max_dev (sup-norm; 0 iff the truth
#'   is constant within every leaf, i.e. the partition refines tau*), n_leaves,
#'   min_leaf_n, and n_cells_star (number of distinct true values = |tau*|, since
#'   tau* is the partition induced by the truth's level sets).
projection_diagnostics <- function(truth, leaf, fit_on = NULL) {
  n <- length(truth)
  if (is.null(fit_on)) fit_on <- rep(TRUE, n)
  na_out <- list(delta = NA_real_, max_dev = NA_real_, n_leaves = NA_integer_,
                 min_leaf_n = NA_integer_, n_cells_star = NA_integer_)
  if (length(leaf) != n) return(na_out)

  ok <- fit_on & !is.na(leaf) & !is.na(truth)
  if (!any(ok)) return(na_out)

  cell  <- tapply(truth[ok], leaf[ok], mean)
  proj  <- unname(cell[leaf])          # NA for a leaf with no fit_on observations
  resid <- truth - proj

  list(
    delta        = sqrt(mean(resid^2, na.rm = TRUE)),
    max_dev      = max(abs(resid), na.rm = TRUE),
    n_leaves     = length(cell),
    min_leaf_n   = as.integer(min(table(leaf[ok]))),
    # Round before uniquing: level sets of the truth are exact by construction in
    # these DGPs, but expit/arithmetic can leave last-bit noise.
    n_cells_star = length(unique(round(truth, 12)))
  )
}

#' Classify a selected partition against the truth it must represent
#'
#' Replaces the binary exact-recovery flag with the taxonomy validity actually
#' turns on. Only "coarsening" is a validity failure.
#'
#' @param d Output of projection_diagnostics().
#' @param tol Sup-norm tolerance for "the truth is constant within every leaf".
#' @return "exact", "over_refinement", "coarsening", or NA.
classify_recovery <- function(d, tol = 1e-8) {
  if (is.na(d$max_dev)) return(NA_character_)
  if (d$max_dev > tol) return("coarsening")
  if (!is.na(d$n_cells_star) && d$n_leaves > d$n_cells_star) return("over_refinement")
  "exact"
}

#' Run one replication of the single-tree ATT estimator; record all endpoints
#'
#' Wraps estimate_att_single_tree in tryCatch. The estimator ERRORS when the
#' Rashomon intersection is empty (structural margin fails); this is a NON-recovery
#' event recorded as converged = FALSE (never dropped -- it is endpoint 1's
#' denominator). Coverage indicators are computed against the oracle theta*_n
#' (endpoint 2) and the nominal true ATT theta_0 (endpoint 3).
#'
#' @param dgp,n,rep_id Cell + replication identity.
#' @param theta_star Oracle working-model ATT for this (dgp, n) (coverage target).
#' @param theta0 True (nominal) ATT, default THETA0 = 0.15.
#' @return One-row tibble: identity + seed, converged, point estimates + diagnostic,
#'   CI bounds (honest/single/crossfit), and six coverage indicators.
run_one_rep <- function(dgp, n, rep_id, theta_star, theta0 = THETA0) {
  seed <- make_seed(dgp, n, rep_id)
  gen  <- get(paste0("generate_dgp_", dgp))
  set.seed(seed)
  d <- gen(n)

  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(
    estimate_att_single_tree(
      d$X, d$A, d$Y, K = 5, outcome_type = "binary",
      discretize_bins = "adaptive", seed = seed
    ),
    error = function(e) e
  )
  elapsed_s <- proc.time()[["elapsed"]] - t0

  base_row <- tibble::tibble(
    dgp = dgp, n = n, rep_id = rep_id, seed = seed,
    status = "run", theta_star = theta_star, theta0 = theta0,
    elapsed_s = elapsed_s
  )

  if (inherits(res, "error")) {
    # Empty intersection (margin fails) or other estimator error -> non-recovery.
    return(dplyr::bind_cols(base_row, tibble::tibble(
      converged = FALSE, error_msg = substr(conditionMessage(res), 1, 160),
      theta_single = NA_real_, theta_crossfit = NA_real_,
      delta = NA_real_, delta_over_se = NA_real_,
      ci_honest_lo = NA_real_, ci_honest_hi = NA_real_,
      ci_single_lo = NA_real_, ci_single_hi = NA_real_,
      ci_cf_lo = NA_real_, ci_cf_hi = NA_real_,
      covers_star_honest = NA, covers_star_single = NA, covers_star_cf = NA,
      covers_theta0_honest = NA, covers_theta0_single = NA, covers_theta0_cf = NA
    ), empty_structure_cols()))
  }

  # Structural endpoints: leaf budget (licenses the no-splitting argument) and the
  # projection error of the selected partition against the known truth.
  pd_e       <- projection_diagnostics(d$e_true,   res$leaf_e)
  pd_m0      <- projection_diagnostics(d$mu0_true, res$leaf_m0, fit_on = d$A == 0)
  pd_m0_unw  <- projection_diagnostics(d$mu0_true, res$leaf_m0)
  delta_prod <- pd_e$delta * pd_m0$delta

  structure_cols <- tibble::tibble(
    # Leaf budget: assump:leaves requires L = o(sqrt n) with balanced leaves.
    n_leaves_e = res$n_leaves_e, n_leaves_m0 = res$n_leaves_m0,
    min_leaf_n_e = res$min_leaf_n_e, min_leaf_n_m0 = res$min_leaf_n_m0,
    leaves_over_sqrtn = max(res$n_leaves_e, res$n_leaves_m0) / sqrt(n),
    # Projection error of each selected partition against the true nuisance.
    delta_e = pd_e$delta, delta_m0 = pd_m0$delta,
    max_dev_e = pd_e$max_dev, max_dev_m0 = pd_m0$max_dev,
    # The bias bound is a PRODUCT; this is the quantity that must be o(n^{-1/2}),
    # so rootn_delta_product is the direct empirical read on the rate condition.
    delta_product = delta_prod, rootn_delta_product = sqrt(n) * delta_prod,
    # Unweighted mu_0 projection: the gap to delta_m0 measures how much the
    # (1 - e_0)-weighted vs unweighted projection distinction actually costs. It is
    # exactly 0 whenever the partition refines tau*_mu (mu_0 constant in the leaf),
    # so a nonzero gap flags the regime where that distinction becomes load-bearing.
    delta_m0_unweighted = pd_m0_unw$delta,
    # Sufficiency: TRUE iff the truth is constant within every selected leaf.
    sufficient_e = pd_e$max_dev <= 1e-8, sufficient_m0 = pd_m0$max_dev <= 1e-8,
    recovery_class_e = classify_recovery(pd_e),
    recovery_class_m0 = classify_recovery(pd_m0),
    n_cells_star_e = pd_e$n_cells_star, n_cells_star_m0 = pd_m0$n_cells_star
  )

  in_ci <- function(v, ci) v >= ci[1] & v <= ci[2]
  dplyr::bind_cols(base_row, tibble::tibble(
    converged = TRUE, error_msg = NA_character_,
    theta_single = res$theta_single, theta_crossfit = res$theta_crossfit,
    delta = res$delta, delta_over_se = res$delta_over_se,
    ci_honest_lo = res$ci_95_honest[1], ci_honest_hi = res$ci_95_honest[2],
    ci_single_lo = res$ci_95_single[1], ci_single_hi = res$ci_95_single[2],
    ci_cf_lo = res$ci_95_crossfit[1], ci_cf_hi = res$ci_95_crossfit[2],
    covers_star_honest   = in_ci(theta_star, res$ci_95_honest),
    covers_star_single   = in_ci(theta_star, res$ci_95_single),
    covers_star_cf       = in_ci(theta_star, res$ci_95_crossfit),
    covers_theta0_honest = in_ci(theta0, res$ci_95_honest),
    covers_theta0_single = in_ci(theta0, res$ci_95_single),
    covers_theta0_cf     = in_ci(theta0, res$ci_95_crossfit)
  ), structure_cols)
}

#' All-NA structural columns, for rows where no tree was fit
#'
#' Kept in ONE place so run_one_rep's error branch and skipped_rep_row cannot
#' drift out of schema parity with the converged branch (the driver's merge warns
#' on column mismatch rather than erroring, so drift would be absorbed silently).
#'
#' @return One-row tibble of NA structural columns.
empty_structure_cols <- function() {
  tibble::tibble(
    n_leaves_e = NA_integer_, n_leaves_m0 = NA_integer_,
    min_leaf_n_e = NA_integer_, min_leaf_n_m0 = NA_integer_,
    leaves_over_sqrtn = NA_real_,
    delta_e = NA_real_, delta_m0 = NA_real_,
    max_dev_e = NA_real_, max_dev_m0 = NA_real_,
    delta_product = NA_real_, rootn_delta_product = NA_real_,
    delta_m0_unweighted = NA_real_,
    sufficient_e = NA, sufficient_m0 = NA,
    recovery_class_e = NA_character_, recovery_class_m0 = NA_character_,
    n_cells_star_e = NA_integer_, n_cells_star_m0 = NA_integer_
  )
}

#' Build a skipped-cell row (ran nothing; distinct from converged = FALSE)
#'
#' Used for BOTH skip reasons, which must stay distinguishable in analysis:
#'   "skipped_infeasible_local" -- cell never attempted (DGP-level guard).
#'   "skipped_oom_local"        -- cell WAS attempted and exceeded the memory cap.
#' Neither is converged = FALSE (an empty intersection on a cell that did run), and
#' neither may be silently dropped: they belong in endpoint 1's denominator.
#'
#' @param dgp,n,rep_id Cell + replication identity.
#' @param theta_star,theta0 Oracle + nominal targets (recorded for completeness).
#' @param status Skip status string.
#' @param msg Explanation recorded in error_msg.
#' @param n_grid Reference n-grid for seed spacing. Pass the FULL study grid explicitly
#'   when a caller has subset N_GRID (e.g. a smoke run), or the recorded seed will not
#'   match the seed the worker actually used for that rep.
#' @return One-row tibble mirroring run_one_rep's schema.
skipped_rep_row <- function(dgp, n, rep_id, theta_star, theta0 = THETA0,
                            status = "skipped_infeasible_local",
                            msg = "locally infeasible (Rashomon OOM); run on cluster",
                            n_grid = N_GRID) {
  tibble::tibble(
    dgp = dgp, n = n, rep_id = rep_id,
    seed = make_seed(dgp, n, rep_id, n_grid = n_grid),
    status = status, theta_star = theta_star, theta0 = theta0,
    elapsed_s = NA_real_,
    converged = NA, error_msg = msg,
    theta_single = NA_real_, theta_crossfit = NA_real_,
    delta = NA_real_, delta_over_se = NA_real_,
    ci_honest_lo = NA_real_, ci_honest_hi = NA_real_,
    ci_single_lo = NA_real_, ci_single_hi = NA_real_,
    ci_cf_lo = NA_real_, ci_cf_hi = NA_real_,
    covers_star_honest = NA, covers_star_single = NA, covers_star_cf = NA,
    covers_theta0_honest = NA, covers_theta0_single = NA, covers_theta0_cf = NA
  ) |>
    dplyr::bind_cols(empty_structure_cols())
}
