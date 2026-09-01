# ============================================================
# run_cell.R
# Study: threshold_superconsistency  (doubletree)
#
# Runs ONE sample-size cell and writes ONE checkpoint. Deliberately one cell per
# PROCESS, invoked from the shell:
#
#     Rscript code/run_cell.R <n> <reps>
#
# Why one process per cell rather than a loop over n inside a single session:
# this codebase has previously exhausted host memory during tree simulations
# (see doubletree/simulations/docs/MEMORY_SAFE_SIMULATIONS.md and MEMORY.md's
# [LEARN:rashomon-memory]). Process-per-cell means (a) nothing accumulates
# across cells -- the OS reclaims everything at exit, (b) a crash costs at most
# one cell, and (c) the staged rollout (smallest cell first, inspect, then step
# up) is enforced by how the script is invoked rather than by discipline inside
# a loop.
#
# MEMORY RULES OBSERVED HERE (all non-negotiable for this study):
#   * NO Rashomon enumeration anywhere. fit_tree() -> optimaltrees(single_tree =
#     TRUE), which is the single-best-tree path; fit_rashomon() is never called.
#     Rashomon set enumeration was the documented cause of every previous
#     memory-exhaustion incident, not ordinary single-tree fitting.
#   * NO parallelism. worker_limit = 1, no mclapply/furrr/future, BLAS and OpenMP
#     thread counts pinned to 1 by the caller.
#   * max_depth = 2 explicitly. Single-tree REGRESSION with continuous covariates
#     explores an exponential space without a depth cap -- optimaltrees itself
#     auto-caps at 2 for this path (treefarms.R ~L783) after a probe showed
#     n=50,p=3 timing out unbounded. Setting it explicitly makes the cap
#     deterministic instead of dependent on a heuristic, and 2 is exactly the
#     depth of the true tree.
#   * gc(full = TRUE) after every replicate, and memory recorded per replicate.
#   * Checkpoint written to disk at the end of the cell.
#
# Sourced-and-callable: run_cell() is a plain function with working defaults, so
# the file can be sourced in a live session and stepped through. The invocation
# at the bottom fires only for a genuine top-level Rscript run.
# ============================================================

suppressPackageStartupMessages({
  library(optimaltrees)
  library(cli)
})

STUDY <- "/Users/dagniel/RAND/rprojects/global-scholars/doubletree/simulations/threshold_superconsistency"
source(file.path(STUDY, "code", "dgp.R"))
source(file.path(STUDY, "code", "stage_b.R"))

# ---- 0 configuration ----------------------------------------------------
# Single visible block; every value has a working default so sourcing this file
# yields a runnable configuration.

SEED_BASE <- 20260831L

#' Stage A discretization resolution: the COARSE default schedule
#'
#' ceiling(n^(1/3)) is the "adaptive" bin schedule optimaltrees shipped as its
#' default before the 2026-08-31 rho change (rho = 1/3 -> 0.45). It is passed
#' EXPLICITLY rather than via discretize_bins = "adaptive" for two reasons:
#'   1. Independence. The rho default change is a separate piece of work; if this
#'      simulation read the default it would silently change meaning when rho
#'      changed, and "Stage A at the current coarse resolution" is precisely the
#'      condition the check is about.
#'   2. Feasibility. At rho = 0.45 and n = 32000 the grid is ~106 thresholds per
#'      coordinate (318 binary features) instead of 32 (93 features), which is a
#'      materially larger search for the solver and exactly the kind of blow-up
#'      the memory rules exist to avoid.
#' The whole point of Stage B is to beat this mesh WITHOUT refining it.
coarse_bins <- function(n) as.integer(ceiling(n^(1 / 3)))

#' Regularization: the theory rate lambda ~ log(n)/n on the mean-loss scale
#'
#' Must be small enough that the OPTIMAL tree contains both true splits. The
#' second split's population gain is constant in n (see
#' population_split_gains()); lambda shrinks, so recovery gets easier with n and
#' n = 2000 is the binding case. Checked explicitly in run_cell().
reg_for_n <- function(n) log(n) / n

# ---- 1 memory instrumentation -------------------------------------------

# phys_footprint probe. Compiled once at startup.
#
# WHY NOT `ps`: on this host `ps` is unavailable to this session entirely
# ("operation not permitted", and the process list cannot be read at all), so any
# RSS-via-subprocess measurement is not an option here. It would be the wrong
# metric anyway: macOS RSS EXCLUDES compressed pages, so a process's RSS can fall
# while its real memory charge climbs -- the exact blind spot behind the
# 2026-08-04 incident. phys_footprint (resident + compressed + IOKit) is what
# jetsam bands on. Probe copied from single_tree_coverage/code/proc_footprint.cpp.
#
# WHY NOT gc() ALONE: gc() measures only R's vector/cons heap. The solver
# allocates natively (new / std::vector) and never touches R's allocator, so the
# R heap can sit flat at tens of MB while the process grows by gigabytes. Both
# numbers are recorded because they answer different questions.
FOOTPRINT_OK <- FALSE
local({
  ok <- tryCatch({
    Rcpp::sourceCpp(file.path(STUDY, "code", "proc_footprint.cpp"))
    TRUE
  }, error = function(e) {
    cli::cli_warn(c(
      "phys_footprint probe failed to compile; memory reporting is DEGRADED to the R heap only.",
      "i" = conditionMessage(e)
    ))
    FALSE
  })
  FOOTPRINT_OK <<- ok
})

#' Process phys_footprint in MB, or NA when the probe is unavailable
proc_footprint_mb <- function() {
  if (!FOOTPRINT_OK) return(c(now = NA_real_, lifetime_max = NA_real_))
  fi <- proc_footprint(as.integer(Sys.getpid()))
  if (fi[["ok"]] != 1) return(c(now = NA_real_, lifetime_max = NA_real_))
  c(now = fi[["footprint"]] / 1024^2, lifetime_max = fi[["lifetime_max"]] / 1024^2)
}

#' R heap usage in MB reported directly by gc()
#'
#' gc() returns megabyte columns for the cons cells and the vector heap; their
#' sum is the R-side footprint. Column names differ across R builds, so the Mb
#' columns are located by name match rather than by position.
gc_mb <- function() {
  g <- gc(verbose = FALSE, full = TRUE)
  mb_cols <- grep("Mb", colnames(g), fixed = TRUE)
  used_col <- mb_cols[[1L]]   # first Mb column is "used"
  sum(g[, used_col])
}

# ---- 2 one replicate ----------------------------------------------------

#' Fit Stage A + Stage B on one replicate and extract the threshold diagnostics
#'
#' @param n Sample size.
#' @param rep Replicate index (recorded, not used for seeding -- the seed is set
#'   once per cell).
#' @return One-row data frame.
one_rep <- function(n, rep) {
  dat <- generate_jump_tree(n)
  X <- dat$X; y <- dat$Y
  n_bins <- coarse_bins(n)
  lambda <- reg_for_n(n)

  t0 <- Sys.time()
  # STAGE A: single best tree from the existing solver. single_tree = TRUE is
  # forced by fit_tree(); no Rashomon set is ever built.
  fit <- fit_tree(
    X, y,
    loss_function       = "squared_error",
    regularization      = lambda,
    discretize_method   = "quantiles",
    discretize_bins     = n_bins,
    max_depth           = 2L,
    worker_limit        = 1L,
    store_training_data = TRUE,
    verbose             = FALSE
  )
  t_stage_a <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (fit@n_trees != 1L) {
    cli::cli_abort("Expected exactly 1 tree from fit_tree(), got {fit@n_trees}.")
  }

  md <- fit@discretization_metadata
  lk <- bin_lookup(md)
  ct_raw <- as_coord_tree(fit@trees[[1L]], lk, md)

  # Collapse any grid-adjacent transition-leaf pairs (prop:transition-leaves)
  # BEFORE Stage B refines anything -- the collapsed node's bracket is what
  # Stage B must search, not the two original, mutually-narrowing grid cuts.
  collapse_out <- collapse_transitions(ct_raw)
  ct <- collapse_out$tree
  n_collapses <- collapse_out$n_collapses

  # STAGE B: refine every threshold inside its own identifying bracket.
  t1 <- Sys.time()
  ref <- refine_tree(ct, X, y)
  t_stage_b <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

  # ---- diagnostics -------------------------------------------------------
  # Root split: expected to be on X1 (its population gain dominates every
  # alternative by an order of magnitude). Recorded, and asserted below only via
  # the topology flag so a failure shows up in the results rather than killing
  # the cell.
  root_coord <- ct$coord
  root_is_x1 <- identical(root_coord, "X1")

  tab <- ref$refined
  # After collapse there should be at most ONE X2 split. More than one is a
  # surprise the pipeline must not silently paper over by picking the biggest.
  x2_rows <- tab[tab$coord == "X2", , drop = FALSE]
  if (nrow(x2_rows) > 1L) {
    cli::cli_abort("Expected at most 1 X2 split post-collapse, got {nrow(x2_rows)}.")
  }
  has_x2 <- nrow(x2_rows) == 1L
  x2_row <- if (has_x2) x2_rows[1L, , drop = FALSE] else NULL

  root_row <- tab[tab$path == "", , drop = FALSE]
  mesh <- 1 / n_bins

  data.frame(
    n = n, rep = rep,
    n_bins = n_bins, lambda = lambda,
    n_leaves = count_leaves(ct),
    n_collapses = n_collapses,
    root_coord = root_coord,
    # Exact, not >=: after collapse a correctly recovered topology has
    # EXACTLY 3 leaves (root split + one child split). >= would silently
    # accept an uncollapsed transition pair (4 leaves) as "ok".
    topology_ok = root_is_x1 && has_x2 && count_leaves(ct) == 3L,

    # split 1 (root, X1) -- the headline diagnostic
    t1_star   = T1_STAR,
    t1_grid   = root_row$grid_cut[[1L]],
    t1_refined = root_row$refined_cut[[1L]],
    err1_stage_a = abs(root_row$grid_cut[[1L]]    - T1_STAR),
    err1_stage_b = abs(root_row$refined_cut[[1L]] - T1_STAR),
    n_cand1 = root_row$n_candidates[[1L]],
    collapsed1 = root_row$collapsed[[1L]],
    mesh_moved1 = abs(root_row$refined_cut[[1L]] - root_row$grid_cut[[1L]]) / mesh,
    # t1_lo/t1_hi are now the REALISED candidate range Stage B actually
    # searched (cand_lo/cand_hi from refine_tree()), not a fixed grid window --
    # the identifying bracket itself (node_bracket()'s lo/hi) can be +/-Inf and
    # is logged separately for inspection, not for this containment check.
    t1_lo = root_row$cand_lo[[1L]],
    t1_hi = root_row$cand_hi[[1L]],
    t1_br_lo = root_row$lo[[1L]],
    t1_br_hi = root_row$hi[[1L]],
    bracket1_contains_truth = isTRUE(T1_STAR > root_row$cand_lo[[1L]] && T1_STAR < root_row$cand_hi[[1L]]),

    # split 2 (X2, inside the left child)
    t2_star    = T2_STAR,
    t2_grid    = if (has_x2) x2_row$grid_cut[[1L]]    else NA_real_,
    t2_refined = if (has_x2) x2_row$refined_cut[[1L]] else NA_real_,
    err2_stage_a = if (has_x2) abs(x2_row$grid_cut[[1L]]    - T2_STAR) else NA_real_,
    err2_stage_b = if (has_x2) abs(x2_row$refined_cut[[1L]] - T2_STAR) else NA_real_,
    collapsed2 = if (has_x2) x2_row$collapsed[[1L]] else NA,
    mesh_moved2 = if (has_x2) abs(x2_row$refined_cut[[1L]] - x2_row$grid_cut[[1L]]) / mesh else NA_real_,
    t2_lo = if (has_x2) x2_row$cand_lo[[1L]] else NA_real_,
    t2_hi = if (has_x2) x2_row$cand_hi[[1L]] else NA_real_,
    bracket2_contains_truth = if (has_x2) {
      isTRUE(T2_STAR > x2_row$cand_lo[[1L]] && T2_STAR < x2_row$cand_hi[[1L]])
    } else NA,
    n_node2 = if (has_x2) x2_row$n_node[[1L]] else NA_integer_,

    sse_stage_a = tree_sse(ct_raw, X, y),
    sse_after_collapse_pre_refine = tree_sse(ct, X, y),
    sse_stage_b = tree_sse(ref$tree, X, y),
    secs_stage_a = t_stage_a,
    secs_stage_b = t_stage_b,
    stringsAsFactors = FALSE
  )
}

# ---- 3 one cell --------------------------------------------------------

#' Run one sample-size cell sequentially and checkpoint it
#'
#' @param n Sample size (default 2000, the smallest cell -- so a bare source +
#'   run_cell() is the micro-test, not the biggest job).
#' @param reps Replicates (default 2, micro-test size).
#' @param out Checkpoint path.
#' @return The cell's results data frame, invisibly.
run_cell <- function(n = 2000L, reps = 2L,
                     out = file.path(STUDY, "results",
                                     sprintf("cell_n%d_r%d.rds", n, reps))) {
  stopifnot(is.numeric(n), length(n) == 1L, n >= 100)
  stopifnot(is.numeric(reps), length(reps) == 1L, reps >= 1L)

  # Seed ONCE per cell process, before any random draw.
  set.seed(SEED_BASE + as.integer(n))

  gains <- population_split_gains()
  lambda <- reg_for_n(n)
  cli::cli_h1("cell n = {n}, reps = {reps}")
  cli::cli_inform(c(
    "*" = "coarse grid: {coarse_bins(n)} bins -> {coarse_bins(n) - 1L} thresholds/coord, mesh ~ {signif(1/coarse_bins(n), 3)}",
    "*" = "lambda = log(n)/n = {signif(lambda, 4)}",
    "*" = "population gain split1 = {signif(gains$gain_split1, 4)}, split2 = {signif(gains$gain_split2, 4)}",
    "*" = "second split is worth including iff gain2 > lambda: {gains$gain_split2 > lambda} (ratio {signif(gains$gain_split2/lambda, 3)})"
  ))
  if (gains$gain_split2 <= lambda) {
    cli::cli_warn(paste0(
      "lambda exceeds the second split's population gain at n = ", n,
      "; the OPTIMAL tree will not contain split 2 and topology recovery must fail."
    ))
  }

  rows <- vector("list", reps)
  mem <- vector("list", reps)
  for (r in seq_len(reps)) {
    rows[[r]] <- one_rep(n, r)
    # Memory bookkeeping AFTER a full gc, once per replicate.
    heap <- gc_mb()
    fp <- proc_footprint_mb()
    mem[[r]] <- data.frame(n = n, rep = r,
                           footprint_mb = fp[["now"]],
                           footprint_lifetime_max_mb = fp[["lifetime_max"]],
                           r_heap_mb = heap)
    cli::cli_inform(paste0(
      "rep ", r, "/", reps,
      "  leaves=", rows[[r]]$n_leaves,
      " topo_ok=", rows[[r]]$topology_ok,
      " collapses=", rows[[r]]$n_collapses,
      " bracket1_ok=", rows[[r]]$bracket1_contains_truth,
      "  |t1hat-t1*| A=", signif(rows[[r]]$err1_stage_a, 3),
      " B=", signif(rows[[r]]$err1_stage_b, 3),
      "  ", round(rows[[r]]$secs_stage_a, 1), "s+", round(rows[[r]]$secs_stage_b, 1), "s",
      "  footprint=", round(fp[["now"]]), "MB (peak ", round(fp[["lifetime_max"]]),
      "MB) heap=", round(heap), "MB"
    ))
  }

  res <- do.call(rbind, rows)
  memdf <- do.call(rbind, mem)
  payload <- list(
    results = res,
    memory = memdf,
    meta = list(
      n = n, reps = reps, seed = SEED_BASE + as.integer(n),
      n_bins = coarse_bins(n), lambda = lambda,
      kappa = KAPPA, sigma = SIGMA, t1_star = T1_STAR, t2_star = T2_STAR,
      max_depth = 2L, use_rashomon = FALSE, worker_limit = 1L,
      loss_function = "squared_error",
      bin_schedule = "explicit ceiling(n^(1/3)) (pre-2026-08-31 coarse default)",
      footprint_probe_ok = FOOTPRINT_OK,
      r_version = R.version.string,
      optimaltrees_version = as.character(utils::packageVersion("optimaltrees")),
      run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
  )
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(payload, out)
  cli::cli_alert_success("wrote {out}")
  cli::cli_inform(paste0(
    "cell summary: mean|t1hat-t1*| Stage A = ", signif(mean(res$err1_stage_a), 4),
    ", Stage B = ", signif(mean(res$err1_stage_b), 4),
    "; peak phys_footprint = ", round(max(memdf$footprint_lifetime_max_mb, na.rm = TRUE)), " MB"
  ))
  invisible(res)
}

# ---- 4 script entry ----------------------------------------------------
# `!interactive()` alone is NOT enough: it is also FALSE when another script
# sources this one under Rscript. sys.nframe() == 0L is the discriminator.
if (!interactive() && sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  n_arg    <- if (length(args) >= 1L) as.integer(args[[1L]]) else 2000L
  reps_arg <- if (length(args) >= 2L) as.integer(args[[2L]]) else 2L
  if (is.na(n_arg) || is.na(reps_arg)) {
    cli::cli_abort("Usage: Rscript run_cell.R <n> <reps>")
  }
  run_cell(n = n_arg, reps = reps_arg)
}
