# ============================================================
# slurm/units.R
# Study: honest_inference_sparsity_failure  (doubletree, S2)
# Spec:  quality_reports/specs/2026-08-21_honest-inference-sparsity-failure.md
#        ("Sample sizes and reps")
#
# THE JOB-ARRAY MAP, AND NOTHING ELSE. Sourced by slurm/run_single_replication.R,
# slurm/print_sizing.R and slurm/combine_results.R, so the array index ->
# unit-of-work mapping is defined EXACTLY ONCE. Two independent copies of this
# arithmetic is the failure that contaminated the 2026-07-10 run of the
# six-approach-arbitration study (its slurm/array.slurm records it): the submitter
# and the worker disagreed about which block a task owned, and every array
# silently recomputed and mislabelled the same first block.
#
# REQUIRES code/common.R to have been sourced first (DGP_IDS, N_GRID,
# reps_for_n(), DIR_RESULTS, cell_path). It defines objects and functions only;
# sourcing it has no side effects.
#
# ------------------------------------------------------------
# WHAT THIS DEPLOYMENT IS FOR
# ------------------------------------------------------------
#
# The spec's grid is 5 DGP variants x n in {500, 2000, 8000}. As of 2026-09-09 the
# five n = 500 cells are complete on disk (results/cell_*_n500_r300.rds, written
# 2026-09-01). The TEN remaining cells -- 5 variants x n in {2000, 8000} -- are what
# this directory exists to run. n = 500 is still in the tables below, because a
# backfill or a re-run of a single n = 500 cell must be submittable through the same
# machinery rather than through a second, divergent path.
#
# REPLICATIONS PER CELL are read from code/common.R's reps_for_n(), NOT restated
# here: R = 300 at n in {500, 2000} and R = 150 at n = 8000. The reduction at
# n = 8000 is the spec's OWN sanctioned exception ("consider dropping the n = 8000
# cell to a smaller R if timing is prohibitive, stating the reduced R explicitly
# rather than silently"), stated at REPS_DEFAULT/REPS_AT_MAX_N in common.R, in
# run_sweep.R's header, in every checkpoint's metadata and in analyze.R's output.
# Overriding HIS_REPS/HIS_REPS_MAX_N changes the array map here automatically,
# because this file calls reps_for_n() rather than hard-coding 300/150.
#
# ------------------------------------------------------------
# WHAT ONE UNIT OF WORK IS, AND WHY IT IS NOT ONE REPLICATION
# ------------------------------------------------------------
#
# One unit = one (dgp, n, CONTIGUOUS BLOCK OF REPLICATIONS). Not one replication,
# and not one whole cell.
#
#   * Not one replication. A SLURM task pays seconds of module-load and R-startup
#     overhead. At n = 2000 one replication costs 3.8 s, so one-rep tasks would
#     spend a comparable amount of wall time starting R as estimating anything, and
#     would submit 2,250 tasks to compute 4 hours of arithmetic.
#
#   * Not one whole cell. Cell-level parallelism floors the wall clock at the
#     slowest cell (28 min at n = 8000) and wastes the array on 10 tasks.
#
# Blocks are sized by COST: reps_per_batch = target_secs / (measured secs per rep at
# this n), so a task is ~TARGET_SECS of work whether it holds 150 replications at
# n = 2000 or 75 at n = 8000.
#
# BATCHING IS EXACTLY REPRODUCIBLE AGAINST THE SEQUENTIAL SWEEP. rep_seed() is a
# pure function of (SEED_MASTER, dgp, n, rep) -- see common.R section 3 -- so
# replication 137 of a cell draws byte-identical data whether it is the 137th
# replication of a sequential run_cell() or the 62nd replication of batch 2 here.
# Splitting a cell across tasks changes nothing about its results, and the five
# n = 500 cells already on disk remain directly comparable to whatever this
# deployment produces. This is the property that makes the whole design safe, and
# it is why the seed scheme must not be changed to anything sequential (a
# `set.seed(1); for (r in ...)` stream would make batch boundaries visible in the
# numbers).
#
# ------------------------------------------------------------
# THE COST MODEL: PER-n MEASUREMENTS, NOT AN EXTRAPOLATION
# ------------------------------------------------------------
#
# head_to_head_comparison/slurm/units.R anchors one measured cost per regime and
# scales it LINEARLY in n. That is deliberately NOT done here, because this study's
# own cost is measurably sublinear in n over the range that matters: 3.83 s/rep at
# n = 2000 against 11.25 s/rep at n = 8000 is a 2.94x cost for a 4x sample size, so
# a linear extrapolation from the n = 2000 anchor would predict 15.3 s/rep at
# n = 8000 and over-reserve by 36 %. Both of the n values this deployment needs were
# measured DIRECTLY, so there is nothing to extrapolate and the model is a lookup.
#
# 98 % of the cost is the ANCHOR, not the method under test: estimate_att_crossfit()
# runs nested CV (K = 5 outer folds x an adaptive lambda search x cv_K = 5 inner
# folds x 2 nuisances) against a single pair of tree fits for estimate_att().
# Measured 2026-09-09 through the installed-library path: at n = 8000, secs_full
# ~ 0.07 and secs_cf ~ 11.2. Any future re-timing should therefore be read as a
# measurement of estimate_att_crossfit(), and a change in its CV defaults is the one
# package change that invalidates this table.
#
# Two mechanisms absorb an under-estimate instead of one number being trusted:
# walltime_for() requests SAFETY_FACTOR x the projected batch cost, and
# run_single_replication.R flushes a partial every FLUSH_EVERY replications, so a
# task killed at the wall limit loses at most FLUSH_EVERY replications and
# resubmitting the same run-id resumes from the flush.
# ============================================================

if (!exists("DGP_IDS", inherits = TRUE) || !exists("reps_for_n", inherits = TRUE)) {
  stop("slurm/units.R requires code/common.R to be sourced first.", call. = FALSE)
}

# ---- 1 waves: one SLURM array per sample size ----------------------------

# ONE ARRAY PER n, not one array for the sweep, and not one per (dgp, n) cell.
#
# Per-n is the right grain for THIS study because n is the only axis the cost and
# the replication count vary along: all five DGP variants are the same size at the
# same n (same p, same leaf_budget, same estimator pair -- the variants differ only
# in eps and in which coordinate the residual lives on, neither of which changes
# the fit cost), while n changes the per-replication cost 3x and the replication
# count from 300 to 150. A single --time/--mem sizing across n would either time
# out the n = 8000 tasks or reserve an hour for n = 2000 tasks that finish in ten
# minutes.
#
# Per-n arrays also make the task-id -> block map offset-free: each array's task ids
# restart at 1 and index directly into that wave's own unit table.
# (six-approach-arbitration carries a GLOBAL unit index plus
# --unit-offset/--max-unit flags to reconcile it, and its own comments record the
# contamination bug that the reconciliation existed to fix. Per-wave tables remove
# the class of bug rather than guarding against it.)

#' Wave id for a sample size, e.g. 2000 -> "n2000"
#'
#' A bare shell-safe token, because it is a job name, a directory name and an
#' exported environment variable, not just an R value.
wave_id <- function(n) paste0("n", as.integer(n))

#' Sample size of a wave id, e.g. "n2000" -> 2000L
wave_n <- function(wave) {
  if (!grepl("^n[0-9]+$", wave)) {
    cli::cli_abort("{.val {wave}} is not a wave id; expected {.val n<integer>}.")
  }
  as.integer(sub("^n", "", wave))
}

# Every wave the design grid admits, in the grid's own order.
WAVE_IDS <- vapply(N_GRID, wave_id, character(1), USE.NAMES = FALSE)

# The waves this deployment exists to run: the cells the spec asks for that are NOT
# already on disk. n500 is deliberately absent from the DEFAULT while remaining a
# valid argument to launch_subset.sh, so a backfill is possible without editing
# code but is never accidental.
WAVES_MISSING <- setdiff(WAVE_IDS, wave_id(500L))

# ---- 2 the measured cost anchors ----------------------------------------

#' Per-replication wall cost, measured, by sample size
#'
#' One replication = one \code{estimate_att()} call PLUS one
#' \code{estimate_att_crossfit()} call on the same draw, which is what
#' \code{run_one_rep()} does.
#'
#' Provenance:
#'   n =  500  1.90 s  -- the study's own 2026-09-01 pilot, as recorded in the note
#'                        above REPS_DEFAULT in code/common.R.
#'   n = 2000  3.83 s  -- measured 2026-09-09, 3 replications of dgp `eps_sev`
#'                        through the installed-library path (3.75, 4.16, 3.57).
#'   n = 8000 11.25 s  -- measured 2026-09-09, same protocol (10.79, 11.36, 11.59).
#'
#' These are DATA, not tuning knobs. Re-measuring supersedes them; edit the numbers
#' here and nothing else changes. The 2026-09-09 figures reproduce the 2026-09-01
#' pilot's 3.8 / 11.1 to within 1-2 %, on a different package-loading path, which is
#' the reason they are trusted enough to size an array from three replications.
COST_ANCHOR_BY_N <- c(`500` = 1.90, `2000` = 3.83, `8000` = 11.25)

#' Target wall-clock cost of ONE array task, in seconds
#'
#' 900 s (15 min) by default -- HALF head_to_head_comparison's 1800 s, on purpose.
#' The whole remaining sweep here is only ~3.9 h of compute, so a 1800 s target
#' would map each of the ten cells onto exactly one task and floor the wall clock at
#' the slowest single cell (28 min). At 900 s every cell splits in two, the array
#' doubles to 20 tasks, and the longest task drops to ~14 min. The cost is one extra
#' R startup per cell (~10 s with both packages, under 2 % of a task).
#'
#' @return Numeric seconds. Override with \code{HIS_TARGET_SECS}.
his_target_secs <- function() {
  v <- Sys.getenv("HIS_TARGET_SECS", "")
  if (!nzchar(v)) return(900)
  x <- suppressWarnings(as.numeric(v))
  if (is.na(x) || x <= 0) {
    cli::cli_abort("HIS_TARGET_SECS = {.val {v}} is not a positive number.")
  }
  x
}

# FLOOR on the number of replications between partial flushes inside a task.
# run_single_replication.R raises it to ceiling(reps / 20) when that is larger, so
# no block writes more than ~20 partials regardless of how many replications it
# holds. At the sizes this deployment produces (150 reps at n = 2000, 75 at
# n = 8000) the floor is what binds: 6 and 3 flushes per block respectively, each
# costing one rewrite of a <=150-row data.frame -- negligible against a 3.8 s
# replication, and it bounds a timeout's loss at 25 replications.
FLUSH_EVERY <- 25L

# Walltime requested per task = SAFETY_FACTOR x projected cost + WALLTIME_FLOOR_SECS,
# rounded up to the next whole minute. SAFETY_FACTOR is 3 rather than a tighter
# number even though the cost model here is a direct measurement rather than an
# extrapolation: the measurement was taken on the dev box, and an O2 compute node's
# per-core throughput is not the dev box's. Over-reserving inside the `short`
# partition is free; a timeout costs a resubmit.
SAFETY_FACTOR <- 3
WALLTIME_FLOOR_SECS <- 900   # 15 min of headroom for module load + R startup

#' Memory request per wave, in GB
#'
#' NOT PROFILED ON O2 -- stated so nobody reads these as measured. They are
#' conservative: both estimator calls run with \code{worker_limit = 1L}, no Rashomon
#' enumeration anywhere, \code{parallel = FALSE} on the cross-fit call and BLAS/OpenMP
#' pinned to one thread (code/common.R section 1), so the live footprint is a handful
#' of n x p matrices plus one tree fit. p = 5 binary covariates and n <= 8000 is a
#' trivially small design matrix; the reason n = 8000 gets more is
#' \code{estimate_att_crossfit()}'s K = 5 nested CV holding several fold fits at once.
#'
#' README_O2.md's monitoring section asks for `sacct -o MaxRSS` on the first wave and
#' for these numbers to be revised from that, rather than left as guesses.
MEM_GB_BY_WAVE <- c(n500 = 4L, n2000 = 6L, n8000 = 8L)

# O2's `short` partition allows up to 12 h; the longest projected task here is
# ~14 min and the longest REQUEST is ~58 min, so no wave needs medium/long. Kept as
# a per-wave vector anyway, because a re-measurement that raises one cost anchor
# should be able to move ONE wave.
PARTITION_BY_WAVE <- c(n500 = "short", n2000 = "short", n8000 = "short")

# ---- 3 the cost model ----------------------------------------------------

#' Projected seconds for one replication at a given sample size
#'
#' A LOOKUP where the size was measured (see \code{COST_ANCHOR_BY_N}); for any other
#' n it falls back to linear scaling from the nearest anchor and says so, rather than
#' silently inventing a number. Linear scaling is known to OVER-state cost above the
#' anchor for this study (the measured 2000 -> 8000 growth is 2.94x for a 4x n), so
#' the fallback errs toward over-reserving, which is the safe direction.
#'
#' @param n Sample size.
#' @return Numeric seconds per replication (both estimator calls).
secs_per_rep <- function(n) {
  n <- as.integer(n)
  key <- as.character(n)
  if (key %in% names(COST_ANCHOR_BY_N)) return(unname(COST_ANCHOR_BY_N[[key]]))
  anchor_n <- as.integer(names(COST_ANCHOR_BY_N))
  nearest <- which.min(abs(anchor_n - n))
  cli::cli_warn(c(
    "No measured cost anchor at n = {n}; scaling linearly from n = {anchor_n[nearest]}.",
    i = "Linear scaling OVER-states this study's cost above its anchor, so the walltime request errs long."
  ))
  unname(COST_ANCHOR_BY_N[[nearest]]) * (n / anchor_n[nearest])
}

# ---- 4 the unit table ----------------------------------------------------

#' The array-index -> unit-of-work map for ONE wave
#'
#' Row \code{i} is what SLURM array task \code{i} of that wave's array must run.
#'
#' Batch boundaries are EVEN: after \code{n_batches} is fixed by the target cost, the
#' block size is recomputed as \code{ceiling(reps / n_batches)}, so the last task is
#' not a stub of 3 replications while its siblings hold 150.
#'
#' Ordered dgp-major (in \code{DGP_IDS} order), then batch ascending -- the same DGP
#' order \code{design_grid()} and \code{run_sweep.R} use, so a task id read out of a
#' log lines up with the sweep's own cell listing.
#'
#' @param wave A wave id, e.g. \code{"n2000"}.
#' @param target_secs Target wall cost per task; see \code{his_target_secs()}.
#' @return A data.frame, one row per array task: \code{task}, \code{wave},
#'   \code{dgp}, \code{eps}, \code{n}, \code{reps_total}, \code{batch},
#'   \code{n_batches}, \code{rep_from}, \code{rep_to}, \code{reps}, \code{est_secs}.
his_unit_table <- function(wave, target_secs = his_target_secs()) {
  if (!wave %in% WAVE_IDS) {
    cli::cli_abort("Unknown wave {.val {wave}}; expected {.val {WAVE_IDS}}.")
  }
  n <- wave_n(wave)
  reps_total <- as.integer(reps_for_n(n))
  if (is.na(reps_total) || reps_total < 1L) {
    cli::cli_abort("reps_for_n({n}) = {reps_total}; check HIS_REPS / HIS_REPS_MAX_N.")
  }
  spr <- secs_per_rep(n)

  rows <- list()
  for (dgp in DGP_IDS) {
    per_batch0 <- max(1L, min(reps_total, as.integer(floor(target_secs / spr))))
    n_batches <- as.integer(ceiling(reps_total / per_batch0))
    per_batch <- as.integer(ceiling(reps_total / n_batches))
    for (b in seq_len(n_batches)) {
      from <- (b - 1L) * per_batch + 1L
      to <- min(b * per_batch, reps_total)
      if (from > to) next   # possible only if rounding produced a spare batch
      rows[[length(rows) + 1L]] <- data.frame(
        wave = wave, dgp = dgp, eps = DGP_SPECS[[dgp]]$eps,
        n = n, reps_total = reps_total,
        batch = b, n_batches = n_batches,
        rep_from = from, rep_to = to, reps = to - from + 1L,
        est_secs = (to - from + 1L) * spr,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out$task <- seq_len(nrow(out))
  rownames(out) <- NULL
  out[, c("task", "wave", "dgp", "eps", "n", "reps_total", "batch", "n_batches",
          "rep_from", "rep_to", "reps", "est_secs")]
}

#' Per-wave SLURM sizing derived from the unit table
#'
#' @param wave A wave id.
#' @param target_secs Target wall cost per task.
#' @return A one-row data.frame: \code{wave}, \code{n}, \code{reps_per_cell},
#'   \code{n_tasks}, \code{max_task_secs}, \code{total_hours}, \code{walltime}
#'   (HH:MM:SS), \code{mem_gb}, \code{partition}.
his_sizing <- function(wave, target_secs = his_target_secs()) {
  ut <- his_unit_table(wave, target_secs)
  max_secs <- max(ut$est_secs)
  data.frame(
    wave = wave,
    n = wave_n(wave),
    reps_per_cell = ut$reps_total[[1L]],
    n_tasks = nrow(ut),
    max_task_secs = max_secs,
    total_hours = sum(ut$est_secs) / 3600,
    walltime = walltime_for(max_secs),
    mem_gb = MEM_GB_BY_WAVE[[wave]],
    partition = PARTITION_BY_WAVE[[wave]],
    stringsAsFactors = FALSE
  )
}

#' Format a SLURM \code{--time} string with the standard safety margin
#'
#' @param secs Projected seconds for the most expensive task.
#' @return "HH:MM:SS".
walltime_for <- function(secs) {
  total <- ceiling((SAFETY_FACTOR * secs + WALLTIME_FLOOR_SECS) / 60) * 60
  sprintf("%02d:%02d:%02d", total %/% 3600L, (total %% 3600L) %/% 60L, total %% 60L)
}

# ---- 5 file naming -------------------------------------------------------

#' Scratch filename for one task's results shard
#'
#' Carries the full cell identity plus the batch index, so a directory listing is
#' readable and \code{combine_results.R} can group shards into cells without a
#' side-car index. \code{reps_total} is in the name because \code{cell_path()}
#' encodes it too, and the combined file must land at exactly the path
#' \code{analyze.R} constructs with \code{cell_path(dgp, n, reps_for_n(n))}.
#'
#' @param dgp,n,reps_total,batch Unit identity.
#' @return A bare filename (no directory).
batch_filename <- function(dgp, n, reps_total, batch) {
  sprintf("batch_%s_n%d_r%d_b%s.rds",
          dgp, as.integer(n), as.integer(reps_total),
          formatC(as.integer(batch), width = 4, flag = "0"))
}

#' Partial-checkpoint filename for a task in flight
#'
#' Same stem as the shard with a \code{part_} prefix, so an interrupted task's
#' leftovers are obvious in a listing and cannot be mistaken for a finished shard.
partial_filename <- function(dgp, n, reps_total, batch) {
  paste0("part_", batch_filename(dgp, n, reps_total, batch))
}
