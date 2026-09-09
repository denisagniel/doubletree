# ============================================================
# slurm/units.R
# Study: head_to_head_comparison  (doubletree, S5)
# Spec:  quality_reports/specs/2026-09-08_head-to-head-comparison.md  (§4, §6)
#
# THE JOB-ARRAY MAP, AND NOTHING ELSE. Sourced by slurm/run_single_replication.R,
# slurm/combine_results.R and (through Rscript) by the launch scripts, so the
# array index -> unit-of-work mapping is defined EXACTLY ONCE. Two independent
# copies of this arithmetic is the failure that contaminated the 2026-07-10 run of
# the six-approach-arbitration study (its slurm/array.slurm records it): the
# submitter and the worker disagreed about which block a task owned, and every
# array silently recomputed and mislabelled the same first block.
#
# REQUIRES common.R to have been sourced first (REGIMES, REGIME_IDS, NSIM_FULL,
# N_GRID_FULL, DIR_RESULTS, cell_path). It defines objects and functions only;
# sourcing it has no side effects.
#
# ------------------------------------------------------------
# WHAT ONE UNIT OF WORK IS, AND WHY IT IS NOT ONE REPLICATION
# ------------------------------------------------------------
#
# One unit = one (regime, dgp, n, CONTIGUOUS BLOCK OF REPLICATIONS). Not one
# replication, and not one whole cell.
#
#   * Not one replication. The cheap cells cost 7 ms/rep (R6 at n = 200); a SLURM
#     task has seconds of module-load and R-startup overhead, so one-rep tasks
#     would spend more wall time starting R than estimating anything, and R6 alone
#     would submit 24,000 tasks to compute 14 minutes of arithmetic.
#
#   * Not one whole cell. README.md §10's projection is that ONE cell -- R2 at
#     n = 8000, nsim = 2000 -- costs 11.6 h by itself, and the whole sweep is
#     ~54 h with 77 % of it in R2 and R5. Cell-level parallelism therefore floors
#     the wall clock at the slowest cell and wastes the array on 37 tasks whose
#     costs span four orders of magnitude.
#
# Blocks are sized by COST, not by count: reps_per_batch = target_secs / (measured
# secs per rep at this n), so a task is ~TARGET_SECS of work whether it holds 86
# replications of R2 at n = 8000 or all 2000 replications of R6 at n = 200.
#
# BATCHING IS EXACTLY REPRODUCIBLE AGAINST THE SEQUENTIAL SWEEP. rep_seed() is a
# pure function of (SEED_MASTER, regime, dgp, n, rep) -- see common.R section 5 --
# so replication 137 of a cell draws byte-identical data whether it is the 137th
# replication of a sequential run_cell() or the 7th replication of batch 3 here.
# Splitting a cell across tasks changes nothing about its results. This is the
# property that makes the whole design safe, and it is why the seed scheme must
# not be changed to anything sequential (a `set.seed(1); for (r in ...)` stream
# would make batch boundaries visible in the numbers).
#
# ------------------------------------------------------------
# THE COST MODEL
# ------------------------------------------------------------
#
# Anchors are the MEASURED per-replication costs (summed over the cell's arms)
# from the 2026-09-09 pilot, read off the projection tables the pilot itself wrote
# to tables/pilot_projection_*.csv, and scaled LINEARLY in n from the n at which
# each was measured -- the same extrapolation run_pilot.R section 4 uses, and the
# same one whose totals README.md §10 quotes (R1 0.94 h, R2 22.8 h, R3 11.6 h,
# R4 0.45 h, R5 18.1 h, R6 0.23 h; 54.3 h in total).
#
# THE MODEL IS KNOWN TO UNDER-STATE THE CROSSFIT-HEAVY REGIMES, and that is
# recorded rather than smoothed over: estimate_att_crossfit()'s nested CV is worse
# than linear in n, so R2's n = 8000 and R5's n = 2000 batches are the ones most
# likely to over-run. Two mechanisms absorb it instead of one estimate being
# trusted: walltime_for() requests SAFETY_FACTOR x the projected batch cost, and
# run_single_replication.R flushes a partial every FLUSH_EVERY replications, so a
# task that is killed at the wall limit loses at most FLUSH_EVERY replications and
# resubmitting the same run-id resumes from the flush.
# ============================================================

if (!exists("REGIME_IDS", inherits = TRUE)) {
  stop("slurm/units.R requires code/common.R to be sourced first.", call. = FALSE)
}

# ---- 1 the measured cost anchors ----------------------------------------

#' Per-replication cost anchors from the 2026-09-09 pilot
#'
#' \code{secs} is the wall-clock cost of ONE replication of one cell of that
#' regime -- summed over every arm the regime runs -- measured at sample size
#' \code{at_n}. Source: \code{tables/pilot_projection_*.csv}, latest run per
#' regime (R5's anchor is the 20260909-071859 run, the one that covers both of
#' R5's DGPs; the others are 20260909-070128).
#'
#' These are DATA, not tuning knobs. Re-running the pilot supersedes them; edit
#' the numbers here and nothing else changes.
COST_ANCHOR <- list(
  R1 = list(secs = 0.904764, at_n = 2000L),
  R2 = list(secs = 5.236222, at_n = 2000L),
  R3 = list(secs = 2.094853, at_n =  500L),
  R4 = list(secs = 0.464356, at_n = 1000L),
  R5 = list(secs = 6.529410, at_n =  500L),
  R6 = list(secs = 0.018556, at_n =  500L)
)

#' Target wall-clock cost of ONE array task, in seconds
#'
#' 1800 s (30 min) by default. The trade-off is explicit: smaller targets shorten
#' the wall clock (more tasks run concurrently) and pay more R-startup overhead
#' per unit of work; larger targets do the reverse. 30 min keeps every task far
#' inside O2's `short` partition limit while holding the array to a few hundred
#' tasks, and keeps the overhead of one R startup (~10 s with both packages) under
#' 1 % of a task.
#'
#' @return Numeric seconds. Override with \code{H2H_TARGET_SECS}.
h2h_target_secs <- function() {
  v <- Sys.getenv("H2H_TARGET_SECS", "")
  if (!nzchar(v)) return(1800)
  x <- suppressWarnings(as.numeric(v))
  if (is.na(x) || x <= 0) {
    cli::cli_abort("H2H_TARGET_SECS = {.val {v}} is not a positive number.")
  }
  x
}

# FLOOR on the number of replications between partial flushes inside a task.
# run_single_replication.R raises it to ceiling(reps / 20) when that is larger, so
# no block writes more than ~20 partials regardless of how many replications it
# holds; see the comment there for the measured cost of a fixed cadence.
FLUSH_EVERY <- 25L

# Walltime requested per task = SAFETY_FACTOR x projected cost + WALLTIME_FLOOR_SECS,
# rounded up to the next whole minute. SAFETY_FACTOR is 3 rather than a tighter
# number for the reason stated in the header: the linear-in-n cost model
# under-states estimate_att_crossfit(), so the tasks most likely to over-run are
# exactly the expensive ones whose loss hurts most.
SAFETY_FACTOR <- 3
WALLTIME_FLOOR_SECS <- 900   # 15 min of headroom for module load + R startup

#' Memory request per regime, in GB
#'
#' NOT PROFILED ON O2 -- stated so nobody reads these as measured. The pilot ran
#' on the dev box without RSS instrumentation. They are conservative: every arm
#' runs with \code{worker_limit = 1L}, no Rashomon enumeration and single-threaded
#' BLAS (common.R section 1), so the live footprint is a handful of n x p matrices
#' plus one tree fit; 4 GB is already ~100x the data. The crossfit regimes (R2,
#' R5) and the large-n regime get more because K = 5 nested CV holds several fold
#' fits at once.
#'
#' README_O2.md's monitoring section asks for `sacct -o MaxRSS` on the first wave
#' and for these numbers to be revised from that, rather than left as guesses.
MEM_GB <- c(R1 = 4L, R2 = 8L, R3 = 6L, R4 = 4L, R5 = 6L, R6 = 4L)

# O2's `short` partition allows up to 12 h; every projected batch here is well
# under 2 h, so no regime needs medium/long. Kept as a per-regime vector anyway,
# because a re-pilot that raises a cost anchor should be able to move ONE regime.
PARTITION <- c(R1 = "short", R2 = "short", R3 = "short",
               R4 = "short", R5 = "short", R6 = "short")

# ---- 2 the cost model ----------------------------------------------------

#' Projected seconds for one replication of one cell
#'
#' @param regime One of \code{REGIME_IDS}.
#' @param n Sample size.
#' @return Numeric seconds per replication (all arms of the regime).
secs_per_rep <- function(regime, n) {
  a <- COST_ANCHOR[[regime]]
  if (is.null(a)) cli::cli_abort("No cost anchor for regime {.val {regime}}.")
  a$secs * (as.numeric(n) / a$at_n)
}

# ---- 3 the unit table ----------------------------------------------------

#' The array-index -> unit-of-work map for ONE regime
#'
#' Row \code{i} is what SLURM array task \code{i} of that regime's array must run.
#' Per-regime rather than global on purpose: each regime gets its own sbatch array
#' with its own \code{--time}/\code{--mem}, and task ids that restart at 1 inside
#' each array need no offset arithmetic to be correct. (six-approach-arbitration
#' carries a global unit index plus \code{--unit-offset}/\code{--max-unit} flags to
#' reconcile it, and its own comments record the contamination bug that the
#' reconciliation existed to fix. Per-regime tables remove the class of bug rather
#' than guarding against it.)
#'
#' Batch boundaries are EVEN: after n_batches is fixed by the target cost, the
#' block size is recomputed as ceiling(reps / n_batches) so the last task is not a
#' stub of 3 replications while its siblings hold 300.
#'
#' @param regime One of \code{REGIME_IDS}.
#' @param target_secs Target wall cost per task; see \code{h2h_target_secs()}.
#' @return A data.frame, one row per array task, ordered by (dgp, n, batch):
#'   \code{task}, \code{regime}, \code{dgp}, \code{n}, \code{reps_total},
#'   \code{batch}, \code{n_batches}, \code{rep_from}, \code{rep_to},
#'   \code{reps}, \code{est_secs}.
h2h_unit_table <- function(regime, target_secs = h2h_target_secs()) {
  if (!regime %in% REGIME_IDS) {
    cli::cli_abort("Unknown regime {.val {regime}}; expected {.val {REGIME_IDS}}.")
  }
  spec <- REGIMES[[regime]]
  reps_total <- NSIM_FULL[[regime]]
  n_grid <- N_GRID_FULL[[regime]]

  rows <- list()
  # dgp-major, then n ascending: the same order design_grid()/run_sweep.R use, so
  # a task id read out of a log lines up with the sweep's own cell listing.
  for (dgp in spec$dgps) {
    for (n in n_grid) {
      spr <- secs_per_rep(regime, n)
      per_batch0 <- max(1L, min(reps_total, as.integer(floor(target_secs / spr))))
      n_batches <- as.integer(ceiling(reps_total / per_batch0))
      per_batch <- as.integer(ceiling(reps_total / n_batches))
      for (b in seq_len(n_batches)) {
        from <- (b - 1L) * per_batch + 1L
        to <- min(b * per_batch, reps_total)
        if (from > to) next   # possible only if rounding produced a spare batch
        rows[[length(rows) + 1L]] <- data.frame(
          regime = regime, dgp = dgp, n = as.integer(n),
          reps_total = as.integer(reps_total),
          batch = b, n_batches = n_batches,
          rep_from = from, rep_to = to, reps = to - from + 1L,
          est_secs = (to - from + 1L) * spr,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  out$task <- seq_len(nrow(out))
  rownames(out) <- NULL
  out[, c("task", "regime", "dgp", "n", "reps_total", "batch", "n_batches",
          "rep_from", "rep_to", "reps", "est_secs")]
}

#' Per-regime SLURM sizing derived from the unit table
#'
#' @param regime One of \code{REGIME_IDS}.
#' @param target_secs Target wall cost per task.
#' @return A one-row data.frame: \code{n_tasks}, \code{max_task_secs},
#'   \code{total_hours}, \code{walltime} (HH:MM:SS), \code{mem_gb},
#'   \code{partition}.
h2h_sizing <- function(regime, target_secs = h2h_target_secs()) {
  ut <- h2h_unit_table(regime, target_secs)
  max_secs <- max(ut$est_secs)
  data.frame(
    regime = regime,
    n_tasks = nrow(ut),
    max_task_secs = max_secs,
    total_hours = sum(ut$est_secs) / 3600,
    walltime = walltime_for(max_secs),
    mem_gb = MEM_GB[[regime]],
    partition = PARTITION[[regime]],
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

# ---- 4 file naming -------------------------------------------------------

#' Scratch filename for one task's results shard
#'
#' Carries the full cell identity plus the batch index, so a directory listing is
#' readable and \code{combine_results.R} can group shards into cells without a
#' side-car index. \code{reps_total} is in the name because \code{cell_path()}
#' encodes it too, and the combined file must land at exactly the path
#' \code{analyze.R::load_results()} globs for.
#'
#' @param regime,dgp,n,reps_total,batch Unit identity.
#' @return A bare filename (no directory).
batch_filename <- function(regime, dgp, n, reps_total, batch) {
  sprintf("batch_%s_%s_n%d_r%d_b%s.rds",
          regime, dgp, as.integer(n), as.integer(reps_total),
          formatC(as.integer(batch), width = 4, flag = "0"))
}

#' Partial-checkpoint filename for a task in flight
#'
#' Same stem as the shard with a \code{part_} prefix, so an interrupted task's
#' leftovers are obvious in a listing and cannot be mistaken for a finished shard.
partial_filename <- function(regime, dgp, n, reps_total, batch) {
  paste0("part_", batch_filename(regime, dgp, n, reps_total, batch))
}
