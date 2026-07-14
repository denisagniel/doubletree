#!/usr/bin/env Rscript
# =============================================================================
# profile_timing.R -- LOCAL profiling to size the cluster job array
# =============================================================================
# Runs a handful of work units locally, measures median per-unit wall time and
# peak memory, then computes job sizing that targets a 1-3 hr wall time per task
# while respecting O2 limits:
#     * <= 1000 tasks per array          (MAX_ARRAY_SIZE)
#     * <= 10000 jobs queued at once      (MAX_CONCURRENT_JOBS)
#
# Writes two artifacts into config/:
#     sizing.json  -- human-readable record (git-tracked or logged)
#     sizing.env   -- KEY=VALUE, sourced by submit.sh and array.slurm (no jq needed)
#
# Run locally BEFORE submitting: profiling on the cluster wastes an allocation.
# Usage:
#   Rscript slurm/profile_timing.R [--n-units N] [--target-hours H] [--study-dir DIR]
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# --- O2 scheduler limits (named constants; edit only if O2 policy changes) ---
MAX_ARRAY_SIZE      <- 1000L    # max tasks in a single --array
MAX_CONCURRENT_JOBS <- 10000L   # max jobs queued across all arrays at once
WALL_MIN_HOURS      <- 1        # target window lower bound
WALL_MAX_HOURS      <- 3        # target window upper bound
# TIME_SAFETY multiplies the projected task time to set --time. 1.5x is too tight
# for this study: per-unit cost is heavy-tailed (Rashomon-set size varies widely at
# the stress DGP), so a task whose units land above the profiled median blows past a
# 1.5x allocation and times out. 3.5x absorbs that tail (reliability > packing; the
# 3 hr WALL_MAX cap still bounds it). Per-method .env walltimes are hand-tuned and
# compensate independently -- this only fixes the stock profiler default.
TIME_SAFETY         <- 3.5      # multiply estimated task time for --time headroom
MEM_SAFETY          <- 1.5      # multiply observed peak mem for --mem headroom
MEM_FLOOR_GB        <- 2L       # never request less than this

option_list <- list(
  make_option("--n-units",      type = "integer",   dest = "n_units_probe", default = 10L,
              help = "How many units to time locally [default %default]"),
  make_option("--target-hours", type = "double",    dest = "target_hours", default = 2,
              help = "Target wall time per task in hours, within [1,3] [default %default]"),
  make_option("--study-dir",    type = "character", dest = "study_dir", default = "..",
              help = "Path to study dir (contains config/, R/) [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

target_hours <- min(max(opt$target_hours, WALL_MIN_HOURS), WALL_MAX_HOURS)
if (target_hours != opt$target_hours) {
  cat(sprintf("NOTE: target-hours clamped to [%g, %g] -> %g\n",
              WALL_MIN_HOURS, WALL_MAX_HOURS, target_hours))
}

# --- Load study code ---------------------------------------------------------
source(file.path(opt$study_dir, "config", "grid.R"))
source(file.path(opt$study_dir, "R", "dgp.R"))
source(file.path(opt$study_dir, "R", "estimators.R"))
source(file.path(opt$study_dir, "R", "run_one.R"))
library(doubletree)
library(optimaltrees)

ut <- unit_table()
total_units <- nrow(ut)
n_probe <- min(opt$n_units_probe, total_units)
cat(sprintf("Profiling %d of %d total work units...\n", n_probe, total_units))

# --- Time each probe unit; track peak memory ---------------------------------
# Capture each result and check it is non-degenerate: a unit whose estimators all
# errored (e.g. the input .rds files are missing on this machine) returns
# instantly with all-NA estimates. Sizing off such runs yields a meaningless
# (tiny) walltime, so we refuse to proceed rather than emit a bogus sizing.env
# (Constitution Section 9: fail loudly, no silent fallback).
# NOTE: this treats an all-NA `estimate` column as a failed unit; run_one.R must
# surface an `estimate` column for the check to apply (see run_one.R contract).
invisible(gc(reset = TRUE))
per_unit_secs <- numeric(n_probe)
n_degenerate  <- 0L
for (i in seq_len(n_probe)) {
  res <- NULL
  t <- system.time(res <- run_one(ut[i, , drop = FALSE]))
  per_unit_secs[i] <- unname(t["elapsed"])
  if (is.null(res) || !("estimate" %in% names(res)) || all(is.na(res$estimate))) {
    n_degenerate <- n_degenerate + 1L
  }
}
if (n_degenerate == n_probe) {
  stop(sprintf(paste0(
    "All %d probe units produced only NA estimates -- every replication failed ",
    "(see the ERROR lines above; commonly the input .rds files are missing on ",
    "this machine).\nRefusing to write a sizing based on failed runs. Fix the ",
    "inputs, then re-run profile_timing.R."), n_probe))
}
if (n_degenerate > 0) {
  cat(sprintf("WARNING: %d/%d probe units returned all-NA estimates; timing may be unrepresentative.\n",
              n_degenerate, n_probe))
}
# gc() reports peak since reset. The LAST column is always "max used (Mb)"
# (column count varies depending on whether a memory limit is set, so index
# from the end rather than a fixed position). Sum Ncells + Vcells rows.
gcinfo <- gc()
peak_mb <- sum(gcinfo[, ncol(gcinfo)])   # last column = "max used (Mb)"
peak_gb <- peak_mb / 1024

med_secs  <- median(per_unit_secs)
mean_secs <- mean(per_unit_secs)
cat(sprintf("Per-unit wall time: median %.3f s, mean %.3f s (n=%d)\n",
            med_secs, mean_secs, n_probe))
cat(sprintf("Peak memory during profiling: %.2f GB\n", peak_gb))

if (med_secs <= 0) stop("Non-positive per-unit time; cannot size jobs.")

# --- Sizing math -------------------------------------------------------------
target_secs <- target_hours * 3600
max_secs    <- WALL_MAX_HOURS * 3600

# Start from the wall-time target, but never pack more than all the work into
# one job (reps_per_job is capped at total_units).
reps_per_job <- max(1L, min(as.integer(floor(target_secs / med_secs)), total_units))
total_tasks  <- as.integer(ceiling(total_units / reps_per_job))

# If we'd exceed the concurrent-jobs cap, pack more reps per job to cut task
# count -- but only up to the 3 hr ceiling.
if (total_tasks > MAX_CONCURRENT_JOBS) {
  reps_needed <- as.integer(ceiling(total_units / MAX_CONCURRENT_JOBS))
  wall_at_needed <- reps_needed * med_secs
  if (wall_at_needed <= max_secs) {
    reps_per_job <- reps_needed
    total_tasks  <- as.integer(ceiling(total_units / reps_per_job))
    cat(sprintf("Packed to %d reps/job to keep total tasks <= %d (wall ~%.2f hr)\n",
                reps_per_job, MAX_CONCURRENT_JOBS, wall_at_needed / 3600))
  } else {
    # Cannot fit under the cap within 3 hr; cap wall at 3 hr and submit in waves.
    reps_per_job <- max(1L, floor(max_secs / med_secs))
    total_tasks  <- as.integer(ceiling(total_units / reps_per_job))
    cat(sprintf("WARNING: cannot fit <= %d jobs within %g hr. Using %d reps/job (~%g hr);\n",
                MAX_CONCURRENT_JOBS, WALL_MAX_HOURS, reps_per_job, WALL_MAX_HOURS))
    cat(sprintf("         total tasks = %d > %d -> submit.sh will submit in WAVES.\n",
                total_tasks, MAX_CONCURRENT_JOBS))
  }
}

# Partition tasks into arrays of <= MAX_ARRAY_SIZE.
n_arrays <- as.integer(ceiling(total_tasks / MAX_ARRAY_SIZE))

# Concurrency cap (%N on --array) so a single user does not swamp the queue and
# so total queued stays within MAX_CONCURRENT_JOBS.
concurrency_cap <- min(MAX_ARRAY_SIZE, MAX_CONCURRENT_JOBS)

est_task_secs <- reps_per_job * med_secs
walltime_secs <- ceiling(est_task_secs * TIME_SAFETY)
# Floor at 10 min so trivially fast tasks still get a sane allocation.
walltime_secs <- max(walltime_secs, 600)

# Format as SLURM D-HH:MM:SS.
fmt_slurm_time <- function(secs) {
  secs <- as.integer(ceiling(secs))
  d <- secs %/% 86400; secs <- secs %% 86400
  h <- secs %/% 3600;  secs <- secs %% 3600
  m <- secs %/% 60;    s <- secs %% 60
  sprintf("%d-%02d:%02d:%02d", d, h, m, s)
}
walltime <- fmt_slurm_time(walltime_secs)

mem_gb <- max(MEM_FLOOR_GB, as.integer(ceiling(peak_gb * MEM_SAFETY)))

cat("\n===== SIZING =====\n")
cat(sprintf("total_units      : %d\n", total_units))
cat(sprintf("reps_per_job     : %d\n", reps_per_job))
cat(sprintf("total_tasks      : %d\n", total_tasks))
cat(sprintf("tasks_per_array  : <= %d (%d array(s))\n", MAX_ARRAY_SIZE, n_arrays))
cat(sprintf("concurrency_cap  : %d\n", concurrency_cap))
cat(sprintf("est wall / task  : %.2f hr\n", est_task_secs / 3600))
cat(sprintf("--time           : %s\n", walltime))
cat(sprintf("--mem            : %dG\n", mem_gb))
cat("==================\n\n")

# --- Write artifacts ---------------------------------------------------------
config_dir <- file.path(opt$study_dir, "config")
dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)

# sizing.env -- sourced by bash (KEY=VALUE, no spaces, no quotes needed).
env_lines <- c(
  sprintf("TOTAL_UNITS=%d", total_units),
  sprintf("REPS_PER_JOB=%d", reps_per_job),
  sprintf("TOTAL_TASKS=%d", total_tasks),
  sprintf("MAX_ARRAY_SIZE=%d", MAX_ARRAY_SIZE),
  sprintf("MAX_CONCURRENT_JOBS=%d", MAX_CONCURRENT_JOBS),
  sprintf("N_ARRAYS=%d", n_arrays),
  sprintf("CONCURRENCY_CAP=%d", concurrency_cap),
  sprintf("WALLTIME=%s", walltime),
  sprintf("MEM_GB=%d", mem_gb)
)
writeLines(env_lines, file.path(config_dir, "sizing.env"))

# sizing.json -- human-readable record.
json <- sprintf(
'{
  "study": "%s",
  "profiled_units": %d,
  "median_secs_per_unit": %.6f,
  "mean_secs_per_unit": %.6f,
  "peak_gb": %.3f,
  "target_hours": %g,
  "total_units": %d,
  "reps_per_job": %d,
  "total_tasks": %d,
  "max_array_size": %d,
  "max_concurrent_jobs": %d,
  "n_arrays": %d,
  "concurrency_cap": %d,
  "walltime": "%s",
  "mem_gb": %d
}',
  STUDY_NAME, n_probe, med_secs, mean_secs, peak_gb, target_hours,
  total_units, reps_per_job, total_tasks, MAX_ARRAY_SIZE, MAX_CONCURRENT_JOBS,
  n_arrays, concurrency_cap, walltime, mem_gb)
writeLines(json, file.path(config_dir, "sizing.json"))

cat(sprintf("Wrote %s and %s\n",
            file.path(config_dir, "sizing.env"),
            file.path(config_dir, "sizing.json")))
