#!/usr/bin/env Rscript
# =============================================================================
# resize_sizing.R -- repackage EXISTING config/sizing_<method>.env for a new
# TOTAL-task target WITHOUT reprofiling.
# =============================================================================
# WHY THIS EXISTS: profiling (profile_per_method.R / profile.slurm) is a long,
# memory-heavy O2 job -- it runs the real estimators at the n=2000/continuous
# Rashomon blow-up to measure per-unit runtime + peak memory. But repackaging the
# SAME workload from many short tasks into ~1000 long tasks does NOT change either
# measurement: the estimator code and DGP are identical, and a task runs its units
# SEQUENTIALLY (run_replication.R), so peak RSS is one unit's peak regardless of how
# many units share a task. Only REPS_PER_JOB (how we pack units) changes -- pure
# arithmetic. So if valid sizing_<method>.env already exist, we can rewrite them for
# a --target-tasks target in seconds on the login node, skipping profiling entirely.
#
# HOW: the profiler stored WALLTIME = reps_per_job * med * TIME_SAFETY, so per-unit
# cost is recovered as med = WALLTIME_secs / (REPS_PER_JOB * TIME_SAFETY). Then the
# cost-proportional allocation (identical to profile_per_method.R --target-tasks)
# gives every method the same expected per-task walltime ~= (total CPU-time)/target,
# which auto-picks the partition. MEM_GB / UNIT_OFFSET / N_UNITS_METHOD are preserved.
#
# WHEN YOU MUST REPROFILE INSTEAD (this script refuses, loudly):
#   - a sizing_<method>.env is missing (nothing to resize from), or
#   - a method's WALLTIME is pinned at the old wall cap (--old-wall-cap-hours): then
#     med is only a LOWER BOUND (the true per-unit cost was clamped away) and resizing
#     would under-set --time. Reprofile that run with profile.slurm.
#   - (your judgement) the estimator code changed since those .env were written --
#     then the stored timing is stale; reprofile.
#
# RUN ON O2 from the study dir (login node is fine -- this is just arithmetic):
#   Rscript slurm/resize_sizing.R --target-tasks 1000
# then:  bash slurm/submit_per_method.sh
# Originals are backed up to config/sizing_<method>.env.bak before overwriting.
# =============================================================================

suppressPackageStartupMessages(library(optparse))

# Must match profile_per_method.R.
METHODS          <- c("full", "crossfit", "doubletree", "dt_averaged",
                      "msplit", "msplit_averaged", "single_tree")
TIME_SAFETY_DEF  <- 3.5
CONCURRENCY_CAP_DEF <- 200L
MAX_CONCURRENT_JOBS <- 10000L
PARTITION_CAPS_H <- c(short = 12, medium = 5 * 24, long = 30 * 24)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--target-tasks", type = "integer", dest = "target_tasks", default = 1000L,
              help = "Total array tasks across ALL methods [default %default]"),
  make_option("--study-dir", type = "character", dest = "study_dir", default = NA_character_,
              help = "Study dir [default: parent of this script's slurm/ dir]"),
  make_option("--time-safety", type = "double", dest = "time_safety", default = TIME_SAFETY_DEF,
              help = "TIME_SAFETY the .env were written with (WALLTIME = reps*med*this) [default %default]"),
  make_option("--old-wall-cap-hours", type = "double", dest = "old_wall_cap_h", default = 3,
              help = "Wall cap of the run that WROTE these .env; a WALLTIME at this cap means med was clamped [default %default]"),
  make_option("--new-wall-cap-hours", type = "double", dest = "new_wall_cap_h", default = 30 * 24,
              help = "Wall cap for the NEW sizing (hours) [default %default = 30d]")
)))

# Resolve study dir from this script's own path (parent of slurm/) if not given.
if (is.na(opt$study_dir)) {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE))
  opt$study_dir <- if (length(file_arg) == 1 && nzchar(file_arg))
    dirname(dirname(normalizePath(file_arg))) else "."
}
config_dir <- file.path(opt$study_dir, "config")
if (!dir.exists(config_dir))
  stop(sprintf("config dir not found under study dir '%s'.", opt$study_dir), call. = FALSE)
if (opt$target_tasks < length(METHODS))
  stop(sprintf("--target-tasks must be >= %d (one task minimum per method).", length(METHODS)),
       call. = FALSE)

# --- Helpers (byte-identical intent to profile_per_method.R) -----------------
fmt_hms <- function(secs) {
  secs <- max(60, ceiling(secs))
  sprintf("%02d:%02d:%02d", secs %/% 3600, (secs %% 3600) %/% 60, secs %% 60)
}
parse_hms <- function(s) {
  p <- as.numeric(strsplit(trimws(s), ":", fixed = TRUE)[[1]])
  if (length(p) != 3 || any(is.na(p)))
    stop(sprintf("cannot parse WALLTIME '%s' (expected HH:MM:SS).", s), call. = FALSE)
  p[1] * 3600 + p[2] * 60 + p[3]
}
pick_partition <- function(secs) {
  hrs <- secs / 3600; caps <- sort(PARTITION_CAPS_H)
  fit <- names(caps)[hrs <= caps]
  if (length(fit)) fit[1] else names(caps)[length(caps)]
}
read_env <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[grepl("=", ln, fixed = TRUE) & !grepl("^\\s*#", ln)]
  kv <- strsplit(ln, "=", fixed = TRUE)
  setNames(vapply(kv, function(x) paste(x[-1], collapse = "="), character(1)),
           vapply(kv, `[`, character(1), 1L))
}

# --- PASS 1: read every method's env, recover per-unit cost ------------------
envs <- list(); prof <- list()
old_cap_secs <- opt$old_wall_cap_h * 3600
for (m in METHODS) {
  path <- file.path(config_dir, sprintf("sizing_%s.env", m))
  if (!file.exists(path))
    stop(sprintf(paste0("Missing %s. Cannot resize what does not exist -- reprofile ",
                        "(sbatch --export=ALL,TARGET_TASKS=%d slurm/profile.slurm)."),
                 path, opt$target_tasks), call. = FALSE)
  e <- read_env(path)
  for (k in c("N_UNITS_METHOD", "REPS_PER_JOB", "WALLTIME", "MEM_GB", "UNIT_OFFSET"))
    if (is.na(e[k]) || !nzchar(e[k]))
      stop(sprintf("%s: missing required key %s.", path, k), call. = FALSE)

  reps_old  <- as.integer(e[["REPS_PER_JOB"]])
  wall_secs <- parse_hms(e[["WALLTIME"]])
  n_units   <- as.integer(e[["N_UNITS_METHOD"]])
  # Recover per-unit median: WALLTIME = reps * med * TIME_SAFETY  =>  med = ...
  med <- wall_secs / (reps_old * opt$time_safety)
  clamped <- wall_secs >= old_cap_secs - 1  # WALLTIME sat at the old cap => med is a floor
  if (clamped)
    stop(sprintf(paste0("%s: WALLTIME %s is at the old wall cap (%.1f h), so the recovered ",
                        "per-unit cost is only a LOWER BOUND -- resizing would under-set --time. ",
                        "Reprofile instead: sbatch --export=ALL,TARGET_TASKS=%d slurm/profile.slurm"),
                 path, e[["WALLTIME"]], opt$old_wall_cap_h, opt$target_tasks), call. = FALSE)
  if (!is.finite(med) || med <= 0)
    stop(sprintf("%s: recovered non-positive per-unit cost (%.3f).", path, med), call. = FALSE)

  envs[[m]] <- e
  prof[[m]] <- list(med = med, n_units = n_units)
  cat(sprintf("[%-15s] recovered med %.2fs/unit (reps_old %d, walltime %s)\n",
              m, med, reps_old, e[["WALLTIME"]]))
}

# --- Cost-proportional task allocation (uniform expected walltime) -----------
cost   <- vapply(prof, function(p) p$n_units * p$med, numeric(1))
budget <- pmax(1L, as.integer(round(opt$target_tasks * cost / sum(cost))))
names(budget) <- names(prof)
new_cap_secs <- opt$new_wall_cap_h * 3600

# --- PASS 2: rewrite each env (backing up the original) ----------------------
summary_rows <- list()
for (m in METHODS) {
  e <- envs[[m]]; p <- prof[[m]]
  reps_new  <- max(1L, as.integer(ceiling(p$n_units / budget[[m]])))
  tasks_new <- as.integer(ceiling(p$n_units / reps_new))
  wall_raw  <- reps_new * p$med * opt$time_safety
  walltime  <- fmt_hms(min(new_cap_secs, wall_raw))
  if (wall_raw > new_cap_secs)
    cat(sprintf(paste0("[%s] WARNING: sized walltime %.1f h exceeds new cap %.1f h; clamped. ",
                       "Raise --target-tasks (more, shorter tasks).\n"),
                m, wall_raw / 3600, opt$new_wall_cap_h))
  partition <- pick_partition(min(new_cap_secs, wall_raw))
  cap       <- min(CONCURRENCY_CAP_DEF, tasks_new)

  # Preserve everything else; update the packing/time/partition fields.
  e[["TOTAL_TASKS"]]    <- as.character(tasks_new)
  e[["REPS_PER_JOB"]]   <- as.character(reps_new)
  e[["CONCURRENCY_CAP"]] <- as.character(cap)
  e[["WALLTIME"]]       <- walltime
  e[["PARTITION"]]      <- partition
  if (is.na(e["MAX_ARRAY_SIZE"]))      e[["MAX_ARRAY_SIZE"]] <- "1000"
  if (is.na(e["MAX_CONCURRENT_JOBS"])) e[["MAX_CONCURRENT_JOBS"]] <- as.character(MAX_CONCURRENT_JOBS)

  # Emit in the same key order profile_per_method.R uses (PARTITION last).
  ord <- c("METHOD", "UNIT_OFFSET", "N_UNITS_METHOD", "TOTAL_TASKS", "REPS_PER_JOB",
           "MAX_ARRAY_SIZE", "MAX_CONCURRENT_JOBS", "CONCURRENCY_CAP", "WALLTIME",
           "MEM_GB", "PARTITION")
  ord <- c(ord[ord %in% names(e)], setdiff(names(e), ord))
  lines <- vapply(ord, function(k) sprintf("%s=%s", k, e[[k]]), character(1))

  path <- file.path(config_dir, sprintf("sizing_%s.env", m))
  file.copy(path, paste0(path, ".bak"), overwrite = TRUE)   # backup original
  writeLines(lines, path)

  cat(sprintf("[%-15s] reps/job %d  tasks %d  --time %s  --mem %sG  -p %s\n",
              m, reps_new, tasks_new, walltime, e[["MEM_GB"]], partition))
  summary_rows[[m]] <- tasks_new
}

grand <- sum(unlist(summary_rows))
cat(sprintf("\nRewrote config/sizing_<method>.env for %d methods (originals -> *.env.bak).\n",
            length(METHODS)))
cat(sprintf("TOTAL array tasks across all methods: %d  (target %d)\n", grand, opt$target_tasks))
if (grand > MAX_CONCURRENT_JOBS)
  cat(sprintf("WARNING: %d tasks exceeds MaxSubmitJobs=%d.\n", grand, MAX_CONCURRENT_JOBS))
cat("Next: bash slurm/submit_per_method.sh   (submits one array per method)\n")
