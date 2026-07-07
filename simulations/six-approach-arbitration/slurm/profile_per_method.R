#!/usr/bin/env Rscript
# =============================================================================
# profile_per_method.R -- per-method WORST-CELL profiling for per-method arrays
# =============================================================================
# WHY THIS EXISTS (not the stock profile_timing.R):
# The grid enumerates units method-major, and each method occupies a contiguous
# 12,000-unit block (full 1-12000, crossfit 12001-24000, ..., single_tree
# 72001-84000). The seven methods differ in cost by ~50x (M-split does M*K=50
# tree fits per unit; full does ~2). A single global sizing off the first N units
# would size everything to the cheapest method (`full`, at n=500) and the M-split
# tasks would time out. So we size EACH METHOD SEPARATELY and submit one array
# per method (see submit_per_method.sh), each with method-appropriate --time/--mem.
#
# To be safe we profile each method at its SLOWEST cell (largest n, hardest DGP:
# n=2000, dgp="continuous"), so the per-job wall time is sized for the worst unit
# that method will encounter, not its average. This over-sizes cheap cells within
# a method slightly (fine -- reliability > packing efficiency), but never times out.
#
# RUN THIS ON O2 (matches run hardware) from the study dir:
#   Rscript slurm/profile_per_method.R --target-hours 2 [--n-units 5]
# Output: config/sizing_<method>.env for each method (sourced by submit_per_method.sh).
# =============================================================================

suppressPackageStartupMessages(library(optparse))

# --- O2 scheduler limits (mirror profile_timing.R) ---------------------------
MAX_ARRAY_SIZE      <- 1000L
MAX_CONCURRENT_JOBS <- 10000L
WALL_MAX_HOURS      <- 3
CONCURRENCY_CAP_DEF <- 200L
TIME_SAFETY         <- 1.5
MEM_SAFETY          <- 1.5

opt <- parse_args(OptionParser(option_list = list(
  make_option("--n-units",      type = "integer",   dest = "n_units_probe", default = 5L,
              help = "Probe units per method at the worst cell [default %default]"),
  make_option("--target-hours", type = "double",    dest = "target_hours", default = 2,
              help = "Target wall time per task, hours [default %default]"),
  make_option("--study-dir",    type = "character", dest = "study_dir", default = NA_character_,
              help = "Study directory [default: auto-detected as the parent of this script's slurm/ dir]")
)))

# Resolve the study dir robustly: explicit --study-dir wins; otherwise derive it
# from this script's own path (parent of slurm/), so it works from any cwd.
if (is.na(opt$study_dir)) {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE))
  if (length(file_arg) == 1 && nzchar(file_arg)) {
    opt$study_dir <- dirname(dirname(normalizePath(file_arg)))
  } else {
    opt$study_dir <- "."   # fallback (e.g. sourced interactively)
  }
}
if (!file.exists(file.path(opt$study_dir, "config", "grid.R"))) {
  stop(sprintf("Cannot find config/grid.R under study dir '%s'. Pass --study-dir explicitly.",
               opt$study_dir), call. = FALSE)
}

source(file.path(opt$study_dir, "config", "grid.R"))
source(file.path(opt$study_dir, "R", "dgp.R"))
source(file.path(opt$study_dir, "R", "estimators.R"))
source(file.path(opt$study_dir, "R", "run_one.R"))
suppressPackageStartupMessages({ library(doubletree); library(optimaltrees) })

ut <- unit_table()
methods <- unique(GRID$method)
# Worst cell = largest n, hardest DGP. Fall back to the method's max-n cell if
# "continuous" is ever removed from the grid.
worst_n   <- max(GRID$n)
worst_dgp <- if ("continuous" %in% GRID$dgp) "continuous" else tail(sort(unique(GRID$dgp)), 1)

target_secs <- opt$target_hours * 3600
max_secs    <- WALL_MAX_HOURS * 3600

fmt_hms <- function(secs) {
  secs <- max(60, ceiling(secs))
  h <- secs %/% 3600; m <- (secs %% 3600) %/% 60; s <- secs %% 60
  sprintf("%02d:%02d:%02d", h, m, s)
}

cat(sprintf("Per-method profiling at worst cell (n=%d, dgp=%s), %d probe units each.\n\n",
            worst_n, worst_dgp, opt$n_units_probe))

summary_rows <- list()
for (meth in methods) {
  # Probe units for this method at the worst cell.
  cell <- ut[ut$method == meth & ut$n == worst_n & ut$dgp == worst_dgp, , drop = FALSE]
  probe <- head(cell, opt$n_units_probe)
  if (nrow(probe) == 0) { cat(sprintf("[%s] no worst-cell units; skipping.\n", meth)); next }

  invisible(gc(reset = TRUE))
  secs <- numeric(nrow(probe)); n_deg <- 0L
  for (i in seq_len(nrow(probe))) {
    res <- NULL
    t <- system.time(res <- tryCatch(run_one(probe[i, , drop = FALSE]), error = function(e) NULL))
    secs[i] <- unname(t["elapsed"])
    if (is.null(res) || !("estimate" %in% names(res)) || all(is.na(res$estimate))) n_deg <- n_deg + 1L
  }
  gcinfo <- gc(); peak_gb <- sum(gcinfo[, ncol(gcinfo)]) / 1024
  med <- median(secs)
  if (n_deg == nrow(probe)) {
    cat(sprintf("[%s] WARNING: all %d probe units degenerate (all-NA). Sizing off wall time only.\n",
                meth, nrow(probe)))
  }
  if (med <= 0) med <- 1  # guard; a method that returns instantly still needs a floor

  # Units in this method's block = nrow(GRID cells for this method) * TOTAL_REPS.
  n_units_method <- sum(ut$method == meth)

  reps_per_job <- max(1L, min(as.integer(floor(target_secs / med)), n_units_method))
  total_tasks  <- as.integer(ceiling(n_units_method / reps_per_job))
  # Respect the concurrent-jobs cap by packing more reps if needed (<= 3 hr).
  if (total_tasks > MAX_CONCURRENT_JOBS) {
    reps_per_job <- as.integer(ceiling(n_units_method / MAX_CONCURRENT_JOBS))
    if (reps_per_job * med > max_secs) reps_per_job <- max(1L, as.integer(floor(max_secs / med)))
    total_tasks <- as.integer(ceiling(n_units_method / reps_per_job))
  }
  walltime <- fmt_hms(min(max_secs, reps_per_job * med * TIME_SAFETY))
  mem_gb   <- max(4L, as.integer(ceiling(peak_gb * MEM_SAFETY)))
  cap      <- min(CONCURRENCY_CAP_DEF, total_tasks)

  # Global unit offset for this method's block (units are 1-based, contiguous).
  unit_offset <- min(ut$unit[ut$method == meth]) - 1L

  env_path <- file.path(opt$study_dir, "config", sprintf("sizing_%s.env", meth))
  writeLines(c(
    sprintf("METHOD=%s", meth),
    sprintf("UNIT_OFFSET=%d", unit_offset),
    sprintf("N_UNITS_METHOD=%d", n_units_method),
    sprintf("TOTAL_TASKS=%d", total_tasks),
    sprintf("REPS_PER_JOB=%d", reps_per_job),
    sprintf("MAX_ARRAY_SIZE=%d", MAX_ARRAY_SIZE),
    sprintf("MAX_CONCURRENT_JOBS=%d", MAX_CONCURRENT_JOBS),
    sprintf("CONCURRENCY_CAP=%d", cap),
    sprintf("WALLTIME=%s", walltime),
    sprintf("MEM_GB=%d", mem_gb)
  ), env_path)

  cat(sprintf("[%-15s] med %.2fs/unit  peak %.2fGB | reps/job %d  tasks %d  --time %s  --mem %dG\n",
              meth, med, peak_gb, reps_per_job, total_tasks, walltime, mem_gb))
  summary_rows[[meth]] <- data.frame(method = meth, med_secs = med, reps_per_job = reps_per_job,
                                     total_tasks = total_tasks, walltime = walltime, mem_gb = mem_gb)
}

cat("\nWrote config/sizing_<method>.env for each method.\n")
cat("Next: bash slurm/submit_per_method.sh   (submits one array per method)\n")
