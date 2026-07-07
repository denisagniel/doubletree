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
# MEMORY GUARD (why this script cannot crash its host):
# The worst cell (n=2000, continuous) is exactly the documented Rashomon blow-up
# (a propensity Rashomon set can explode to ~123k trees at n>=1000; M-split does
# M*K=50 tree fits/unit). Left unguarded this OOMs the whole process -- on a laptop
# it restarted the machine; on O2 it would let SLURM OS-kill the entire salloc and
# lose ALL seven methods' timing. So every probe unit is run in a KILLABLE callr
# subprocess whose resident memory (RSS, incl. C++/Rcpp -- where GOSDT actually
# lives, and which gc() cannot see) is polled via `ps` and killed if it exceeds
# --mem-cap-gb. A runaway unit dies ALONE, recorded as degenerate, and its death
# is itself a sizing signal (that cell needs a big --mem array, or exclusion).
# Set --mem-cap-gb BELOW your salloc --mem so the cap fires before SLURM does.
#
# RUN THIS ON O2 (matches run hardware) from the study dir, inside a generously
# sized interactive allocation, e.g.:
#   salloc -p interactive -c 4 --mem=48G -t 0-04:00
#   Rscript slurm/profile_per_method.R --target-hours 2 [--n-units 5] [--mem-cap-gb 24]
# Output: config/sizing_<method>.env for each method (sourced by submit_per_method.sh).
# =============================================================================

suppressPackageStartupMessages(library(optparse))

# Memory guard depends on callr (killable subprocess) + ps (portable RSS polling,
# macOS + Linux). Fail loudly if absent rather than silently profiling unguarded.
for (pkg in c("callr", "ps")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(paste0(
      "Package '%s' is required for the profiling memory guard but is not ",
      "installed. Install it (install.packages('%s')) before profiling."),
      pkg, pkg), call. = FALSE)
  }
}

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
  make_option("--mem-cap-gb",   type = "double",    dest = "mem_cap_gb", default = 24,
              help = paste0("Per-unit resident-memory ceiling in GB; a probe unit ",
                            "exceeding this is killed and recorded as degenerate. ",
                            "Set BELOW your salloc --mem. [default %default]")),
  make_option("--poll-secs",    type = "double",    dest = "poll_secs", default = 0.5,
              help = "How often to poll child RSS, seconds [default %default]"),
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

# --- Memory-capped unit runner ------------------------------------------------
# Run one probe unit in a killable callr subprocess, polling its resident memory
# (RSS, via ps -- captures C++/Rcpp, unlike gc()) and killing it if it exceeds
# `cap_gb`. Returns elapsed seconds, whether the result is degenerate (all-NA or
# errored/killed), the peak RSS observed (GB), and the outcome reason. Never
# throws for a unit-level failure: a blow-up is DATA, not a script error.
run_unit_capped <- function(unit_row, study_dir, cap_gb, poll_secs, libpath) {
  cap_bytes <- cap_gb * 1024^3
  study_dir_abs <- normalizePath(study_dir)

  proc <- callr::r_bg(
    function(unit_row, study_dir) {
      source(file.path(study_dir, "config", "grid.R"))
      source(file.path(study_dir, "R", "dgp.R"))
      source(file.path(study_dir, "R", "estimators.R"))
      source(file.path(study_dir, "R", "run_one.R"))
      suppressPackageStartupMessages({
        library(doubletree); library(optimaltrees)
      })
      run_one(unit_row)
    },
    args = list(unit_row = unit_row, study_dir = study_dir_abs),
    libpath = libpath,
    supervise = TRUE   # child is killed if the parent (this profiler) dies
  )

  t0 <- proc.time()[["elapsed"]]
  peak_bytes <- 0
  killed <- FALSE
  ps_handle <- tryCatch(ps::ps_handle(proc$get_pid()), error = function(e) NULL)

  while (proc$is_alive()) {
    if (!is.null(ps_handle)) {
      rss <- tryCatch(ps::ps_memory_info(ps_handle)[["rss"]], error = function(e) NA_real_)
      if (is.finite(rss)) {
        peak_bytes <- max(peak_bytes, rss)
        if (rss > cap_bytes) {
          proc$kill_tree()   # kill child AND any GOSDT worker descendants
          killed <- TRUE
          break
        }
      }
    }
    Sys.sleep(poll_secs)
  }
  proc$wait(timeout = 5000)
  elapsed <- proc.time()[["elapsed"]] - t0

  # Resolve outcome. Killed => degenerate by construction. Otherwise inspect the
  # result: an errored subprocess or an all-NA `estimate` is a failed unit.
  degenerate <- TRUE
  reason <- if (killed) "killed (mem cap)" else "unknown"
  if (!killed) {
    res <- tryCatch(proc$get_result(), error = function(e) e)
    if (inherits(res, "error") || inherits(res, "condition")) {
      reason <- sprintf("error: %s", conditionMessage(res))
    } else if (is.null(res) || !("estimate" %in% names(res)) || all(is.na(res$estimate))) {
      reason <- "degenerate (all-NA estimate)"
    } else {
      degenerate <- FALSE
      reason <- "ok"
    }
  }

  list(elapsed = elapsed, degenerate = degenerate,
       peak_gb = peak_bytes / 1024^3, killed = killed, reason = reason)
}

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

cat(sprintf(paste0("Per-method profiling at worst cell (n=%d, dgp=%s), %d probe units each.\n",
                   "Per-unit memory cap: %.1f GB (units exceeding it are killed & recorded).\n\n"),
            worst_n, worst_dgp, opt$n_units_probe, opt$mem_cap_gb))

summary_rows <- list()
for (meth in methods) {
  # Probe units for this method at the worst cell.
  cell <- ut[ut$method == meth & ut$n == worst_n & ut$dgp == worst_dgp, , drop = FALSE]
  probe <- head(cell, opt$n_units_probe)
  if (nrow(probe) == 0) { cat(sprintf("[%s] no worst-cell units; skipping.\n", meth)); next }

  secs <- numeric(nrow(probe)); unit_peak_gb <- numeric(nrow(probe))
  n_deg <- 0L; n_killed <- 0L
  for (i in seq_len(nrow(probe))) {
    r <- run_unit_capped(probe[i, , drop = FALSE], opt$study_dir,
                         opt$mem_cap_gb, opt$poll_secs, .libPaths())
    secs[i] <- r$elapsed
    unit_peak_gb[i] <- r$peak_gb
    if (r$degenerate) n_deg <- n_deg + 1L
    if (r$killed)     n_killed <- n_killed + 1L
  }
  peak_gb <- max(unit_peak_gb, 0)   # measured RSS incl. C++/Rcpp (not gc-based)
  med <- median(secs)
  if (n_deg == nrow(probe)) {
    cat(sprintf("[%s] WARNING: all %d probe units degenerate (all-NA/killed). Sizing off wall time only.\n",
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
  cap      <- min(CONCURRENCY_CAP_DEF, total_tasks)

  # Memory sizing. If any unit was killed, the true peak EXCEEDS the cap -- the
  # measured peak_gb is only a lower bound, so size --mem above the cap and warn
  # loudly (Constitution S9: no silent under-sizing). Otherwise use measured peak.
  if (n_killed > 0) {
    mem_gb <- max(4L, as.integer(ceiling(opt$mem_cap_gb * MEM_SAFETY)))
    cat(sprintf(paste0("[%s] WARNING: %d/%d probe units hit the %.1f GB cap and were KILLED. ",
                       "True memory need EXCEEDS the cap; --mem set to %dG is a FLOOR, not a ",
                       "measurement. Re-profile with a higher --mem-cap-gb (in a bigger salloc) ",
                       "or consider excluding this cell.\n"),
                meth, n_killed, nrow(probe), opt$mem_cap_gb, mem_gb))
  } else {
    mem_gb <- max(4L, as.integer(ceiling(peak_gb * MEM_SAFETY)))
  }

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

  peak_note <- if (n_killed > 0) sprintf(">%.2fGB", opt$mem_cap_gb) else sprintf("%.2fGB", peak_gb)
  cat(sprintf("[%-15s] med %.2fs/unit  peak %s  killed %d/%d | reps/job %d  tasks %d  --time %s  --mem %dG\n",
              meth, med, peak_note, n_killed, nrow(probe), reps_per_job, total_tasks, walltime, mem_gb))
  summary_rows[[meth]] <- data.frame(method = meth, med_secs = med, peak_gb = peak_gb,
                                     n_killed = n_killed, reps_per_job = reps_per_job,
                                     total_tasks = total_tasks, walltime = walltime, mem_gb = mem_gb)
}

cat("\nWrote config/sizing_<method>.env for each method.\n")
cat("Next: bash slurm/submit_per_method.sh   (submits one array per method)\n")
