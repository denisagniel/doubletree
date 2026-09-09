#!/usr/bin/env Rscript
# ============================================================
# slurm/print_sizing.R
# Study: honest_inference_sparsity_failure  (doubletree, S2)
#
# Run:  Rscript simulations/honest_inference_sparsity_failure/slurm/print_sizing.R \
#         [--format table|env|units] [--waves n2000,n8000]
#
# Prints the job-array sizing that slurm/units.R derives from the measured
# per-replication costs. Three formats:
#
#   table  (default) human-readable: tasks, walltime, memory, projected hours per
#          wave, plus the totals README_O2.md quotes.
#   env    one shell-parseable line per wave, consumed by launch_subset.sh:
#          WAVE N_TASKS WALLTIME MEM_GB PARTITION
#   units  the full per-task unit table (every array index and the block it owns).
#
# WHY THE LAUNCHERS CALL THIS instead of hard-coding numbers: the array size and the
# worker's idea of which block a task owns are two views of the SAME table. If they
# are computed separately they can drift, and a drifted map does not fail -- it
# silently recomputes and mislabels somebody else's block (see the 2026-07-10
# contamination bug recorded in six-approach-arbitration/slurm/array.slurm).
#
# Reads no results and writes no files.
# ============================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--format", type = "character", dest = "format", default = "table",
              help = "One of table, env, units [default %default]"),
  make_option("--waves", type = "character", dest = "waves", default = "",
              help = "Comma-separated wave ids [default: the missing waves, n2000,n8000]"),
  make_option("--all-waves", action = "store_true", dest = "all_waves", default = FALSE,
              help = "Include every wave the design grid admits, n500 included"),
  make_option("--pkg-root", type = "character", dest = "pkg_root", default = "",
              help = "doubletree package root [default: three levels above this script]")
)

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (!length(hit)) return(NA_character_)
  normalizePath(dirname(sub("^--file=", "", hit[[1L]])), mustWork = FALSE)
}

if (!interactive() && sys.nframe() == 0L) {

  opt <- parse_args(OptionParser(option_list = option_list))

  pkg_root <- if (nzchar(opt$pkg_root)) {
    opt$pkg_root
  } else {
    d <- script_dir()
    if (is.na(d)) getwd() else normalizePath(file.path(d, "..", "..", ".."), mustWork = FALSE)
  }
  if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    stop("--pkg-root has no DESCRIPTION: ", pkg_root, call. = FALSE)
  }

  setwd(pkg_root)
  # Sizing needs only the grid objects, but common.R loads packages on the way to
  # defining them; take the installed path so this works on the cluster unchanged.
  Sys.setenv(HIS_USE_INSTALLED = "1")
  suppressMessages(
    source(file.path("simulations", "honest_inference_sparsity_failure", "code",
                     "common.R"))
  )
  source(file.path("simulations", "honest_inference_sparsity_failure", "slurm",
                   "units.R"))

  waves <- if (nzchar(opt$waves)) {
    trimws(strsplit(opt$waves, ",", fixed = TRUE)[[1L]])
  } else if (opt$all_waves) {
    WAVE_IDS
  } else {
    WAVES_MISSING
  }
  unknown <- setdiff(waves, WAVE_IDS)
  if (length(unknown)) {
    stop("Unknown wave(s): ", paste(unknown, collapse = ", "),
         "; expected ", paste(WAVE_IDS, collapse = ", "), call. = FALSE)
  }

  target <- his_target_secs()
  sizing <- do.call(rbind, lapply(waves, his_sizing, target_secs = target))

  if (identical(opt$format, "env")) {
    # Deliberately bare and positional: the launcher reads it with `read`, so any
    # decoration here would become a parsing bug there.
    for (i in seq_len(nrow(sizing))) {
      cat(sprintf("%s %d %s %d %s\n", sizing$wave[i], sizing$n_tasks[i],
                  sizing$walltime[i], sizing$mem_gb[i], sizing$partition[i]))
    }
  } else if (identical(opt$format, "units")) {
    ut <- do.call(rbind, lapply(waves, his_unit_table, target_secs = target))
    ut$est_secs <- round(ut$est_secs, 1)
    print(ut, row.names = FALSE)
  } else {
    cat(sprintf("study : honest_inference_sparsity_failure (S2)\n"))
    cat(sprintf("waves : %s   (of %s; n500 is already on disk)\n",
                paste(waves, collapse = ", "), paste(WAVE_IDS, collapse = ", ")))
    cat(sprintf("dgps  : %s\n", paste(DGP_IDS, collapse = ", ")))
    cat(sprintf("reps  : R = %d at n < %d, R = %d at n = %d (HIS_REPS / HIS_REPS_MAX_N)\n",
                REPS_DEFAULT, max(N_GRID), REPS_AT_MAX_N, max(N_GRID)))
    cat(sprintf("cost  : measured s/rep -- %s\n",
                paste(sprintf("n=%s: %.2f", names(COST_ANCHOR_BY_N),
                              COST_ANCHOR_BY_N), collapse = "  ")))
    cat(sprintf("target seconds per array task: %.0f  (override with HIS_TARGET_SECS)\n",
                target))
    cat(sprintf("walltime requested = %g x projected + %g s, rounded up to the minute\n\n",
                SAFETY_FACTOR, WALLTIME_FLOOR_SECS))
    show <- sizing
    show$max_task_secs <- round(show$max_task_secs, 1)
    show$total_hours <- round(show$total_hours, 2)
    print(show, row.names = FALSE)
    cat(sprintf("\ncells: %d    total array tasks: %d    total compute: %.2f h    longest task: %.1f min\n",
                length(waves) * length(DGP_IDS), sum(sizing$n_tasks),
                sum(sizing$total_hours), max(sizing$max_task_secs) / 60))
    cat(sprintf("wall clock if every task runs concurrently: ~%.0f min (the longest task)\n",
                max(sizing$max_task_secs) / 60))
  }
}
