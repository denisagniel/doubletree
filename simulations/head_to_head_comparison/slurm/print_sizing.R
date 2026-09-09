#!/usr/bin/env Rscript
# ============================================================
# slurm/print_sizing.R
# Study: head_to_head_comparison  (doubletree, S5)
#
# Run:  Rscript simulations/head_to_head_comparison/slurm/print_sizing.R [--format table|env|units]
#
# Prints the job-array sizing that slurm/units.R derives from the pilot's measured
# per-replication costs. Three formats:
#
#   table  (default) human-readable: tasks, walltime, memory, projected hours per
#          regime, plus the totals README_O2.md quotes.
#   env    one shell-parseable line per regime, consumed by launch_all_simulations.sh
#          and launch_subset.sh:  REGIME N_TASKS WALLTIME MEM_GB PARTITION
#   units  the full per-task unit table (every array index and the block it owns).
#
# WHY THE LAUNCHERS CALL THIS instead of hard-coding numbers: the array size and
# the worker's idea of which block a task owns are two views of the SAME table. If
# they are computed separately they can drift, and a drifted map does not fail --
# it silently recomputes and mislabels somebody else's block (see the 2026-07-10
# contamination bug recorded in six-approach-arbitration/slurm/array.slurm).
#
# Reads no results and writes no files.
# ============================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--format", type = "character", dest = "format", default = "table",
              help = "One of table, env, units [default %default]"),
  make_option("--regimes", type = "character", dest = "regimes", default = "",
              help = "Comma-separated regime ids [default: all]"),
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
  Sys.setenv(H2H_USE_INSTALLED = "1")
  suppressMessages(
    source(file.path("simulations", "head_to_head_comparison", "code", "common.R"))
  )
  source(file.path("simulations", "head_to_head_comparison", "slurm", "units.R"))

  regimes <- if (nzchar(opt$regimes)) {
    trimws(strsplit(opt$regimes, ",", fixed = TRUE)[[1L]])
  } else {
    REGIME_IDS
  }
  unknown <- setdiff(regimes, REGIME_IDS)
  if (length(unknown)) {
    stop("Unknown regime(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  target <- h2h_target_secs()
  sizing <- do.call(rbind, lapply(regimes, h2h_sizing, target_secs = target))

  if (identical(opt$format, "env")) {
    # Deliberately bare and positional: the launchers read it with `read`, so any
    # decoration here would become a parsing bug there.
    for (i in seq_len(nrow(sizing))) {
      cat(sprintf("%s %d %s %d %s\n", sizing$regime[i], sizing$n_tasks[i],
                  sizing$walltime[i], sizing$mem_gb[i], sizing$partition[i]))
    }
  } else if (identical(opt$format, "units")) {
    ut <- do.call(rbind, lapply(regimes, h2h_unit_table, target_secs = target))
    ut$est_secs <- round(ut$est_secs, 1)
    print(ut, row.names = FALSE)
  } else {
    cat(sprintf("target seconds per array task: %.0f  (override with H2H_TARGET_SECS)\n",
                target))
    cat(sprintf("walltime requested = %g x projected + %g s, rounded up to the minute\n\n",
                SAFETY_FACTOR, WALLTIME_FLOOR_SECS))
    show <- sizing
    show$max_task_secs <- round(show$max_task_secs, 1)
    show$total_hours <- round(show$total_hours, 2)
    print(show, row.names = FALSE)
    cat(sprintf("\ntotal array tasks: %d    total compute: %.1f h    longest task: %.1f min\n",
                sum(sizing$n_tasks), sum(sizing$total_hours),
                max(sizing$max_task_secs) / 60))
    cat(sprintf("wall clock if every task runs concurrently: ~%.0f min (the longest task)\n",
                max(sizing$max_task_secs) / 60))
  }
}
