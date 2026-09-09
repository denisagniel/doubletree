#!/usr/bin/env Rscript
# =============================================================================
# purge_failed_tasks.R -- delete poisoned task/unit result files so a resubmit
#                          regenerates them
# =============================================================================
# A "poisoned" result is one where the replication failed: the row carries a
# non-NA `error_msg` (M4) or its `estimate` is all-NA. Because tasks are
# idempotent (run_replication.R / array.slurm skip a unit whose partial exists
# and skip a task whose task_*.rds exists), a poisoned file would otherwise be
# kept forever on `--resume`. This tool finds them and deletes BOTH the assembled
# task_*.rds AND the offending unit_*.rds partials, so re-running the array
# recomputes exactly the failed units.
#
# Reports first, deletes only on --apply (or when you confirm interactively).
#
# Usage:
#   Rscript slurm/purge_failed_tasks.R --scratch-dir DIR            # report only
#   Rscript slurm/purge_failed_tasks.R --scratch-dir DIR --apply    # delete
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--scratch-dir", type = "character", dest = "scratch_dir",
              help = "Run-specific scratch dir containing task_*.rds and partials/"),
  make_option("--apply",       action = "store_true", dest = "apply", default = FALSE,
              help = "Actually delete poisoned files (default: report only)")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$scratch_dir))
if (!dir.exists(opt$scratch_dir)) stop(sprintf("scratch dir not found: %s", opt$scratch_dir))

is_poisoned <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(TRUE)
  bad_err <- ("error_msg" %in% names(df)) && any(!is.na(df$error_msg))
  bad_est <- ("estimate"  %in% names(df)) && all(is.na(df$estimate))
  isTRUE(bad_err) || isTRUE(bad_est)
}

# --- Scan assembled task files -----------------------------------------------
task_files <- list.files(opt$scratch_dir, pattern = "^task_[0-9]+\\.rds$", full.names = TRUE)
poison_tasks <- character(0)
poison_units <- integer(0)
for (f in task_files) {
  df <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(df)) { poison_tasks <- c(poison_tasks, f); next }
  if ("error_msg" %in% names(df)) {
    bad <- df$unit[!is.na(df$error_msg)]
    if (length(bad)) poison_units <- c(poison_units, bad)
  }
  if ("estimate" %in% names(df) && "unit" %in% names(df)) {
    bad <- df$unit[is.na(df$estimate)]
    if (length(bad)) poison_units <- c(poison_units, bad)
  }
  if (is_poisoned(df)) poison_tasks <- c(poison_tasks, f)
}

# --- Scan per-unit partials (source of truth for resume) ---------------------
partials_dir <- file.path(opt$scratch_dir, "partials")
poison_partials <- character(0)
if (dir.exists(partials_dir)) {
  pfiles <- list.files(partials_dir, pattern = "^unit_[0-9]+\\.rds$", full.names = TRUE)
  for (f in pfiles) {
    df <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is_poisoned(df)) {
      poison_partials <- c(poison_partials, f)
      if (!is.null(df) && "unit" %in% names(df)) poison_units <- c(poison_units, df$unit)
    }
  }
}

poison_units <- sort(unique(poison_units))
# Any task file whose block includes a poisoned unit must also be dropped so it
# is re-assembled after the units are recomputed.
poison_tasks <- unique(poison_tasks)

cat(sprintf("Scanned %d task file(s) and %d partial(s).\n",
            length(task_files),
            if (dir.exists(partials_dir)) length(list.files(partials_dir, pattern = "^unit_")) else 0L))
cat(sprintf("Poisoned: %d unit(s), %d task file(s), %d partial(s).\n",
            length(poison_units), length(poison_tasks), length(poison_partials)))
if (length(poison_units)) {
  cat("  poisoned units:", paste(head(poison_units, 50), collapse = ", "),
      if (length(poison_units) > 50) "..." else "", "\n")
}

if (length(poison_tasks) == 0 && length(poison_partials) == 0) {
  cat("Nothing to purge.\n")
  quit(save = "no", status = 0)
}

if (!isTRUE(opt$apply)) {
  cat("\nReport only. Re-run with --apply to delete the above, then resubmit to regenerate.\n")
  quit(save = "no", status = 0)
}

removed <- 0L
for (f in c(poison_tasks, poison_partials)) {
  if (file.exists(f) && unlink(f) == 0L) removed <- removed + 1L
}
cat(sprintf("Deleted %d file(s). Resubmit (bash slurm/submit.sh or the array) to regenerate.\n", removed))
