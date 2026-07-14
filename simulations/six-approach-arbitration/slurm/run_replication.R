#!/usr/bin/env Rscript
# =============================================================================
# run_replication.R -- run ONE array task's worth of work units
# =============================================================================
# Invoked by array.slurm once per SLURM array index. Each task runs a contiguous
# block of `reps_per_job` work units, then writes ONE result file to the run's
# scratch directory. Output is idempotent: if the task's file already exists it
# exits immediately (cheap resume after partial failure).
#
# Crash safety (M3): each unit's result is flushed to partials/unit_NNNNNNN.rds
# (keyed by ABSOLUTE unit id) the moment it finishes. A wall-time TIMEOUT (SLURM
# SIGTERM->SIGKILL, see array.slurm's trap) loses at most the one in-flight unit;
# resubmitting the SAME run-id skips every completed unit (via the partials) and
# finishes the remainder cheaply. The task_*.rds file is assembled from the block's
# partials at the end (that is what combine.R reads). Because partials are keyed by
# ABSOLUTE unit id, the per-method UNIT_OFFSET slicing below is preserved exactly.
#
# No silent all-NA (M4): if run_one() errors on a unit, the persisted row carries
# the error MESSAGE in a character `error_msg` column (NA on success) instead of a
# silent all-NA row that looks like a legitimate non-converged replication.
#
# Deliberately uses library() (module R has no devtools) and explicit dest= on
# every optparse option (guards against hyphenated-arg parsing bugs).
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--task-id",      type = "integer", dest = "task_id",
              help = "SLURM array task id (1-based)"),
  make_option("--reps-per-job", type = "integer", dest = "reps_per_job",
              help = "Number of work units this task should run"),
  make_option("--study-dir",    type = "character", dest = "study_dir",
              help = "Absolute path to the study directory (contains config/, R/)"),
  make_option("--scratch-dir",  type = "character", dest = "scratch_dir",
              help = "Run-specific scratch dir for per-task result files"),
  # Per-method sizing (submit_per_method.sh): task ids restart at 1 within each
  # method, so a task's units are offset into that method's contiguous block and
  # capped at the block end. Defaults (0, Inf) reproduce the original global
  # behaviour, so the stock submit.sh keeps working unchanged.
  make_option("--unit-offset",  type = "integer", dest = "unit_offset", default = 0L,
              help = "Add this to the local unit index (start of this method's block) [default %default]"),
  make_option("--max-unit",     type = "integer", dest = "max_unit", default = NA_integer_,
              help = "Last global unit this task may run (method block end) [default: last unit]")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$task_id), !is.null(opt$reps_per_job),
          !is.null(opt$study_dir), !is.null(opt$scratch_dir))

# --- Load study code (order matters) -----------------------------------------
source(file.path(opt$study_dir, "config", "grid.R"))
source(file.path(opt$study_dir, "R", "dgp.R"))
source(file.path(opt$study_dir, "R", "estimators.R"))
source(file.path(opt$study_dir, "R", "run_one.R"))
library(doubletree)  # e.g. library(mypackage); blank if no project package
library(optimaltrees)          # e.g. library(ranger); library(glmnet)

# --- Determine this task's block of units ------------------------------------
# Global-unit indexing. With --unit-offset O and --max-unit U, task t covers
# global units (O + (t-1)*reps + 1) .. min(O + t*reps, U). Defaults O=0, U=last
# reproduce the original single-block behaviour.
ut <- unit_table()
max_unit <- if (is.na(opt$max_unit)) nrow(ut) else opt$max_unit
start <- opt$unit_offset + (opt$task_id - 1L) * opt$reps_per_job + 1L
end   <- min(opt$unit_offset + opt$task_id * opt$reps_per_job, max_unit)
if (start > max_unit) {
  # Task index beyond this block's last unit (padding in final array) -- skip.
  cat(sprintf("[task %d] no units (start %d > block end %d); exiting.\n",
              opt$task_id, start, max_unit))
  quit(save = "no", status = 0)
}
block <- ut[start:end, , drop = FALSE]

# --- Drop infeasible cells ----------------------------------------------------
# is_feasible() (config/grid.R) excludes (method,n,dgp) cells whose Rashomon set
# exceeds feasible memory at the stress DGP (see INFEASIBLE_CELLS). Unit numbering
# is unchanged; we simply do not run these units. If the whole block is excluded,
# write nothing and exit 0 (so the array task is "done", not failed).
n_before <- nrow(block)
block <- block[is_feasible(block), , drop = FALSE]
n_excluded <- n_before - nrow(block)
if (n_excluded > 0) {
  cat(sprintf("[task %d] excluded %d infeasible unit(s) (INFEASIBLE_CELLS).\n",
              opt$task_id, n_excluded))
}
if (nrow(block) == 0) {
  cat(sprintf("[task %d] all %d units infeasible; nothing to run; exiting.\n",
              opt$task_id, n_before))
  quit(save = "no", status = 0)
}

# --- Idempotent skip: whole task already assembled ---------------------------
out_file <- file.path(opt$scratch_dir,
                      sprintf("task_%s.rds", formatC(opt$task_id, width = 6, flag = "0")))
if (file.exists(out_file)) {
  cat(sprintf("[task %d] result already exists (%s); skipping.\n",
              opt$task_id, out_file))
  quit(save = "no", status = 0)
}

# --- Per-unit checkpoint dir (M3) --------------------------------------------
# partials/ lives in this task's scratch dir (the per-method subdir under
# submit_per_method.sh), keyed by ABSOLUTE unit id so it is unambiguous even
# across methods and resubmits.
partials_dir <- file.path(opt$scratch_dir, "partials")
dir.create(partials_dir, recursive = TRUE, showWarnings = FALSE)
partial_path <- function(u) file.path(partials_dir,
                                      sprintf("unit_%s.rds", formatC(u, width = 7, flag = "0")))

# --- rbind that tolerates differing columns ----------------------------------
# Error rows (M4) carry only the identity columns + error_msg, so they lack the
# full success schema; bind by name, filling absent columns with NA.
.rbind_fill <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(NULL)
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)))
  }
  cols <- unique(unlist(lapply(rows, names)))
  filled <- lapply(rows, function(r) {
    for (c in setdiff(cols, names(r))) r[[c]] <- NA
    r[, cols, drop = FALSE]
  })
  do.call(rbind, filled)
}

# --- Run one unit, persisting an error MESSAGE rather than a silent all-NA row -
# (M4) A thrown error in run_one() (estimator/DGP failure) would otherwise kill
# the whole block; instead capture it, record the message, and keep going so one
# bad unit never loses its block-mates' completed work.
run_unit <- function(unit_row) {
  err <- NA_character_
  res <- tryCatch(run_one(unit_row),
                  error = function(e) { err <<- conditionMessage(e); NULL })
  if (is.null(res)) {
    # Minimal identity row + the failure text; combine.R fills the rest with NA.
    res <- data.frame(unit = unit_row$unit, config_id = unit_row$config_id,
                      rep_id = unit_row$rep_id, estimate = NA_real_,
                      stringsAsFactors = FALSE)
  }
  # A success row has no error_msg yet; an error row has err. Never overwrite a
  # non-NA error_msg that run_one() might itself have set.
  if (!("error_msg" %in% names(res))) res$error_msg <- err
  else if (is.na(res$error_msg[1]))   res$error_msg <- err
  res
}

# --- Run the block, checkpointing each unit (M3) -----------------------------
cat(sprintf("[task %d] running units %d-%d (%d units)\n",
            opt$task_id, start, end, nrow(block)))
t0 <- Sys.time()
n_done_resumed <- 0L
n_failed <- 0L
for (i in seq_len(nrow(block))) {
  u  <- block$unit[i]
  pf <- partial_path(u)
  if (file.exists(pf)) { n_done_resumed <- n_done_resumed + 1L; next }  # resume skip
  row <- run_unit(block[i, , drop = FALSE])
  if (!is.na(row$error_msg[1])) {
    n_failed <- n_failed + 1L
    cat(sprintf("[task %d] unit %d ERRORED: %s\n", opt$task_id, u, row$error_msg[1]))
  }
  tmp <- paste0(pf, ".tmp")
  saveRDS(row, tmp)
  file.rename(tmp, pf)                                  # atomic per-unit flush
}
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (n_done_resumed > 0)
  cat(sprintf("[task %d] resumed: %d unit(s) already done, skipped.\n",
              opt$task_id, n_done_resumed))
cat(sprintf("[task %d] block done in %.1f s (%d failed unit(s))\n",
            opt$task_id, elapsed, n_failed))

# --- Assemble the task file from this block's partials -----------------------
rows <- lapply(block$unit, function(u) {
  pf <- partial_path(u)
  if (file.exists(pf)) readRDS(pf) else NULL
})
result <- .rbind_fill(rows)
if (is.null(result) || nrow(result) < nrow(block)) {
  stop(sprintf("[task %d] only %s/%d unit partials present; not assembling task file.",
               opt$task_id, if (is.null(result)) 0 else nrow(result), nrow(block)))
}
result <- result[order(result$unit), , drop = FALSE]

# --- Write the task file atomically to scratch (tmp then rename) --------------
dir.create(opt$scratch_dir, recursive = TRUE, showWarnings = FALSE)
tmp <- paste0(out_file, ".tmp")
saveRDS(result, tmp)
file.rename(tmp, out_file)
cat(sprintf("[task %d] wrote %s (%d rows)\n", opt$task_id, out_file, nrow(result)))
