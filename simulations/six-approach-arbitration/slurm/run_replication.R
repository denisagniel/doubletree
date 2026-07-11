#!/usr/bin/env Rscript
# =============================================================================
# run_replication.R -- run ONE array task's worth of work units
# =============================================================================
# Invoked by array.slurm once per SLURM array index. Each task runs a contiguous
# block of `reps_per_job` work units, then writes ONE result file to the run's
# scratch directory. Output is idempotent: if the task's file already exists it
# exits immediately (cheap resume after partial failure).
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

# --- Idempotent skip ----------------------------------------------------------
out_file <- file.path(opt$scratch_dir,
                      sprintf("task_%s.rds", formatC(opt$task_id, width = 6, flag = "0")))
if (file.exists(out_file)) {
  cat(sprintf("[task %d] result already exists (%s); skipping.\n",
              opt$task_id, out_file))
  quit(save = "no", status = 0)
}

# --- Run the block ------------------------------------------------------------
cat(sprintf("[task %d] running units %d-%d (%d units)\n",
            opt$task_id, start, end, nrow(block)))
t0 <- Sys.time()
rows <- vector("list", nrow(block))
for (i in seq_len(nrow(block))) {
  rows[[i]] <- run_one(block[i, , drop = FALSE])
}
result <- do.call(rbind, rows)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("[task %d] done in %.1f s (%.2f s/unit)\n",
            opt$task_id, elapsed, elapsed / nrow(block)))

# --- Write atomically to scratch (tmp then rename) ---------------------------
dir.create(opt$scratch_dir, recursive = TRUE, showWarnings = FALSE)
tmp <- paste0(out_file, ".tmp")
saveRDS(result, tmp)
file.rename(tmp, out_file)
cat(sprintf("[task %d] wrote %s (%d rows)\n", opt$task_id, out_file, nrow(result)))
