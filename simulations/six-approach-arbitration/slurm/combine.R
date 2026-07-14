#!/usr/bin/env Rscript
# =============================================================================
# combine.R -- aggregate per-task scratch results into ONE final file in home
# =============================================================================
# Reads every task_*.rds from the run's scratch dir, verifies completeness
# against the expected task count, streams them into a single data frame, and
# writes the ONLY home-directory artifact of the run: results/<run-id>.rds.
#
# Loud-failure guarantees (Constitution Section 9, no silent fallbacks):
#   * Errors if the scratch dir does not match the current grid hash (stale code).
#   * Reports every missing task id; refuses to write unless --allow-partial.
#
# Usage:
#   Rscript slurm/combine.R --run-id RID --scratch-dir DIR --study-dir DIR \
#           [--allow-partial]
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--run-id",       type = "character", dest = "run_id"),
  make_option("--scratch-dir",  type = "character", dest = "scratch_dir",
              help = "Run-specific scratch dir containing task_*.rds"),
  make_option("--study-dir",    type = "character", dest = "study_dir", default = ".."),
  make_option("--allow-partial", action = "store_true", dest = "allow_partial",
              default = FALSE, help = "Write result even if some tasks are missing")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$run_id), !is.null(opt$scratch_dir))

source(file.path(opt$study_dir, "config", "grid.R"))

# --- Stale-code guard: compare grid hash recorded at submit time -------------
# Fingerprint the study's SOURCE FILES (not just the GRID object) so any change
# to the design OR the DGP/estimator/run code invalidates an in-flight run.
# Files are hashed in a fixed order; must stay byte-identical to submit.sh.
CODE_FILES <- c("config/grid.R", "R/dgp.R", "R/estimators.R", "R/run_one.R")

code_hash <- function(study_dir) {
  # Bitwise-free polynomial rolling hash mod the Mersenne prime 2^31-1. All
  # intermediates stay < 2^53 so double arithmetic is exact (no bitwXor, which
  # overflows past 2^31 in base R; no digest dependency).
  MOD <- 2147483647          # 2^31 - 1
  h <- 0
  for (f in CODE_FILES) {
    bytes <- readBin(file.path(study_dir, f), what = "raw", n = file.size(file.path(study_dir, f)))
    for (b in as.integer(bytes)) h <- (h * 257 + b) %% MOD
  }
  sprintf("%.0f", h)
}

recorded_hash_file <- file.path(opt$scratch_dir, "GRID_HASH")
current_hash <- code_hash(opt$study_dir)
if (file.exists(recorded_hash_file)) {
  recorded <- trimws(readLines(recorded_hash_file, warn = FALSE)[1])
  if (!identical(recorded, current_hash)) {
    stop(sprintf(paste0(
      "STALE RESULTS: scratch dir was created with code hash %s but the study ",
      "source files now hash to %s.\nThe simulation code (grid/dgp/estimators/",
      "run_one) changed after this run started. Re-profile and re-submit, or ",
      "clean this run-id with clean.sh."),
      recorded, current_hash))
  }
} else {
  warning("No GRID_HASH found in scratch dir; cannot verify code freshness.")
}

# --- Discover task files -----------------------------------------------------
# Two layouts are supported:
#   (a) global (submit.sh):        <scratch>/task_NNNNNN.rds
#   (b) per-method (submit_per_method.sh): <scratch>/<method>/task_NNNNNN.rds
# recursive=TRUE finds both; completeness is checked per layout below.
files <- list.files(opt$scratch_dir, pattern = "^task_[0-9]+\\.rds$",
                    full.names = TRUE, recursive = TRUE)
if (length(files) == 0) stop(sprintf("No task_*.rds files in %s", opt$scratch_dir))

per_method_envs <- list.files(file.path(opt$study_dir, "config"),
                              pattern = "^sizing_.*\\.env$", full.names = TRUE)

get_total_tasks <- function(env_file) {
  kv <- readLines(env_file, warn = FALSE)
  tt <- grep("^TOTAL_TASKS=", kv, value = TRUE)
  if (length(tt)) as.integer(sub("^TOTAL_TASKS=", "", tt[1])) else NA_integer_
}

if (length(per_method_envs) > 0) {
  # Per-method: check each method's subdir against its own TOTAL_TASKS.
  total_missing <- 0L
  for (env_file in per_method_envs) {
    meth <- sub("^sizing_(.*)\\.env$", "\\1", basename(env_file))
    exp_m <- get_total_tasks(env_file)
    mdir <- file.path(opt$scratch_dir, meth)
    mfiles <- list.files(mdir, pattern = "^task_[0-9]+\\.rds$", full.names = FALSE)
    found_m <- as.integer(sub("^task_0*([0-9]+)\\.rds$", "\\1", mfiles))
    miss_m <- if (is.na(exp_m)) integer(0) else setdiff(seq_len(exp_m), found_m)
    if (length(miss_m) > 0) {
      total_missing <- total_missing + length(miss_m)
      cat(sprintf("  [%s] MISSING %d/%d tasks: %s\n", meth, length(miss_m), exp_m,
                  paste(head(miss_m, 20), collapse = ", ")))
    }
  }
  if (total_missing > 0) {
    msg <- sprintf("MISSING %d tasks across methods (see above).", total_missing)
    if (opt$allow_partial) warning(msg, "\nProceeding with --allow-partial.")
    else stop(msg, "\nRefusing to write partial result. Re-run missing tasks or pass --allow-partial.")
  }
} else {
  # Global layout: single TOTAL_TASKS from sizing.env.
  found_ids <- as.integer(sub("^task_0*([0-9]+)\\.rds$", "\\1", basename(files)))
  expected_tasks <- NA_integer_
  sizing_env <- file.path(opt$study_dir, "config", "sizing.env")
  if (file.exists(sizing_env)) expected_tasks <- get_total_tasks(sizing_env)
  if (!is.na(expected_tasks)) {
    missing <- setdiff(seq_len(expected_tasks), found_ids)
    if (length(missing) > 0) {
      msg <- sprintf("MISSING %d/%d tasks: %s", length(missing), expected_tasks,
                     paste(head(missing, 50), collapse = ", "))
      if (opt$allow_partial) warning(msg, "\nProceeding with --allow-partial.")
      else stop(msg, "\nRefusing to write partial result. Re-run missing tasks or pass --allow-partial.")
    }
  }
}

# --- Stream + bind (newest file LAST so the newest row wins de-dup) -----------
# Read task files in ASCENDING mtime order so that when two files re-emit the
# same (unit, method) -- e.g. a timeout resume/backfill rewrote a task -- the
# NEWEST file's rows come last and win the `fromLast=TRUE` de-dup below.
cat(sprintf("Combining %d task files...\n", length(files)))
files <- files[order(file.info(files)$mtime)]
parts <- lapply(files, readRDS)
# bind_fill: error rows (M4) carry only identity cols + error_msg, so per-file
# schemas can differ; bind by name, filling absent columns with NA.
if (requireNamespace("data.table", quietly = TRUE)) {
  result <- as.data.frame(data.table::rbindlist(parts, use.names = TRUE, fill = TRUE))
} else {
  cat("NOTE: data.table not installed; using do.call(rbind) (needs matching cols).\n")
  cols  <- unique(unlist(lapply(parts, names)))
  parts <- lapply(parts, function(d) { for (c in setdiff(cols, names(d))) d[[c]] <- NA; d[, cols, drop = FALSE] })
  result <- do.call(rbind, parts)
}

# --- De-duplicate by (unit, method) keeping the NEWEST -------------------------
# A `unit` is the canonical (config, rep) id; because `method` is a GRID dimension
# here, unit is already globally unique across methods, so (unit, method) keys on a
# superset that can never over-collapse legitimate rows -- but it is the robust,
# self-documenting key if the schema ever grows multiple rows per unit. Task files
# overlap in scratch: per-method subdirs restart task ids at 1, and timeout-driven
# resumes/backfills re-emit a unit's file; the recursive glob above sweeps them all.
# Left unchecked this rbinds the SAME (unit, method) multiple times (observed 3.75x
# in an interim run), silently inflating N and corrupting every downstream coverage/
# bias summary. We keep the NEWEST occurrence (fromLast=TRUE, after mtime-ascending
# read): a resume that RE-RAN a unit -- e.g. after a code fix within the same run --
# is authoritative over the stale earlier copy, not a coin-flip "keep first". Warn
# loudly (Constitution S9: no silent absorption) so an unexpected overlap surfaces.
key <- if ("method" %in% names(result)) {
  paste(result$unit, result$method, sep = "\r")
} else {
  as.character(result$unit)
}
n_raw <- nrow(result)
dup   <- duplicated(key, fromLast = TRUE)   # TRUE on all but the LAST (newest) copy
result <- result[!dup, , drop = FALSE]
result <- result[order(result$unit), , drop = FALSE]
n_dup <- sum(dup)
if (n_dup > 0) {
  warning(sprintf(paste0(
    "De-duplicated %d duplicate (unit, method) row(s) (%.2fx inflation): %d raw rows -> ",
    "%d unique. Overlapping scratch task files (e.g. timeout resumes/backfills). Kept ",
    "the NEWEST file per (unit, method). Inspect scratch if the factor is unexpected."),
    n_dup, n_raw / nrow(result), n_raw, nrow(result)), call. = FALSE)
}
cat(sprintf("Combined %d unique (unit, method) rows (from %d raw; expected %d units).\n",
            nrow(result), n_raw, n_units()))

# --- Surface failed replications (M4: no silent all-NA) ----------------------
if ("error_msg" %in% names(result)) {
  n_err <- sum(!is.na(result$error_msg))
  if (n_err > 0) {
    cat(sprintf("WARNING: %d/%d row(s) (%.1f%%) carry an error_msg (failed replications).\n",
                n_err, nrow(result), 100 * n_err / nrow(result)))
    ex <- head(unique(result$error_msg[!is.na(result$error_msg)]), 3)
    cat("  example error(s):\n"); for (e in ex) cat(sprintf("    - %s\n", e))
  } else {
    cat("All units succeeded (no error_msg set).\n")
  }
}

# --- Attach provenance metadata ----------------------------------------------
attr(result, "run_id")    <- opt$run_id
attr(result, "grid_hash") <- current_hash
attr(result, "combined_at") <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

# --- Write the ONE home artifact ---------------------------------------------
results_dir <- file.path(opt$study_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(results_dir, sprintf("%s.rds", opt$run_id))
saveRDS(result, out)
cat(sprintf("Wrote FINAL result to home: %s\n", out))
cat("Scratch per-task files can now be cleaned with clean.sh.\n")
