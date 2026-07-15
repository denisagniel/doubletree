#!/usr/bin/env Rscript
# =============================================================================
# peek_stale.R -- READ-ONLY interim look at a run whose code hash went STALE
# =============================================================================
# combine.R refuses to aggregate a run once the study source files change
# (Constitution S9 stale-code guard). That guard is correct: those rows were
# produced by the OLD estimator code and must NOT become the official
# results/<run-id>.rds. This script lets you LOOK at the partial results
# anyway, with three hard guarantees so the interim data can never be mistaken
# for a real combine:
#
#   1. It writes to results/INTERIM_STALE_<run-id>.rds  (never <run-id>.rds).
#   2. It stamps stale=TRUE + the recorded vs current hash onto the object.
#   3. It prints a loud banner and per-method coverage/bias so you can eyeball
#      the arbitration BEFORE re-profiling and re-submitting on current code.
#
# It reuses combine.R's mtime-ordered read + (unit, method) de-dup so the N is
# not inflated by overlapping scratch files (submit resumes/backfills).
#
# Usage (on O2, from the study dir):
#   Rscript slurm/peek_stale.R --run-id RID --scratch-dir DIR --study-dir DIR
# =============================================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--run-id",      type = "character", dest = "run_id"),
  make_option("--scratch-dir", type = "character", dest = "scratch_dir"),
  make_option("--study-dir",   type = "character", dest = "study_dir", default = "..")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$run_id), !is.null(opt$scratch_dir))

# --- Read recorded vs current hash for the banner (do NOT stop on mismatch) ---
recorded <- NA_character_
rh_file <- file.path(opt$scratch_dir, "GRID_HASH")
if (file.exists(rh_file)) recorded <- trimws(readLines(rh_file, warn = FALSE)[1])

# --- Discover task files (both layouts), newest LAST for correct de-dup -------
files <- list.files(opt$scratch_dir, pattern = "^task_[0-9]+\\.rds$",
                    full.names = TRUE, recursive = TRUE)
if (length(files) == 0) stop(sprintf("No task_*.rds files in %s", opt$scratch_dir))
files <- files[order(file.info(files)$mtime)]
cat(sprintf("Reading %d task files (stale-code interim look)...\n", length(files)))
parts <- lapply(files, readRDS)

if (requireNamespace("data.table", quietly = TRUE)) {
  result <- as.data.frame(data.table::rbindlist(parts, use.names = TRUE, fill = TRUE))
} else {
  cols  <- unique(unlist(lapply(parts, names)))
  parts <- lapply(parts, function(d) { for (c in setdiff(cols, names(d))) d[[c]] <- NA; d[, cols, drop = FALSE] })
  result <- do.call(rbind, parts)
}

# --- De-dup by (unit, method) keeping NEWEST (mirror combine.R:157-164) --------
key <- if ("method" %in% names(result)) paste(result$unit, result$method, sep = "\r") else as.character(result$unit)
n_raw <- nrow(result)
dup   <- duplicated(key, fromLast = TRUE)
result <- result[!dup, , drop = FALSE]
result <- result[order(result$unit), , drop = FALSE]
if (sum(dup) > 0) {
  cat(sprintf("De-duplicated %d overlapping row(s) (%.2fx): %d raw -> %d unique.\n",
              sum(dup), n_raw / nrow(result), n_raw, nrow(result)))
}

# --- Loud stale banner --------------------------------------------------------
cat("\n=============================================================\n")
cat(" INTERIM STALE-CODE LOOK -- NOT AN OFFICIAL COMBINE\n")
cat(sprintf(" run-id        : %s\n", opt$run_id))
cat(sprintf(" recorded hash : %s  (code the run was submitted with)\n", recorded))
cat(" NOTE: current source hashes differently; these rows are from the\n")
cat("       OLD estimator code. Re-profile + re-submit for real results.\n")
cat("=============================================================\n\n")

# --- Per-method coverage / bias eyeball (goal-(i) and goal-(ii)) --------------
if (all(c("method", "covered", "error") %in% names(result))) {
  agg <- aggregate(
    cbind(n_rows = rep(1, nrow(result)),
          coverage = result$covered,
          bias = result$error,
          cover_cf = result$covered_crossfit) ~ method,
    data = transform(result, cover_cf = result$covered_crossfit),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  # n_rows via table (the formula mean above turns it into a proportion)
  counts <- as.data.frame(table(result$method), stringsAsFactors = FALSE)
  names(counts) <- c("method", "n_rows")
  agg$n_rows <- counts$n_rows[match(agg$method, counts$method)]
  cat("Per-method summary (interim, stale code):\n")
  print(agg, row.names = FALSE)
}

# --- Write quarantined artifact (NEVER results/<run-id>.rds) -------------------
attr(result, "run_id")        <- opt$run_id
attr(result, "stale")         <- TRUE
attr(result, "recorded_hash") <- recorded
attr(result, "combined_at")   <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

results_dir <- file.path(opt$study_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(results_dir, sprintf("INTERIM_STALE_%s.rds", opt$run_id))
saveRDS(result, out)
cat(sprintf("\nWrote quarantined interim file: %s\n", out))
cat("scp THAT single file back to inspect locally. Do NOT treat it as final.\n")
