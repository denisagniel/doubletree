#!/usr/bin/env Rscript
# ============================================================
# slurm/combine_results.R
# Study: head_to_head_comparison  (doubletree, S5)
#
# Run ON O2, after every array has finished:
#   Rscript simulations/head_to_head_comparison/slurm/combine_results.R \
#     --scratch-dir /n/scratch/.../head_to_head_comparison/<RUN_ID> \
#     --pkg-root /path/to/doubletree --run-id <RUN_ID>
#
# (launch_subset.sh writes this exact command, with the paths already filled in,
# into the study's MANIFEST.md.)
#
# ------------------------------------------------------------
# WHAT IT PRODUCES, AND WHY THAT SHAPE
# ------------------------------------------------------------
#
# The scratch tree holds one shard per array task: a BLOCK of replications of one
# (regime, dgp, n) cell. code/analyze.R's load_results() globs
# results/cell_*.rds and expects each file to be ONE CELL's payload --
# list(results, dgp_summary, meta) -- exactly as run_cell() writes it.
#
# So this script does not invent a "combined results" format. It reassembles the
# shards of each cell into the SAME payload run_cell() would have written, at the
# SAME path cell_path() would have chosen. Every downstream consumer
# (load_results(), summarise_cell(), rmse_ratio_table(), run_sweep.R's tables)
# then works on cluster output with no changes and no cluster-specific branch. A
# cluster-only "combined.rds" would have forked the analysis path in two, and the
# fork is where the two stop agreeing.
#
# ------------------------------------------------------------
# COMPLETENESS IS CHECKED, NOT ASSUMED
# ------------------------------------------------------------
#
# A cell is written only when its shards hold replications 1..reps_total with no
# gaps. A cell missing shards is REPORTED with the array task ids to resubmit and
# then skipped, because a checkpoint named `..._r2000.rds` that holds 1,400
# replications is worse than a missing one: every MC standard error computed from
# it is right for 1,400 and reported as if the study had run 2,000.
#
# --allow-incomplete writes short cells anyway, but names them for the number of
# replications actually present, so the filename cannot lie.
#
# ------------------------------------------------------------
# THE PILOT-CHECKPOINT COLLISION (a real hazard in this study's results/ dir)
# ------------------------------------------------------------
#
# load_results() globs EVERY cell_*.rds in results/ and rbinds them. The 2026-09-09
# pilot left `cell_R1_A_n500_r30.rds` and friends in that directory. Those files do
# not collide by NAME with the sweep's `cell_R1_A_n500_r1000.rds` -- but they DO
# collide in load_results()'s output, which would then hold replications 1..30
# twice for every cell the pilot also ran, at pilot scale, silently doubling their
# weight in every average. This script refuses to write into a results/ directory
# that still holds pilot-scale checkpoints, and --archive-pilot moves them aside.
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
})

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (!length(hit)) return(NA_character_)
  normalizePath(dirname(sub("^--file=", "", hit[[1L]])), mustWork = FALSE)
}

default_pkg_root <- function() {
  d <- script_dir()
  if (is.na(d)) return(getwd())
  normalizePath(file.path(d, "..", "..", ".."), mustWork = FALSE)
}

option_list <- list(
  make_option("--scratch-dir", type = "character", dest = "scratch_dir", default = "",
              help = "Run's scratch dir (holds one subdir of shards per regime) [required]"),
  make_option("--pkg-root", type = "character", dest = "pkg_root",
              default = default_pkg_root(),
              help = "doubletree package root [default: three levels above this script]"),
  make_option("--regimes", type = "character", dest = "regimes", default = "",
              help = "Comma-separated regime ids to combine [default: whatever the shards contain]"),
  make_option("--run-id", type = "character", dest = "run_id", default = NA_character_,
              help = "Run id, stamped into each cell's metadata"),
  make_option("--out-dir", type = "character", dest = "out_dir", default = "",
              help = "Where cell_*.rds go [default: the study's results/ dir]"),
  make_option("--target-secs", type = "double", dest = "target_secs", default = NA_real_,
              help = "Target seconds per task used at LAUNCH; needed to reproduce the expected shard set"),
  make_option("--allow-incomplete", action = "store_true", dest = "allow_incomplete",
              default = FALSE,
              help = "Write cells that are missing shards, named for the replications actually present"),
  make_option("--archive-pilot", action = "store_true", dest = "archive_pilot",
              default = FALSE,
              help = "Move pilot-scale cell_*.rds out of the output dir instead of aborting"),
  make_option("--dry-run", action = "store_true", dest = "dry_run", default = FALSE,
              help = "Report what would be written, write nothing")
)

if (!interactive() && sys.nframe() == 0L) {

  opt <- parse_args(OptionParser(option_list = option_list))

  # ---- 1 assert the inputs -----------------------------------------------
  if (!nzchar(opt$scratch_dir)) stop("--scratch-dir is required.", call. = FALSE)
  if (!dir.exists(opt$scratch_dir)) {
    stop("--scratch-dir does not exist: ", opt$scratch_dir, call. = FALSE)
  }
  if (!file.exists(file.path(opt$pkg_root, "DESCRIPTION"))) {
    stop("--pkg-root has no DESCRIPTION: ", opt$pkg_root, call. = FALSE)
  }

  setwd(opt$pkg_root)
  Sys.setenv(H2H_USE_INSTALLED = "1")
  if (!is.na(opt$target_secs)) {
    Sys.setenv(H2H_TARGET_SECS = format(opt$target_secs, scientific = FALSE))
  }
  source(file.path("simulations", "head_to_head_comparison", "code", "common.R"))
  source(file.path("simulations", "head_to_head_comparison", "slurm", "units.R"))

  out_dir <- if (nzchar(opt$out_dir)) opt$out_dir else DIR_RESULTS
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cli::cli_h1("Combine shards -> cell checkpoints")
  cli::cli_inform(c(
    "*" = "scratch: {.path {opt$scratch_dir}}",
    "*" = "output : {.path {out_dir}}",
    "*" = "run id : {opt$run_id}",
    "*" = "target secs/task: {h2h_target_secs()} (must match the launch value)"
  ))

  # ---- 2 read every shard -------------------------------------------------
  # [.]rds rather than \\.rds: this file is read by Rscript from disk, so the
  # backslash form would be fine here, but the two spellings are kept identical
  # across the slurm/ scripts so the pattern can be grepped for as one string.
  files <- list.files(opt$scratch_dir, pattern = "^batch_.*[.]rds$",
                      full.names = TRUE, recursive = TRUE)
  if (!length(files)) {
    cli::cli_abort(c(
      "No shards ({.file batch_*.rds}) under {.path {opt$scratch_dir}}.",
      i = "Have the arrays finished? {.code squeue -u $USER} / {.code bash slurm/check_progress.sh}"
    ))
  }
  cli::cli_alert_info("{length(files)} shard file(s) found.")

  shards <- lapply(files, function(f) {
    p <- readRDS(f)
    if (!all(c("unit", "results") %in% names(p))) {
      cli::cli_abort("Shard {.file {basename(f)}} is not a task shard (no unit/results).")
    }
    p$file <- f
    p
  })

  keep_regimes <- if (nzchar(opt$regimes)) {
    trimws(strsplit(opt$regimes, ",", fixed = TRUE)[[1L]])
  } else {
    unique(vapply(shards, function(p) p$unit$regime, character(1)))
  }
  unknown <- setdiff(keep_regimes, REGIME_IDS)
  if (length(unknown)) cli::cli_abort("Unknown regime(s) {.val {unknown}}.")
  shards <- shards[vapply(shards, function(p) p$unit$regime %in% keep_regimes, logical(1))]

  # ---- 3 pilot-collision guard -------------------------------------------
  # Any existing cell file for a (regime, dgp, n) this run also produced, at a
  # DIFFERENT reps_total, would be rbound alongside the new one by load_results().
  cell_key <- function(regime, dgp, n) paste(regime, dgp, n, sep = "|")
  new_keys <- unique(vapply(shards,
                            function(p) cell_key(p$unit$regime, p$unit$dgp, p$unit$n),
                            character(1)))
  new_reps <- vapply(shards, function(p) as.integer(p$unit$reps_total), integer(1))
  existing <- list.files(out_dir, pattern = "^cell_.*[.]rds$", full.names = TRUE)
  stale <- character(0)
  for (f in existing) {
    m <- regmatches(basename(f),
                    regexec("^cell_(.+)_(n[0-9]+)_r([0-9]+)[.]rds$", basename(f)))[[1L]]
    if (length(m) != 4L) next
    n_val <- as.integer(sub("^n", "", m[[3L]]))
    reps_val <- as.integer(m[[4L]])
    # m[[2]] is "<regime>_<dgp>"; the regime is the first underscore-delimited field.
    head_parts <- strsplit(m[[2L]], "_", fixed = TRUE)[[1L]]
    regime <- head_parts[[1L]]
    dgp <- paste(head_parts[-1L], collapse = "_")
    if (cell_key(regime, dgp, n_val) %in% new_keys && !(reps_val %in% new_reps)) {
      stale <- c(stale, f)
    }
  }
  if (length(stale)) {
    if (opt$archive_pilot) {
      arch <- file.path(out_dir, sprintf("archive_pilot_%s", format(Sys.time(), "%Y%m%d-%H%M%S")))
      dir.create(arch, recursive = TRUE, showWarnings = FALSE)
      ok <- file.rename(stale, file.path(arch, basename(stale)))
      if (!all(ok)) cli::cli_abort("Could not move every pilot checkpoint into {.path {arch}}.")
      cli::cli_alert_success("Moved {length(stale)} pilot-scale checkpoint(s) to {.path {arch}}.")
    } else {
      cli::cli_abort(c(
        "{length(stale)} pilot-scale checkpoint(s) in {.path {out_dir}} cover cells this run also produced.",
        x = "{.fun load_results} globs every {.file cell_*.rds} and would rbind BOTH, double-counting the pilot's replications at pilot scale.",
        i = "Re-run with {.code --archive-pilot} to move them aside, or move them yourself:",
        " " = "mkdir -p {out_dir}/archive_pilot && mv {paste(basename(stale), collapse = ' ')} {out_dir}/archive_pilot/"
      ))
    }
  }

  # ---- 4 group shards into cells and check completeness -------------------
  key_of <- function(p) paste(p$unit$regime, p$unit$dgp, p$unit$n, p$unit$reps_total,
                              sep = "|")
  groups <- split(shards, vapply(shards, key_of, character(1)))

  report <- list()
  written <- 0L
  for (k in names(groups)) {
    g <- groups[[k]]
    u0 <- g[[1L]]$unit
    regime <- u0$regime; dgp_id <- u0$dgp
    n <- as.integer(u0$n); reps_total <- as.integer(u0$reps_total)

    res <- do.call(rbind, lapply(g, `[[`, "results"))
    res <- res[order(res$rep, match(res$arm, REGIMES[[regime]]$arms)), , drop = FALSE]
    rownames(res) <- NULL

    have <- sort(unique(res$rep))
    missing_reps <- setdiff(seq_len(reps_total), have)
    # Which array tasks own the missing replications -- the resubmit instruction.
    ut <- h2h_unit_table(regime)
    missing_tasks <- integer(0)
    if (length(missing_reps)) {
      rows <- which(ut$dgp == dgp_id & ut$n == n)
      for (i in rows) {
        if (any(missing_reps >= ut$rep_from[i] & missing_reps <= ut$rep_to[i])) {
          missing_tasks <- c(missing_tasks, ut$task[i])
        }
      }
    }

    complete <- length(missing_reps) == 0L
    report[[length(report) + 1L]] <- data.frame(
      regime = regime, dgp = dgp_id, n = n, reps_target = reps_total,
      reps_present = length(have), shards = length(g),
      shards_expected = sum(ut$dgp == dgp_id & ut$n == n),
      complete = complete,
      missing_tasks = if (length(missing_tasks)) {
        paste(sort(unique(missing_tasks)), collapse = ",")
      } else NA_character_,
      stringsAsFactors = FALSE
    )

    if (!complete && !opt$allow_incomplete) {
      cli::cli_alert_danger(
        "{regime}/{dgp_id}/n={n}: {length(have)}/{reps_total} reps -- SKIPPED. Resubmit task(s) {paste(sort(unique(missing_tasks)), collapse = ', ')}."
      )
      next
    }

    reps_written <- if (complete) reps_total else length(have)
    out <- cell_path(regime, dgp_id, n, reps_written)
    if (nzchar(opt$out_dir)) out <- file.path(out_dir, basename(out))

    dgp <- make_dgp(dgp_id)
    meta <- cell_metadata(regime, dgp, n, reps_written, res)
    # Cluster provenance, alongside the metadata run_cell() already stamps. A
    # results file that is argued over should be able to name the run, the hosts
    # and the shards it came from without anyone consulting a shell history.
    meta$cluster <- list(
      run_id = opt$run_id,
      scratch_dir = opt$scratch_dir,
      target_secs = h2h_target_secs(),
      n_shards = length(g),
      shard_files = basename(vapply(g, `[[`, character(1), "file")),
      hosts = unique(vapply(g, function(p) p$runtime$host, character(1))),
      slurm_job_ids = unique(vapply(g, function(p) p$runtime$slurm_job, character(1))),
      elapsed_secs_sum = sum(vapply(g, function(p) p$runtime$elapsed_secs, numeric(1))),
      complete = complete,
      reps_target = reps_total
    )

    payload <- list(results = res, dgp_summary = g[[1L]]$dgp_summary, meta = meta)

    if (opt$dry_run) {
      cli::cli_alert_info("would write {.file {basename(out)}} ({nrow(res)} rows, {length(have)} reps)")
    } else {
      tmp <- paste0(out, ".tmp")
      saveRDS(payload, tmp)
      file.rename(tmp, out)
      written <- written + 1L
      cli::cli_alert_success(
        "{basename(out)} -- {length(have)}/{reps_total} reps x {length(unique(res$arm))} arms = {nrow(res)} rows, {meta$n_failed_reps} failed rep(s)"
      )
    }
  }

  # ---- 5 the report -------------------------------------------------------
  rep_df <- do.call(rbind, report)
  rep_df <- rep_df[order(match(rep_df$regime, REGIME_IDS), rep_df$dgp, rep_df$n), ,
                   drop = FALSE]
  rownames(rep_df) <- NULL
  cli::cli_h2("Cell completeness")
  print(rep_df, row.names = FALSE)

  if (!opt$dry_run) {
    stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
    csv <- file.path(DIR_TABLES, sprintf("combine_report_%s.csv", stamp))
    utils::write.csv(rep_df, csv, row.names = FALSE)
    cli::cli_alert_success("Wrote {written} cell checkpoint(s); report at {.file {basename(csv)}}.")
  }

  n_incomplete <- sum(!rep_df$complete)
  if (n_incomplete) {
    cli::cli_warn(c(
      "{n_incomplete} cell(s) are incomplete.",
      i = "Resubmit with the SAME run id -- finished tasks exit immediately and interrupted ones resume from their last flush:",
      " " = "RUN_ID={opt$run_id} bash slurm/launch_subset.sh {paste(unique(rep_df$regime[!rep_df$complete]), collapse = ' ')}"
    ))
  } else {
    cli::cli_alert_success("Every cell is complete. Next: the tables in code/run_sweep.R section 3, or code/analyze.R directly.")
  }
  cli::cli_h1("Read spec §0 before quoting any number. This is not a leaderboard.")
}
