#!/usr/bin/env Rscript
# ============================================================
# slurm/combine_results.R
# Study: honest_inference_sparsity_failure  (doubletree, S2)
#
# Run ON O2, after every array has finished:
#   Rscript --no-init-file simulations/honest_inference_sparsity_failure/slurm/combine_results.R \
#     --scratch-dir /n/scratch/.../honest_inference_sparsity_failure/<RUN_ID> \
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
# (dgp, n) cell. code/analyze.R does NOT glob -- it constructs the exact path
# cell_path(dgp, n, reps_for_n(n)) for all 15 cells of design_grid() and aborts
# naming any that is absent. Each of those files must be ONE CELL's payload,
# list(results, dgp_summary, meta), exactly as run_cell() writes it.
#
# So this script does not invent a "combined results" format. It reassembles the
# shards of each cell into the SAME payload run_cell() would have written, at the
# SAME path cell_path() would have chosen. analyze.R, run_sweep.R and every table
# then work on cluster output with no changes and no cluster-specific branch. A
# cluster-only "combined.rds" would have forked the analysis path in two, and the
# fork is where the two stop agreeing.
#
# The five n = 500 cells already on disk are NOT touched and NOT re-read: they are
# already in exactly this format, written sequentially on 2026-09-01, and
# rep_seed()'s purity means the cells this script writes are directly comparable to
# them (slurm/units.R, "batching is exactly reproducible").
#
# ------------------------------------------------------------
# COMPLETENESS IS CHECKED, NOT ASSUMED
# ------------------------------------------------------------
#
# A cell is written only when its shards hold replications 1..reps_total with no
# gaps. A cell missing shards is REPORTED with the array task ids to resubmit and
# then skipped, because a checkpoint named `..._r300.rds` that holds 210
# replications is worse than a missing one: every MC standard error computed from it
# is right for 210 and reported as if the study had run 300, and analyze.R -- which
# looks the file up by that exact name -- would accept it silently.
#
# --allow-incomplete writes short cells anyway, but names them for the number of
# replications actually present, so the filename cannot lie. Note that analyze.R
# will then NOT find them (it looks for r300 / r150) and will abort naming them,
# which is the intended outcome: a short cell is inspectable but not analysable.
#
# ------------------------------------------------------------
# THE REPS-MISMATCH GUARD
# ------------------------------------------------------------
#
# reps_for_n() is driven by HIS_REPS / HIS_REPS_MAX_N. If the arrays were launched
# under one value and this script runs under another, the shards' own reps_total
# (recorded in each unit at launch time) would disagree with what cell_path() now
# computes, and the assembled cell would be written under a filename describing a
# replication count it does not have. That is checked per cell and aborts.
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
              help = "Run's scratch dir (holds one subdir of shards per wave) [required]"),
  make_option("--pkg-root", type = "character", dest = "pkg_root",
              default = default_pkg_root(),
              help = "doubletree package root [default: three levels above this script]"),
  make_option("--waves", type = "character", dest = "waves", default = "",
              help = "Comma-separated wave ids to combine [default: whatever the shards contain]"),
  make_option("--run-id", type = "character", dest = "run_id", default = NA_character_,
              help = "Run id, stamped into each cell's metadata"),
  make_option("--out-dir", type = "character", dest = "out_dir", default = "",
              help = "Where cell_*.rds go [default: the study's results/ dir]"),
  make_option("--target-secs", type = "double", dest = "target_secs", default = NA_real_,
              help = "Target seconds per task used at LAUNCH; needed to reproduce the expected shard set"),
  make_option("--allow-incomplete", action = "store_true", dest = "allow_incomplete",
              default = FALSE,
              help = "Write cells that are missing shards, named for the replications actually present"),
  make_option("--overwrite", action = "store_true", dest = "overwrite", default = FALSE,
              help = "Overwrite an existing cell checkpoint instead of refusing"),
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
  Sys.setenv(HIS_USE_INSTALLED = "1")
  if (!is.na(opt$target_secs)) {
    Sys.setenv(HIS_TARGET_SECS = format(opt$target_secs, scientific = FALSE))
  }
  source(file.path("simulations", "honest_inference_sparsity_failure", "code",
                   "common.R"))
  source(file.path("simulations", "honest_inference_sparsity_failure", "slurm",
                   "units.R"))

  out_dir <- if (nzchar(opt$out_dir)) opt$out_dir else DIR_RESULTS
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cli::cli_h1("Combine shards -> cell checkpoints (S2)")
  cli::cli_inform(c(
    "*" = "scratch: {.path {opt$scratch_dir}}",
    "*" = "output : {.path {out_dir}}",
    "*" = "run id : {opt$run_id}",
    "*" = "reps   : HIS_REPS={REPS_DEFAULT} HIS_REPS_MAX_N={REPS_AT_MAX_N} (must match the launch values)",
    "*" = "target secs/task: {his_target_secs()} (must match the launch value)"
  ))

  # ---- 2 read every shard -------------------------------------------------
  # [.]rds rather than \\.rds: the two spellings are kept identical across the
  # slurm/ scripts so the pattern can be grepped for as one string, and the bracket
  # form survives `Rscript -e`'s own backslash processing.
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

  keep_waves <- if (nzchar(opt$waves)) {
    trimws(strsplit(opt$waves, ",", fixed = TRUE)[[1L]])
  } else {
    unique(vapply(shards, function(p) p$unit$wave, character(1)))
  }
  unknown <- setdiff(keep_waves, WAVE_IDS)
  if (length(unknown)) cli::cli_abort("Unknown wave(s) {.val {unknown}}.")
  shards <- shards[vapply(shards, function(p) p$unit$wave %in% keep_waves, logical(1))]

  # ---- 3 group shards into cells and check completeness -------------------
  key_of <- function(p) paste(p$unit$wave, p$unit$dgp, p$unit$n, p$unit$reps_total,
                              sep = "|")
  groups <- split(shards, vapply(shards, key_of, character(1)))

  report <- list()
  written <- 0L
  for (k in names(groups)) {
    g <- groups[[k]]
    u0 <- g[[1L]]$unit
    wave <- u0$wave; dgp_id <- u0$dgp
    n <- as.integer(u0$n); reps_total <- as.integer(u0$reps_total)

    # The reps-mismatch guard (see the header). reps_for_n() is what cell_path()
    # and analyze.R use NOW; reps_total is what the array was built with.
    reps_now <- as.integer(reps_for_n(n))
    if (!identical(reps_now, reps_total)) {
      cli::cli_abort(c(
        "Shards for {dgp_id}/n={n} were produced with reps_total = {reps_total}, but {.fun reps_for_n}({n}) = {reps_now} in this session.",
        x = "The assembled cell would be written as {.file {basename(cell_path(dgp_id, n, reps_now))}} while holding {reps_total} replications.",
        i = "Re-run this script with the SAME HIS_REPS / HIS_REPS_MAX_N the arrays were launched under (see the run's MANIFEST.md)."
      ))
    }

    res <- do.call(rbind, lapply(g, `[[`, "results"))
    res <- res[order(res$rep), , drop = FALSE]
    rownames(res) <- NULL

    # A replication appearing twice means two shards claim the same block -- a map
    # drift, which must never be averaged over.
    if (anyDuplicated(res$rep)) {
      dup <- sort(unique(res$rep[duplicated(res$rep)]))
      cli::cli_abort(c(
        "{dgp_id}/n={n}: replication(s) {.val {dup}} appear in more than one shard.",
        i = "Two shards claim the same block: the array map drifted between launch and now (different {.envvar HIS_TARGET_SECS}?). Do not combine; re-launch with a fresh run id."
      ))
    }

    have <- sort(unique(res$rep))
    missing_reps <- setdiff(seq_len(reps_total), have)
    # Which array tasks own the missing replications -- the resubmit instruction.
    ut <- his_unit_table(wave)
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
      wave = wave, dgp = dgp_id, n = n, reps_target = reps_total,
      reps_present = length(have), shards = length(g),
      shards_expected = sum(ut$dgp == dgp_id & ut$n == n),
      n_failed = sum(!is.na(res$error_message)),
      complete = complete,
      missing_tasks = if (length(missing_tasks)) {
        paste(sort(unique(missing_tasks)), collapse = ",")
      } else NA_character_,
      stringsAsFactors = FALSE
    )

    if (!complete && !opt$allow_incomplete) {
      cli::cli_alert_danger(
        "{dgp_id}/n={n}: {length(have)}/{reps_total} reps -- SKIPPED. Resubmit task(s) {paste(sort(unique(missing_tasks)), collapse = ', ')} of wave {wave}."
      )
      next
    }

    reps_written <- if (complete) reps_total else length(have)
    out <- cell_path(dgp_id, n, reps_written)
    if (nzchar(opt$out_dir)) out <- file.path(out_dir, basename(out))

    if (file.exists(out) && !opt$overwrite) {
      cli::cli_alert_warning(
        "{basename(out)} already exists -- NOT overwritten. Pass {.code --overwrite} if that is intended."
      )
      next
    }

    dgp <- make_dgp(dgp_id)
    lambda_n <- log(n) / n
    meta <- cell_metadata(dgp, n, reps_written, lambda_n, res)
    # Cluster provenance, alongside the metadata run_cell() already stamps. A results
    # file that is argued over should be able to name the run, the hosts and the
    # shards it came from without anyone consulting a shell history.
    meta$cluster <- list(
      run_id = opt$run_id,
      wave = wave,
      scratch_dir = opt$scratch_dir,
      target_secs = his_target_secs(),
      n_shards = length(g),
      shard_files = basename(vapply(g, `[[`, character(1), "file")),
      hosts = unique(vapply(g, function(p) p$runtime$host, character(1))),
      slurm_job_ids = unique(vapply(g, function(p) p$runtime$slurm_job, character(1))),
      pkg_load_paths = unique(vapply(g, function(p) p$runtime$pkg_load_path, character(1))),
      elapsed_secs_sum = sum(vapply(g, function(p) p$runtime$elapsed_secs, numeric(1))),
      complete = complete,
      reps_target = reps_total
    )
    # Every shard must have taken the installed path. A mixture would mean some
    # replications came from uninstalled source -- see run_single_replication.R's
    # assert_installed_load().
    if (!identical(meta$cluster$pkg_load_paths, "installed")) {
      cli::cli_abort(c(
        "{dgp_id}/n={n}: shards report package load path(s) {.val {meta$cluster$pkg_load_paths}}.",
        x = "At least one shard did not use the installed packages, so its replications are not reproducible from {.path {opt$pkg_root}}'s installed code.",
        i = "Re-run those tasks with {.code Rscript --no-init-file} (see run_simulations.slurm)."
      ))
    }

    payload <- list(results = res, dgp_summary = g[[1L]]$dgp_summary, meta = meta)

    if (opt$dry_run) {
      cli::cli_alert_info("would write {.file {basename(out)}} ({nrow(res)} reps)")
    } else {
      tmp <- paste0(out, ".tmp")
      saveRDS(payload, tmp)
      file.rename(tmp, out)
      written <- written + 1L
      cli::cli_alert_success(
        "{basename(out)} -- {length(have)}/{reps_total} reps, {meta$n_failures} failure(s), cov_anch={signif(mean(res$covered_anchor, na.rm = TRUE), 3)} cov_hon={signif(mean(res$covered_honest, na.rm = TRUE), 3)} cov_cf={signif(mean(res$covered_cf, na.rm = TRUE), 3)}"
      )
    }
  }

  # ---- 4 the report -------------------------------------------------------
  rep_df <- do.call(rbind, report)
  rep_df <- rep_df[order(rep_df$n, match(rep_df$dgp, DGP_IDS)), , drop = FALSE]
  rownames(rep_df) <- NULL
  cli::cli_h2("Cell completeness")
  print(rep_df, row.names = FALSE)

  if (!opt$dry_run) {
    stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
    csv <- file.path(DIR_TABLES, sprintf("combine_report_%s.csv", stamp))
    utils::write.csv(rep_df, csv, row.names = FALSE)
    cli::cli_alert_success("Wrote {written} cell checkpoint(s); report at {.file {basename(csv)}}.")
  }

  # ---- 5 is the WHOLE design grid now on disk? ----------------------------
  # The point of this deployment is to complete the spec's grid, so the last thing
  # printed is whether it IS complete -- checked against design_grid() and
  # cell_path(), which is precisely what analyze.R will do next.
  grid <- design_grid()
  grid$path <- vapply(seq_len(nrow(grid)),
                      function(i) cell_path(grid$dgp[[i]], grid$n[[i]], grid$reps[[i]]),
                      character(1))
  grid$on_disk <- file.exists(grid$path)
  cli::cli_h2("Spec grid ({nrow(grid)} cells = {length(DGP_IDS)} variants x {length(N_GRID)} n)")
  print(grid[, c("dgp", "eps", "n", "reps", "on_disk")], row.names = FALSE)

  n_incomplete <- sum(!rep_df$complete)
  if (n_incomplete) {
    cli::cli_warn(c(
      "{n_incomplete} cell(s) are incomplete.",
      i = "Resubmit with the SAME run id -- finished tasks exit immediately and interrupted ones resume from their last flush:",
      " " = "RUN_ID={opt$run_id} bash slurm/launch_subset.sh {paste(unique(rep_df$wave[!rep_df$complete]), collapse = ' ')}"
    ))
  } else if (all(grid$on_disk)) {
    cli::cli_alert_success(
      "Every cell of the spec grid is on disk. Next: {.code Rscript --no-init-file simulations/honest_inference_sparsity_failure/code/analyze.R}"
    )
  } else {
    cli::cli_alert_warning(
      "This run's cells are complete, but {sum(!grid$on_disk)} cell(s) of the spec grid are still absent: {paste(basename(grid$path[!grid$on_disk]), collapse = ', ')}"
    )
  }
  cli::cli_h1("analyze.R executes the spec's Decision rule: if ci_95_anchor coverage is below 1-alpha at ANY eps, escalate to Oracle before interpreting anything else.")
}
