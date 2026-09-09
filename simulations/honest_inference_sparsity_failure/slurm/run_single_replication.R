#!/usr/bin/env Rscript
# ============================================================
# slurm/run_single_replication.R
# Study: honest_inference_sparsity_failure  (doubletree, S2)
# Spec:  quality_reports/specs/2026-08-21_honest-inference-sparsity-failure.md
#
# Run:  Rscript simulations/honest_inference_sparsity_failure/slurm/run_single_replication.R \
#         --wave n8000 --task-id 7 --pkg-root /path/to/doubletree --scratch-dir /path/to/scratch
#
# ONE SLURM ARRAY TASK'S WORTH OF WORK: one contiguous block of replications of one
# (dgp, n) cell, as mapped by slurm/units.R. Writes ONE shard to the run's scratch
# directory; slurm/combine_results.R reassembles shards into the cell_*.rds
# checkpoints that code/analyze.R already reads.
#
# The name is the skill's standard name for this file and is kept for consistency
# with the other studies' slurm/ directories, but the unit is a BLOCK of
# replications, not a single replication -- slurm/units.R's header explains why.
#
# ------------------------------------------------------------
# THE PACKAGE-LOADING PATH (the thing that makes this a CLUSTER script)
# ------------------------------------------------------------
#
# code/common.R's own default is pkgload::load_all() against this source tree AND
# against the sibling `../optimaltrees` source tree. That is a DEV-BOX pattern and
# it does not survive a SLURM job: module R has no working pkgload/devtools
# dev-load, and `../optimaltrees` is a relative path that is only correct when both
# repositories happen to be checked out side by side. This script therefore sets
#
#     HIS_USE_INSTALLED = "1"
#
# BEFORE sourcing common.R, which switches its section-0 block to
# `library(optimaltrees); library(doubletree)` against packages that were
# R CMD INSTALLed on the cluster (README_O2.md, prerequisites). Nothing else about
# common.R differs between the two paths: every package call in it is either
# `doubletree::`-qualified or reached through asNamespace(). The path that actually
# ran is stamped into each cell's metadata as `pkg_load_path`.
#
# AND THE PATH IS VERIFIED, NOT ASSUMED -- see assert_installed_load() below. The
# doubletree repository ships a tracked `.Rprofile` that runs
# devtools::load_all("../optimaltrees") on startup, so on a side-by-side checkout an
# R process can end up with a DEV-LOADED optimaltrees even though this script asked
# for the installed one, and every result would then have been produced by
# uninstalled source. That is exactly the silent-wrong-answer failure the gate
# exists to prevent, so it is checked and it aborts.
#
# common.R does much more than load packages -- it defines the DGP dial, the design
# grid, the seed scheme, the three intervals, the fidelity diagnostic and
# run_one_rep()/run_cell() -- so it is sourced in full here rather than being
# partially copied. Copying the study logic into the SLURM script would create a
# second definition of the experiment that could drift from the one the 2026-09-01
# pilot validated and the five n = 500 cells were produced under.
#
# ------------------------------------------------------------
# CRASH SAFETY AND IDEMPOTENCE
# ------------------------------------------------------------
#
#  * If this task's shard already exists, exit 0 immediately. Resubmitting a run-id
#    therefore costs nothing for finished tasks.
#  * Every FLUSH_EVERY replications (slurm/units.R, default 25) the block's
#    completed rows are written atomically to part_<shard>.rds. A wall-time
#    TIMEOUT (SLURM SIGTERM -> SIGKILL) loses at most FLUSH_EVERY replications,
#    and resubmitting the same run-id resumes from the flush rather than from zero.
#  * on_error = "record": a failed replication persists ONE ROW carrying the
#    condition message in `error_message`, per common.R's rep_row_template().
#    Never a plausible-looking NA estimate with no explanation, and never a short
#    row that silently drops out of an rbind.
#
# Reproducibility note: rep_seed() is a pure function of
# (SEED_MASTER, dgp, n, rep), so replication r of a cell is byte-identical whether
# it runs here as part of batch 2 or in a sequential run_sweep(). Batching is
# invisible in the results, and the n = 8000 cells this produces are directly
# comparable to the n = 500 cells produced sequentially on 2026-09-01.
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
})

# ---- 0 thread pinning (must precede any BLAS-touching load) --------------
# Identical to run_sweep.R's child_env. Not cosmetic: it is one leg of this study's
# "no parallelism" memory rule (simulations/docs/MEMORY_SAFE_SIMULATIONS.md and
# MEMORY.md's [LEARN:rashomon-memory]); the others are worker_limit = 1L on both
# estimator calls and parallel = FALSE on the cross-fit call, both in common.R. On a
# cluster it additionally means one task cannot oversubscribe its single allocated
# CPU.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
           NUMEXPR_NUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1")

#' Directory holding this script, from Rscript's own --file argument
#'
#' Used only to DEFAULT --pkg-root, so that an interactive/manual invocation with
#' no flags still works (r-interactive-entry M1). SLURM always passes --pkg-root
#' explicitly.
#'
#' @return Absolute path, or NA_character_ when the script path is unavailable.
script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (!length(hit)) return(NA_character_)
  normalizePath(dirname(sub("^--file=", "", hit[[1L]])), mustWork = FALSE)
}

# <pkg root>/simulations/honest_inference_sparsity_failure/slurm -> up three levels.
default_pkg_root <- function() {
  d <- script_dir()
  if (is.na(d)) return(getwd())
  normalizePath(file.path(d, "..", "..", ".."), mustWork = FALSE)
}

#' Abort unless both packages were loaded from a real installed library
#'
#' The doubletree repository ships a TRACKED `.Rprofile` that calls
#' \code{devtools::load_all("../optimaltrees")}. R reads \code{./.Rprofile} as the
#' user profile whenever \code{R_PROFILE_USER} is unset, so a plain
#' \code{Rscript} launched from the package root -- which is what SLURM does --
#' silently dev-loads optimaltrees from source BEFORE this script runs a single
#' line. common.R's \code{isNamespaceLoaded("optimaltrees")} guard would then find
#' it already present and \code{library(optimaltrees)} would be a no-op on the dev
#' namespace.
#'
#' The consequence is not an error: it is a completed sweep whose numbers came from
#' uninstalled source, which is precisely what
#' \code{HIS_USE_INSTALLED = "1"} exists to rule out. So the load path is verified
#' by asking where each namespace actually came from.
#'
#' run_simulations.slurm and quick_test.sh additionally pass
#' \code{--no-init-file} so the profile never runs; this check is the backstop for
#' every other way of invoking the script.
#'
#' @return Invisibly, a named character vector of the resolved package directories.
assert_installed_load <- function() {
  libs <- normalizePath(.libPaths(), mustWork = FALSE)
  where <- vapply(c("optimaltrees", "doubletree"), function(p) {
    d <- getNamespaceInfo(asNamespace(p), "path")
    normalizePath(d, mustWork = FALSE)
  }, character(1))
  dev <- where[!dirname(where) %in% libs]
  if (length(dev)) {
    stop(sprintf(paste0(
      "package(s) %s were DEV-LOADED from source, not from an installed library:\n",
      "  %s\n",
      "This is almost certainly the tracked .Rprofile at the doubletree package root\n",
      "running devtools::load_all(\"../optimaltrees\") at startup. Re-invoke with\n",
      "  Rscript --no-init-file ...      (what run_simulations.slurm does)\n",
      "or export R_PROFILE_USER=/dev/null. Results from a dev-loaded package are not\n",
      "what --pkg-root's installed packages would produce, so this is fatal, not a warning."),
      paste(names(dev), collapse = ", "),
      paste(sprintf("%s -> %s", names(dev), dev), collapse = "\n  ")),
      call. = FALSE)
  }
  invisible(where)
}

option_list <- list(
  make_option("--wave", type = "character", dest = "wave", default = "n2000",
              help = "Wave id (one array per sample size), e.g. n2000 or n8000 [default %default]"),
  make_option("--task-id", type = "integer", dest = "task_id", default = 1L,
              help = "SLURM array task id within this wave's array, 1-based [default %default]"),
  make_option("--pkg-root", type = "character", dest = "pkg_root",
              default = default_pkg_root(),
              help = "Absolute path to the doubletree package root [default: three levels above this script]"),
  make_option("--scratch-dir", type = "character", dest = "scratch_dir",
              default = file.path(tempdir(), "his_scratch"),
              help = "Directory for this run's result shards [default %default]"),
  make_option("--target-secs", type = "double", dest = "target_secs",
              default = NA_real_,
              help = "Target wall cost per task, in seconds; MUST match the launcher's value or the array map differs [default: HIS_TARGET_SECS, else 900]"),
  make_option("--reps-per-job", type = "integer", dest = "reps_per_job",
              default = NA_integer_,
              help = "Override the cost-derived block size with a fixed number of replications (quick_test.sh uses this)"),
  make_option("--flush-every", type = "integer", dest = "flush_every",
              default = NA_integer_,
              help = "Replications between partial checkpoints [default: FLUSH_EVERY in slurm/units.R]"),
  make_option("--force", action = "store_true", dest = "force", default = FALSE,
              help = "Recompute even if this task's shard already exists")
)

if (!interactive() && sys.nframe() == 0L) {

  opt <- parse_args(OptionParser(option_list = option_list))

  # ---- 1 assert the inputs, once, right here (r-interactive-entry M2) ----
  if (!dir.exists(opt$pkg_root)) {
    stop("--pkg-root does not exist: ", opt$pkg_root, call. = FALSE)
  }
  if (!file.exists(file.path(opt$pkg_root, "DESCRIPTION"))) {
    stop("--pkg-root has no DESCRIPTION; it must be the doubletree package root: ",
         opt$pkg_root, call. = FALSE)
  }
  if (is.na(opt$task_id) || opt$task_id < 1L) {
    stop("--task-id must be a positive integer; got ", opt$task_id, call. = FALSE)
  }

  # ---- 2 load the study, on the INSTALLED-package path ------------------
  # common.R derives every study path from getwd() and requires it to be the
  # package root (it checks DESCRIPTION's Package field), so setwd() first.
  setwd(opt$pkg_root)
  Sys.setenv(HIS_USE_INSTALLED = "1")
  if (!is.na(opt$target_secs)) {
    Sys.setenv(HIS_TARGET_SECS = format(opt$target_secs, scientific = FALSE))
  }

  source(file.path("simulations", "honest_inference_sparsity_failure", "code",
                   "common.R"))
  source(file.path("simulations", "honest_inference_sparsity_failure", "slurm",
                   "units.R"))

  if (!isTRUE(USE_INSTALLED_PKGS)) {
    stop("common.R took the SOURCE loading path despite HIS_USE_INSTALLED = 1; ",
         "its section-0 gate has changed.", call. = FALSE)
  }
  pkg_dirs <- assert_installed_load()

  # ---- 3 resolve this task's unit ---------------------------------------
  target_secs <- his_target_secs()
  ut <- his_unit_table(opt$wave, target_secs = target_secs)
  if (opt$task_id > nrow(ut)) {
    # Array padding beyond the wave's task count. Not an error: exit 0 so SLURM
    # records the task as done rather than failed.
    cat(sprintf("[%s task %d] beyond this wave's %d task(s); nothing to do.\n",
                opt$wave, opt$task_id, nrow(ut)))
    quit(save = "no", status = 0)
  }
  u <- ut[opt$task_id, ]

  # A fixed block size (quick_test.sh) overrides the cost-derived one: run the
  # FIRST `reps_per_job` replications of this unit and label the shard as such, so
  # a smoke-test shard can never be mistaken for a production one -- the shard's
  # `unit$reps` will not match the block the production map assigns, and
  # combine_results.R's completeness check refuses to write a cell from it.
  if (!is.na(opt$reps_per_job)) {
    u$rep_to <- min(u$rep_from + opt$reps_per_job - 1L, u$reps_total)
    u$reps <- u$rep_to - u$rep_from + 1L
  }

  # Flush cadence. Each flush rewrites the whole block-so-far, so a FIXED cadence
  # costs O(reps^2 / cadence) bytes written. Capping the block at ~20 flushes keeps
  # that overhead negligible while still bounding a timeout's loss at 5 % of a
  # task; at this study's block sizes (150 reps at n = 2000, 75 at n = 8000) the
  # FLUSH_EVERY floor is what binds, giving 6 and 3 flushes respectively. An
  # explicit --flush-every always wins, for tests that want per-rep flushes.
  flush_every <- if (is.na(opt$flush_every)) {
    max(FLUSH_EVERY, as.integer(ceiling(u$reps / 20L)))
  } else {
    opt$flush_every
  }
  if (flush_every < 1L) stop("--flush-every must be >= 1.", call. = FALSE)

  dir.create(opt$scratch_dir, recursive = TRUE, showWarnings = FALSE)
  shard <- file.path(opt$scratch_dir,
                     batch_filename(u$dgp, u$n, u$reps_total, u$batch))
  partial <- file.path(opt$scratch_dir,
                       partial_filename(u$dgp, u$n, u$reps_total, u$batch))

  cat(sprintf(
    paste0("[%s task %d/%d] dgp=%s eps=%.2f n=%d reps %d-%d of %d (batch %d/%d), ",
           "projected %.1f s, target %.0f s, flush every %d\n"),
    u$wave, opt$task_id, nrow(ut), u$dgp, u$eps, u$n, u$rep_from, u$rep_to,
    u$reps_total, u$batch, u$n_batches, u$est_secs, target_secs, flush_every))
  cat(sprintf("[%s task %d] host=%s R=%s doubletree=%s (%s) optimaltrees=%s (%s)\n",
              u$wave, opt$task_id, Sys.info()[["nodename"]], getRversion(),
              utils::packageVersion("doubletree"), pkg_dirs[["doubletree"]],
              utils::packageVersion("optimaltrees"), pkg_dirs[["optimaltrees"]]))

  # ---- 4 idempotent skip -------------------------------------------------
  if (file.exists(shard) && !opt$force) {
    cat(sprintf("[%s task %d] shard already exists (%s); skipping.\n",
                u$wave, opt$task_id, basename(shard)))
    quit(save = "no", status = 0)
  }

  # ---- 5 resume from the last partial flush -----------------------------
  # The partial holds the rows of the replications already finished IN THIS UNIT.
  # It is trusted only if it belongs to this exact unit -- a stale partial from a
  # differently-sized array (different --target-secs, hence different batch
  # boundaries) must not be silently absorbed into this one's results.
  done_rows <- list()
  next_rep <- u$rep_from
  if (file.exists(partial) && !opt$force) {
    p <- readRDS(partial)
    same_unit <- identical(p$unit$wave, u$wave) &&
      identical(p$unit$dgp, u$dgp) &&
      identical(as.integer(p$unit$n), as.integer(u$n)) &&
      identical(as.integer(p$unit$rep_from), as.integer(u$rep_from)) &&
      identical(as.integer(p$unit$rep_to), as.integer(u$rep_to))
    if (same_unit) {
      done_rows <- list(p$results)
      next_rep <- max(p$results$rep) + 1L
      cat(sprintf("[%s task %d] resuming: %d replication(s) already flushed; next rep = %d.\n",
                  u$wave, opt$task_id, nrow(p$results), next_rep))
    } else {
      cat(sprintf("[%s task %d] partial %s belongs to a different unit; ignoring it.\n",
                  u$wave, opt$task_id, basename(partial)))
    }
  }

  # ---- 6 run the block ---------------------------------------------------
  dgp <- make_dgp(u$dgp)
  lambda_n <- log(u$n) / u$n
  t0 <- Sys.time()

  #' Write an atomic checkpoint (tmp then rename) so a kill mid-write cannot
  #' leave a truncated .rds that a later resume would read as valid.
  write_atomic <- function(obj, path) {
    tmp <- paste0(path, ".tmp")
    saveRDS(obj, tmp)
    file.rename(tmp, path)
  }

  n_since_flush <- 0L
  if (next_rep <= u$rep_to) {
    for (r in seq.int(next_rep, u$rep_to)) {
      done_rows[[length(done_rows) + 1L]] <-
        run_one_rep(dgp, u$n, r, on_error = "record")
      n_since_flush <- n_since_flush + 1L
      if (n_since_flush >= flush_every || r == u$rep_to) {
        acc <- do.call(rbind, done_rows)
        done_rows <- list(acc)   # collapse so the next rbind stays cheap
        write_atomic(list(unit = u, results = acc), partial)
        n_since_flush <- 0L
        gc(full = TRUE, verbose = FALSE)
        cat(sprintf("[%s task %d] rep %d/%d flushed (%.1f s elapsed, cov_anch=%.3f cov_hon=%.3f cov_cf=%.3f)\n",
                    u$wave, opt$task_id, r, u$rep_to,
                    as.numeric(difftime(Sys.time(), t0, units = "secs")),
                    mean(acc$covered_anchor, na.rm = TRUE),
                    mean(acc$covered_honest, na.rm = TRUE),
                    mean(acc$covered_cf, na.rm = TRUE)))
      }
    }
  } else {
    cat(sprintf("[%s task %d] every replication already flushed; assembling shard.\n",
                u$wave, opt$task_id))
  }

  res <- do.call(rbind, done_rows)
  res <- res[order(res$rep), , drop = FALSE]
  rownames(res) <- NULL

  # A shard that does not hold every replication of its block is a bug, not a
  # partial success: combine_results.R counts replications per cell and would
  # otherwise report a silently short cell as complete.
  got <- sort(unique(res$rep))
  want <- seq.int(u$rep_from, u$rep_to)
  if (!identical(as.integer(got), as.integer(want))) {
    stop(sprintf("[%s task %d] shard would hold reps {%s} but the block is %d-%d.",
                 u$wave, opt$task_id, paste(got, collapse = ","),
                 u$rep_from, u$rep_to), call. = FALSE)
  }

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  n_fail <- sum(!is.na(res$error_message))

  # cell_metadata() is common.R's own provenance stamp, computed here on the BLOCK
  # so a shard is self-describing; combine_results.R recomputes it on the assembled
  # cell, which is the copy analyze.R reads. Both come from the same function, so a
  # shard and its cell cannot disagree about the design.
  write_atomic(list(
    unit = u,
    results = res,
    dgp_summary = dgp_population_summary_dial(dgp),
    meta = cell_metadata(dgp, u$n, u$reps, lambda_n, res),
    runtime = list(
      host = Sys.info()[["nodename"]],
      slurm_job = Sys.getenv("SLURM_JOB_ID", NA_character_),
      slurm_array_task = Sys.getenv("SLURM_ARRAY_TASK_ID", NA_character_),
      target_secs = target_secs,
      est_secs = u$est_secs,
      elapsed_secs = elapsed,
      secs_in_estimators = sum(res$secs_full + res$secs_cf, na.rm = TRUE),
      n_failed_reps = n_fail,
      r_version = R.version.string,
      pkg_load_path = if (USE_INSTALLED_PKGS) "installed" else "source",
      doubletree_dir = pkg_dirs[["doubletree"]],
      optimaltrees_dir = pkg_dirs[["optimaltrees"]],
      doubletree_version = as.character(utils::packageVersion("doubletree")),
      optimaltrees_version = as.character(utils::packageVersion("optimaltrees")),
      run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
  ), shard)
  unlink(partial)

  cat(sprintf("[%s task %d] wrote %s -- %d reps, %d failed, %.1f s (projected %.1f s, %.2f s/rep)\n",
              u$wave, opt$task_id, basename(shard), u$reps, n_fail,
              elapsed, u$est_secs, elapsed / max(1L, u$reps)))

  # Distinct error texts once per shard: a per-replication failure rate is
  # recoverable from the rows, but the MESSAGE is what says whether to fix or
  # resubmit.
  em <- unique(stats::na.omit(res$error_message))
  if (length(em)) {
    cat(sprintf("[%s task %d] %d distinct error message(s):\n",
                u$wave, opt$task_id, length(em)))
    for (m in em) cat("  x ", substr(m, 1, 300), "\n", sep = "")
  }
  # Distinct warning texts once per shard, for the same reason.
  wm <- unique(stats::na.omit(res$warning_messages))
  if (length(wm)) {
    cat(sprintf("[%s task %d] %d muffled warning(s), %d distinct text(s):\n",
                u$wave, opt$task_id, sum(res$n_warnings, na.rm = TRUE), length(wm)))
    for (m in wm) cat("  * ", substr(m, 1, 200), "\n", sep = "")
  }
}
