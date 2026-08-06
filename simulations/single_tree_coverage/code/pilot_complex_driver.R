# ============================================================
# pilot_complex_driver.R  (TEMP pilot -- delete after analysis)
# Study: 2026-08-04_single-tree-coverage  (doubletree)
#
# Drives the R=100 complex-DGP pilot as a sequence of FRESH R subprocesses,
# one per (n, rep-chunk). Rationale: the native TreeFARMS/Rcpp RSS leak is
# unreclaimable by R gc() (heap stayed flat at 66.5 MB while RSS grew to a
# SIGKILL at call ~125). Ending the worker process after each chunk returns
# that memory to the OS. Each chunk checkpoints to disk, so an OOM kill costs
# at most CHUNK reps of one n, and the driver is RESUMABLE (skips done chunks).
#
# Chunk size shrinks with n (high-n cells cost more RSS): n<=2000 -> 25 reps,
# n=5000 -> 10, n=10000 -> 1 (ONE rep per subprocess).
#
# 2026-08-05 CORRECTION: chunk=5 at n=10000 was NOT safe. That is ~50 leaky native
# calls per process (5 reps x K=5 folds x 2 nuisances); the 2026-08-04 run was
# SIGKILLed (status 137) on n=10000 chunks r1-5, r11-15, r26-30, and then the OOM
# killer escalated to the DRIVER (chunks 1-28 logged, chunk 29 never logged at all)
# and took the 16 GB host down hard. Two fixes: (a) n=10000 runs 1 rep per
# subprocess, capping a worker at 10 calls; (b) the driver enforces its own RSS
# budget per child so an over-budget worker dies ALONE and the driver survives.
#
# Run (from repo root):
#   Rscript doubletree/simulations/single_tree_coverage/code/pilot_complex_driver.R
# Outputs: results/chunks_complexR100/chunk_n<>_r<lo>-<hi>.rds (checkpoints)
#          results/results_<runid>_complexR100.rds (merged; analyze.R-ready)
# ============================================================

suppressPackageStartupMessages({ library(purrr); library(readr); library(fs); library(dplyr) })

STUDY_DIR <- fs::path("doubletree", "simulations", "single_tree_coverage")
source(fs::path(STUDY_DIR, "code", "oracle_theta_star.R"))   # THETA0, N_GRID, DGP_NAMES, load_oracle
source(fs::path(STUDY_DIR, "code", "harness.R"))             # make_seed, skipped_rep_row
source(fs::path(STUDY_DIR, "code", "memguard.R"))            # mg_* watchdog + host gate
WORKER <- fs::path(STUDY_DIR, "code", "pilot_complex_chunk_worker.R")

# The FULL study grid, captured BEFORE any test-time override. make_seed() spaces seeds
# by position in the n-grid, so the skipped-rep rows built at merge time must be seeded
# against the full grid or their recorded seeds would not match what the workers used.
N_GRID_ALL <- N_GRID

# GS_TAG / GS_N_GRID exist so a smoke run can exercise this driver END-TO-END without
# touching live checkpoints or the merged results (GS_TAG=smoke GS_N_GRID=500 GS_REPS=2).
TAG <- Sys.getenv("GS_TAG", unset = "complexR100")
ng  <- Sys.getenv("GS_N_GRID", unset = "")
if (nzchar(ng)) N_GRID <- as.integer(strsplit(ng, ",", fixed = TRUE)[[1]])

# GS_RETRY_OOM=1 clears the permanent .oom blocks and retries those reps. Use only after
# raising MG$hard_cap or moving to a larger host -- otherwise they fail identically.
RETRY_OOM <- identical(Sys.getenv("GS_RETRY_OOM", unset = "0"), "1")

CHUNK_DIR <- fs::path(STUDY_DIR, "results", paste0("chunks_", TAG))
LOG_DIR <- fs::path(STUDY_DIR, "logs")
fs::dir_create(CHUNK_DIR)

RUNID <- format(Sys.time(), "%Y%m%d_%H%M%S")
CHILD_LOG_DIR <- fs::path(LOG_DIR, paste0("chunks_", RUNID))
fs::dir_create(CHILD_LOG_DIR)
mg_init(STUDY_DIR, fs::path(LOG_DIR, paste0("events_", RUNID, ".tsv")))

REPS <- as.integer(Sys.getenv("GS_REPS", unset = "100"))
chunk_for_n <- function(n) if (n <= 2000) 25L else if (n <= 5000) 10L else 1L

# Build the chunk plan: (n, rep_lo, rep_hi, outfile)
plan <- purrr::map(N_GRID, function(n) {
  cs <- chunk_for_n(n)
  los <- seq(1L, REPS, by = cs)
  purrr::map(los, function(lo) {
    hi <- min(lo + cs - 1L, REPS)
    tibble::tibble(n = n, rep_lo = lo, rep_hi = hi,
                   out = fs::path(CHUNK_DIR, sprintf("chunk_n%d_r%d-%d.rds", n, lo, hi)))
  }) |> purrr::list_rbind()
}) |> purrr::list_rbind()

cli::cli_inform("Chunk plan: {nrow(plan)} chunks, REPS={REPS}, dgp=complex, n={.val {N_GRID}}.")
cli::cli_inform("Memory guard: cap {round(.mg$hard_cap/GB, 2)} GB/worker{if (.mg$degraded) ' (DEGRADED: RSS fallback)' else ''}; events -> {.path {fs::path(LOG_DIR, paste0('events_', RUNID, '.tsv'))}}")

mg_event("DRIVER_START", runid = RUNID, pid = Sys.getpid(), chunks = nrow(plan),
         reps = REPS, hard_cap_gb = round(.mg$hard_cap / GB, 2),
         degraded = .mg$degraded)

t0 <- Sys.time()
consec_abort <- 0L
tally <- character(0)
for (i in seq_len(nrow(plan))) {
  p <- plan[i, ]
  label <- sprintf("n%d_r%d-%d", p$n, p$rep_lo, p$rep_hi)

  # RESUMABLE: sentinel, or a legacy .rds that actually deserialises (integrity test).
  if (mg_chunk_done(p$out)) {
    cli::cli_inform("  [{i}/{nrow(plan)}] skip (done) {fs::path_file(p$out)}")
    next
  }

  # PERMANENTLY BLOCKED: this chunk already exceeded the memory cap. Retrying costs a
  # full worker launch and fails identically, so without this guard the run can never
  # finish. The reps are NOT lost -- the merge below emits an explicit
  # status = "skipped_oom_local" row for each one.
  if (mg_chunk_blocked(p$out)) {
    if (!RETRY_OOM) {
      cli::cli_inform("  [{i}/{nrow(plan)}] skip (oom-blocked) {fs::path_file(p$out)} -- GS_RETRY_OOM=1 to retry")
      next
    }
    fs::file_delete(paste0(p$out, ".oom"))   # forced retry: clear the block first
    mg_event("OOM_BLOCK_CLEARED", chunk = label)
  }

  # Never launch onto an already-critical host.
  if (!mg_wait_for_host()) {
    mg_event("DRIVER_BAIL", chunk = label, reason = "host_never_cleared")
    cli::cli_warn("Host stayed under critical memory pressure; stopping before {label}.")
    break
  }

  res <- mg_run_guarded(
    worker   = WORKER,
    args     = c(as.character(p$n), as.character(p$rep_lo), as.character(p$rep_hi), p$out),
    out_log  = fs::path(CHILD_LOG_DIR, sprintf("%03d_%s.out", i, label)),
    err_log  = fs::path(CHILD_LOG_DIR, sprintf("%03d_%s.err", i, label)),
    sentinel = paste0(p$out, ".done"),
    label    = label
  )
  tally <- c(tally, res$status)
  el <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

  if (identical(res$status, "CHUNK_OK")) {
    consec_abort <- 0L
    cli::cli_inform("  [{i}/{nrow(plan)}] {.strong done} {label} peak={round(res$peak/GB, 2)}GB ({el} min elapsed)")
  } else {
    consec_abort <- if (isTRUE(res$fatal)) consec_abort + 1L else 0L
    cli::cli_warn("  [{i}/{nrow(plan)}] {.strong {res$status}} {label} exit={res$exit} peak={round(res$peak/GB, 2)}GB")
    err_f <- fs::path(CHILD_LOG_DIR, sprintf("%03d_%s.err", i, label))
    if (fs::file_exists(err_f)) {
      tl <- tail(readLines(err_f, warn = FALSE), 4)
      if (length(tl)) cat(paste0("      | ", tl, collapse = "\n"), "\n", sep = "")
    }
    # Repeated host-pressure aborts mean the HOST is the problem, not this cell.
    # Keep going after a single trip, but never grind a stressed machine into the ground.
    if (consec_abort >= MG$bail_aborts) {
      mg_event("DRIVER_BAIL", reason = "consecutive_host_aborts", n = consec_abort)
      cli::cli_warn("{consec_abort} consecutive host-pressure aborts; stopping to protect the host.")
      break
    }
  }
}

if (length(tally)) {
  tb <- sort(table(tally), decreasing = TRUE)
  cli::cli_inform(c("Chunk outcomes this run:",
                    stats::setNames(paste0(names(tb), ": ", as.integer(tb)), rep("*", length(tb)))))
  if ("OOM_KILLED_BY_OS" %in% names(tb)) {
    cli::cli_warn(c("!" = "{tb[['OOM_KILLED_BY_OS']]} chunk{?s} were killed by the OS with NO watchdog trip.",
                    "i" = "The guard was too slow: lower MG$hard_cap or MG$poll_s in memguard.R."))
  }
}

# --- merge all checkpoints into one analyze.R-ready results file -------------
chunk_files <- fs::dir_ls(CHUNK_DIR, regexp = "chunk_n\\d+_r\\d+-\\d+\\.rds$")
if (length(chunk_files) == 0) stop("No chunk checkpoints produced.", call. = FALSE)
results_raw <- purrr::map(chunk_files, readr::read_rds) |> purrr::list_rbind()

# Dropping n=10000 from chunk=5 to chunk=1 changes checkpoint FILENAMES but leaves
# the old files on disk, so the glob above can read the same rep twice (e.g. both
# chunk_n10000_r6-10.rds and chunk_n10000_r6-6.rds). Per-rep seeds are deterministic
# in (dgp, n, rep_id) via make_seed(), so duplicate rows are IDENTICAL and dedup is
# safe -- but warn, so a real schema/seed mismatch is never silently absorbed.
results <- results_raw |>
  dplyr::distinct(dgp, n, rep_id, .keep_all = TRUE) |>
  dplyr::arrange(n, rep_id)
n_dup <- nrow(results_raw) - nrow(results)
if (n_dup > 0) {
  cli::cli_warn(c(
    "Dropped {n_dup} duplicate rep row{?s} while merging.",
    "i" = "Expected after a chunk-size change (stale coarse checkpoints still on disk)."
  ))
}

# --- account for EVERY planned rep, so the denominator stays honest -----------
# A SIGKILLed worker writes no file, so its reps would simply be ABSENT from the frame
# rather than visible as failures. That is NOT missing-at-random: the reps that blow the
# memory cap are the ones with the largest Rashomon sets -- flattest objective, smallest
# structural margin -- i.e. exactly the reps most likely to fail structure recovery.
# Dropping them silently biases coverage UP and recovery UP. So every planned
# (n, rep_id) with no data gets an explicit row, with the two reasons distinguishable.
expected <- purrr::pmap(plan, function(n, rep_lo, rep_hi, out) {
  tibble::tibble(n = n, rep_id = seq.int(rep_lo, rep_hi), oom = mg_chunk_blocked(out))
}) |> purrr::list_rbind() |> dplyr::distinct(n, rep_id, .keep_all = TRUE)

missing <- dplyr::anti_join(expected, results, by = c("n", "rep_id"))
if (nrow(missing) > 0) {
  oracle_tbl <- load_oracle(STUDY_DIR)
  ts_for_n <- function(n) {
    v <- oracle_tbl$theta_star[oracle_tbl$dgp == "complex" & oracle_tbl$n == n]
    if (length(v) == 1L) v else NA_real_
  }
  filler <- purrr::pmap(missing, function(n, rep_id, oom) {
    skipped_rep_row(
      "complex", n, rep_id, ts_for_n(n),
      status = if (oom) "skipped_oom_local" else "missing_not_run",
      msg = if (oom) {
        sprintf("exceeded the %.1f GB per-worker memory cap; not evaluable on this host",
                .mg$hard_cap / GB)
      } else {
        "planned but never completed (run interrupted before this chunk)"
      },
      n_grid = N_GRID_ALL
    )
  }) |> purrr::list_rbind()
  results <- dplyr::bind_rows(results, filler) |> dplyr::arrange(n, rep_id)
  cli::cli_warn(c(
    "!" = "{nrow(missing)} of {nrow(expected)} planned rep{?s} produced no data; recorded explicitly.",
    "i" = "oom-capped: {sum(missing$oom)} | never run: {sum(!missing$oom)}.",
    "i" = "These rows carry status != {.val run} and must STAY in the denominator -- OOM failures track large Rashomon sets, so they are not missing-at-random."
  ))
}

out_path <- fs::path(STUDY_DIR, "results", paste0("results_", RUNID, "_", TAG, ".rds"))
attr(results, "run_meta") <- list(runid = RUNID, reps = REPS, dgp = "complex",
  sequential = TRUE, chunked_subprocess = TRUE, n_rows = nrow(results),
  reps_planned = nrow(expected),
  chunk_sizes = vapply(N_GRID, chunk_for_n, integer(1)),
  mem_hard_cap_gb = round(.mg$hard_cap / GB, 2), mem_metric_degraded = .mg$degraded,
  chunk_status = if (length(tally)) as.list(table(tally)) else list(),
  # Rep-level denominator: analyze.R must treat status != "run" as non-recovery, not
  # as absent data. Recorded here so the split is auditable from the saved object.
  rep_status = as.list(table(results$status)),
  elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins")), timestamp = Sys.time())
readr::write_rds(results, out_path)

mg_event("DRIVER_DONE", runid = RUNID, rows = nrow(results),
         chunks = length(chunk_files),
         elapsed_min = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))

cli::cli_inform(c(
  "v" = "Merged {length(chunk_files)} chunks -> {nrow(results)} rows in {round(as.numeric(difftime(Sys.time(), t0, units='mins')),1)} min.",
  "i" = "Saved {.path {out_path}}.",
  "i" = "Next: Rscript {fs::path(STUDY_DIR, 'code', 'analyze.R')} {out_path}"
))
