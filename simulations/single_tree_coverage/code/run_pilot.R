# ============================================================
# run_pilot.R
# Study: 2026-08-04_single-tree-coverage  (doubletree)
# Purpose: R = 10 end-to-end smoke test of the coverage harness at n in {500, 2000}
#   over all four DGPs. Verifies that `converged`/`delta` populate, the oracle
#   theta*_n loads, and coverage is computed -- BEFORE scaling to the full study.
#   The `continuous` DGP is locally infeasible (Rashomon OOM on 16 GB; see harness.R)
#   so it is recorded as skipped, not run. Also prints a full-run cost estimate.
# Inputs : results/theta_star_oracle.rds (built by oracle_theta_star.R).
# Outputs: results/pilot_<runid>.rds (per-rep rows + a small summary attribute).
# Run    : from repo root ->
#   Rscript doubletree/simulations/single_tree_coverage/code/run_pilot.R
# Run time: ~5-10 min (60 estimator calls, serial).
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(fs)
})
suppressMessages(devtools::load_all("doubletree", quiet = TRUE))

STUDY_DIR <- fs::path("doubletree", "simulations", "single_tree_coverage")
# oracle_theta_star.R supplies THETA0, N_GRID, DGP_NAMES, load_oracle() and sources
# dgps.R; its compute-and-print main block is guarded so source() does NOT recompute.
source(fs::path(STUDY_DIR, "code", "oracle_theta_star.R"))
source(fs::path(STUDY_DIR, "code", "harness.R"))

set.seed(BASE_SEED)                       # study base seed (per-rep seeds via make_seed)

PILOT_REPS <- 10L
PILOT_N    <- c(500L, 2000L)

# --- oracle target (cached) --------------------------------------------------
oracle_tbl <- load_oracle(STUDY_DIR)

# --- build the pilot grid ----------------------------------------------------
grid <- tidyr::expand_grid(dgp = DGP_NAMES, n = PILOT_N, rep_id = seq_len(PILOT_REPS)) |>
  dplyr::left_join(dplyr::select(oracle_tbl, dgp, n, theta_star), by = c("dgp", "n")) |>
  dplyr::mutate(feasible = purrr::map2_lgl(dgp, n, is_feasible_local))

run_grid  <- dplyr::filter(grid, feasible)
skip_grid <- dplyr::filter(grid, !feasible)

cli::cli_inform(c(
  "i" = "Pilot: {nrow(run_grid)} estimator calls to run, {nrow(skip_grid)} skipped (locally infeasible).",
  "i" = "Feasible DGPs: {.val {sort(unique(run_grid$dgp))}}; skipped: {.val {sort(unique(skip_grid$dgp))}}."
))

# --- run ---------------------------------------------------------------------
t0 <- Sys.time()
run_rows <- run_grid |>
  dplyr::select(dgp, n, rep_id, theta_star) |>
  purrr::pmap(run_one_rep, .progress = "pilot") |>
  purrr::list_rbind()
elapsed_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

skip_rows <- if (nrow(skip_grid) > 0) {
  skip_grid |>
    dplyr::select(dgp, n, rep_id, theta_star) |>
    purrr::pmap(skipped_rep_row) |>
    purrr::list_rbind()
} else {
  NULL
}

results <- dplyr::bind_rows(run_rows, skip_rows)

# --- save --------------------------------------------------------------------
runid <- format(Sys.time(), "%Y%m%d_%H%M%S")
fs::dir_create(fs::path(STUDY_DIR, "results"))
out_path <- fs::path(STUDY_DIR, "results", paste0("pilot_", runid, ".rds"))
readr::write_rds(results, out_path)

# --- per-cell recovery + coverage (feasible cells only) ----------------------
cell_summary <- run_rows |>
  dplyr::summarise(
    reps               = dplyr::n(),
    recovery_rate      = mean(converged),
    n_converged        = sum(converged),
    cov_star_honest    = mean(covers_star_honest[converged]),
    cov_theta0_honest  = mean(covers_theta0_honest[converged]),
    mean_delta         = mean(delta[converged]),
    mean_theta_single  = mean(theta_single[converged]),
    .by = c(dgp, n)
  ) |>
  dplyr::arrange(dgp, n)

cat("\n=== PILOT per-cell summary (feasible cells; R =", PILOT_REPS, ") ===\n")
cell_summary |>
  dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, 4))) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\n=== theta*_n vs theta_0 at pilot cells (oracle) ===\n")
oracle_tbl |>
  dplyr::filter(n %in% PILOT_N) |>
  dplyr::select(dgp, n, theta_star, theta0, theta0_realized, gap_vs_realized) |>
  dplyr::mutate(dplyr::across(c(theta_star, theta0_realized, gap_vs_realized),
                              \(x) round(x, 4))) |>
  as.data.frame() |>
  print(row.names = FALSE)

# --- full-run cost estimate --------------------------------------------------
# Per-call seconds by (dgp, n) measured in the pilot (n in {500, 2000}); the pilot
# cannot measure n in {5000, 10000}, so those rows use the guarded-probe timings
# (2026-08-04: simple/moderate n=10000 ~26s; complex n=5000 ~15s, n=10000 ~27s).
# Extrapolate to R = 500 over the FEASIBLE grid (3 binary DGPs x 5 n). Serial
# estimate; divide by worker count for the furrr path in run_coverage.R.
pilot_time_by_cell <- run_rows |>
  dplyr::summarise(mean_sec = mean(elapsed_s), .by = c(dgp, n))

cat("\n=== Per-call timing (pilot cells) ===\n")
pilot_time_by_cell |>
  dplyr::mutate(mean_sec = round(mean_sec, 1)) |>
  as.data.frame() |>
  print(row.names = FALSE)

# Probe-measured s/call at the larger n the pilot does not reach (feasible DGPs).
probe_sec <- tibble::tribble(
  ~n,      ~sec_binary,
  1000L,    6,
  5000L,   15,
  10000L,  27
)
mean_pilot_500  <- mean(pilot_time_by_cell$mean_sec[pilot_time_by_cell$n == 500L])
mean_pilot_2000 <- mean(pilot_time_by_cell$mean_sec[pilot_time_by_cell$n == 2000L])
sec_by_n <- tibble::tibble(
  n   = N_GRID,
  sec = c(mean_pilot_500, probe_sec$sec_binary[probe_sec$n == 1000L],
          mean_pilot_2000, probe_sec$sec_binary[probe_sec$n == 5000L],
          probe_sec$sec_binary[probe_sec$n == 10000L])
)

FULL_REPS       <- 500L
feasible_dgps   <- setdiff(DGP_NAMES, LOCAL_INFEASIBLE_DGPS)
n_feasible_dgps <- length(feasible_dgps)
serial_secs     <- sum(sec_by_n$sec) * n_feasible_dgps * FULL_REPS
full_calls      <- n_feasible_dgps * length(N_GRID) * FULL_REPS
workers         <- max(1L, parallel::detectCores() - 1L)

cat(sprintf("\n=== Timing / full-run estimate ===\n"))
cat(sprintf("Pilot: %d calls in %.1f s (mean %.2f s/call).\n",
            nrow(run_rows), elapsed_s, elapsed_s / nrow(run_rows)))
cat(sprintf("Feasible full run: %d DGPs x %d n x %d reps = %s calls.\n",
            n_feasible_dgps, length(N_GRID), FULL_REPS, format(full_calls, big.mark = ",")))
cat(sprintf("Estimated SERIAL wall time: %.1f h.\n", serial_secs / 3600))
cat(sprintf("Estimated PARALLEL wall time on %d workers: ~%.1f h.\n",
            workers, serial_secs / 3600 / workers))
cat(sprintf("(continuous DGP excluded: locally infeasible; cluster required.)\n"))
cat(sprintf("Saved: %s\n", out_path))

