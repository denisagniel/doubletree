# ============================================================
# run_coverage.R
# Study: 2026-08-04_single-tree-coverage  (doubletree)
# Purpose: Full coverage + structure-recovery study. R = 500 replications per
#   (DGP, n) cell over all four DGPs and n in {500, 1000, 2000, 5000, 10000},
#   driving estimate_att_single_tree and recording recovery (converged), theta*_n
#   coverage, and theta_0 coverage. Parallelised with furrr; each worker self-loads
#   the SOURCE package (pkgload::load_all) so results match the pilot exactly.
#
# FEASIBILITY (READ THIS BEFORE RUNNING LOCALLY):
#   The `continuous` DGP is LOCALLY INFEASIBLE at every n: continuous covariates
#   produce many quantile-threshold binary features, so the propensity Rashomon set
#   explodes (~1e5 trees, >10 GB) and OOMs a 16 GB host (guarded probe 2026-08-04;
#   sister study sacct 2026-08-03: 29-32 GB). Those cells are SKIPPED locally
#   (status = "skipped_infeasible_local") and must be run on a high-memory cluster
#   (see setup-cluster-simulations). The binary DGPs (simple/moderate/complex) run
#   locally to n = 10000 (<0.3 GB/call).
#
#   Estimated wall time for the FEASIBLE run (3 DGPs x 5 n x 500): ~7 h serial,
#   ~1 h on 9 workers (see run_pilot.R for the measured basis). Override REPS /
#   WORKERS below (or via env GS_REPS / GS_WORKERS) to shrink a trial run.
#
# Inputs : results/theta_star_oracle.rds (built by oracle_theta_star.R).
# Outputs: results/results_<runid>.rds (all per-rep rows incl. skipped cells).
# Run    : from repo root ->
#   Rscript doubletree/simulations/single_tree_coverage/code/run_coverage.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(fs)
  library(future)
  library(furrr)
  library(progressr)
})

# TWO PACKAGE-LOADING PATHS, ONE GATE (added 2026-09-09 for SLURM deployment):
#
#   DOUBLETREE_USE_INSTALLED unset/"0"  DEV path (default, unchanged behaviour):
#       devtools::load_all() on the source tree. Correct on the dev box.
#
#   DOUBLETREE_USE_INSTALLED="1"        CLUSTER path: library() against packages
#       that were R CMD INSTALLed on the cluster. Module R does not carry a working
#       pkgload/devtools dev-load.
USE_INSTALLED_PKGS <- Sys.getenv("DOUBLETREE_USE_INSTALLED", "0") == "1"
if (USE_INSTALLED_PKGS) {
  suppressPackageStartupMessages({
    library(doubletree)
  })
} else {
  suppressMessages(devtools::load_all("doubletree", quiet = TRUE))
}

REPO_ROOT <- normalizePath(".")     # scripts run from the repo root
STUDY_DIR <- fs::path("doubletree", "simulations", "single_tree_coverage")
source(fs::path(STUDY_DIR, "code", "oracle_theta_star.R"))   # THETA0, N_GRID, load_oracle
source(fs::path(STUDY_DIR, "code", "harness.R"))             # run_one_rep, is_feasible_local, ...

set.seed(BASE_SEED)                 # study base seed (per-rep seeds via make_seed)

# --- run configuration (override via env for trial runs) ---------------------
REPS         <- as.integer(Sys.getenv("GS_REPS", unset = "500"))
default_wrk  <- as.character(max(1L, parallel::detectCores() - 1L))
WORKERS      <- as.integer(Sys.getenv("GS_WORKERS", unset = default_wrk))
# Set GS_RUN_INFEASIBLE=1 ONLY on a high-memory host/cluster to attempt continuous.
RUN_INFEASIBLE <- identical(Sys.getenv("GS_RUN_INFEASIBLE", unset = "0"), "1")

# --- oracle target (cached) --------------------------------------------------
oracle_tbl <- load_oracle(STUDY_DIR)

# --- build the full grid -----------------------------------------------------
grid <- tidyr::expand_grid(dgp = DGP_NAMES, n = N_GRID, rep_id = seq_len(REPS)) |>
  dplyr::left_join(dplyr::select(oracle_tbl, dgp, n, theta_star), by = c("dgp", "n")) |>
  dplyr::mutate(feasible = RUN_INFEASIBLE | purrr::map2_lgl(dgp, n, is_feasible_local))

run_grid  <- dplyr::filter(grid, feasible)
skip_grid <- dplyr::filter(grid, !feasible)

cli::cli_inform(c(
  "i" = "Full run: REPS={REPS}, WORKERS={WORKERS}.",
  "i" = "{nrow(run_grid)} calls to run; {nrow(skip_grid)} skipped (locally infeasible: {.val {sort(unique(skip_grid$dgp))}})."
))

# --- worker: self-load SOURCE package + study code, then run one rep ---------
# Multisession workers are fresh R sessions; they must load the SOURCE package
# (the installed build is stale) exactly once, then reuse it across tasks.
# Uses the same DOUBLETREE_USE_INSTALLED gate as the main script.
worker_run_one <- function(dgp, n, rep_id, theta_star) {
  if (!isTRUE(getOption("gs_stc_loaded"))) {
    use_installed <- Sys.getenv("DOUBLETREE_USE_INSTALLED", "0") == "1"
    if (use_installed) {
      suppressPackageStartupMessages({
        library(doubletree)
      })
    } else {
      suppressMessages(pkgload::load_all(file.path(REPO_ROOT, "doubletree"), quiet = TRUE))
    }
    source(file.path(REPO_ROOT, "doubletree", "simulations", "single_tree_coverage",
                     "code", "oracle_theta_star.R"))
    source(file.path(REPO_ROOT, "doubletree", "simulations", "single_tree_coverage",
                     "code", "harness.R"))
    options(gs_stc_loaded = TRUE)
  }
  run_one_rep(dgp, n, rep_id, theta_star)
}

# --- run (parallel over feasible cells) --------------------------------------
future::plan(future::multisession, workers = WORKERS)
on.exit(future::plan(future::sequential), add = TRUE)
progressr::handlers(global = TRUE)
progressr::handlers("cli")

t0 <- Sys.time()
run_rows <- progressr::with_progress({
  p <- progressr::progressor(steps = nrow(run_grid))
  run_grid |>
    dplyr::select(dgp, n, rep_id, theta_star) |>
    furrr::future_pmap(
      function(dgp, n, rep_id, theta_star) {
        out <- worker_run_one(dgp, n, rep_id, theta_star)
        p()
        out
      },
      .options = furrr::furrr_options(seed = TRUE,
                                      globals = c("REPO_ROOT", "worker_run_one", "p"))
    ) |>
    purrr::list_rbind()
})
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

# --- save (with run metadata) ------------------------------------------------
runid <- format(Sys.time(), "%Y%m%d_%H%M%S")
attr(results, "run_meta") <- list(
  runid = runid, reps = REPS, workers = WORKERS, n_grid = N_GRID,
  dgp_names = DGP_NAMES, run_infeasible = RUN_INFEASIBLE,
  base_seed = BASE_SEED, elapsed_s = elapsed_s,
  n_run = nrow(run_rows), n_skipped = nrow(skip_grid),
  timestamp = Sys.time()
)
fs::dir_create(fs::path(STUDY_DIR, "results"))
out_path <- fs::path(STUDY_DIR, "results", paste0("results_", runid, ".rds"))
readr::write_rds(results, out_path)

cli::cli_inform(c(
  "v" = "Ran {nrow(run_rows)} calls in {round(elapsed_s/3600, 2)} h on {WORKERS} workers.",
  "i" = "Saved {.path {out_path}}.",
  "i" = "Next: Rscript {fs::path(STUDY_DIR, 'code', 'analyze.R')}"
))
