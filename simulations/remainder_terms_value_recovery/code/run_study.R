# ---------------------------------------------------------------------------
# Driver: empirical size of the remainder terms T_a and T_epsilon in the
# value-recovery (L2-closeness) route to EIF-ATT inference.
#
# Usage (from the study directory):
#   Rscript code/run_study.R                 # full study, defaults below
#   RTVR_REPS=5 Rscript code/run_study.R     # pilot
#
# Every parameter has a working default, so sourcing this file interactively
# yields a runnable configuration.
#
# Run time: dominated by the exact optimal-tree fits at the largest n, Lbar
# cells (a single depth-6 outcome fit at n = 20000, Lbar = 20 is ~45 s).
#
# Lbar = 30 is DELIBERATELY DROPPED from the default grid (2026-09-08): the
# 2026-09-08 run crashed at Lbar = 30 calibration (n = 2000, the cheapest
# cell) with "Model limit exceeded" -- Configuration::model_limit (10,000)
# was hit by tied-optimal-score trees during single-tree extraction, not by
# a Rashomon-set request. See
# quality_reports/2026-09-08_optimaltrees-p5-feasibility-status.md and
# Oracle's tie-explosion analysis. Re-add Lbar = 30 via RTVR_L once the
# extractor's deterministic-single-optimum fix (models.hpp) lands and the
# runtime curve from Lbar in {5,10,20} says the architecture is still worth
# extending that far.
# ---------------------------------------------------------------------------

library(optimaltrees)
library(future)
library(furrr)
library(readr)
library(fs)
library(dplyr)

# ---- configuration (env-overridable, all with working defaults) ------------
REPS <- as.integer(Sys.getenv("RTVR_REPS", "500"))
WORKERS <- as.integer(Sys.getenv("RTVR_WORKERS", "8"))
N_GRID <- as.integer(strsplit(Sys.getenv("RTVR_N", "2000,5000,10000,20000"), ",")[[1]])
L_GRID <- as.integer(strsplit(Sys.getenv("RTVR_L", "5,10,20"), ",")[[1]])  # 30 dropped, see header note
FIT_TIME_LIMIT <- as.numeric(Sys.getenv("RTVR_FIT_TIME_LIMIT", "600"))
SEED <- as.integer(Sys.getenv("RTVR_SEED", "20260908"))
RESULT_DIR <- Sys.getenv("RTVR_RESULT_DIR", "results")

stopifnot(REPS >= 1L, WORKERS >= 1L, length(N_GRID) >= 1L, length(L_GRID) >= 1L,
          all(is.finite(N_GRID)), all(is.finite(L_GRID)))

source(fs::path("code", "dgp.R"))
source(fs::path("code", "terms.R"))
source(fs::path("code", "bounds.R"))

fs::dir_create(RESULT_DIR)
fs::dir_create(fs::path(RESULT_DIR, "cells"))

cat(sprintf(paste0("config | reps=%d workers=%d n={%s} L={%s}\n",
                   "        | p=%d thresholds/coord=%d binary features=%d grid cells=%d\n"),
            REPS, WORKERS, paste(N_GRID, collapse = ","), paste(L_GRID, collapse = ","),
            P_COVARIATES, N_THRESH, P_COVARIATES * N_THRESH, N_GRID_CELL))

TRUTHS <- setNames(lapply(L_GRID, make_truth), as.character(L_GRID))
CELLS <- expand.grid(n = N_GRID, L = L_GRID, KEEP.OUT.ATTRS = FALSE)
CELLS$depth <- vapply(as.character(CELLS$L), function(k) TRUTHS[[k]]$depth_mu, numeric(1))


# ---- step 1: calibrate lambda once per cell -------------------------------
# lambda is a per-cell nuisance of the fitting algorithm, not an object of
# study; calibrating once and reusing keeps the main loop to ~2 fits per
# replication instead of the ~25 a per-replication bisection would need.
calibration_path <- fs::path(RESULT_DIR, "lambda_calibration.rds")

calibrate_cells <- function() {
  out <- vector("list", nrow(CELLS))
  for (i in seq_len(nrow(CELLS))) {
    n <- CELLS$n[i]; L <- CELLS$L[i]
    truth <- TRUTHS[[as.character(L)]]
    set.seed(SEED + 1000L * i)
    dat <- simulate_data(n, truth)
    ctrl <- dat$A == 0L
    t0 <- Sys.time()
    cal_e <- calibrate_lambda(as.data.frame(dat$Xbin), dat$A, "log_loss",
                              L, truth$depth_e, time_limit = FIT_TIME_LIMIT)
    cal_mu <- calibrate_lambda(as.data.frame(dat$Xbin[ctrl, , drop = FALSE]),
                               dat$Y[ctrl], "squared_error",
                               L, truth$depth_mu, time_limit = FIT_TIME_LIMIT)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("calibrate | n=%6d L=%2d depth=%d | lam_e=%.3e (%2d leaves) lam_mu=%.3e (%2d leaves) | %.0fs\n",
                n, L, truth$depth_mu, cal_e$lambda, cal_e$n_leaves,
                cal_mu$lambda, cal_mu$n_leaves, el))
    out[[i]] <- data.frame(n = n, L = L, depth = truth$depth_mu,
                           lambda_e = cal_e$lambda, leaves_e_cal = cal_e$n_leaves,
                           binding_e = cal_e$binding,
                           lambda_mu = cal_mu$lambda, leaves_mu_cal = cal_mu$n_leaves,
                           binding_mu = cal_mu$binding, calib_secs = el)
  }
  do.call(rbind, out)
}

if (fs::file_exists(calibration_path)) {
  CALIB <- readr::read_rds(calibration_path)
  cat("calibration | reusing", calibration_path, "\n")
} else {
  CALIB <- calibrate_cells()
  readr::write_rds(CALIB, calibration_path)
}
stopifnot(nrow(CALIB) == nrow(CELLS))
if (!all(CALIB$binding_e & CALIB$binding_mu)) {
  cat("NOTE | leaf budget not attainable in these cells (data do not support Lbar leaves):\n")
  print(CALIB[!(CALIB$binding_e & CALIB$binding_mu), c("n", "L", "leaves_e_cal", "leaves_mu_cal")])
}


# ---- step 2: replications, cell by cell, resumable ------------------------
future::plan(future::multisession, workers = WORKERS)
on.exit(future::plan(future::sequential), add = TRUE)

run_cell <- function(i) {
  n <- CELLS$n[i]; L <- CELLS$L[i]
  truth <- TRUTHS[[as.character(L)]]
  cal <- CALIB[CALIB$n == n & CALIB$L == L, ]
  stopifnot(nrow(cal) == 1L)
  path <- fs::path(RESULT_DIR, "cells", sprintf("n%06d_L%02d.rds", n, L))
  if (fs::file_exists(path)) {
    cached <- readr::read_rds(path)
    if (nrow(cached) >= REPS) {
      cat(sprintf("cell      | n=%6d L=%2d | cached (%d reps)\n", n, L, nrow(cached)))
      return(cached)
    }
  }
  t0 <- Sys.time()
  res <- furrr::future_map_dfr(
    seq_len(REPS),
    function(r) run_replication(n, truth, cal$lambda_e, cal$lambda_mu,
                                time_limit = FIT_TIME_LIMIT),
    .options = furrr::furrr_options(seed = SEED + 7L * i,
                                    globals = TRUE,
                                    packages = c("optimaltrees")),
    .progress = FALSE
  )
  res$rep <- seq_len(nrow(res))
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  readr::write_rds(res, path)
  cat(sprintf("cell      | n=%6d L=%2d | %d reps in %.0fs (%.2fs/rep) | max identity gap %.1e\n",
              n, L, nrow(res), el, el / nrow(res), max(res$identity_gap)))
  res
}

# cheap cells first, so a truncated run still yields a usable picture
order_by_cost <- order(CELLS$depth, CELLS$n)
results <- do.call(rbind, lapply(order_by_cost, run_cell))

# The exact refit identity T = T_eps - T_a is algebra, not an approximation: if
# it fails, the leaf values are not within-leaf control means and every reported
# number is suspect.
stopifnot("refit identity T = T_eps - T_a violated" =
            max(results$identity_gap) < IDENTITY_TOL)

readr::write_rds(results, fs::path(RESULT_DIR, "replications.rds"))
readr::write_rds(bound_table(CELLS, n_feat = P_COVARIATES * N_THRESH),
                 fs::path(RESULT_DIR, "bounds.rds"))
cat(sprintf("done | %d replications | max identity gap %.1e\n",
            nrow(results), max(results$identity_gap)))
