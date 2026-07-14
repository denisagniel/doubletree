# =============================================================================
# verify_msplit_parity.R -- confirm M-split ~ cross-fitting on the ARBITRATION
# DGPs (strengthened propensity, 2026-07-08), using the live harness.
# =============================================================================
# Re-verification of an earlier finding that was measured on the now-deprecated
# six_approach_comparison/ weak-propensity DGPs. Here we drive the arbitration
# study's own run_one() so the DGPs, estimators (real package calls), depth caps,
# and per-(config,rep) seeds are exactly the study's.
#
# Scope (local): methods crossfit & msplit, all 4 DGPs, n=1000, REPS reps.
# The full grid (7 methods x 3 n x 4 dgp x 1000 reps) is the cluster job.
#
# Run from the six-approach-arbitration/ directory.
# =============================================================================

suppressPackageStartupMessages({
  library(optimaltrees)
  devtools::load_all("../..", quiet = TRUE)  # local doubletree R/ (repo root is doubletree/)
})

source("config/grid.R")
source("R/dgp.R")
source("R/estimators.R")
source("R/run_one.R")

N       <- 1000L
METHODS <- c("crossfit", "msplit")
DGPS    <- c("simple", "moderate", "complex", "continuous")
# Binary DGPs are cheap -> 50 reps (decisive for the parity question). The
# continuous STRESS DGP is ~10-15x slower per unit (depth-capped GOSDT on
# continuous X), so use fewer reps: enough for a coverage/SD signal without the
# multi-hour cost. This is the stress regime, not the core parity comparison.
reps_for <- function(dg) if (dg == "continuous") 15L else 50L

# Build unit rows straight from the study's own enumeration so seeds match the
# real run. We take the first REPS rep_ids for each (method, dgp) cell at n=N.
ut <- unit_table()
rows <- list()
for (meth in METHODS) {
  for (dg in DGPS) {
    reps <- reps_for(dg)
    cell <- ut[ut$n == N & ut$dgp == dg & ut$method == meth, ]
    cell <- cell[order(cell$rep_id), ][seq_len(reps), ]
    if (!all(is_feasible(cell))) {
      cat(sprintf("SKIP %-8s %-10s (infeasible cell)\n", meth, dg)); next
    }
    res <- vector("list", nrow(cell))
    for (i in seq_len(nrow(cell))) {
      res[[i]] <- tryCatch(run_one(cell[i, ]),
                           error = function(e) NULL)
    }
    res <- do.call(rbind, res[!vapply(res, is.null, logical(1))])
    ok  <- res[!is.na(res$estimate), ]
    summ <- data.frame(
      method = meth, dgp = dg, n = N,
      reps_ok = nrow(ok), reps_tried = reps,
      bias = mean(ok$error), sd_theta = sd(ok$estimate),
      coverage = mean(ok$covered), stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1]] <- summ
    cat(sprintf("%-8s %-10s bias=%+.4f SD=%.4f cov=%.2f (ok %d/%d)\n",
                meth, dg, summ$bias, summ$sd_theta, summ$coverage,
                summ$reps_ok, reps)); flush.console()
  }
}

out <- do.call(rbind, rows)
readr::write_csv(out, "dev-scripts/msplit_parity_arbitration.csv")
cat("\nWrote dev-scripts/msplit_parity_arbitration.csv\n")

# Quick verdict: msplit vs crossfit SD ratio per DGP (want ~1).
cat("\n=== SD ratio msplit/crossfit per DGP ===\n")
for (dg in DGPS) {
  cf <- out[out$method == "crossfit" & out$dgp == dg, ]
  ms <- out[out$method == "msplit"   & out$dgp == dg, ]
  if (nrow(cf) && nrow(ms)) {
    cat(sprintf("  %-10s %.2f  (cf cov %.2f / ms cov %.2f)\n",
                dg, ms$sd_theta / cf$sd_theta, cf$coverage, ms$coverage))
  }
}
