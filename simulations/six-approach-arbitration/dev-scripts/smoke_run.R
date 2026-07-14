# =============================================================================
# smoke_run.R -- fast LOCAL sanity run for the arbitration study (not the cluster)
# =============================================================================
# Runs a tiny grid (small n, few reps) covering all seven estimators PLUS a couple
# of escalate=TRUE cells, and prints a per-cell summary: coverage of the reported CI,
# cross-fit twin coverage, mean |delta|, mean realized rashomon_c, and convergence.
#
# Purpose: verify end-to-end wiring before submitting the full SLURM study --
#   (1) all 7 methods produce finite estimates,
#   (2) the averaged methods (dt_averaged, msplit_averaged) return the honest bias-aware
#       CI + cross-fit twin fields,
#   (3) escalate=TRUE cells report rashomon_c >= 1 (escalation is LIVE).
#
# Run time: a few minutes on a laptop (K, M, n and reps are intentionally small).
# Usage:  Rscript simulations/six-approach-arbitration/dev-scripts/smoke_run.R
# =============================================================================

suppressPackageStartupMessages({
  library(doubletree)
})

# Resolve the study root robustly: from the --file arg (Rscript), else common CWDs.
.find_study_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if (length(fa) == 1) {
    here <- dirname(normalizePath(sub("^--file=", "", fa)))
    cand <- normalizePath(file.path(here, ".."), mustWork = FALSE)
    if (dir.exists(file.path(cand, "R"))) return(cand)
  }
  for (cand in c("simulations/six-approach-arbitration",
                 "doubletree/simulations/six-approach-arbitration", ".")) {
    if (dir.exists(file.path(cand, "R"))) return(normalizePath(cand))
  }
  stop("Could not locate the six-approach-arbitration study root.", call. = FALSE)
}
study_root <- .find_study_root()

source(file.path(study_root, "config", "grid.R"))
source(file.path(study_root, "R", "dgp.R"))
source(file.path(study_root, "R", "estimators.R"))
source(file.path(study_root, "R", "run_one.R"))

set.seed(20260714)

# Small, fast smoke grid: one small n, one easy + one moderate DGP, all methods, plus
# escalate=TRUE for the intersection methods. Not for inference -- wiring check only.
SMOKE <- rbind(
  expand.grid(n = 300L, dgp = c("simple", "moderate"),
              method = c("full", "crossfit", "doubletree", "dt_averaged",
                         "msplit", "msplit_averaged", "single_tree"),
              escalate = FALSE, stringsAsFactors = FALSE),
  expand.grid(n = 300L, dgp = "moderate",
              method = c("doubletree", "dt_averaged", "single_tree"),
              escalate = TRUE, stringsAsFactors = FALSE)
)
SMOKE$config_id <- seq_len(nrow(SMOKE))

N_REPS <- 5L

rows <- list()
for (cid in SMOKE$config_id) {
  cfg <- SMOKE[cid, , drop = FALSE]
  for (rep_id in seq_len(N_REPS)) {
    unit_row <- cfg
    unit_row$rep_id <- rep_id
    unit_row$unit <- (cid - 1L) * N_REPS + rep_id
    unit_row$seed <- .unit_seed(cid, rep_id, total_reps = N_REPS, base_seed = 20260714L)
    rows[[length(rows) + 1L]] <- tryCatch(
      run_one(unit_row),
      error = function(e) {
        message(sprintf("  [cell %d rep %d] ERROR: %s", cid, rep_id, conditionMessage(e)))
        NULL
      }
    )
  }
  cat(sprintf("cell %2d/%2d done: %-16s dgp=%-9s escalate=%s\n",
              cid, nrow(SMOKE), cfg$method, cfg$dgp, cfg$escalate))
}

res <- do.call(rbind, rows)

# Per-cell summary.
agg <- function(df) {
  data.frame(
    method    = df$method[1],
    dgp       = df$dgp[1],
    escalate  = df$escalate[1],
    reps      = nrow(df),
    n_finite  = sum(is.finite(df$estimate)),
    coverage  = mean(df$covered, na.rm = TRUE),
    cov_twin  = mean(df$covered_crossfit, na.rm = TRUE),
    mean_abs_delta = mean(abs(df$delta), na.rm = TRUE),
    mean_c_e  = mean(df$rashomon_c_e, na.rm = TRUE),
    conv_rate = mean(df$converged, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
# Group by cell (method + dgp + escalate) and summarize.
key <- paste(res$method, res$dgp, res$escalate, sep = "|")
by_cell <- do.call(rbind, lapply(split(res, key), agg))
rownames(by_cell) <- NULL

cat("\n================ SMOKE SUMMARY (", N_REPS, "reps/cell ) ================\n", sep = "")
print(by_cell, digits = 3, row.names = FALSE)

# Sanity assertions (soft -- warn, don't stop, since a few reps can be unlucky).
esc <- res[res$escalate & is.finite(res$rashomon_c_e), ]
if (nrow(esc) > 0 && all(esc$rashomon_c_e >= 1 - 1e-8)) {
  cat("\nOK: escalate=TRUE cells report rashomon_c_e >= 1 (escalation is live).\n")
} else if (nrow(esc) > 0) {
  cat("\nNOTE: some escalate cells have rashomon_c_e < 1 (m0-scale reporting or fallback).\n")
}
avg <- res[res$method %in% c("dt_averaged", "msplit_averaged") & is.finite(res$delta), ]
if (nrow(avg) > 0) {
  cat("OK: averaged methods returned finite delta (cross-fit twin wired).\n")
}
cat("\nSmoke run complete.\n")
