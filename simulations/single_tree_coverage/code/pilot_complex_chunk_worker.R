# ============================================================
# pilot_complex_chunk_worker.R  (TEMP pilot -- delete after analysis)
# Study: 2026-08-04_single-tree-coverage  (doubletree)
#
# Runs ONE chunk: dgp=complex, a single n, rep_id in [REP_LO, REP_HI].
# Invoked as a FRESH R subprocess by pilot_complex_driver.R. Exiting the
# process after each chunk is what releases the native (TreeFARMS/Rcpp) RSS
# that R's gc() cannot reclaim -- the confirmed cause of the OOM at call ~125
# (R heap stayed flat at 66.5 MB while process RSS grew until SIGKILL).
#
# Writes ONE checkpoint per chunk so an OOM kill loses at most one chunk.
# Output schema == run_coverage.R (reuses run_one_rep) -> analyze.R-compatible.
#
# Args (positional): <n> <rep_lo> <rep_hi> <chunk_out.rds>
# Run: invoked by the driver; not called directly.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(readr); library(fs)
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

STUDY_DIR <- fs::path("doubletree", "simulations", "single_tree_coverage")
source(fs::path(STUDY_DIR, "code", "oracle_theta_star.R"))   # THETA0, N_GRID, load_oracle
source(fs::path(STUDY_DIR, "code", "harness.R"))             # run_one_rep, make_seed, BASE_SEED

args   <- commandArgs(trailingOnly = TRUE)
n_val  <- as.integer(args[[1]])
rep_lo <- as.integer(args[[2]])
rep_hi <- as.integer(args[[3]])
out    <- args[[4]]

set.seed(BASE_SEED)                       # per-rep seeds are deterministic via make_seed
oracle_tbl <- load_oracle(STUDY_DIR)
ts <- oracle_tbl$theta_star[oracle_tbl$dgp == "complex" & oracle_tbl$n == n_val]

reps <- rep_lo:rep_hi
rows <- vector("list", length(reps))
t0 <- Sys.time()
for (j in seq_along(reps)) {
  rows[[j]] <- run_one_rep("complex", n_val, reps[j], ts)
}
res <- purrr::list_rbind(rows)
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

readr::write_rds(res, out)

# Completion sentinel, published by ATOMIC RENAME and only after the .rds above is
# fully written and closed. The driver's resume test keys on THIS, never on the mere
# existence of the .rds -- a write_rds truncated by a SIGKILL would otherwise be
# mistaken for a finished chunk and silently poison the merged results.
done_tmp <- paste0(out, ".done.tmp")
cat(sprintf("n=%d reps=%d-%d rows=%d\n", n_val, rep_lo, rep_hi, nrow(res)), file = done_tmp)
file.rename(done_tmp, paste0(out, ".done"))   # atomic within a filesystem (APFS)

cat(sprintf("CHUNK_OK n=%d reps=%d-%d (%d rows) %.1fs -> %s\n",
            n_val, rep_lo, rep_hi, nrow(res), el, out))
