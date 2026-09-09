#!/bin/bash
# =============================================================================
# slurm/quick_test.sh -- 30-second smoke test of the cluster path, no SLURM
# =============================================================================
# Run LOCALLY (before pushing) and again ON O2 (after installing the packages),
# from anywhere:
#
#   bash slurm/quick_test.sh
#   R_LIBS=/path/to/lib bash slurm/quick_test.sh     # test against a chosen library
#
# WHAT IT ACTUALLY VERIFIES -- it is not a results check, it is a WIRING check:
#
#   1. run_single_replication.R loads the study through `library(doubletree)` /
#      `library(optimaltrees)` (H2H_USE_INSTALLED=1) rather than through
#      code/common.R's default pkgload::load_all() dev chain. If the packages are
#      not installed, this fails HERE -- on a laptop in 10 seconds -- instead of
#      failing 134 times on the cluster.
#   2. The array-index -> block map in slurm/units.R resolves, and the shard lands
#      at the filename slurm/combine_results.R expects.
#   3. A real replication of a real cell executes end to end: DGP draw, every arm
#      of the regime, the derived anchor arm, the result row schema.
#
# It runs R6 (the cheapest regime, ~7 ms/rep at n = 200) with --reps-per-job 2, so
# it costs seconds and cannot be confused with production output: shards go to a
# throwaway directory under the study's own results/ tree, NOT to scratch, and
# combine_results.R is never pointed at them.
#
# Exit status 0 = the cluster path works. Anything else = do not deploy.
# =============================================================================

set -euo pipefail

SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY_DIR="$(dirname "${SLURM_DIR}")"
PKG_ROOT="$(cd "${STUDY_DIR}/../.." && pwd)"

TEST_DIR="${STUDY_DIR}/results/quick_test"
REGIME="${REGIME:-R6}"
TASK_ID="${TASK_ID:-1}"
REPS="${REPS:-2}"

echo "=========================================================="
echo " quick_test.sh -- head_to_head_comparison cluster path"
echo "=========================================================="
echo " pkg root : ${PKG_ROOT}"
echo " regime   : ${REGIME}   task id: ${TASK_ID}   reps: ${REPS}"
echo " shard dir: ${TEST_DIR}"
echo " R_LIBS   : ${R_LIBS:-<default>}"
echo

rm -rf "${TEST_DIR}"
mkdir -p "${TEST_DIR}"

# --- 1 are the packages INSTALLED (not dev-loaded)? --------------------------
# The single most common cluster failure for this study, so it is checked first
# and reported as itself rather than as a downstream "object not found".
echo "--- 1. installed-package check (library(), not pkgload) ---"
Rscript -e '
  pkgs <- c("optimaltrees", "doubletree")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop(sprintf(paste0("package(s) not installed: %s\n",
                        "  local: R CMD INSTALL optimaltrees && R CMD INSTALL doubletree\n",
                        "  on O2: see slurm/README_O2.md, prerequisites"),
                 paste(miss, collapse = ", ")), call. = FALSE)
  }
  for (p in pkgs) cat(sprintf("  %s %s\n", p, utils::packageVersion(p)))
'
echo

# --- 2 does the sizing table resolve? ---------------------------------------
echo "--- 2. job-array sizing (slurm/units.R) ---"
Rscript "${SLURM_DIR}/print_sizing.R" --pkg-root "${PKG_ROOT}"
echo

# --- 3 run one real block, tiny ---------------------------------------------
echo "--- 3. one real block: ${REGIME} task ${TASK_ID}, ${REPS} replication(s) ---"
Rscript "${SLURM_DIR}/run_single_replication.R" \
  --regime      "${REGIME}" \
  --task-id     "${TASK_ID}" \
  --pkg-root    "${PKG_ROOT}" \
  --scratch-dir "${TEST_DIR}" \
  --reps-per-job "${REPS}" \
  --flush-every 1
echo

# --- 4 inspect the shard ----------------------------------------------------
# The shard is opened and its contract checked, not just its existence: a file
# that exists but holds zero usable rows is the failure this step exists to catch.
# NOTE: the regex below uses [.] rather than \\. on purpose -- `Rscript -e`
# performs its own backslash processing, so a backslash escape written here
# reaches R as an invalid single backslash.
echo "--- 4. shard contract ---"
Rscript -e '
  dir <- commandArgs(trailingOnly = TRUE)[[1]]
  f <- list.files(dir, pattern = "^batch_.*[.]rds$", full.names = TRUE)
  if (length(f) != 1L) stop(sprintf("expected 1 shard in %s, found %d", dir, length(f)))
  p <- readRDS(f)
  stopifnot(all(c("unit", "results", "dgp_summary", "runtime") %in% names(p)))
  res <- p$results
  cat(sprintf("  shard        : %s\n", basename(f)))
  cat(sprintf("  unit         : %s / %s / n = %d, reps %d-%d\n",
              p$unit$regime, p$unit$dgp, p$unit$n, p$unit$rep_from, p$unit$rep_to))
  cat(sprintf("  rows         : %d  (%d rep(s) x %d arm(s))\n",
              nrow(res), length(unique(res$rep)), length(unique(res$arm))))
  cat(sprintf("  arms         : %s\n", paste(unique(res$arm), collapse = ", ")))
  cat(sprintf("  load path    : %s (doubletree %s)\n",
              p$runtime$pkg_load_path, p$runtime$doubletree_version))
  cat(sprintf("  elapsed      : %.2f s (projected %.2f s for the full block)\n",
              p$runtime$elapsed_secs, p$unit$est_secs))
  if (!identical(p$runtime$pkg_load_path, "installed")) {
    stop("shard reports pkg_load_path other than installed; the cluster path did NOT take effect.")
  }
  fail <- res[!is.na(res$error_message), , drop = FALSE]
  if (nrow(fail)) {
    cat("  ERRORS:\n"); for (m in unique(fail$error_message)) cat("    x ", substr(m, 1, 300), "\n", sep = "")
    stop("replication(s) failed; fix before deploying.")
  }
  if (!all(is.finite(res$theta), is.finite(res$sigma))) {
    stop("non-finite theta/sigma in the shard.")
  }
  cat(sprintf("  theta range  : [%.5g, %.5g]  (theta0 = %.5g)\n",
              min(res$theta), max(res$theta), res$theta0[[1]]))
  cat(sprintf("  covered      : %d/%d\n", sum(res$covered), nrow(res)))
' "${TEST_DIR}"

echo
echo "=========================================================="
echo " PASS -- the library()-based cluster path works end to end."
echo " Shards left in ${TEST_DIR} (throwaway; not read by combine_results.R)."
echo "=========================================================="
