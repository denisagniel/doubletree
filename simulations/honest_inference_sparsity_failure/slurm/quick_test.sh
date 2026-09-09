#!/bin/bash
# =============================================================================
# slurm/quick_test.sh -- ~1-minute smoke test of the cluster path, no SLURM
# =============================================================================
# Run LOCALLY (before pushing) and again ON O2 (after installing the packages),
# from anywhere:
#
#   bash slurm/quick_test.sh
#   REPS=1 bash slurm/quick_test.sh                  # faster
#   WAVES="n8000" bash slurm/quick_test.sh           # one wave only
#   R_LIBS=/path/to/lib bash slurm/quick_test.sh     # test against a chosen library
#
# WHAT IT ACTUALLY VERIFIES -- it is not a results check, it is a WIRING check:
#
#   1. run_single_replication.R loads the study through `library(doubletree)` /
#      `library(optimaltrees)` (HIS_USE_INSTALLED=1) rather than through
#      code/common.R's default pkgload::load_all() dev chain, AND the namespaces
#      really did come from an installed library (assert_installed_load()). If the
#      packages are not installed, this fails HERE -- on a laptop in seconds --
#      instead of failing 20 times on the cluster.
#   2. The array-index -> block map in slurm/units.R resolves, and the shard lands
#      at the filename slurm/combine_results.R expects.
#   3. A real replication of a real PRODUCTION cell executes end to end at BOTH
#      missing sample sizes: DGP dial construction with exact enumeration, the draw,
#      estimate_att(), estimate_att_crossfit(), all three intervals, the fidelity
#      diagnostic, and the deterministic anchor-containment assertion inside
#      run_one_rep().
#
# BOTH n = 2000 AND n = 8000 are exercised, not just the cheap one. This study's
# cost and its failure surface are both dominated by estimate_att_crossfit()'s
# nested CV, whose behaviour at n = 8000 is what the deployment actually depends on;
# a smoke test that only ran the cheap wave would verify the wiring and none of the
# risk. The price is ~25 s of extra wall time.
#
# Shards go to a throwaway directory under the study's own results/ tree, NOT to
# scratch, and combine_results.R is never pointed at them. They cannot be mistaken
# for production output: a --reps-per-job shard holds fewer replications than the
# production map assigns to that block, so combine_results.R's completeness check
# would refuse to write a cell from it, and the throwaway dir is deleted and
# recreated on every run.
#
# Exit status 0 = the cluster path works. Anything else = do not deploy.
# =============================================================================

set -euo pipefail

SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY_DIR="$(dirname "${SLURM_DIR}")"
PKG_ROOT="$(cd "${STUDY_DIR}/../.." && pwd)"

TEST_DIR="${STUDY_DIR}/results/quick_test"
WAVES="${WAVES:-n2000 n8000}"
TASK_ID="${TASK_ID:-1}"
REPS="${REPS:-2}"

echo "=========================================================="
echo " quick_test.sh -- honest_inference_sparsity_failure (S2)"
echo "=========================================================="
echo " pkg root : ${PKG_ROOT}"
echo " waves    : ${WAVES}   task id: ${TASK_ID}   reps/wave: ${REPS}"
echo " shard dir: ${TEST_DIR}"
echo " R_LIBS   : ${R_LIBS:-<default>}"
echo

rm -rf "${TEST_DIR}"
mkdir -p "${TEST_DIR}"

# --- 1 are the packages INSTALLED (not dev-loaded)? --------------------------
# The single most common cluster failure for this study, so it is checked first and
# reported as itself rather than as a downstream "object not found". optimaltrees
# >= 0.4.1 is checked by VERSION, not just presence: this study was blocked on
# exactly that constraint from 2026-09-04 to 2026-09-09.
#
# --no-init-file everywhere below, for the reason run_simulations.slurm documents:
# the doubletree repository ships a tracked .Rprofile that dev-loads
# ../optimaltrees, which would make this check report the SOURCE version.
echo "--- 1. installed-package check (library(), not pkgload) ---"
Rscript --no-init-file -e '
  need <- c(optimaltrees = "0.4.1", doubletree = "0.0.0.9000")
  for (p in names(need)) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf(paste0("package %s is not installed\n",
                          "  local: R CMD INSTALL optimaltrees && R CMD INSTALL doubletree\n",
                          "  on O2: see slurm/README_O2.md, prerequisites"), p),
           call. = FALSE)
    }
    v <- utils::packageVersion(p)
    if (v < need[[p]]) stop(sprintf("%s %s installed, need >= %s", p, v, need[[p]]), call. = FALSE)
    cat(sprintf("  %s %s  (%s)\n", p, v, dirname(find.package(p))))
  }
'
echo

# --- 2 does the sizing table resolve? ---------------------------------------
echo "--- 2. job-array sizing (slurm/units.R) ---"
Rscript --no-init-file "${SLURM_DIR}/print_sizing.R" --pkg-root "${PKG_ROOT}"
echo

# --- 3 run one real block per wave, tiny ------------------------------------
for WAVE in ${WAVES}; do
  echo "--- 3.${WAVE} one real block: ${WAVE} task ${TASK_ID}, ${REPS} replication(s) ---"
  Rscript --no-init-file "${SLURM_DIR}/run_single_replication.R" \
    --wave        "${WAVE}" \
    --task-id     "${TASK_ID}" \
    --pkg-root    "${PKG_ROOT}" \
    --scratch-dir "${TEST_DIR}/${WAVE}" \
    --reps-per-job "${REPS}" \
    --flush-every 1
  echo
done

# --- 4 inspect the shards ---------------------------------------------------
# The shards are opened and their contract checked, not just their existence: a file
# that exists but holds zero usable rows is the failure this step exists to catch.
# NOTE: the regex below uses [.] rather than \\. on purpose -- `Rscript -e` performs
# its own backslash processing, so a backslash escape written here reaches R as an
# invalid single backslash.
echo "--- 4. shard contract ---"
Rscript --no-init-file -e '
  dir <- commandArgs(trailingOnly = TRUE)[[1]]
  f <- list.files(dir, pattern = "^batch_.*[.]rds$", full.names = TRUE, recursive = TRUE)
  if (!length(f)) stop(sprintf("no shard found under %s", dir))
  for (path in sort(f)) {
    p <- readRDS(path)
    stopifnot(all(c("unit", "results", "dgp_summary", "meta", "runtime") %in% names(p)))
    res <- p$results
    cat(sprintf("  shard      : %s\n", basename(path)))
    cat(sprintf("  unit       : %s / %s (eps = %.2f) / n = %d, reps %d-%d of %d\n",
                p$unit$wave, p$unit$dgp, p$unit$eps, p$unit$n,
                p$unit$rep_from, p$unit$rep_to, p$unit$reps_total))
    cat(sprintf("  rows       : %d\n", nrow(res)))
    cat(sprintf("  load path  : %s (doubletree %s from %s)\n",
                p$runtime$pkg_load_path, p$runtime$doubletree_version,
                p$runtime$doubletree_dir))
    cat(sprintf("  elapsed    : %.2f s for %d rep(s) = %.2f s/rep (block projected %.1f s)\n",
                p$runtime$elapsed_secs, nrow(res),
                p$runtime$elapsed_secs / nrow(res), p$unit$est_secs))
    if (!identical(p$runtime$pkg_load_path, "installed")) {
      stop("shard reports pkg_load_path other than installed; the cluster path did NOT take effect.")
    }
    fail <- res[!is.na(res$error_message), , drop = FALSE]
    if (nrow(fail)) {
      cat("  ERRORS:\n"); for (m in unique(fail$error_message)) cat("    x ", substr(m, 1, 400), "\n", sep = "")
      stop("replication(s) failed; fix before deploying.")
    }
    # Every quantity the spec asks for must be present and finite -- an NA here
    # would mean a silently short row rather than a loud failure.
    need <- c("theta_full", "sigma_full", "theta_cf", "sigma_cf", "delta",
              "ci_anchor_lo", "ci_anchor_hi", "ci_honest_lo", "ci_honest_hi",
              "fid_se_hat", "fid_stat")
    bad <- need[!vapply(need, function(k) all(is.finite(res[[k]])), logical(1))]
    if (length(bad)) stop(sprintf("non-finite %s in the shard.", paste(bad, collapse = ", ")))
    if (any(is.na(res$covered), is.na(res$covered_anchor), is.na(res$covered_honest),
            is.na(res$covered_cf), is.na(res$fid_reject))) {
      stop("NA in a coverage or rejection indicator.")
    }
    cat(sprintf("  theta0     : %.5g   theta_full: [%.5g, %.5g]   theta_cf: [%.5g, %.5g]\n",
                res$theta0[[1]], min(res$theta_full), max(res$theta_full),
                min(res$theta_cf), max(res$theta_cf)))
    cat(sprintf("  covered    : wald %d/%d  anchor %d/%d  honest %d/%d  crossfit %d/%d   fid_reject %d/%d\n",
                sum(res$covered), nrow(res), sum(res$covered_anchor), nrow(res),
                sum(res$covered_honest), nrow(res), sum(res$covered_cf), nrow(res),
                sum(res$fid_reject), nrow(res)))
    cat(sprintf("  widths     : anchor %.4g  honest %.4g  wald %.4g\n",
                mean(res$ci_anchor_width), mean(res$ci_honest_width),
                mean(res$ci_width)))
    cat("\n")
  }
  cat(sprintf("  %d shard(s) checked.\n", length(f)))
' "${TEST_DIR}"

echo
echo "=========================================================="
echo " PASS -- the library()-based cluster path works end to end"
echo " at every wave tested (${WAVES})."
echo " Shards left in ${TEST_DIR} (throwaway; not read by combine_results.R)."
echo "=========================================================="
