#!/bin/bash
# =============================================================================
# slurm/launch_all_simulations.sh -- submit every MISSING wave (n = 2000 and 8000)
# =============================================================================
# Run ON O2:
#   bash slurm/launch_all_simulations.sh
#
# Deliberately a one-line wrapper around launch_subset.sh rather than a second copy
# of the submit logic: two launchers that drift apart is how a sweep ends up with
# half its arrays sized from a stale table.
#
# WHAT "ALL" MEANS HERE. The spec's grid is 5 DGP variants x n in {500, 2000, 8000},
# R = 300 (R = 150 at n = 8000 by the spec's own sanctioned exception). The five
# n = 500 cells were completed sequentially on 2026-09-01 and are on disk as
# results/cell_*_n500_r300.rds. This launcher therefore submits the TEN remaining
# cells -- 5 variants x n in {2000, 8000}, 2,250 replications -- and NOT n = 500. To
# re-run or backfill an n = 500 cell, ask for it explicitly:
#
#   bash slurm/launch_subset.sh n500
#
# Budget for the missing waves (from the 2026-09-09 measured costs, via
# slurm/units.R; see README_O2.md "Job array size and budget"):
#   20 array tasks across 2 arrays, 3.94 h of compute, longest task ~14 min,
#   wall clock ~15 min if the cluster runs them concurrently.
#   60 % of the compute is the n = 8000 wave; 98 % of every task is
#   estimate_att_crossfit()'s nested CV, not the estimator under test.
#
# All of launch_subset.sh's environment gates apply (RUN_ID to resume,
# HIS_TARGET_SECS to re-batch, HIS_REPS/HIS_REPS_MAX_N, CONCURRENCY_CAP,
# DRY_RUN=1 to inspect first). Start with DRY_RUN=1 if this is the first submit
# after a code change.
# =============================================================================

set -euo pipefail

SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "${SLURM_DIR}/launch_subset.sh" n2000 n8000
