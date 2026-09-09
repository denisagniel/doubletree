#!/bin/bash
# =============================================================================
# slurm/launch_all_simulations.sh -- submit the FULL sweep (all six regimes)
# =============================================================================
# Run ON O2:
#   bash slurm/launch_all_simulations.sh
#
# Deliberately a one-line wrapper around launch_subset.sh rather than a second
# copy of the submit logic: two launchers that drift apart is how a sweep ends up
# with half its arrays sized from a stale table.
#
# Budget for the full sweep (from the 2026-09-09 pilot's measured costs, via
# slurm/units.R; see README_O2.md §"Job array size and budget"):
#   134 array tasks across 6 arrays, 54.2 h of compute, longest task ~29 min,
#   wall clock ~30-45 min if the cluster runs them concurrently.
#   77 % of the compute is R2 (22.8 h) + R5 (18.1 h), both estimate_att_crossfit()-heavy.
#
# All of launch_subset.sh's environment gates apply (RUN_ID to resume,
# H2H_TARGET_SECS to re-batch, CONCURRENCY_CAP, DRY_RUN=1 to inspect first).
# Start with DRY_RUN=1 if this is the first submit after a code change.
# =============================================================================

set -euo pipefail

SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "${SLURM_DIR}/launch_subset.sh" R1 R2 R3 R4 R5 R6
