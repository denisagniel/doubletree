#!/bin/bash
# =============================================================================
# profile.slurm -- submittable batch job for per-method WORST-CELL profiling
# =============================================================================
# Runs slurm/profile_per_method.R as a BATCH job (not interactively). Profiling
# hammers the worst cell (n=2000, dgp=continuous) across all 7 methods -- the
# documented Rashomon/M-split blow-up -- so it is long and memory-heavy and does
# NOT belong in an interactive salloc. Submit and walk away:
#
#   sbatch slurm/profile.slurm
#
# then, once config/sizing_<method>.env exist for all 7 methods:
#   bash slurm/submit_per_method.sh
#
# GATE (same as submit_per_method.sh): the agent cannot `git push`, so REINSTALL
# both packages on this node first, else profiling sizes STALE code:
#   R CMD INSTALL doubletree optimaltrees   # from the repo root
#
# MEMORY: --mem below is the SLURM allocation; --mem-cap-gb (passed to the
# profiler) is the PER-UNIT ceiling and MUST be strictly below --mem so the
# in-process guard kills a runaway unit BEFORE SLURM OS-kills the whole job (which
# would lose every method's timing). 24 GB cap inside a 48 GB allocation leaves
# headroom for the R session + a mid-kill spike. Override both together if needed.
# =============================================================================
#SBATCH --job-name=arb-profile
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=0-06:00
#SBATCH --output=logs/profile_%j.out
#SBATCH --error=logs/profile_%j.err

set -euo pipefail

# --- Tunables (override at submit time, e.g. sbatch --export=ALL,TARGET_HOURS=3 ...) ---
TARGET_HOURS="${TARGET_HOURS:-2}"     # target wall time per PRODUCTION task (sizing input)
MEM_CAP_GB="${MEM_CAP_GB:-24}"        # per-unit kill ceiling; MUST be < --mem (48G)
N_UNITS_PROBE="${N_UNITS_PROBE:-5}"   # probe units per method at the worst cell

# --- Resolve study dir --------------------------------------------------------
# Under `sbatch`, SLURM runs a SPOOLED COPY of this script (e.g.
# /var/spool/slurmd/job<id>/slurm_script), so ${BASH_SOURCE[0]} does NOT point
# into the repo -- deriving STUDY_DIR from it yields /var/spool/slurmd and the
# mkdir below fails "Permission denied", aborting under set -e before any sizing
# runs. SLURM sets SLURM_SUBMIT_DIR to the dir `sbatch` was invoked from (the
# study dir, per the header) and cd's into it, so prefer that. Fall back to
# BASH_SOURCE for an interactive `bash slurm/profile.slurm`.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  STUDY_DIR="${SLURM_SUBMIT_DIR}"
  SLURM_DIR="${STUDY_DIR}/slurm"
else
  SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  STUDY_DIR="$(dirname "${SLURM_DIR}")"
fi

# Fail loudly if we did not land in the study dir (e.g. sbatch'd from elsewhere):
# the profiler needs config/grid.R and slurm/profile_per_method.R to be here.
if [[ ! -f "${SLURM_DIR}/profile_per_method.R" || ! -f "${STUDY_DIR}/config/grid.R" ]]; then
  echo "ERROR: study dir '${STUDY_DIR}' is missing slurm/profile_per_method.R or config/grid.R." >&2
  echo "       Submit from the study dir:  cd <study-dir> && sbatch slurm/profile.slurm" >&2
  exit 1
fi

mkdir -p "${STUDY_DIR}/logs"

# --- Module environment (module R has no devtools; profiler uses base + callr/ps) ---
module load gcc/14.2.0 2>/dev/null || module load gcc || true
module load R/4.4.2   2>/dev/null || module load R

# --- Freshness gate: fail loudly if project packages are missing/stale-uninstalled ---
Rscript -e '
  pkgs <- c("doubletree", "optimaltrees")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop(sprintf("package(s) NOT installed on this node: %s -- R CMD INSTALL them first (agent cannot git push).",
                 paste(miss, collapse = ", ")))
  for (p in pkgs) cat(sprintf("preflight: %s %s OK\n", p, utils::packageVersion(p)))
'

echo "[$(date '+%F %T')] host=$(hostname) job=${SLURM_JOB_ID:-<none>} profiling: target=${TARGET_HOURS}h mem-cap=${MEM_CAP_GB}G n-units=${N_UNITS_PROBE}"

# --- Profile every method at its worst cell; writes config/sizing_<method>.env ---
# --mem-cap-gb (24) < --mem (48G): the per-unit guard fires before SLURM.
Rscript "${SLURM_DIR}/profile_per_method.R" \
  --study-dir    "${STUDY_DIR}" \
  --target-hours "${TARGET_HOURS}" \
  --mem-cap-gb   "${MEM_CAP_GB}" \
  --n-units      "${N_UNITS_PROBE}"

echo "[$(date '+%F %T')] profiling finished with status $?; sizing files:"
ls -1 "${STUDY_DIR}/config"/sizing_*.env 2>/dev/null || echo "  (none written -- check the log above)"
