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
# UNATTENDED ONE-SHOT (profile THEN auto-submit ~1000 long tasks):
#   sbatch --export=ALL,TARGET_TASKS=1000,AUTO_SUBMIT=1 slurm/profile.slurm
# This sizes for ~1000 tasks and, ONLY if all 7 sizing envs come out valid,
# chain-launches submit_per_method.sh from within this job. If profiling is
# truncated (timeout/OOM leaves <7 envs), it submits NOTHING and exits non-zero.
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
# `medium` (not `short`): profiling runs all 7 methods SERIALLY, and the heavy
# ones (doubletree/dt_averaged, the n=2000 continuous Rashomon blow-up) burn hours
# at the worst cell. On `short` (12h cap) the 6h run TIMED OUT after only 3 of 7
# methods (job 47933203), leaving later methods' sizing_<method>.env unwritten.
# `medium` (5-day cap) removes the 12h ceiling so profiling can't be truncated by
# the cap, and schedules sooner than `long`; the 24h wall below is a generous
# ceiling for a run-once sizing job. If runtime is still a concern, cut probe cost
# with --export=ALL,N_UNITS_PROBE=3.
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=1-00:00
#SBATCH --output=logs/profile_%j.out
#SBATCH --error=logs/profile_%j.err

set -euo pipefail

# --- Tunables (override at submit time, e.g. sbatch --export=ALL,TARGET_HOURS=3 ...) ---
TARGET_HOURS="${TARGET_HOURS:-2}"     # target wall time per PRODUCTION task (sizing input)
# TARGET_TASKS: if set, size for this TOTAL number of production array tasks across
# all 7 methods (cost-proportional -> uniform per-task walltime, partition auto-picked)
# instead of a per-task wall-time target. Use to stay under the 10k MaxSubmitJobs cap
# and submit ONE shot, unattended, e.g.  sbatch --export=ALL,TARGET_TASKS=1000 slurm/profile.slurm
TARGET_TASKS="${TARGET_TASKS:-}"      # empty => classic target-hours sizing
MEM_CAP_GB="${MEM_CAP_GB:-24}"        # per-unit kill ceiling; MUST be < --mem (48G)
N_UNITS_PROBE="${N_UNITS_PROBE:-5}"   # probe units per method at the worst cell
# AUTO_SUBMIT=1: after profiling, if AND ONLY IF all 7 sizing_<method>.env are valid,
# chain-launch the production arrays via submit_per_method.sh from within this job --
# so a single `sbatch` profiles then submits, unattended. If any method's sizing is
# missing/invalid (e.g. profiling timed out after N<7 methods, as job 47933203 did),
# NOTHING is submitted and this job exits non-zero (better than launching a partial,
# broken study). Requires the cluster to permit sbatch-from-within-a-job (O2 does).
AUTO_SUBMIT="${AUTO_SUBMIT:-0}"

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

# Sizing mode: TARGET_TASKS (if set) wins over TARGET_HOURS.
SIZING_ARGS=(--target-hours "${TARGET_HOURS}")
if [[ -n "${TARGET_TASKS}" ]]; then
  SIZING_ARGS=(--target-tasks "${TARGET_TASKS}")
fi

echo "[$(date '+%F %T')] host=$(hostname) job=${SLURM_JOB_ID:-<none>} profiling: ${SIZING_ARGS[*]} mem-cap=${MEM_CAP_GB}G n-units=${N_UNITS_PROBE}"

# --- Profile every method at its worst cell; writes config/sizing_<method>.env ---
# --mem-cap-gb (24) < --mem (48G): the per-unit guard fires before SLURM.
# Do not let a profiler non-zero exit abort the script under `set -e` before we can
# report it; capture the status explicitly.
prof_status=0
Rscript "${SLURM_DIR}/profile_per_method.R" \
  --study-dir    "${STUDY_DIR}" \
  "${SIZING_ARGS[@]}" \
  --mem-cap-gb   "${MEM_CAP_GB}" \
  --n-units      "${N_UNITS_PROBE}" || prof_status=$?

echo "[$(date '+%F %T')] profiling finished with status ${prof_status}; sizing files:"
ls -1 "${STUDY_DIR}/config"/sizing_*.env 2>/dev/null || echo "  (none written -- check the log above)"

# --- Optional gated auto-submit ----------------------------------------------
# Chain the production submission ONLY if profiling succeeded AND all 7 methods have
# a valid sizing env (positive TOTAL_TASKS + non-empty WALLTIME/PARTITION). A partial
# set means profiling was truncated (timeout/OOM); submitting it would launch a broken
# study, so we refuse and exit non-zero instead.
if [[ "${AUTO_SUBMIT}" == "1" ]]; then
  METHODS=(full crossfit doubletree dt_averaged msplit msplit_averaged single_tree)
  if (( prof_status != 0 )); then
    echo "AUTO_SUBMIT: profiling exited ${prof_status}; NOT submitting." >&2
    exit "${prof_status}"
  fi
  ok=1
  for m in "${METHODS[@]}"; do
    envf="${STUDY_DIR}/config/sizing_${m}.env"
    if [[ ! -f "${envf}" ]]; then
      echo "AUTO_SUBMIT: missing ${envf}; profiling incomplete. NOT submitting." >&2; ok=0; continue
    fi
    # Source in a subshell to read the values without polluting this shell.
    tt=$(bash -c "set -a; source '${envf}'; printf '%s' \"\${TOTAL_TASKS:-}\"")
    wt=$(bash -c "set -a; source '${envf}'; printf '%s' \"\${WALLTIME:-}\"")
    pt=$(bash -c "set -a; source '${envf}'; printf '%s' \"\${PARTITION:-}\"")
    if [[ -z "${tt}" || "${tt}" -lt 1 || -z "${wt}" || -z "${pt}" ]]; then
      echo "AUTO_SUBMIT: ${envf} invalid (TOTAL_TASKS='${tt}' WALLTIME='${wt}' PARTITION='${pt}'). NOT submitting." >&2
      ok=0
    fi
  done
  if (( ok != 1 )); then
    echo "AUTO_SUBMIT: sizing incomplete/invalid across the 7 methods; submitted NOTHING." >&2
    exit 1
  fi
  echo "[$(date '+%F %T')] AUTO_SUBMIT: all 7 sizing envs valid; launching submit_per_method.sh"
  bash "${SLURM_DIR}/submit_per_method.sh"
fi
