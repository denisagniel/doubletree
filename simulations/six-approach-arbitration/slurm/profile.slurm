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
#   TARGET_TASKS=1000 AUTO_SUBMIT=1 sbatch --export=ALL slurm/profile.slurm
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
# with N_UNITS_PROBE=3 and a plain --export=ALL.
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=1-00:00
#SBATCH --output=logs/profile_%j.out
#SBATCH --error=logs/profile_%j.err

# --- Module environment -------------------------------------------------------
# Backfilled 2026-09-23. `module` is a bash FUNCTION that Lmod `export -f`s into the
# environment, so it is INHERITED, not found on PATH: present when a human submits from
# an interactive shell, ABSENT under `ssh host 'bash -s'`. With strict mode already on,
# the bare `module load` below then aborted the task with no useful output.
#
# Bootstrap MUST precede `set -euo pipefail`: /etc/profile.d/* scripts reference unset
# variables and return non-zero, both fatal under `set -eu`.
#
# (Separately: the 2026-09-18 mass cancellation of this study was NOT this bug -- it was
# `sbatch --export=ALL,KEY=value`, fixed in the launcher. See O2_SSH_GOTCHAS.md §12.)
set +eu
if ! command -v module >/dev/null 2>&1; then
  for profile_script in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh /etc/profile; do
    [ -r "${profile_script}" ] && . "${profile_script}" && break
  done
fi
set -euo pipefail

# Distinguishable exit codes: 78 module env, 79 Rscript missing.
die() { echo "[$(date '+%F %T')] FATAL(${2}): ${1}" >&2; exit "${2}"; }
command -v module >/dev/null 2>&1 || die "module unavailable after profile bootstrap" 78
module purge 2>/dev/null || true
module load gcc/14.2.0 || die "module load gcc/14.2.0" 78
module load R/4.4.2 || die "module load R/4.4.2" 78
command -v Rscript >/dev/null 2>&1 || die "Rscript not on PATH after module load" 79
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
