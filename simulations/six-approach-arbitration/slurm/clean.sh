#!/bin/bash
# =============================================================================
# clean.sh -- remove stale scratch runs and orphaned logs
# =============================================================================
# Solves "stale results/logs hang around after I update simulation code".
# Lists every run-id present in scratch, marks which have a combined result in
# home results/, and removes superseded scratch runs on confirmation. NEVER
# touches home results/ (final artifacts) unless you explicitly ask.
#
# Usage:
#   bash slurm/clean.sh                 # interactive: list + confirm scratch cleanup
#   bash slurm/clean.sh --keep-latest   # remove all scratch runs EXCEPT the newest
#   bash slurm/clean.sh --combined-only # remove only scratch runs already combined to home
#   bash slurm/clean.sh --run-id RID    # remove one specific scratch run
#   bash slurm/clean.sh --dry-run       # show what would be removed, do nothing
# =============================================================================

set -euo pipefail

HMS_ID="dma12"
PROJECT_NAME="global-scholars"
STUDY_NAME="six-approach-arbitration"

SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDY_DIR="$(dirname "${SLURM_DIR}")"
SCRATCH_ROOT="/n/scratch/users/${HMS_ID:0:1}/${HMS_ID}/${PROJECT_NAME}/${STUDY_NAME}"
RESULTS_DIR="${STUDY_DIR}/results"

MODE="interactive"
TARGET_RUN=""
DRY_RUN=0
for arg in "$@"; do
  case "${arg}" in
    --keep-latest)   MODE="keep-latest" ;;
    --combined-only) MODE="combined-only" ;;
    --dry-run)       DRY_RUN=1 ;;
    --run-id) MODE="one" ;;
    *) if [[ "${MODE}" == "one" && -z "${TARGET_RUN}" ]]; then TARGET_RUN="${arg}"; fi ;;
  esac
done

if [[ ! -d "${SCRATCH_ROOT}" ]]; then
  echo "Nothing to clean: ${SCRATCH_ROOT} does not exist."
  exit 0
fi

# Enumerate scratch runs, newest first.
mapfile -t RUNS < <(ls -1t "${SCRATCH_ROOT}" 2>/dev/null || true)
if (( ${#RUNS[@]} == 0 )); then
  echo "No scratch runs under ${SCRATCH_ROOT}."
  exit 0
fi

LATEST="${RUNS[0]}"
has_combined() { [[ -f "${RESULTS_DIR}/$1.rds" ]]; }

echo "Scratch runs under ${SCRATCH_ROOT}:"
for r in "${RUNS[@]}"; do
  tag=""
  [[ "${r}" == "${LATEST}" ]] && tag="${tag} [latest]"
  has_combined "${r}" && tag="${tag} [combined->home]"
  sz=$(du -sh "${SCRATCH_ROOT}/${r}" 2>/dev/null | cut -f1)
  printf "  %-28s %6s%s\n" "${r}" "${sz:-?}" "${tag}"
done
echo

# Decide which runs to remove.
TO_REMOVE=()
case "${MODE}" in
  one)
    [[ -z "${TARGET_RUN}" ]] && { echo "ERROR: --run-id needs a value." >&2; exit 1; }
    TO_REMOVE=("${TARGET_RUN}") ;;
  keep-latest)
    for r in "${RUNS[@]}"; do [[ "${r}" != "${LATEST}" ]] && TO_REMOVE+=("${r}"); done ;;
  combined-only)
    for r in "${RUNS[@]}"; do has_combined "${r}" && TO_REMOVE+=("${r}"); done ;;
  interactive)
    # Default: offer to remove every run that has already been combined to home.
    for r in "${RUNS[@]}"; do has_combined "${r}" && TO_REMOVE+=("${r}"); done
    if (( ${#TO_REMOVE[@]} == 0 )); then
      echo "No scratch runs have been combined to home yet; nothing suggested for removal."
      echo "Use --keep-latest or --run-id RID to force cleanup."
      exit 0
    fi ;;
esac

if (( ${#TO_REMOVE[@]} == 0 )); then
  echo "Nothing matched for removal."
  exit 0
fi

echo "Will remove ${#TO_REMOVE[@]} scratch run(s):"
for r in "${TO_REMOVE[@]}"; do echo "  ${SCRATCH_ROOT}/${r}"; done
echo "(home results/ are NOT touched.)"

if (( DRY_RUN == 1 )); then
  echo "--dry-run: no changes made."
  exit 0
fi

read -r -p "Proceed? [y/N] " ans
if [[ "${ans}" != "y" && "${ans}" != "Y" ]]; then
  echo "Aborted."
  exit 0
fi

for r in "${TO_REMOVE[@]}"; do
  rm -rf "${SCRATCH_ROOT:?}/${r}"
  echo "removed ${r}"
  # If we removed the run that logs/latest pointed at, drop the dangling link.
  if [[ -L "${STUDY_DIR}/logs/latest" && ! -e "${STUDY_DIR}/logs/latest" ]]; then
    rm -f "${STUDY_DIR}/logs/latest"
    echo "cleared dangling ./logs/latest"
  fi
done
echo "Done."
