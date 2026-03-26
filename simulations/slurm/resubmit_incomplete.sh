#!/bin/bash

# Resubmit incomplete DML-ATT simulation jobs
#
# Based on current completion status, resubmit only configurations
# that don't have 1000 replications yet.

set -e

cd "$(dirname "$0")/.."
mkdir -p logs

echo "==========================================="
echo "Resubmitting Incomplete Simulations"
echo "==========================================="
echo ""

# Configurations that need resubmission (< 1000 reps completed)
# Format: DGP N METHOD N_BATCHES

declare -a INCOMPLETE=(
  # dgp2 - missing tree/rashomon/partial forest
  "dgp2 1600 rashomon 4"
  "dgp2 1600 tree 4"
  "dgp2 1600 forest 4"
  "dgp2 800 rashomon 2"
  "dgp2 800 tree 2"
  "dgp2 800 forest 2"
  "dgp2 400 rashomon 1"
  "dgp2 400 forest 1"

  # dgp3 - missing most tree/rashomon
  "dgp3 1600 rashomon 4"
  "dgp3 1600 tree 4"
  "dgp3 800 rashomon 2"
  "dgp3 800 tree 2"
  "dgp3 400 rashomon 1"
  "dgp3 400 tree 1"
)

TOTAL_JOBS=0

for config in "${INCOMPLETE[@]}"; do
  read -r DGP N METHOD N_BATCHES <<< "$config"

  echo "Submitting: ${DGP} n=${N} ${METHOD} (${N_BATCHES} batches)"

  sbatch --export=DGP=${DGP},N=${N},METHOD=${METHOD} \
         --array=1-${N_BATCHES} \
         slurm/run_dml_batch.slurm

  TOTAL_JOBS=$((TOTAL_JOBS + N_BATCHES))
  sleep 0.2
done

echo ""
echo "==========================================="
echo "Resubmission complete!"
echo "Total jobs submitted: ${TOTAL_JOBS}"
echo ""
echo "Monitor progress:"
echo "  squeue -u \$USER"
echo "  find /n/scratch/users/d/dma12/global-scholars/results/o2_primary/ -name '*.rds' | wc -l"
echo "==========================================="
