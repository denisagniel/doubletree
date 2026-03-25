#!/bin/bash

# Launch forest + linear simulations only (tree/rashomon have S7 issues)
#
# Grid: 3 DGPs × 3 sample sizes × 2 methods = 18 configurations
# Each configuration: 500 replications (array job)
# Total: 9,000 simulations
#
# Estimated time: 15-30 minutes (with parallelization)

set -e

mkdir -p logs

JOB_IDS=()

echo "=========================================="
echo "Launching DML-ATT Simulations (Forest + Linear Only)"
echo "=========================================="
echo ""
echo "Configuration:"
echo "  DGPs: dgp1 (binary), dgp2 (continuous), dgp3 (moderate)"
echo "  Sample sizes: 400, 800, 1600"
echo "  Methods: forest, linear (tree/rashomon have S7 issues)"
echo "  Replications per config: 500"
echo "  Total simulations: 9,000"
echo ""

CONFIG_NUM=0

for DGP in dgp1 dgp2 dgp3; do
  for N in 400 800 1600; do
    for METHOD in forest linear; do
      CONFIG_NUM=$((CONFIG_NUM + 1))

      echo "[$CONFIG_NUM/18] Submitting: DGP=$DGP, N=$N, METHOD=$METHOD"

      JOB_OUTPUT=$(sbatch --export=DGP=$DGP,N=$N,METHOD=$METHOD \
                          slurm/run_dml_simulations.slurm)

      JOB_ID=$(echo $JOB_OUTPUT | awk '{print $4}')
      JOB_IDS+=($JOB_ID)

      echo "  Job ID: $JOB_ID (500 tasks)"

      sleep 0.5
    done
  done
done

echo ""
echo "=========================================="
echo "All jobs submitted!"
echo "=========================================="
echo ""
echo "Job IDs: ${JOB_IDS[@]}"
echo ""
echo "Monitor: bash slurm/check_progress.sh"
echo "Expected completion: 15-30 minutes"
echo ""
