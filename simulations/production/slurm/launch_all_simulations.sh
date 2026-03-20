#!/bin/bash

# Launch all 36 DML-ATT simulation configurations on O2
#
# Grid: 3 DGPs × 3 sample sizes × 4 methods = 36 configurations
# Each configuration: 500 replications (array job)
# Total: 18,000 simulations
#
# Estimated time: 30-60 minutes (with parallelization)

set -e  # Exit on error

# Create logs directory
mkdir -p logs

# Track job IDs for dependency management
JOB_IDS=()

echo "=========================================="
echo "Launching DML-ATT Primary Simulations"
echo "=========================================="
echo ""
echo "Configuration:"
echo "  DGPs: dgp1 (binary), dgp2 (continuous), dgp3 (moderate)"
echo "  Sample sizes: 400, 800, 1600"
echo "  Methods: tree, rashomon, forest, linear"
echo "  Replications per config: 500"
echo "  Total simulations: 18,000"
echo ""
echo "O2 Settings:"
echo "  Memory: 6G per task"
echo "  Time limit: 1 hour per task"
echo "  Partition: short"
echo ""

# Counter for launched jobs
CONFIG_NUM=0

# Loop over all configurations
for DGP in dgp1 dgp2 dgp3; do
  for N in 400 800 1600; do
    for METHOD in tree rashomon forest linear; do
      CONFIG_NUM=$((CONFIG_NUM + 1))

      echo "[$CONFIG_NUM/36] Submitting: DGP=$DGP, N=$N, METHOD=$METHOD"

      # Submit job array
      JOB_OUTPUT=$(sbatch --export=DGP=$DGP,N=$N,METHOD=$METHOD \
                          slurm/run_dml_simulations.slurm)

      # Extract job ID
      JOB_ID=$(echo $JOB_OUTPUT | awk '{print $4}')
      JOB_IDS+=($JOB_ID)

      echo "  Job ID: $JOB_ID (500 tasks)"

      # Brief pause to avoid overwhelming scheduler
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
echo "Monitor progress:"
echo "  squeue -u \$USER"
echo "  watch -n 10 'squeue -u \$USER | grep dml_sim'"
echo ""
echo "Check status:"
echo "  bash slurm/check_progress.sh"
echo ""
echo "Cancel all:"
echo "  scancel ${JOB_IDS[@]}"
echo ""
echo "Expected completion: 30-60 minutes"
echo ""
