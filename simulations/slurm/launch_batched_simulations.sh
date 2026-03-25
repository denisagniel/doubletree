#!/bin/bash

# Launch DML-ATT Batched Simulations on O2
#
# Submits 36 configurations × variable batches = 420 total array jobs
#
# Configurations:
#   3 DGPs: dgp1, dgp2, dgp3
#   3 Sample sizes: 400, 800, 1600
#   4 Methods: tree, rashomon, forest, linear
#
# Batching strategy (target: 5 minutes per batch):
#   N=400:  100 reps/batch × 5 batches = 500 reps (60 jobs)
#   N=800:  50 reps/batch × 10 batches = 500 reps (120 jobs)
#   N=1600: 25 reps/batch × 20 batches = 500 reps (240 jobs)
#   Total: 420 array jobs
#
# Total replications: 36 configs × 500 reps = 18,000 replications

set -e  # Exit on error

# Navigate to simulations directory
cd "$(dirname "$0")/.."

# Create logs directory if it doesn't exist
mkdir -p logs

echo "==========================================="
echo "DML-ATT Batched Simulations Launcher"
echo "==========================================="
echo ""
echo "Configuration:"
echo "  - 36 total configurations (3 DGPs × 3 sample sizes × 4 methods)"
echo "  - 500 replications per configuration"
echo "  - 420 total batched array jobs"
echo "  - Target: 5 minutes per batch"
echo ""

# Counter for submitted jobs
TOTAL_JOBS=0

# Loop over all configurations
for DGP in dgp1 dgp2 dgp3; do
  for N in 400 800 1600; do
    for METHOD in tree rashomon forest linear; do

      # Calculate number of batches based on sample size
      if [ "$N" -eq 400 ]; then
        N_BATCHES=5
      elif [ "$N" -eq 800 ]; then
        N_BATCHES=10
      elif [ "$N" -eq 1600 ]; then
        N_BATCHES=20
      fi

      echo "Submitting: ${DGP} n=${N} ${METHOD} (${N_BATCHES} batches)"

      # Submit array job
      sbatch --export=DGP=${DGP},N=${N},METHOD=${METHOD} \
             --array=1-${N_BATCHES} \
             slurm/run_dml_batch.slurm

      # Increment counter
      TOTAL_JOBS=$((TOTAL_JOBS + N_BATCHES))

      # Small delay to avoid overwhelming scheduler
      sleep 0.2

    done
  done
done

echo ""
echo "==========================================="
echo "Submission complete!"
echo "Total batched jobs submitted: ${TOTAL_JOBS}"
echo ""
echo "Monitor progress:"
echo "  squeue -u \$USER"
echo "  tail -f logs/dml_batch_*.out"
echo ""
echo "Check results:"
echo "  ls -lh /n/scratch/users/\${USER:0:1}/\${USER}/global-scholars/results/o2_primary/"
echo "  find /n/scratch/users/\${USER:0:1}/\${USER}/global-scholars/results/o2_primary/ -name '*.rds' | wc -l"
echo "==========================================="
