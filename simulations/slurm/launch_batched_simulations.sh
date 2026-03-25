#!/bin/bash

# Launch DML-ATT Batched Simulations on O2
#
# Submits 36 configurations × variable batches = 84 total array jobs
#
# Configurations:
#   3 DGPs: dgp1, dgp2, dgp3
#   3 Sample sizes: 400, 800, 1600
#   4 Methods: tree, rashomon, forest, linear
#
# Batching strategy (target: 5-10 minutes per batch, 1000 reps per config):
#   N=400:  1000 reps/batch × 1 batch = 1000 reps (12 jobs)
#   N=800:  500 reps/batch × 2 batches = 1000 reps (24 jobs)
#   N=1600: 250 reps/batch × 4 batches = 1000 reps (48 jobs)
#   Total: 84 array jobs
#
# Total replications: 36 configs × 1000 reps = 36,000 replications

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
echo "  - 1000 replications per configuration"
echo "  - 84 total batched array jobs"
echo "  - Target: 5-10 minutes per batch"
echo ""

# Counter for submitted jobs
TOTAL_JOBS=0

# Loop over all configurations
for DGP in dgp1 dgp2 dgp3; do
  for N in 400 800 1600; do
    for METHOD in tree rashomon forest linear; do

      # Calculate number of batches based on sample size
      # 1000 reps per config: N=400 (1 batch), N=800 (2 batches), N=1600 (4 batches)
      if [ "$N" -eq 400 ]; then
        N_BATCHES=1
      elif [ "$N" -eq 800 ]; then
        N_BATCHES=2
      elif [ "$N" -eq 1600 ]; then
        N_BATCHES=4
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
