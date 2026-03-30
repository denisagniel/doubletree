#!/bin/bash

# Resubmit DGP5 Forest and Linear Jobs (Fixed for Continuous Outcomes)
#
# After fixing baseline methods to support continuous outcomes, this script
# resubmits ONLY the DGP5 forest and linear jobs that failed with
# "Y must be binary (0/1)" error.
#
# Total: 2 methods × 3 sample sizes × variable batches = 14 array jobs
# Total replications: 6,000

set -e

cd "$(dirname "$0")/.."
mkdir -p logs

echo "==========================================="
echo "DGP5 Forest/Linear Resubmission"
echo "==========================================="
echo ""
echo "Configuration:"
echo "  - 1 DGP: dgp5 (continuous outcome)"
echo "  - 3 Sample sizes: 400, 800, 1600"
echo "  - 2 Methods: forest, linear (FIXED)"
echo "  - 1000 replications per configuration"
echo "  - 14 total batched array jobs"
echo ""
echo "Fix applied:"
echo "  - Forest: Now uses ranger regression mode"
echo "  - Linear: Now uses lm() instead of glm()"
echo ""

TOTAL_JOBS=0

# Loop over DGP5 only
for N in 400 800 1600; do
  for METHOD in forest linear; do

    # Calculate number of batches based on sample size
    if [ "$N" -eq 400 ]; then
      N_BATCHES=1
    elif [ "$N" -eq 800 ]; then
      N_BATCHES=2
    elif [ "$N" -eq 1600 ]; then
      N_BATCHES=4
    fi

    echo "Submitting: dgp5 n=${N} ${METHOD} (${N_BATCHES} batches)"

    sbatch --export=DGP=dgp5,N=${N},METHOD=${METHOD} \
           --array=1-${N_BATCHES} \
           slurm/run_dml_batch.slurm

    TOTAL_JOBS=$((TOTAL_JOBS + N_BATCHES))
    sleep 0.2

  done
done

echo ""
echo "==========================================="
echo "Submission complete!"
echo "Total jobs submitted: ${TOTAL_JOBS}"
echo ""
echo "Monitor progress:"
echo "  squeue -u \$USER"
echo "  tail -f logs/dml_batch_*.out"
echo ""
echo "Check results:"
echo "  find /n/scratch/users/\${USER:0:1}/\${USER}/global-scholars/results/o2_primary/ -name 'dgp5*.rds' | wc -l"
echo "  # Should reach 6,000 when complete (3000 forest + 3000 linear)"
echo "==========================================="
