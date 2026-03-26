#!/bin/bash

# Launch Phase 1 DGPs Only (DGP4-6)
#
# Submits only the new Phase 1 configurations:
# - DGP4: Continuous features, binary outcome
# - DGP5: Continuous features, continuous outcome
# - DGP6: Mixed features (2 binary + 2 continuous)
#
# Total: 36 configurations × variable batches = 84 array jobs
# Total replications: 36,000

set -e

cd "$(dirname "$0")/.."
mkdir -p logs

echo "==========================================="
echo "Phase 1 DGPs Launcher (DGP4-6)"
echo "==========================================="
echo ""
echo "Configuration:"
echo "  - 3 new DGPs (continuous/mixed features)"
echo "  - 3 Sample sizes: 400, 800, 1600"
echo "  - 4 Methods: tree, rashomon, forest, linear"
echo "  - 1000 replications per configuration"
echo "  - 84 total batched array jobs"
echo ""

TOTAL_JOBS=0

# Loop over Phase 1 DGPs only
for DGP in dgp4 dgp5 dgp6; do
  for N in 400 800 1600; do
    for METHOD in tree rashomon forest linear; do

      # Calculate number of batches based on sample size
      if [ "$N" -eq 400 ]; then
        N_BATCHES=1
      elif [ "$N" -eq 800 ]; then
        N_BATCHES=2
      elif [ "$N" -eq 1600 ]; then
        N_BATCHES=4
      fi

      echo "Submitting: ${DGP} n=${N} ${METHOD} (${N_BATCHES} batches)"

      sbatch --export=DGP=${DGP},N=${N},METHOD=${METHOD} \
             --array=1-${N_BATCHES} \
             slurm/run_dml_batch.slurm

      TOTAL_JOBS=$((TOTAL_JOBS + N_BATCHES))
      sleep 0.2

    done
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
echo "  find /n/scratch/users/\${USER:0:1}/\${USER}/global-scholars/results/o2_primary/ -name 'dgp[456]*.rds' | wc -l"
echo "  # Should reach 36,000 when complete"
echo "==========================================="
