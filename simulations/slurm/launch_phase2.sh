#!/bin/bash

# Launch Phase 2 DGPs (DGP7-9)
#
# Phase 2: Cases where tree outperforms linear
# - DGP7: Deep 3-way interaction
# - DGP8: Double nonlinearity (sin/cos in both e and m0)
# - DGP9: Weak overlap stress test
#
# Total: 36 configurations × variable batches = 84 array jobs
# Total replications: 36,000

set -e

cd "$(dirname "$0")/.."
mkdir -p logs

echo "==========================================="
echo "Phase 2 DGPs Launcher (DGP7-9)"
echo "==========================================="
echo ""
echo "Configuration:"
echo "  - 3 DGPs demonstrating tree advantages:"
echo "    * DGP7: Deep 3-way interaction"
echo "    * DGP8: Double nonlinearity (sin/cos)"
echo "    * DGP9: Weak overlap stress test"
echo "  - 3 Sample sizes: 400, 800, 1600"
echo "  - 4 Methods: tree, rashomon, forest, linear"
echo "  - 1000 replications per configuration"
echo "  - 84 total batched array jobs"
echo ""
echo "Expected outcomes:"
echo "  - DGP7: Tree beats linear (deep interactions)"
echo "  - DGP8: Tree beats linear (double nonlinearity)"
echo "  - DGP9: Both maintain coverage (weak overlap)"
echo ""

TOTAL_JOBS=0

# Loop over Phase 2 DGPs
for DGP in dgp7 dgp8 dgp9; do
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
echo "  find /n/scratch/users/\${USER:0:1}/\${USER}/global-scholars/results/o2_primary/ -name 'dgp[789]*.rds' | wc -l"
echo "  # Should reach 36,000 when complete"
echo ""
echo "Combined Phase 1 + Phase 2:"
echo "  find /n/scratch/users/\${USER:0:1}/\${USER}/global-scholars/results/o2_primary/ -name 'dgp*.rds' | wc -l"
echo "  # Should reach 108,000 when all complete (72k original + 36k Phase 2)"
echo "==========================================="
