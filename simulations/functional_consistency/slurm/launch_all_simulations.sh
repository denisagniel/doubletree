#!/bin/bash
# Launch all functional consistency simulations on O2 cluster

# Parameter grid
N_VALUES=(200 400 800 1600 3200)
DGP_VALUES=("simple" "complex" "sparse")
METHOD_VALUES=("standard_msplit" "averaged_tree" "pattern_aggregation")
K_VALUES=(2 3 5)
M=10

# Total configs: 5 n × 3 dgp × 3 method × 3 K = 135 configs
# Each config: 50 array tasks × 10 reps = 500 reps
# Total: 67,500 replications

echo "=========================================="
echo "Launching Functional Consistency Simulations"
echo "=========================================="
echo "Total configurations: 135"
echo "Replications per config: 500"
echo "Total replications: 67,500"
echo "=========================================="
echo ""

JOB_COUNT=0

for N in "${N_VALUES[@]}"; do
  for DGP in "${DGP_VALUES[@]}"; do
    for METHOD in "${METHOD_VALUES[@]}"; do
      for K in "${K_VALUES[@]}"; do
        JOB_COUNT=$((JOB_COUNT + 1))

        echo "[${JOB_COUNT}/135] Submitting: n=${N}, dgp=${DGP}, method=${METHOD}, K=${K}, M=${M}"

        JOB_ID=$(sbatch --parsable slurm/run_simulations.slurm ${N} ${DGP} ${METHOD} ${K} ${M})

        echo "  Job ID: ${JOB_ID}"
      done
    done
  done
done

echo ""
echo "=========================================="
echo "All jobs submitted!"
echo "Total jobs: ${JOB_COUNT}"
echo "=========================================="
echo ""
echo "Monitor progress with:"
echo "  bash slurm/check_progress.sh"
echo ""
echo "After completion, combine results with:"
echo "  Rscript slurm/combine_results.R"
