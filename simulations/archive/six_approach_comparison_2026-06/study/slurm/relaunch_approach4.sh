#!/bin/bash
#SBATCH --job-name=approach4_rerun
#SBATCH --output=logs/approach4_rerun_%a.out
#SBATCH --error=logs/approach4_rerun_%a.err
#SBATCH --array=13-24
#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

# Rerun approach 4 (doubletree_averaged) jobs after fixing S7 type error
# These are jobs 13-24 from the fast_approaches array
#
# Original mapping from run_fast_approaches.sh:
# 36 jobs = 3 approaches × 4 DGPs × 3 n
# approach ∈ {1, 4, 6}
# Jobs 13-24 correspond to APPROACH_IDX=2 → APPROACH=4 (doubletree_averaged)

module load gcc/14.2.0 R/4.4.2

cd $SLURM_SUBMIT_DIR

JOB_ID=$SLURM_ARRAY_TASK_ID

# Decode parameters (same logic as run_fast_approaches.sh)
APPROACH_IDX=$(( (JOB_ID - 1) / 12 + 1 ))
if [ $APPROACH_IDX -eq 1 ]; then
  APPROACH=1
elif [ $APPROACH_IDX -eq 2 ]; then
  APPROACH=4
else
  APPROACH=6
fi

# Decode DGP and n within group of 12
REMAINDER=$(( (JOB_ID - 1) % 12 ))
DGP=$(( REMAINDER / 3 + 1 ))
N_IDX=$(( REMAINDER % 3 + 1 ))

# Map n_idx to actual n
if [ $N_IDX -eq 1 ]; then
  N=500
elif [ $N_IDX -eq 2 ]; then
  N=1000
else
  N=2000
fi

echo "=========================================="
echo "Approach 4 Rerun - Job $JOB_ID"
echo "=========================================="
echo "Approach: $APPROACH (doubletree_averaged)"
echo "DGP: $DGP"
echo "Sample size: $N"
echo "Replications: 500"
echo "Start time: $(date)"
echo "=========================================="
echo ""

# Run simulation (overwrite old failed results)
Rscript code/run_single_replication.R \
  --approach $APPROACH \
  --dgp $DGP \
  --n $N \
  --reps 500 \
  --output results/raw/fast_approach_${JOB_ID}.rds

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Job $JOB_ID complete"
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "=========================================="

exit $EXIT_CODE
