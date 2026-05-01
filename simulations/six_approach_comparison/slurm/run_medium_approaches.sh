#!/bin/bash
#SBATCH --job-name=medium_approach
#SBATCH --output=logs/medium_%a.out
#SBATCH --error=logs/medium_%a.err
#SBATCH --array=1-24
#SBATCH --time=02:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=shared

# Array 2: Medium approaches (ii, iii)
# 24 jobs = 2 approaches × 4 DGPs × 3 n
# Each job: 500 replications

module load gcc/9.2.0 R/4.2.1

cd $SLURM_SUBMIT_DIR

JOB_ID=$SLURM_ARRAY_TASK_ID

# Mapping: 24 jobs = 2 approaches × 4 DGPs × 3 n
# approach ∈ {2, 3} (crossfit-separate, doubletree)
# dgp ∈ {1, 2, 3, 4}
# n_idx ∈ {1, 2, 3}

# Decode approach (groups of 12)
APPROACH_IDX=$(( (JOB_ID - 1) / 12 + 1 ))
if [ $APPROACH_IDX -eq 1 ]; then
  APPROACH=2
else
  APPROACH=3
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
echo "Medium Approach Job $JOB_ID"
echo "=========================================="
echo "Approach: $APPROACH"
echo "DGP: $DGP"
echo "Sample size: $N"
echo "Replications: 500"
echo "Start time: $(date)"
echo "=========================================="
echo ""

# Run simulation
Rscript code/run_single_replication.R \
  --approach $APPROACH \
  --dgp $DGP \
  --n $N \
  --reps 500 \
  --output results/raw/medium_approach_${JOB_ID}.rds

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Job $JOB_ID complete"
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "=========================================="

exit $EXIT_CODE
