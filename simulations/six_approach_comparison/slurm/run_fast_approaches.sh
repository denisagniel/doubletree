#!/bin/bash
#SBATCH --job-name=fast_approach
#SBATCH --output=logs/fast_%a.out
#SBATCH --error=logs/fast_%a.err
#SBATCH --array=1-36
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=shared

# Array 1: Fast approaches (i, iv, vi)
# 36 jobs = 3 approaches × 4 DGPs × 3 n
# Each job: 500 replications

# Load R module (adjust for your cluster)
module load gcc/9.2.0 R/4.2.1

# Change to simulation directory
cd $SLURM_SUBMIT_DIR

# Parse job array ID to get parameters
JOB_ID=$SLURM_ARRAY_TASK_ID

# Mapping: 36 jobs = 3 approaches × 4 DGPs × 3 n
# approach ∈ {1, 4, 6} (full-sample, doubletree-singlefit, msplit-singlefit)
# dgp ∈ {1, 2, 3, 4} (simple, moderate, complex, continuous)
# n_idx ∈ {1, 2, 3} (500, 1000, 2000)

# Decode approach (groups of 12)
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
echo "Fast Approach Job $JOB_ID"
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
  --output results/raw/fast_approach_${JOB_ID}.rds

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Job $JOB_ID complete"
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "=========================================="

exit $EXIT_CODE
