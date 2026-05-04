#!/bin/bash
#SBATCH --job-name=msplit_approach
#SBATCH --output=logs/msplit_%a.out
#SBATCH --error=logs/msplit_%a.err
#SBATCH --array=1-60
#SBATCH --time=02:30:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

# Array 3: M-split approach (v)
# 60 jobs = 4 DGPs × 3 n × 5 batches
# Each job: 100 replications (total 500 per DGP×n)

module load gcc/9.2.0 R/4.2.1

cd $SLURM_SUBMIT_DIR

JOB_ID=$SLURM_ARRAY_TASK_ID

# Mapping: 60 jobs = 4 DGPs × 3 n × 5 batches
APPROACH=5  # M-split

# Decode DGP (groups of 15)
DGP=$(( (JOB_ID - 1) / 15 + 1 ))

# Decode n and batch within group of 15
REMAINDER=$(( (JOB_ID - 1) % 15 ))
N_IDX=$(( REMAINDER / 5 + 1 ))
BATCH=$(( REMAINDER % 5 + 1 ))

# Map n_idx to actual n
if [ $N_IDX -eq 1 ]; then
  N=500
elif [ $N_IDX -eq 2 ]; then
  N=1000
else
  N=2000
fi

# Determine rep range for this batch
REP_START=$(( (BATCH - 1) * 100 + 1 ))
REP_END=$(( BATCH * 100 ))

echo "=========================================="
echo "M-Split Approach Job $JOB_ID"
echo "=========================================="
echo "Approach: $APPROACH (M-split)"
echo "DGP: $DGP"
echo "Sample size: $N"
echo "Batch: $BATCH / 5"
echo "Replications: $REP_START to $REP_END"
echo "Start time: $(date)"
echo "=========================================="
echo ""

# Run simulation
Rscript code/run_single_replication.R \
  --approach $APPROACH \
  --dgp $DGP \
  --n $N \
  --rep_start $REP_START \
  --rep_end $REP_END \
  --output results/raw/msplit_approach_${JOB_ID}.rds

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Job $JOB_ID complete"
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "=========================================="

exit $EXIT_CODE
