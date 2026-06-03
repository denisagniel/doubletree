#!/bin/bash
#SBATCH --job-name=msplit_approach
#SBATCH --output=logs/msplit_%a.out
#SBATCH --error=logs/msplit_%a.err
#SBATCH --array=1-600
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

# Array 3: M-split approach (v)
# 600 jobs = 4 DGPs × 3 n × 50 batches
# Each job: 10 replications (total 500 per DGP×n)
#
# Updated 2026-05-29: Increased memory 12G→16G and time 2.5h→4h
# Updated 2026-06-03: Increased 4h→6h; DGP3 n=500 observed at ~6h
# Updated 2026-06-03: 50 batches of 10 reps; 6h wall time (conservative).
#   DGP3 n=1000: ~352 sec/rep × 10 = 58 min; n=2000 extrapolated ~93 min.

module load gcc/14.2.0 R/4.4.2

cd $SLURM_SUBMIT_DIR

JOB_ID=$SLURM_ARRAY_TASK_ID

# Mapping: 600 jobs = 4 DGPs × 3 n × 50 batches
APPROACH=5  # M-split

# Decode DGP (groups of 150)
DGP=$(( (JOB_ID - 1) / 150 + 1 ))

# Decode n and batch within group of 150
REMAINDER=$(( (JOB_ID - 1) % 150 ))
N_IDX=$(( REMAINDER / 50 + 1 ))
BATCH=$(( REMAINDER % 50 + 1 ))

# Map n_idx to actual n
if [ $N_IDX -eq 1 ]; then
  N=500
elif [ $N_IDX -eq 2 ]; then
  N=1000
else
  N=2000
fi

# Determine rep range for this batch
REP_START=$(( (BATCH - 1) * 10 + 1 ))
REP_END=$(( BATCH * 10 ))

echo "=========================================="
echo "M-Split Approach Job $JOB_ID"
echo "=========================================="
echo "Approach: $APPROACH (M-split)"
echo "DGP: $DGP"
echo "Sample size: $N"
echo "Batch: $BATCH / 50"
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
