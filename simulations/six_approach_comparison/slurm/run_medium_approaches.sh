#!/bin/bash
#SBATCH --job-name=medium_approach
#SBATCH --output=logs/medium_%a.out
#SBATCH --error=logs/medium_%a.err
#SBATCH --array=1-480
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

# Array 2: Medium approaches (ii, iii) — batched
# 480 jobs = 2 approaches × 4 DGPs × 3 n × 20 batches
# Each job: 25 replications (500 total per approach×DGP×n combination)
#
# Updated 2026-05-29: Increased memory 8G→16G and time 2.5h→4h
# Updated 2026-06-03: Split into 5 batches of 100 reps each (was 1×500).
#   CV-based approach 2 at n=2000/DGP3 was OOM-killed at 6h with 500 reps.
# Updated 2026-06-03: Doubled to 20 batches of 25 reps; targeting 1-2h/job.
#   Approach 2 n=2000 DGP3 at 100 reps/job was 4.5h (cluster estimate).

module load gcc/14.2.0 R/4.4.2

cd $SLURM_SUBMIT_DIR

JOB_ID=$SLURM_ARRAY_TASK_ID

# Mapping: 480 jobs = 2 approaches × 4 DGPs × 3 n × 20 batches
# approach ∈ {2, 3} (crossfit-separate, doubletree)
# dgp ∈ {1, 2, 3, 4}
# n_idx ∈ {1, 2, 3} → n ∈ {500, 1000, 2000}
# batch ∈ {1..20} → reps 1-25, 26-50, ..., 476-500

# Decode approach (groups of 240)
APPROACH_IDX=$(( (JOB_ID - 1) / 240 + 1 ))
if [ $APPROACH_IDX -eq 1 ]; then
  APPROACH=2
else
  APPROACH=3
fi

# Decode DGP, n, batch within group of 240
REMAINDER=$(( (JOB_ID - 1) % 240 ))
DGP=$(( REMAINDER / 60 + 1 ))
SUBREMAINDER=$(( REMAINDER % 60 ))
N_IDX=$(( SUBREMAINDER / 20 + 1 ))
BATCH=$(( SUBREMAINDER % 20 + 1 ))

# Map n_idx to actual n
if [ $N_IDX -eq 1 ]; then
  N=500
elif [ $N_IDX -eq 2 ]; then
  N=1000
else
  N=2000
fi

# Determine rep range for this batch
REP_START=$(( (BATCH - 1) * 25 + 1 ))
REP_END=$(( BATCH * 25 ))

echo "=========================================="
echo "Medium Approach Job $JOB_ID"
echo "=========================================="
echo "Approach: $APPROACH"
echo "DGP: $DGP"
echo "Sample size: $N"
echo "Batch: $BATCH / 20"
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
  --output results/raw/medium_approach_${JOB_ID}.rds

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Job $JOB_ID complete"
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "=========================================="

exit $EXIT_CODE
