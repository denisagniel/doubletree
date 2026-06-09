#!/bin/bash
#SBATCH --job-name=approach3
#SBATCH --output=logs/approach3_%a.out
#SBATCH --error=logs/approach3_%a.err
#SBATCH --array=1-600
# NOTE: Full range is 1-1200; submit in 2 chunks via launch_all.sh
# (SLURM max array size = 1000)
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

# Approach 3: Doubletree (Rashomon intersection, cross-fitted)
# 1200 jobs = 12 configs × 100 batches; 10 reps per batch = 1000 reps per config
# 12 configs = 4 DGPs × 3 n values (500, 1000, 2000)
#
# auto_tune_intersecting=TRUE: starts at eps_n=2*sqrt(log(n)/n), increases if needed.
# Hard failure (exit=1) if intersection still fails after auto-tuning.
#
# Module versions authoritative source: .claude/skills/setup-cluster-simulations/SKILL.md
module load gcc/14.2.0 R/4.4.2

cd $SLURM_SUBMIT_DIR

JOB_ID=$SLURM_ARRAY_TASK_ID
APPROACH=3
N_BATCHES=100
REPS_PER_BATCH=10

CONFIG_IDX=$(( (JOB_ID - 1) / N_BATCHES + 1 ))  # 1..12
BATCH=$(( (JOB_ID - 1) % N_BATCHES + 1 ))        # 1..100

DGP=$(( (CONFIG_IDX - 1) / 3 + 1 ))
N_IDX=$(( (CONFIG_IDX - 1) % 3 + 1 ))

if [ $N_IDX -eq 1 ]; then
  N=500
elif [ $N_IDX -eq 2 ]; then
  N=1000
else
  N=2000
fi

REP_START=$(( (BATCH - 1) * REPS_PER_BATCH + 1 ))
REP_END=$(( BATCH * REPS_PER_BATCH ))

echo "=========================================="
echo "Approach $APPROACH Job $JOB_ID"
echo "=========================================="
echo "Config: $CONFIG_IDX (DGP=$DGP, N=$N)"
echo "Batch: $BATCH / $N_BATCHES"
echo "Replications: $REP_START to $REP_END"
echo "Start time: $(date)"
echo "=========================================="
echo ""

Rscript code/run_single_replication.R \
  --approach $APPROACH \
  --dgp $DGP \
  --n $N \
  --rep_start $REP_START \
  --rep_end $REP_END \
  --output results/raw/approach${APPROACH}_${JOB_ID}.rds

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Job $JOB_ID complete | exit=$EXIT_CODE | $(date)"
echo "=========================================="

exit $EXIT_CODE
