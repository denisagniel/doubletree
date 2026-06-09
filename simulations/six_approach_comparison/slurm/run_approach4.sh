#!/bin/bash
#SBATCH --job-name=approach4
#SBATCH --output=logs/approach4_%a.out
#SBATCH --error=logs/approach4_%a.err
#SBATCH --array=1-900
# DGPs 1-3 only (9 configs × 100 batches = 900 jobs); 1 chunk fits SLURM max array.
# DGP 4 (continuous) runs separately via run_approach4_dgp4.sh (1 rep/batch, 8h wall).
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

# Approach 4: Doubletree Averaged (Rashomon intersection + averaged leaf values)
# 900 jobs = 9 configs × 100 batches; 10 reps per batch = 1000 reps per config
# 9 configs = DGPs 1-3 × n values (500, 1000, 2000)
#
# DGP4 (continuous) uses a separate script because auto-tune failures on continuous
# features exhaust all tiers (~50 GOSDT calls × ~160s/call = ~8000s/rep), making
# 10-rep batches exceed the 6h wall limit.
#
# auto_tune_intersecting=TRUE (via estimate_att_doubletree_averaged default)
#
# Module versions authoritative source: .claude/skills/setup-cluster-simulations/SKILL.md
module load gcc/14.2.0 R/4.4.2

cd $SLURM_SUBMIT_DIR

JOB_ID=$SLURM_ARRAY_TASK_ID
APPROACH=4
N_BATCHES=100
REPS_PER_BATCH=10

CONFIG_IDX=$(( (JOB_ID - 1) / N_BATCHES + 1 ))  # 1..9
BATCH=$(( (JOB_ID - 1) % N_BATCHES + 1 ))        # 1..100

DGP=$(( (CONFIG_IDX - 1) / 3 + 1 ))   # 1..3 (DGPs 1-3 only)
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
