#!/bin/bash
#SBATCH --job-name=approach5
#SBATCH --output=logs/approach5_%a.out
#SBATCH --error=logs/approach5_%a.err
#SBATCH --array=1-1000
# NOTE: Full range is 1-2400; submit in 3 chunks via launch_all.sh:
#   chunk 1: --array=1-1000, chunk 2: --array=1001-2000, chunk 3: --array=2001-2400
# (SLURM max array size = 1000)
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

# Approach 5: M-split Doubletree (modal structure, cross-fitted predictions)
# 2400 jobs = 12 configs × 200 batches; 5 reps per batch = 1000 reps per config
# 12 configs = 4 DGPs × 3 n values (500, 1000, 2000)
#
# Module versions authoritative source: .claude/skills/setup-cluster-simulations/SKILL.md
module load gcc/14.2.0 R/4.4.2

cd $SLURM_SUBMIT_DIR

JOB_ID=$SLURM_ARRAY_TASK_ID
APPROACH=5
N_BATCHES=200
REPS_PER_BATCH=5

CONFIG_IDX=$(( (JOB_ID - 1) / N_BATCHES + 1 ))  # 1..12
BATCH=$(( (JOB_ID - 1) % N_BATCHES + 1 ))        # 1..200

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
