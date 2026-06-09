#!/bin/bash
#SBATCH --job-name=approach4_dgp4
#SBATCH --output=logs/approach4_dgp4_%a.out
#SBATCH --error=logs/approach4_dgp4_%a.err
#SBATCH --array=1-1000
# NOTE: Full range is 1-3000; submit in 3 chunks via launch_all.sh:
#   chunk 1: --array=1-1000, chunk 2: --array=1001-2000, chunk 3: --array=2001-3000
# (SLURM max array size = 1000)
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

# Approach 4: Doubletree Averaged -- DGP 4 (continuous features) ONLY
# 3000 jobs = 3 n-values × 1000 batches; 1 rep per batch = 1000 reps per config
#
# WHY SEPARATE: On the continuous DGP, auto-tune exhausts all tiers before
# failing (~50 GOSDT calls × ~160s/call × 2 nuisances = ~16000s/rep at n=500
# locally; ~4000-8000s on cluster). 10-rep batches would exceed the 6h wall.
# With 1 rep/batch and 8h wall, even the n=2000 worst case (~6h on cluster)
# completes cleanly and is logged as a hard-stop error.
#
# Output files named approach4_dgp4_${JOB_ID}.rds to avoid collision with
# the DGPs-1-3 script which uses approach4_${JOB_ID}.rds.
#
# Module versions authoritative source: .claude/skills/setup-cluster-simulations/SKILL.md
module load gcc/14.2.0 R/4.4.2

cd $SLURM_SUBMIT_DIR

JOB_ID=$SLURM_ARRAY_TASK_ID
APPROACH=4
DGP=4
N_BATCHES=1000
REPS_PER_BATCH=1

# CONFIG_IDX: 1 = n=500, 2 = n=1000, 3 = n=2000
CONFIG_IDX=$(( (JOB_ID - 1) / N_BATCHES + 1 ))  # 1..3
BATCH=$(( (JOB_ID - 1) % N_BATCHES + 1 ))        # 1..1000

if [ $CONFIG_IDX -eq 1 ]; then
  N=500
elif [ $CONFIG_IDX -eq 2 ]; then
  N=1000
else
  N=2000
fi

REP_START=$BATCH
REP_END=$BATCH

echo "=========================================="
echo "Approach $APPROACH (DGP4) Job $JOB_ID"
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
  --output results/raw/approach${APPROACH}_dgp4_${JOB_ID}.rds

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Job $JOB_ID complete | exit=$EXIT_CODE | $(date)"
echo "=========================================="

exit $EXIT_CODE
