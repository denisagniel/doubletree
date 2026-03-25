#!/bin/bash
#SBATCH --job-name=dml_batch
#SBATCH --output=logs/dml_batch_%A_%a.out
#SBATCH --error=logs/dml_batch_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=00:15:00
#SBATCH --partition=short

# DML-ATT Batched Simulations for O2 Cluster
#
# Usage:
#   sbatch --export=DGP=dgp1,N=400,METHOD=tree --array=1 slurm/run_dml_batch.slurm
#
# Required environment variables:
#   DGP     : dgp1, dgp2, or dgp3
#   N       : sample size (400, 800, or 1600)
#   METHOD  : tree, rashomon, forest, or linear
#
# Array job size depends on sample size (1000 reps total per config):
#   N=400:  --array=1     (1000 reps per batch, 1 batch, ~5 min)
#   N=800:  --array=1-2   (500 reps per batch, 2 batches, ~5 min each)
#   N=1600: --array=1-4   (250 reps per batch, 4 batches, ~5 min each)
#
# Target: 5-10 minutes per batch (15 minute limit provides safety margin)
#
# O2 Configuration:
# - partition=short: 12 hour limit, faster scheduling
# - mem=6G: Conservative for batched DML simulations
# - cpus-per-task=1: No inner parallelization (array handles parallelism)
# - time=00:15:00: 5-minute target with safety margin

# Validate environment variables
if [ -z "$DGP" ]; then
    echo "ERROR: DGP environment variable not set"
    echo "Usage: sbatch --export=DGP=dgp1,N=400,METHOD=tree --array=1-5 slurm/run_dml_batch.slurm"
    exit 1
fi

if [ -z "$N" ]; then
    echo "ERROR: N environment variable not set"
    exit 1
fi

if [ -z "$METHOD" ]; then
    echo "ERROR: METHOD environment variable not set"
    exit 1
fi

# Calculate batch size based on sample size
# Target: 1000 reps per config, 5-10 minutes per batch
if [ "$N" -eq 400 ]; then
    BATCH_SIZE=1000  # 1 batch per config
elif [ "$N" -eq 800 ]; then
    BATCH_SIZE=500   # 2 batches per config
elif [ "$N" -eq 1600 ]; then
    BATCH_SIZE=250   # 4 batches per config
else
    echo "ERROR: Invalid sample size N=${N} (must be 400, 800, or 1600)"
    exit 1
fi

# Calculate starting replication for this batch
# N=400, Batch 1: reps 1-1000
# N=800, Batch 1: reps 1-500, Batch 2: reps 501-1000
# N=1600, Batch 1: reps 1-250, Batch 2: reps 251-500, etc.
BATCH_START=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH_SIZE + 1 ))

echo "Configuration: DGP=${DGP}, N=${N}, METHOD=${METHOD}"
echo "Batch: ${SLURM_ARRAY_TASK_ID}, Reps: ${BATCH_START}-$(( BATCH_START + BATCH_SIZE - 1 ))"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Node: ${HOSTNAME}"

# Load gcc first (required for R on O2)
module load gcc/14.2.0

# Load R module (try multiple versions)
for R_VERSION in R/4.4.2 R/4.4.1 R/4.3.1 R/4.2.2 R; do
    if module load $R_VERSION 2>/dev/null; then
        echo "Loaded $R_VERSION"
        break
    fi
done

# Verify R loaded successfully
if ! command -v R &> /dev/null; then
    echo "ERROR: R module failed to load"
    echo "Available versions: module spider R"
    exit 1
fi

echo "R version: $(R --version | head -1)"

# Set R library path to user library
export R_LIBS_USER="${HOME}/R/library"

# Navigate to project directory
cd $SLURM_SUBMIT_DIR

# Set output directory to scratch (not home directory!)
# This prevents filling up home directory quota
OUTPUT_DIR="/n/scratch/users/${USER:0:1}/${USER}/global-scholars/results/o2_primary"
mkdir -p "$OUTPUT_DIR"

echo "Output directory: $OUTPUT_DIR"
echo "Starting batch at $(date)"

# Run batch of replications
Rscript run_batch_replications.R \
  --dgp "$DGP" \
  --sample-size "$N" \
  --method "$METHOD" \
  --batch-start "$BATCH_START" \
  --batch-size "$BATCH_SIZE" \
  --output-dir "$OUTPUT_DIR" \
  --tau 0.10 \
  --k-folds 5 \
  --seed-offset 10000 \
  --worker-limit 1

# Capture exit code
EXIT_CODE=$?

echo "Completed at $(date) with exit code $EXIT_CODE"

# Exit with script's exit code
exit $EXIT_CODE
