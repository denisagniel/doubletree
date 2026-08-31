#!/bin/bash
#SBATCH --job-name=adaptive_cv
#SBATCH --output=logs/adaptive_cv_%A_%a.out
#SBATCH --error=logs/adaptive_cv_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --partition=short

# Adaptive CV Validation Simulation for O2 Cluster
#
# Usage:
#   sbatch --export=N=2000 --array=1-20 adaptive_cv_batch.slurm
#
# Required environment variable:
#   N : sample size (1000 or 2000)
#
# Array job size: 20 batches of 25 reps each = 500 total reps per config
#
# Resources:
# - mem=32G: Adaptive CV is very memory-intensive (nested CV + tree fitting)
# - batch_size=25: Reduced from 50 to further reduce memory pressure
# - time=01:00:00: Adaptive CV takes longer (nested iterations)
# - partition=short: 12 hour limit, faster scheduling

# Module versions - see .claude/skills/setup-cluster-simulations/SKILL.md
module load gcc/14.2.0
module load R/4.4.2

# Validate environment variables
if [ -z "$N" ]; then
    echo "ERROR: N environment variable not set"
    echo "Usage: sbatch --export=N=2000 --array=1-20 adaptive_cv_batch.slurm"
    exit 1
fi

# Batch configuration
BATCH_SIZE=25  # 25 reps per batch (reduced for memory)
BATCH_START=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH_SIZE + 1 ))
SEED_OFFSET=50000  # Different from main simulations

echo "================================================================"
echo "Adaptive CV Validation Batch"
echo "================================================================"
echo "Configuration: Complex DGP, N=${N}"
echo "Batch: ${SLURM_ARRAY_TASK_ID}, Reps: ${BATCH_START}-$(( BATCH_START + BATCH_SIZE - 1 ))"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Node: ${HOSTNAME}"
echo ""

# Verify R loaded successfully
if ! command -v R &> /dev/null; then
    echo "ERROR: R module failed to load"
    exit 1
fi

echo "R version: $(R --version | head -1)"

# Set R library path to user library
export R_LIBS_USER="${HOME}/R/library"

# Navigate to simulation directory
cd $SLURM_SUBMIT_DIR

# Set output directory to scratch
OUTPUT_DIR="/n/scratch/users/${USER:0:1}/${USER}/global-scholars/adaptive_cv_validation"
mkdir -p "$OUTPUT_DIR"

echo "Output directory: $OUTPUT_DIR"
echo "Starting batch at $(date)"
echo ""

# Run batch of replications
Rscript adaptive_cv_worker.R \
  "$N" \
  "$BATCH_START" \
  "$BATCH_SIZE" \
  "$OUTPUT_DIR" \
  "$SEED_OFFSET"

# Capture exit code
EXIT_CODE=$?

echo ""
echo "Completed at $(date) with exit code $EXIT_CODE"

# Exit with script's exit code
exit $EXIT_CODE
