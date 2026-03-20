#!/bin/bash
#SBATCH --job-name=dml_sim
#SBATCH --output=logs/dml_%A_%a.out
#SBATCH --error=logs/dml_%A_%a.err
#SBATCH --array=1-500
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --partition=short

# DML-ATT Primary Simulations for O2 Cluster
#
# Usage:
#   sbatch --export=DGP=dgp1,N=400,METHOD=tree slurm/run_dml_simulations.slurm
#
# Required environment variables:
#   DGP     : dgp1, dgp2, or dgp3
#   N       : sample size (400, 800, or 1600)
#   METHOD  : tree, rashomon, forest, or linear
#
# Array job: 500 replications (1-500)
#
# O2 Configuration:
# - partition=short: 12 hour limit, faster scheduling
# - mem=6G: Conservative for DML simulations
# - cpus-per-task=1: No inner parallelization (array handles parallelism)
# - time=01:00:00: Most replications complete in 10-30 minutes

# Validate environment variables
if [ -z "$DGP" ]; then
    echo "ERROR: DGP environment variable not set"
    echo "Usage: sbatch --export=DGP=dgp1,N=400,METHOD=tree slurm/run_dml_simulations.slurm"
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

echo "Configuration: DGP=${DGP}, N=${N}, METHOD=${METHOD}"
echo "Replication: ${SLURM_ARRAY_TASK_ID}"
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
echo "Starting simulation at $(date)"

# Run single replication
Rscript run_single_replication.R \
  --dgp "$DGP" \
  --sample-size "$N" \
  --method "$METHOD" \
  --replication $SLURM_ARRAY_TASK_ID \
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
