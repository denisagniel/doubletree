#!/bin/bash

# Submit all simulation arrays
# SLURM limits arrays to 1000 jobs max, so we split into chunks
# Command-line --array overrides #SBATCH --array in the script

echo "=============================================="
echo "Submitting Six-Approach Comparison Study"
echo "=============================================="
echo ""

cd "$(dirname "$0")/.."

# Fast approaches: 36 jobs (under 1000, no split needed)
echo "Submitting fast approaches (36 jobs)..."
sbatch slurm/run_fast_approaches.sh
echo ""

# Medium approaches: 1200 jobs = 2 chunks of 600
echo "Submitting medium approaches (1200 jobs in 2 chunks)..."
sbatch --array=1-600 slurm/run_medium_approaches.sh
sbatch --array=601-1200 slurm/run_medium_approaches.sh
echo ""

# M-split approach: 1200 jobs = 2 chunks of 600
echo "Submitting msplit approach (1200 jobs in 2 chunks)..."
sbatch --array=1-600 slurm/run_msplit_approach.sh
sbatch --array=601-1200 slurm/run_msplit_approach.sh
echo ""

echo "=============================================="
echo "Total submitted: 2436 jobs"
echo "  Fast:   36  (1 submission)"
echo "  Medium: 1200 (2 submissions of 600)"
echo "  Msplit: 1200 (2 submissions of 600)"
echo "=============================================="
echo ""
echo "Monitor with: bash slurm/check_progress.sh"
