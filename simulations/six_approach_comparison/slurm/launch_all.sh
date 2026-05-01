#!/bin/bash

# Master launch script for six-approach comparison study
# Submits all 120 jobs across 3 arrays

echo "=============================================="
echo "Six-Approach Comparison Study"
echo "=============================================="
echo ""
echo "Total jobs: 120"
echo "  Array 1 (Fast):   36 jobs (approaches i, iv, vi)"
echo "  Array 2 (Medium): 24 jobs (approaches ii, iii)"
echo "  Array 3 (M-split): 60 jobs (approach v)"
echo ""
echo "DGPs: 4 (simple, moderate, complex, continuous)"
echo "Sample sizes: 3 (500, 1000, 2000)"
echo "Reps per setting: 500"
echo "Total estimations: 36,000"
echo ""
echo "Expected wall time:"
echo "  120 cores: ~2 hours"
echo "  60 cores: ~4 hours"
echo "  40 cores: ~6 hours"
echo ""

# Create output directories
mkdir -p results/raw
mkdir -p logs
mkdir -p results/combined
mkdir -p results/plots

# Check if we're in the right directory
if [ ! -f "code/run_single_replication.R" ]; then
  echo "Error: Must run from six_approach_comparison/ directory"
  exit 1
fi

# Submit job arrays
echo "Submitting Array 1: Fast approaches (i, iv, vi) - 36 jobs..."
FAST_JOB=$(sbatch --parsable slurm/run_fast_approaches.sh)
if [ $? -eq 0 ]; then
  echo "  Job ID: $FAST_JOB"
else
  echo "  ERROR: Failed to submit fast_approaches"
  exit 1
fi

echo ""
echo "Submitting Array 2: Medium approaches (ii, iii) - 24 jobs..."
MEDIUM_JOB=$(sbatch --parsable slurm/run_medium_approaches.sh)
if [ $? -eq 0 ]; then
  echo "  Job ID: $MEDIUM_JOB"
else
  echo "  ERROR: Failed to submit medium_approaches"
  exit 1
fi

echo ""
echo "Submitting Array 3: M-split approach (v) - 60 jobs..."
MSPLIT_JOB=$(sbatch --parsable slurm/run_msplit_approach.sh)
if [ $? -eq 0 ]; then
  echo "  Job ID: $MSPLIT_JOB"
else
  echo "  ERROR: Failed to submit msplit_approach"
  exit 1
fi

echo ""
echo "=============================================="
echo "All jobs submitted successfully!"
echo "=============================================="
echo ""
echo "Job IDs:"
echo "  Fast approaches:   $FAST_JOB"
echo "  Medium approaches: $MEDIUM_JOB"
echo "  M-split approach:  $MSPLIT_JOB"
echo ""
echo "Monitor progress with:"
echo "  bash slurm/check_progress.sh"
echo ""
echo "Or manually:"
echo "  squeue -u \$USER"
echo "  sacct -j $FAST_JOB,$MEDIUM_JOB,$MSPLIT_JOB --format=JobID,JobName,State,Elapsed"
echo ""
echo "When all jobs complete, combine results:"
echo "  Rscript code/combine_results.R"
echo ""
