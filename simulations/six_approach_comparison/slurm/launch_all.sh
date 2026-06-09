#!/bin/bash

# Master launch script for six-approach comparison study
# Submits all 6624 jobs across 6 approaches (11 sbatch calls due to chunking)
#
# Total: 6 approaches × 12 configs × 1000 reps = 72,000 replications

echo "=============================================="
echo "Six-Approach Comparison Study"
echo "=============================================="
echo ""
echo "Total jobs:  6624"
echo "  Approach 1 (fullsample):         24 jobs (1 array)"
echo "  Approach 2 (crossfit):         1200 jobs (2 chunks)"
echo "  Approach 3 (doubletree):       1200 jobs (2 chunks)"
echo "  Approach 4 (dt_averaged):      1200 jobs (2 chunks)"
echo "  Approach 5 (msplit):           2400 jobs (3 chunks)"
echo "  Approach 6 (msplit_averaged):   600 jobs (1 array)"
echo ""
echo "DGPs: 4 (simple, moderate, complex, continuous)"
echo "Sample sizes: 3 (500, 1000, 2000)"
echo "Reps per config: 1000"
echo "Total replications: 72,000"
echo ""

# Create output directories
mkdir -p results/raw logs results/combined results/plots

# Check we're in the right directory
if [ ! -f "code/run_single_replication.R" ]; then
  echo "Error: Must run from six_approach_comparison/ directory"
  exit 1
fi

submit_chunk() {
  local script="$1"
  local array_range="$2"
  local label="$3"
  echo "  Submitting $label (array $array_range)..."
  JOB=$(sbatch --parsable --array="$array_range" "$script")
  if [ $? -eq 0 ]; then
    echo "    Job ID: $JOB"
    echo "$JOB"
  else
    echo "    ERROR: Failed to submit $label"
    exit 1
  fi
}

echo "Approach 1 (fullsample) — 24 jobs..."
A1=$(submit_chunk slurm/run_approach1.sh "1-24" "approach1")

echo ""
echo "Approach 2 (crossfit) — 1200 jobs in 2 chunks..."
A2_1=$(submit_chunk slurm/run_approach2.sh "1-600"     "approach2 chunk1")
A2_2=$(submit_chunk slurm/run_approach2.sh "601-1200"  "approach2 chunk2")

echo ""
echo "Approach 3 (doubletree) — 1200 jobs in 2 chunks..."
A3_1=$(submit_chunk slurm/run_approach3.sh "1-600"     "approach3 chunk1")
A3_2=$(submit_chunk slurm/run_approach3.sh "601-1200"  "approach3 chunk2")

echo ""
echo "Approach 4 (dt_averaged) — 1200 jobs in 2 chunks..."
A4_1=$(submit_chunk slurm/run_approach4.sh "1-600"     "approach4 chunk1")
A4_2=$(submit_chunk slurm/run_approach4.sh "601-1200"  "approach4 chunk2")

echo ""
echo "Approach 5 (msplit) — 2400 jobs in 3 chunks..."
A5_1=$(submit_chunk slurm/run_approach5.sh "1-1000"    "approach5 chunk1")
A5_2=$(submit_chunk slurm/run_approach5.sh "1001-2000" "approach5 chunk2")
A5_3=$(submit_chunk slurm/run_approach5.sh "2001-2400" "approach5 chunk3")

echo ""
echo "Approach 6 (msplit_averaged) — 600 jobs..."
A6=$(submit_chunk slurm/run_approach6.sh "1-600" "approach6")

echo ""
echo "=============================================="
echo "All jobs submitted successfully!"
echo "=============================================="
echo ""
echo "Monitor progress:"
echo "  squeue -u \$USER"
echo "  bash slurm/check_progress.sh"
echo ""
echo "When complete, combine results:"
echo "  Rscript code/combine_results.R"
echo ""
echo "Expected output files: 6624 .rds files in results/raw/"
echo "Naming: approach{1-6}_{job_id}.rds"
echo ""
