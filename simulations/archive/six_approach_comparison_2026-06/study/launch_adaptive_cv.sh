#!/bin/bash

# Launch Adaptive CV Validation Simulation on O2
#
# Tests impact of adaptive CV on bias and coverage in complex DGP
# Compares to baseline results where:
#   - n=1000: Bias -0.020, Coverage 85%
#   - n=2000: Bias -0.020, Coverage 88%

set -e

cd "$(dirname "$0")"
mkdir -p logs

echo "==========================================="
echo "Adaptive CV Validation Simulation"
echo "==========================================="
echo ""
echo "Configuration:"
echo "  - DGP: complex (where bias/coverage issues observed)"
echo "  - Sample sizes: 1000, 2000"
echo "  - Method: estimate_att with adaptive CV"
echo "  - Replications: 500 per configuration"
echo "  - Batches: 20 batches of 25 reps each"
echo "  - Resources: 32GB memory per job"
echo ""
echo "Baseline (standard CV):"
echo "  n=1000: Bias -0.020, Coverage 85%"
echo "  n=2000: Bias -0.020, Coverage 88%"
echo ""
echo "Target (adaptive CV):"
echo "  Bias: < 0.010"
echo "  Coverage: 94-96%"
echo ""

TOTAL_JOBS=0

for N in 1000 2000; do
  echo "Submitting: n=${N} (20 batches)"

  sbatch --export=N=${N} --array=1-20 adaptive_cv_batch.slurm

  TOTAL_JOBS=$((TOTAL_JOBS + 20))
  sleep 0.2
done

echo ""
echo "==========================================="
echo "Submission complete!"
echo "Total jobs submitted: ${TOTAL_JOBS}"
echo "Total replications: 1000 (500 per configuration)"
echo ""
echo "Monitor progress:"
echo "  squeue -u \$USER | grep adaptive_cv"
echo "  tail -f logs/adaptive_cv_*.out"
echo ""
echo "Check results:"
echo "  ls -lh /n/scratch/users/\${USER:0:1}/\${USER}/global-scholars/adaptive_cv_validation/"
echo "  # Should have 40 files when complete (2 configs × 20 batches)"
echo ""
echo "Analyze results:"
echo "  cd /n/scratch/users/\${USER:0:1}/\${USER}/global-scholars/adaptive_cv_validation/"
echo "  Rscript ../../doubletree/simulations/six_approach_comparison/analyze_adaptive_cv.R"
echo "==========================================="
