#!/bin/bash

# Check progress of six-approach comparison study

echo "=============================================="
echo "Six-Approach Comparison Study Progress"
echo "=============================================="
echo ""

# Count completed jobs by checking for output files
FAST_COMPLETE=$(ls results/raw/fast_approach_*.rds 2>/dev/null | wc -l)
MEDIUM_COMPLETE=$(ls results/raw/medium_approach_*.rds 2>/dev/null | wc -l)
MSPLIT_COMPLETE=$(ls results/raw/msplit_approach_*.rds 2>/dev/null | wc -l)
TOTAL_COMPLETE=$((FAST_COMPLETE + MEDIUM_COMPLETE + MSPLIT_COMPLETE))
TOTAL_EXPECTED=816  # 36 fast + 480 medium (20x25 reps) + 300 msplit (25x20 reps)

echo "Completed jobs: $TOTAL_COMPLETE / $TOTAL_EXPECTED"
echo ""
echo "By array:"
echo "  Fast approaches (i, iv, vi):   $FAST_COMPLETE / 36"
echo "  Medium approaches (ii, iii):   $MEDIUM_COMPLETE / 480  (batched: 20x25 reps each)"
echo "  M-split approach (v):          $MSPLIT_COMPLETE / 300  (batched: 25x20 reps each)"
echo ""

# Calculate percentage
if [ $TOTAL_COMPLETE -gt 0 ]; then
  PERCENT=$((100 * TOTAL_COMPLETE / TOTAL_EXPECTED))
  echo "Progress: $PERCENT% complete"
  echo ""
fi

# Check for errors in log files
ERROR_COUNT=$(grep -l "Error\|ERROR\|FATAL" logs/*.err 2>/dev/null | wc -l)
if [ $ERROR_COUNT -gt 0 ]; then
  echo "⚠️  WARNING: $ERROR_COUNT jobs have errors in stderr"
  echo "Check logs/*.err files for details"
  echo ""
  echo "Example errors:"
  grep -h "Error\|ERROR" logs/*.err 2>/dev/null | head -5
else
  echo "✓ No errors detected in stderr files"
fi
echo ""

# Check for failed jobs (exit code != 0)
FAILED_JOBS=$(grep -l "Exit code: [^0]" logs/*.out 2>/dev/null | wc -l)
if [ $FAILED_JOBS -gt 0 ]; then
  echo "⚠️  WARNING: $FAILED_JOBS jobs failed (non-zero exit code)"
  echo "Check logs/*.out files for details"
else
  echo "✓ No failed jobs detected"
fi
echo ""

# Check running jobs
RUNNING=$(squeue -u $USER -h -n fast_approach,medium_approach,msplit_approach 2>/dev/null | wc -l)
if [ $RUNNING -gt 0 ]; then
  echo "Currently running: $RUNNING jobs"
  echo ""
  squeue -u $USER -n fast_approach,medium_approach,msplit_approach
else
  echo "No jobs currently running"
fi
echo ""

# Final status
if [ $TOTAL_COMPLETE -eq $TOTAL_EXPECTED ]; then
  echo "=============================================="
  echo "✓ All $TOTAL_EXPECTED jobs complete!"
  echo "=============================================="
  echo ""
  echo "Next steps:"
  echo "  1. Check for errors above"
  echo "  2. Combine results:"
  echo "       Rscript code/combine_results.R"
  echo "  3. Analyze results:"
  echo "       Rscript code/analyze_results.R"
  echo ""
elif [ $TOTAL_COMPLETE -gt 0 ] && [ $RUNNING -eq 0 ]; then
  echo "=============================================="
  echo "⚠️  WARNING: Jobs completed but count < $TOTAL_EXPECTED"
  echo "=============================================="
  echo ""
  echo "Expected $TOTAL_EXPECTED, found $TOTAL_COMPLETE"
  echo "Some jobs may have failed. Check logs above."
  echo ""
elif [ $RUNNING -gt 0 ]; then
  # Estimate time remaining
  if [ $TOTAL_COMPLETE -gt 5 ]; then
    echo "Estimated completion: Check back in 30 minutes"
  fi
fi
