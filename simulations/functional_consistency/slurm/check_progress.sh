#!/bin/bash
# Check progress of functional consistency simulations

echo "=========================================="
echo "Functional Consistency Simulation Progress"
echo "=========================================="
echo ""

# Count jobs by status
PENDING=$(squeue -u $USER -n fc_sim -t PENDING -h | wc -l)
RUNNING=$(squeue -u $USER -n fc_sim -t RUNNING -h | wc -l)
TOTAL_JOBS=135

# Count completed results
RESULT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/results"
if [ -d "${RESULT_DIR}" ]; then
  COMPLETED=$(find ${RESULT_DIR} -name "*.rds" | wc -l)
  EXPECTED_FILES=$((TOTAL_JOBS * 50 * 10))  # 135 configs × 50 tasks × 10 reps
  PROGRESS=$((100 * COMPLETED / EXPECTED_FILES))
else
  COMPLETED=0
  EXPECTED_FILES=67500
  PROGRESS=0
fi

echo "Job Status:"
echo "  Running:  ${RUNNING}"
echo "  Pending:  ${PENDING}"
echo "  Total:    ${TOTAL_JOBS} configurations"
echo ""
echo "Results:"
echo "  Completed: ${COMPLETED} / ${EXPECTED_FILES} replications (${PROGRESS}%)"
echo ""

if [ $RUNNING -eq 0 ] && [ $PENDING -eq 0 ]; then
  echo "=========================================="
  echo "All jobs complete!"
  echo "=========================================="
  echo ""
  echo "Combine results with:"
  echo "  Rscript slurm/combine_results.R"
else
  # Estimate time remaining
  if [ $PROGRESS -gt 0 ]; then
    echo "Estimated time remaining: TBD (monitor manually)"
  fi
fi
