#!/bin/bash

# Check progress of DML-ATT simulations on O2
#
# Usage:
#   bash slurm/check_progress.sh

# Output directory on scratch
OUTPUT_DIR="/n/scratch/users/${USER:0:1}/${USER}/global-scholars/results/o2_primary"

echo "=========================================="
echo "DML-ATT Simulation Progress"
echo "=========================================="
echo ""

# Check if output directory exists
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Output directory does not exist: $OUTPUT_DIR"
    echo "No simulations have started yet."
    exit 0
fi

# Count completed replications by configuration
echo "Completed replications by configuration:"
echo ""
printf "%-10s %-10s %-10s %-10s\n" "DGP" "N" "METHOD" "COUNT"
echo "----------------------------------------"

TOTAL_COMPLETED=0

for DGP in dgp1 dgp2 dgp3; do
  for N in 400 800 1600; do
    for METHOD in tree rashomon forest linear; do
      COUNT=$(find "$OUTPUT_DIR" -name "${DGP}_n${N}_${METHOD}_rep*.rds" 2>/dev/null | wc -l)
      printf "%-10s %-10s %-10s %-10s\n" "$DGP" "$N" "$METHOD" "$COUNT/500"
      TOTAL_COMPLETED=$((TOTAL_COMPLETED + COUNT))
    done
  done
done

echo "----------------------------------------"
echo "Total completed: $TOTAL_COMPLETED / 18000"
echo ""

# Calculate percentage
PERCENT=$(awk "BEGIN {printf \"%.1f\", ($TOTAL_COMPLETED / 18000) * 100}")
echo "Progress: $PERCENT%"
echo ""

# Check running jobs
echo "Running SLURM jobs:"
RUNNING=$(squeue -u $USER -n dml_sim | wc -l)
RUNNING=$((RUNNING - 1))  # Subtract header line

if [ $RUNNING -gt 0 ]; then
    echo "  $RUNNING tasks running"
    echo ""
    echo "Job details:"
    squeue -u $USER -n dml_sim -o "%.18i %.9P %.30j %.8T %.10M %.6D"
else
    echo "  No running jobs"
fi

echo ""

# Check for failed jobs (exit code != 0 in log files)
FAILED=$(grep -l "exit code 1" logs/dml_*.err 2>/dev/null | wc -l)
if [ $FAILED -gt 0 ]; then
    echo "WARNING: $FAILED failed tasks detected in log files"
    echo "Check logs/dml_*.err for details"
    echo ""
fi

# Estimate time remaining
if [ $TOTAL_COMPLETED -gt 0 ] && [ $RUNNING -gt 0 ]; then
    # Rough estimate: assume average 2 minutes per replication
    REMAINING=$((18000 - TOTAL_COMPLETED))
    EST_MINUTES=$(awk "BEGIN {printf \"%.0f\", ($REMAINING / $RUNNING) * 2}")
    echo "Estimated time remaining: ~$EST_MINUTES minutes"
    echo "(assuming $RUNNING parallel tasks, 2 min/task)"
fi

echo ""
echo "=========================================="
