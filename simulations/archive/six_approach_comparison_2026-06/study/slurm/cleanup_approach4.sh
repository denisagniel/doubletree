#!/bin/bash
# Clean up approach 4 (doubletree_averaged) results and logs before rerunning
# Jobs 13-24 from fast_approaches array

cd "$(dirname "$0")/.."

echo "Cleaning up approach 4 (jobs 13-24) logs and results..."

# Remove old logs for approach 4 jobs
for i in {13..24}; do
  rm -f logs/fast_${i}.out logs/fast_${i}.err
  rm -f logs/approach4_rerun_${i}.out logs/approach4_rerun_${i}.err
done

# Remove old results for approach 4 jobs
for i in {13..24}; do
  rm -f results/raw/fast_approach_${i}.rds
done

echo "Cleanup complete."
echo ""
echo "To relaunch approach 4 jobs:"
echo "  cd $(pwd)"
echo "  sbatch slurm/relaunch_approach4.sh"
