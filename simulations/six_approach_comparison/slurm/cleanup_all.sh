#!/bin/bash
# Clean up all simulation logs and results before full rerun

cd "$(dirname "$0")/.."

echo "WARNING: This will delete all simulation logs and results!"
read -p "Are you sure? (yes/no): " confirm

if [ "$confirm" != "yes" ]; then
  echo "Cleanup cancelled."
  exit 0
fi

echo "Cleaning up all logs and results..."

# Remove all logs
rm -f logs/*.out logs/*.err

# Remove all results
rm -f results/raw/*.rds

# Remove combined results
rm -f results/combined/*.rds

echo "Cleanup complete."
echo ""
echo "To relaunch all simulations:"
echo "  cd $(pwd)"
echo "  sbatch slurm/run_fast_approaches.sh"
echo "  sbatch slurm/run_medium_approaches.sh"
echo "  sbatch slurm/run_msplit_approach.sh"
echo ""
echo "Or use the launch_all.sh script if available."
