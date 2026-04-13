#!/bin/bash
#SBATCH --job-name=fc_sim
#SBATCH --output=logs/fc_sim_%A_%a.out
#SBATCH --error=logs/fc_sim_%A_%a.err
#SBATCH --array=1-50
#SBATCH --time=2:00:00
#SBATCH --mem=4G
#SBATCH --partition=short

# Functional Consistency Simulation - O2 Cluster
# Array job: each task runs 10 replications for one configuration

# Load modules
module load gcc/9.2.0
module load R/4.2.1

# Configuration from command line
N=$1
DGP=$2
METHOD=$3
K=$4
M=${5:-10}  # Default M=10

# Task-specific settings
TASK_ID=$(printf "%03d" $SLURM_ARRAY_TASK_ID)
START_SEED=$(( ($SLURM_ARRAY_TASK_ID - 1) * 10 + 1 ))
END_SEED=$(( $START_SEED + 9 ))

# Directories
SIM_DIR="${HOME}/global-scholars/doubletree/simulations/functional_consistency"
SCRATCH_DIR="/n/scratch/users/${USER:0:1}/${USER}/fc_sim"
RESULT_DIR="${SIM_DIR}/results"

# Create directories
mkdir -p ${SCRATCH_DIR}
mkdir -p ${RESULT_DIR}
mkdir -p ${SIM_DIR}/slurm/logs

# CD to simulation directory
cd ${SIM_DIR}

echo "=========================================="
echo "Functional Consistency Simulation"
echo "=========================================="
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Configuration: n=${N}, dgp=${DGP}, method=${METHOD}, K=${K}, M=${M}"
echo "Seeds: ${START_SEED} to ${END_SEED}"
echo "Scratch: ${SCRATCH_DIR}"
echo "=========================================="

# Run replications
for SEED in $(seq ${START_SEED} ${END_SEED}); do
  OUTPUT_FILE="${SCRATCH_DIR}/batch_${TASK_ID}_seed_${SEED}.rds"

  echo "Running seed ${SEED}..."

  Rscript slurm/run_single_replication.R \
    --n ${N} \
    --dgp ${DGP} \
    --method ${METHOD} \
    --K ${K} \
    --M ${M} \
    --seed ${SEED} \
    --output ${OUTPUT_FILE}

  if [ $? -ne 0 ]; then
    echo "ERROR: Replication failed for seed ${SEED}"
    exit 1
  fi
done

# Copy results to permanent storage
echo "Copying results to ${RESULT_DIR}..."
BATCH_PATTERN="n${N}_dgp${DGP}_method${METHOD}_K${K}_M${M}"
cp ${SCRATCH_DIR}/batch_${TASK_ID}_*.rds ${RESULT_DIR}/${BATCH_PATTERN}_batch${TASK_ID}/

echo "Task complete!"
