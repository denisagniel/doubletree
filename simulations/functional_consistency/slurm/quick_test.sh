#!/bin/bash
# Quick test of simulation infrastructure
# Runs 3 replications of smallest configuration locally

echo "=========================================="
echo "Quick Test: Functional Consistency Simulation"
echo "=========================================="
echo "Running 3 replications locally..."
echo "Configuration: n=200, dgp=simple, method=standard_msplit, K=2, M=10"
echo ""

# Setup
SIM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd ${SIM_DIR}

TEST_DIR="test_output"
mkdir -p ${TEST_DIR}

# Test parameters
N=200
DGP="simple"
METHOD="standard_msplit"
K=2
M=10

# Run 3 replications
for SEED in 1 2 3; do
  echo "Running replication ${SEED}/3..."

  Rscript slurm/run_single_replication.R \
    --n ${N} \
    --dgp ${DGP} \
    --method ${METHOD} \
    --K ${K} \
    --M ${M} \
    --seed ${SEED} \
    --output ${TEST_DIR}/test_seed_${SEED}.rds

  if [ $? -ne 0 ]; then
    echo "ERROR: Test failed!"
    exit 1
  fi
done

echo ""
echo "=========================================="
echo "Test complete!"
echo "=========================================="
echo ""
echo "Results saved to: ${TEST_DIR}/"
echo ""
echo "Checking results..."
Rscript -e "
results <- lapply(1:3, function(i) readRDS('${TEST_DIR}/test_seed_\$i.rds'))
df <- do.call(rbind, results)
print(df)
cat('\n')
cat('All replications successful!\n')
cat('Mean ATT estimate:', mean(df\$att_est), '\n')
cat('Mean coverage:', mean(df\$coverage), '\n')
cat('Max functional consistency diff:', max(df\$max_diff_e, df\$max_diff_m0), '\n')
"

echo ""
echo "If test passed, ready to deploy to O2!"
