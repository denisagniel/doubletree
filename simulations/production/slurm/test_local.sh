#!/bin/bash

# Test O2 infrastructure locally before deploying
#
# This script simulates the O2 environment and tests the single-replication
# script to catch errors before submitting to the cluster.

set -e

echo "=========================================="
echo "Local Test of O2 Simulation Infrastructure"
echo "=========================================="
echo ""

# Check if we're in the right directory
if [ ! -f "run_single_replication.R" ]; then
    echo "ERROR: Must run from doubletree/simulations/production/ directory"
    exit 1
fi

# Check if R is available
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found in PATH"
    exit 1
fi

echo "R version: $(R --version | head -1)"
echo ""

# Create temporary output directory
OUTPUT_DIR="/tmp/o2_test_$(date +%s)"
mkdir -p "$OUTPUT_DIR"

echo "Test output directory: $OUTPUT_DIR"
echo ""

# Test configurations
CONFIGS=(
  "dgp1 400 tree"
  "dgp2 800 rashomon"
  "dgp3 1600 forest"
)

echo "Running 3 test replications (one per DGP)..."
echo ""

TEST_RESULTS=()

for CONFIG in "${CONFIGS[@]}"; do
  read -r DGP N METHOD <<< "$CONFIG"

  echo "Testing: DGP=$DGP, N=$N, METHOD=$METHOD"

  # Run single replication
  if Rscript run_single_replication.R \
      --dgp "$DGP" \
      --sample-size "$N" \
      --method "$METHOD" \
      --replication 1 \
      --output-dir "$OUTPUT_DIR" \
      --tau 0.10 \
      --k-folds 5 \
      --seed-offset 10000 \
      --worker-limit 1 \
      2>&1 | tee "${OUTPUT_DIR}/test_${DGP}_${METHOD}.log"; then

    echo "  ✓ SUCCESS"
    TEST_RESULTS+=("PASS")

    # Check if output file was created
    EXPECTED_FILE="${OUTPUT_DIR}/${DGP}_n${N}_${METHOD}_rep0001.rds"
    if [ -f "$EXPECTED_FILE" ]; then
      echo "  ✓ Output file created: $(basename $EXPECTED_FILE)"

      # Check file size
      SIZE=$(stat -f%z "$EXPECTED_FILE" 2>/dev/null || stat -c%s "$EXPECTED_FILE" 2>/dev/null)
      echo "  ✓ File size: $SIZE bytes"
    else
      echo "  ✗ WARNING: Expected output file not found"
      TEST_RESULTS+=("FAIL")
    fi

  else
    echo "  ✗ FAILED (exit code: $?)"
    TEST_RESULTS+=("FAIL")
  fi

  echo ""
done

# Summary
echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo ""

PASS_COUNT=0
FAIL_COUNT=0

for RESULT in "${TEST_RESULTS[@]}"; do
  if [ "$RESULT" = "PASS" ]; then
    PASS_COUNT=$((PASS_COUNT + 1))
  else
    FAIL_COUNT=$((FAIL_COUNT + 1))
  fi
done

echo "Passed: $PASS_COUNT / ${#CONFIGS[@]}"
echo "Failed: $FAIL_COUNT / ${#CONFIGS[@]}"
echo ""

if [ $FAIL_COUNT -eq 0 ]; then
  echo "✓ All tests passed!"
  echo ""
  echo "Next steps:"
  echo "  1. Review test logs in: $OUTPUT_DIR"
  echo "  2. Transfer code to O2: rsync -avz doubletree/ o2:path/to/doubletree/"
  echo "  3. Run on O2: bash slurm/launch_all_simulations.sh"
  echo ""
  echo "Keep test output? (y/n)"
  read -r KEEP
  if [ "$KEEP" != "y" ]; then
    rm -rf "$OUTPUT_DIR"
    echo "Test output cleaned up."
  else
    echo "Test output preserved: $OUTPUT_DIR"
  fi
else
  echo "✗ Some tests failed. Check logs in: $OUTPUT_DIR"
  echo ""
  echo "Common issues:"
  echo "  - Missing R packages (install dplyr, optparse, ranger)"
  echo "  - Incorrect paths in run_single_replication.R"
  echo "  - optimaltrees not installed"
  echo ""
  exit 1
fi

echo ""
echo "=========================================="
