#!/bin/bash

# Quick test: Just run forest method (fastest, no optimaltrees S7 issues)

set -e

echo "Quick Test: Forest method only (fastest)"
echo ""

OUTPUT_DIR="/tmp/o2_quick_test_$(date +%s)"
mkdir -p "$OUTPUT_DIR"

echo "Testing: DGP=dgp1, N=400, METHOD=forest"

if Rscript run_single_replication.R \
    --dgp "dgp1" \
    --sample-size "400" \
    --method "forest" \
    --replication 1 \
    --output-dir "$OUTPUT_DIR" \
    --tau 0.10 \
    --k-folds 5 \
    --seed-offset 10000 \
    --worker-limit 1; then

  echo ""
  echo "✓ Test PASSED"
  echo ""
  echo "Output file created:"
  ls -lh "$OUTPUT_DIR"/*.rds
  echo ""
  echo "Ready to deploy to O2!"

  rm -rf "$OUTPUT_DIR"
else
  echo ""
  echo "✗ Test FAILED"
  exit 1
fi
