# Batch Simulation System - README

## Overview

Memory-safe batch processing system for Rashomon-DML simulations. Splits 6,400 runs into 16 manageable batches.

## Configuration

- **Total configs:** 64 (4 DGPs × 4 sample sizes × 4 tolerances)
- **Replications per config:** 100
- **Total runs:** 6,400
- **Batch organization:** 16 batches of 4 configs each

## Batch Structure

Each batch:
- **1 DGP** (dgp1, dgp2, dgp3, or dgp4)
- **1 sample size** (200, 400, 800, or 1600)
- **4 tolerances** (0.01, 0.05, 0.1, 0.2)
- **Result:** 400 simulation runs per batch

## Memory & Performance

Per batch (estimated):
- **Runs:** 400
- **Time:** 5-10 minutes
- **Peak memory:** 3-5 GB
- **Output:** 4 result files (~2-8 MB each)

Full run (all 16 batches):
- **Total time:** ~1.5-2.5 hours
- **Peak memory:** 3-5 GB (never accumulates)
- **Total output:** ~250 MB

## Usage

### Option 1: Run All Batches Sequentially (Recommended)

```bash
cd doubletree
Rscript simulations/run_all_batches.R
```

This runs all 16 batches one after another. Each batch runs in a fresh R process, preventing memory accumulation.

### Option 2: Run Individual Batches

```bash
# Run specific batch
Rscript simulations/batches/batch_01_dgp1_n200.R

# Run first 4 batches (all DGPs at n=200)
for i in 01 02 03 04; do
  Rscript simulations/batches/batch_${i}_*.R
done
```

Useful if:
- You want to test with a small subset first
- A batch failed and needs rerunning
- You want to run batches in parallel (on multi-core system)

### Option 3: Test Run (Single Batch)

```bash
# Test with batch 1 only (fastest: n=200)
Rscript simulations/batches/batch_01_dgp1_n200.R

# Should take ~5 minutes, use <3GB memory
```

## After Batches Complete

### 1. Combine Results

```bash
Rscript simulations/combine_batch_results.R
```

This merges all 64 config results into:
- `simulation_summary.rds` (R data frame)
- `simulation_summary.csv` (for viewing)

### 2. Generate Figures and Tables

```bash
Rscript simulations/analyze_results.R
```

Creates:
- `figures/figure1_intersection_existence.pdf`
- `figures/figure2_bias_tradeoff.pdf`
- `figures/figure3a_coverage.pdf`
- `figures/figure3b_ci_width.pdf`
- `figures/figure4_rashomon_overhead.pdf`
- `figures/table1_summary.csv`

## Batch Files Location

- **Batch scripts:** `simulations/batches/batch_*.R`
- **Results:** `simulations/results_extended/result_*.rds`
- **Combined:** `simulations/results_extended/simulation_summary.csv`
- **Figures:** `simulations/figures/*.pdf`

## Troubleshooting

### Batch Failed

Check which configs are missing:

```bash
ls simulations/results_extended/result_*.rds | wc -l
# Should be 64 if all batches completed
```

Rerun specific failed batch:

```bash
Rscript simulations/batches/batch_XX_dgpY_nZZZ.R
```

### Memory Issues

If a batch uses >10GB memory:
1. Check for other R processes: `ps aux | grep "[Rr]"`
2. Kill old R processes if needed
3. Reduce replications in batch script (change `n_reps <- 100` to `50`)

### Time Estimates Wrong

Larger sample sizes (n=1600) take longer than small ones (n=200). Adjust expectations:
- Batches 1-4 (n=200): ~5 min each
- Batches 5-8 (n=400): ~6 min each
- Batches 9-12 (n=800): ~8 min each
- Batches 13-16 (n=1600): ~10 min each

## Progress Monitoring

Watch progress in real-time:

```bash
# Count completed configs
ls simulations/results_extended/result_*.rds 2>/dev/null | wc -l

# Monitor memory during batch run
watch -n 5 'ps aux | grep "[Rr]script" | grep -v grep'
```

## Quick Test Before Full Run

Test infrastructure with batch 1:

```bash
Rscript simulations/batches/batch_01_dgp1_n200.R
```

If successful (4 result files created in ~5 min with <3GB memory), proceed with full run.

## Parallel Execution (Advanced)

If you have a multi-core system and want faster completion, run 2-4 batches in parallel:

```bash
# Terminal 1
Rscript simulations/batches/batch_01_dgp1_n200.R &

# Terminal 2
Rscript simulations/batches/batch_02_dgp2_n200.R &

# Terminal 3
Rscript simulations/batches/batch_03_dgp3_n200.R &

# Wait for all
wait
```

Monitor total memory usage to ensure you don't exceed system limits.
