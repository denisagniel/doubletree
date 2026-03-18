# Batch Simulation System

**Created:** 2026-03-13
**Reason:** Full `run_primary.R` hangs in background mode due to mclapply issues on macOS

**Solution:** Run simulations in 3 smaller batches (one per DGP)

---

## Overview

**Total simulations:** 18,000 (3 DGPs × 4 methods × 3 n × 500 reps)

**Batch structure:**
- **Batch 1 (DGP1):** 6,000 simulations (Binary outcome)
- **Batch 2 (DGP2):** 6,000 simulations (Continuous outcome)
- **Batch 3 (DGP3):** 6,000 simulations (Moderate complexity)

**Runtime:** ~4-5 hours total (sequential batches)

---

## Quick Start

### Option A: Run All Batches Automatically
```bash
./run_all_batches.sh
```

This runs all 3 batches sequentially. Logs saved to `logs/batch_run_TIMESTAMP/`

### Option B: Run Batches Individually
```bash
# Run one at a time
Rscript run_dgp1_batch.R > logs/dgp1.log 2>&1
Rscript run_dgp2_batch.R > logs/dgp2.log 2>&1
Rscript run_dgp3_batch.R > logs/dgp3.log 2>&1

# After all complete, combine results
Rscript combine_batch_results.R
```

### Option C: Run in Background
```bash
# Run all batches in background
nohup ./run_all_batches.sh > batch_run.log 2>&1 &

# Monitor progress
./check_batch_progress.sh

# Watch log
tail -f logs/batch_run_*/dgp1.log
```

---

## Monitoring

### Check Progress
```bash
./check_batch_progress.sh
```

Shows:
- Active processes
- Latest logs
- Completion status
- Results files

### Manual Monitoring
```bash
# Check running processes
ps aux | grep "run_dgp.*batch.R"

# View latest log
tail -50 logs/batch_run_*/dgp1.log

# Check results
ls -lh results/dgp*_batch_*/
```

---

## Output Files

### Individual Batch Results
```
results/dgp1_batch_YYYY-MM-DD/
├── dgp1_results.rds      # Full replication data (6,000 rows)
└── dgp1_summary.csv      # Summary statistics

results/dgp2_batch_YYYY-MM-DD/
├── dgp2_results.rds
└── dgp2_summary.csv

results/dgp3_batch_YYYY-MM-DD/
├── dgp3_results.rds
└── dgp3_summary.csv
```

### Combined Results (after running combine_batch_results.R)
```
results/primary_YYYY-MM-DD/
├── simulation_results.rds   # All 18,000 simulations
└── summary_stats.csv        # Aggregated by DGP × method × n
```

---

## Workflow

1. **Run batches:**
   ```bash
   ./run_all_batches.sh
   ```

2. **Monitor progress:**
   ```bash
   ./check_batch_progress.sh
   ```

3. **Combine results:**
   ```bash
   Rscript combine_batch_results.R
   ```

4. **Generate tables/figures:**
   ```bash
   Rscript analyze_manuscript.R
   ```

---

## Troubleshooting

### Batch hangs or dies
```bash
# Check process
ps aux | grep run_dgp

# Check memory
top -pid <PID>

# Kill if stuck
kill <PID>

# Restart that specific batch
Rscript run_dgpX_batch.R > logs/dgpX_restart.log 2>&1 &
```

### Verify individual batch works
```bash
# Test with mini script
Rscript test_dgp1_mini.R
```

### Check for errors
```bash
# Look for errors in logs
grep -i "error\|Error" logs/batch_run_*/dgp*.log
```

---

## Advantages Over Full Script

✅ **Smaller chunks:** Easier to debug and recover from failures
✅ **Better monitoring:** Progress visible per batch
✅ **Fault tolerance:** One batch failing doesn't lose all work
✅ **Flexible:** Can run batches in parallel on different machines
✅ **No mclapply:** Avoids macOS forking issues

---

## Estimated Runtimes

Based on hash-optimized code:

| Batch | Simulations | Runtime |
|-------|-------------|---------|
| DGP1 (Binary) | 6,000 | ~1.3 hours |
| DGP2 (Continuous) | 6,000 | ~0.5 hours |
| DGP3 (Moderate) | 6,000 | ~0.7 hours |
| **Total** | **18,000** | **~2.5 hours** |

Conservative estimate: 4-5 hours (with buffer)

---

## Next Steps After Completion

1. ✓ Combine batch results
2. ✓ Verify 18,000 rows total
3. ✓ Check convergence rates (should be >95%)
4. Generate manuscript outputs:
   - Table 1 (primary results)
   - Table 2 (stress tests - separate script)
   - Figures 1-2
5. Three-way fidelity check (paper ↔ code ↔ results)

---

## Bug Fix Applied

**Problem:** doubletree hung at n=800 with Rashomon
**Solution:** Hash-based tree comparison (commit 216cc03)
**Verification:** All DGPs complete at n=800 in <0.2 seconds ✓

---

## Files

- `run_dgp1_batch.R` - DGP1 simulations
- `run_dgp2_batch.R` - DGP2 simulations
- `run_dgp3_batch.R` - DGP3 simulations
- `run_all_batches.sh` - Master script
- `check_batch_progress.sh` - Progress monitor
- `combine_batch_results.R` - Merge batch outputs
- `test_dgp1_mini.R` - Quick validation test
