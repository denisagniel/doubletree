# Production Simulations Launched

**Date:** 2026-03-13
**Configuration:** All 3 DGPs (comprehensive)
**Status:** Running in background

---

## Configuration

**Simulation Grid:**
- DGPs: 3 (binary, continuous, moderate)
- Methods: 4 (tree, rashomon, forest, linear)
- Sample sizes: 3 (n = 400, 800, 1600)
- Replications: 500 per config
- **Total: 18,000 simulations**

**Parallelization:**
- Cores: 4 (auto-detected)
- Background: nohup with log file

**Output:**
- Directory: `results/primary_2026-03-13/`
- Results: `simulation_results.rds`
- Summary: `summary_stats.csv`

---

## Estimated Runtime

**With hash optimization fix:**
- DGP1 (binary): ~0.12 sec/rep
- DGP2 (continuous): ~0.03 sec/rep
- DGP3 (moderate): ~0.04 sec/rep

**Conservative estimate:** 4-8 hours with 4 cores

**Previous estimate (60h) was based on pre-fix measurements**

---

## Monitoring Commands

### Check Progress
```bash
# View live output (last 50 lines)
tail -50 production_run_*.log

# Count completed simulations
grep -c "Simulations completed" production_run_*.log

# Check memory usage
ps aux | grep Rscript

# Monitor results directory
ls -lh results/primary_2026-03-13/
```

### Check Status Script
```bash
# Use built-in monitoring
./check_status.sh
```

### If Stuck/Hanging
```bash
# Check if process is running
ps aux | grep run_primary.R

# Kill if needed (use PID from above)
kill <PID>
```

---

## Verification After Completion

1. **Check log for errors:**
   ```bash
   grep -i error production_run_*.log
   grep -i warning production_run_*.log
   ```

2. **Verify output files:**
   ```bash
   ls -lh results/primary_2026-03-13/
   ```

3. **Check completion:**
   ```r
   results <- readRDS("results/primary_2026-03-13/simulation_results.rds")
   nrow(results)  # Should be 18,000
   ```

---

## What Was Fixed

**Bug:** doubletree hung at n=800 with Rashomon due to JSON parsing bottleneck in structure intersection

**Fix:**
- Hash-based comparison (digest::xxhash64) instead of JSON string matching
- 100-1000x speedup for large Rashomon sets
- Verified working: all DGPs complete at n=800 in <0.2 seconds

**Commits:**
- optimaltrees (216cc03): Hash-based intersection
- doubletree (564c31f): Pass discretization parameters

---

## Expected Completion

**Started:** 2026-03-13 evening
**Expected:** 4-8 hours (overnight run)
**Check:** Tomorrow morning

---

## Next Steps After Completion

1. Verify results (see above)
2. Run stress simulations (`run_stress.R`)
3. Generate manuscript outputs (`analyze_manuscript.R`)
4. Create Tables 1-2, Figures 1-2
5. Three-way fidelity check (paper ↔ code ↔ results)
