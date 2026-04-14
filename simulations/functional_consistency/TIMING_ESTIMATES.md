# Timing Estimates for O2 Deployment

**Updated:** 2026-04-13 (after increasing time limit to 2 hours)

---

## Per Replication Estimates

Each replication involves:
1. Generate DGP (n observations, binary covariates)
2. Run estimation method

### By Sample Size

| n | Standard M-split | Averaged Tree | Pattern Aggregation |
|---|-----------------|---------------|---------------------|
| 200 | 30-45 sec | 45-60 sec | 30-45 sec |
| 400 | 1-1.5 min | 1.5-2 min | 1-1.5 min |
| 800 | 2-3 min | 3-5 min | 2-3 min |
| 1600 | 4-6 min | 6-9 min | 4-6 min |
| 3200 | 6-10 min | 10-15 min | 6-10 min |

**Averaged tree takes longer** because it:
- Runs standard M-split first (to get modal structure)
- Then does M×K additional refits (M=10, K=2-5 → 20-50 refits)
- Predicts on full dataset from each refit

**Pattern aggregation is similar to standard M-split** because it:
- Runs standard M-split
- Just aggregates existing predictions (minimal extra work)

### By K

Higher K means:
- More folds per split
- But each fold trains on more data
- Net effect: K=5 ~30% slower than K=2

---

## Per Job Estimates

Each SLURM job runs **10 replications** for one configuration.

### Fast Jobs (n=200)
- Standard M-split: **5-8 minutes**
- Averaged tree: **8-10 minutes**
- Pattern aggregation: **5-8 minutes**

### Medium Jobs (n=800)
- Standard M-split: **20-30 minutes**
- Averaged tree: **30-50 minutes**
- Pattern aggregation: **20-30 minutes**

### Slow Jobs (n=3200, K=5)
- Standard M-split: **60-100 minutes**
- Averaged tree: **100-150 minutes** (longest!)
- Pattern aggregation: **60-100 minutes**

**Time limit set:** 2 hours per job
- Allows slowest jobs to complete
- Most jobs will finish much faster

---

## Total Runtime Estimate

**135 configurations:**
- 27 configs per n value (3 methods × 3 K × 3 DGP)
- Each config: 50 array tasks

**If all jobs run in parallel:**
- Limited by slowest job: ~2 hours
- But O2 has node limits, so jobs queue

**Realistic estimate with queuing:**
- Assume ~20-30 jobs can run simultaneously
- Wave 1 (fast jobs, n=200): 10-20 minutes
- Wave 2 (medium, n=400-800): 30-60 minutes
- Wave 3 (large, n=1600): 1-2 hours
- Wave 4 (largest, n=3200): 2-3 hours
- **Total wall time: 4-6 hours**

**Conservative estimate:** 8-12 hours (accounting for queue delays)

---

## Partition Choice

**Current:** `short` partition
- Time limit: typically 12 hours
- Our jobs: up to 2 hours each
- ✓ Appropriate

**If jobs timeout:**
- Switch to `medium` or `long` partition
- Edit `#SBATCH --partition=` in `run_simulations.slurm`

---

## Monitoring

Check progress with:
```bash
bash check_progress.sh
```

See running/pending jobs:
```bash
squeue -u $USER -n fc_sim
```

See job time remaining:
```bash
squeue -u $USER -n fc_sim -o "%.18i %.9P %.30j %.8u %.8T %.10M %.10l"
```

---

## If Jobs Timeout

**Symptoms:**
- Jobs show "TIMEOUT" in squeue
- Missing result files for specific configs

**Solution:**
1. Identify which configs timed out (check which result files are missing)
2. Re-run just those configs:
   ```bash
   # Edit launch_all_simulations.sh to only include failed configs
   # Or submit manually:
   sbatch run_simulations.slurm 3200 sparse averaged_tree 5 10
   ```

3. Or increase time limit further:
   ```bash
   # Edit run_simulations.slurm
   #SBATCH --time=4:00:00  # 4 hours
   ```

---

## Expected Completion

**Launch:** Now (when you run deployment commands)
**First results:** ~10-20 minutes (fast jobs complete)
**50% complete:** ~2-3 hours
**All complete:** ~4-8 hours (could be up to 12 with heavy queuing)

**Best practice:** Check progress every 1-2 hours, be patient!

---

## Bottleneck Analysis

**Slowest operations:**
1. **Averaged tree with n=3200, K=5**
   - 50 refits × 10 reps = 500 tree fits per job
   - Each fit: ~2,560 training observations (4/5 × 3200)
   - Total: ~100-150 minutes per job

2. **Tree fitting with large n**
   - Even with regularization, large trees take time
   - Depth computation, split finding, pruning

3. **M×K predictions on full dataset**
   - Averaged tree predicts on all n=3200 observations, M×K=50 times
   - Pattern aggregation: just averages (fast)

**Optimization opportunities** (if jobs are too slow):
- Reduce M from 10 to 5 (cuts averaged tree work in half)
- Use smaller K for large n
- Increase regularization (smaller trees)
- But we want realistic settings, so current design is reasonable

---

## Storage Requirements

**Per replication:** ~10KB (one data frame row saved as RDS)
**Per job:** 10 reps × ~10KB = ~100KB
**Total intermediate:** 135 configs × 50 jobs × 100KB = ~675MB
**Combined results:** ~70MB (single data frame)

✓ Well within O2 limits

---

## Summary

**Most likely scenario:**
- 4-6 hours total wall time
- Jobs complete in waves (fast → slow)
- No timeouts with 2-hour limit
- Check progress periodically

**Worst case:**
- Heavy queuing → 8-12 hours
- Some n=3200 jobs timeout → re-run those
- Still completable within a day

**Best case:**
- Light queue, many parallel nodes
- All complete in 2-3 hours
- Combine results immediately
