# Pre-Flight Checklist for Production Simulations

**Use this before running any simulation with N_REPS > 10**

Run date: ________________

---

## 1. Code Review

- [ ] Script sources `simulation_helpers.R`
- [ ] All package loading wrapped in `suppressMessages({})`
- [ ] Worker function uses `suppress_all()` for ALL fitting code
- [ ] All model calls have `verbose = FALSE`
- [ ] No `sink()` calls to file paths (only `/dev/null` or `NUL`)
- [ ] No explicit log file creation (`*.log`, `*_output.txt`, etc.)
- [ ] Progress uses `progress_msg()` with `every >= 50`
- [ ] Memory monitoring uses `monitor_memory()` if long-running
- [ ] Results saved with `safe_save()` for atomic writes

## 2. Test Run (N_REPS = 3)

```r
# Set to test mode
N_REPS <- 3
N_CORES <- 1  # Single-threaded for testing

# Run script
source("run_beta_study.R")  # Or whichever script
```

- [ ] Completes without errors
- [ ] Results file created and loadable
- [ ] **No .log files created** (`ls -lh *.log`)
- [ ] Console output is minimal (< 20 lines)
- [ ] Results look reasonable

## 3. File System Check

Before starting:
```bash
# Check current directory size
du -sh .

# Check for existing large files
find . -name "*.log" -size +1M
find . -name "*.Rout" -size +1M
find . -name "*output*.txt" -size +1M

# Expected: no files found
```

- [ ] Current directory < 100 MB
- [ ] No large log files exist
- [ ] No temp files from previous runs

## 4. Resource Availability

```r
# In R console:
parallel::detectCores()  # Should be >= 2
system("vm_stat | head -10")  # Check free memory

# Expected: "Pages free" > 100k (~400 MB)
```

- [ ] At least 2 CPU cores available
- [ ] At least 8 GB RAM free
- [ ] Not running other heavy processes

## 5. Configuration Sanity Check

```r
# Verify these are set correctly:
N_REPS      # Should be 500 for production
N_CORES     # Should be 2-4 (NOT > 4 to avoid memory issues)
SEED_OFFSET # Unique per study to avoid overlap
```

- [ ] `N_REPS` is appropriate (500 for production, 50 for pilot)
- [ ] `N_CORES` ≤ 4 (more causes memory issues)
- [ ] `SEED_OFFSET` is unique and documented
- [ ] Output directory path is correct

## 6. Expected Runtime

Estimate:
```r
# Rough formula:
estimated_hours = (N_REPS * n_configs * 3_seconds) / (3600 * N_CORES)

# For run_beta_study.R:
# 500 reps × 27 configs × 3 sec / (3600 × 2 cores) ≈ 11 hours

# For run_primary.R:
# 500 reps × 36 configs × 3 sec / (3600 × 2 cores) ≈ 15 hours
```

- [ ] Estimated runtime is acceptable
- [ ] Can leave running overnight if needed
- [ ] Won't interfere with other scheduled work

## 7. Monitoring Plan

During run:
```bash
# In separate terminal, check every hour:
ls -lh results/*.rds     # Should see file growing
du -sh .                 # Directory size should grow slowly
find . -name "*.log"     # Should return NOTHING

# If any .log files appear: STOP THE SIMULATION
```

- [ ] Know how to check progress (`ls -lh results/`)
- [ ] Know how to check for log bloat (`find . -name "*.log"`)
- [ ] Know how to kill R process if needed (`pkill -9 R`)

## 8. Backup Plan

- [ ] Results directory is in git (if appropriate) or backed up
- [ ] Can resume from partial results if interrupted
- [ ] Batch scripts available if memory becomes an issue

## 9. Post-Run Verification

After completion:
```r
# In R:
results <- readRDS("results/primary_2026-03-12/simulation_results.rds")
nrow(results)  # Should equal N_REPS × n_configs
mean(results$converged)  # Should be > 0.95

# Check for log bloat:
check_large_files(".", min_mb = 10)
```

- [ ] Results file is complete
- [ ] Convergence rate > 95%
- [ ] No large unexpected files created
- [ ] Summary statistics look reasonable

---

## Emergency: Log Bloat Detected

**If you find .log files > 100 MB during a run:**

1. **Stop immediately**: `pkill R` (from terminal)
2. **Check damage**: `find . -name "*.log" -ls`
3. **Remove logs**: `rm *.log` (after confirming they're not needed)
4. **Review script**: Check for missing `suppress_all()` calls
5. **Fix and test**: Run N_REPS=3 test first
6. **Resume**: Restart from last saved batch if using batch processing

---

## Sign-Off

Before starting production run with N_REPS ≥ 100:

**Reviewed by:** ________________

**Date:** ________________

**All checks passed:** [ ] YES

**Notes:**

---

**Remember:** It's always faster to spend 5 minutes checking than to spend 10 hours re-running after a crash.
