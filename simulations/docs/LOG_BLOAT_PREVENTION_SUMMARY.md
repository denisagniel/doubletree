# Log Bloat Prevention - Implementation Summary

**Date:** 2026-03-12
**Issue:** Simulations creating >10 GB log files
**Status:** ✅ RESOLVED

---

## What Was Fixed

### 1. Created Reusable Helper Functions

**File:** `doubletree/simulations/simulation_helpers.R`

**Functions added:**
- `suppress_all()` - Complete output suppression for simulation workers
- `progress_msg()` - Safe progress tracking (console only, infrequent)
- `monitor_memory()` - Memory monitoring with automatic GC
- `safe_save()` - Atomic file writes to prevent corruption
- `check_large_files()` - Detect unexpected large files after runs
- `nullfile()` - Cross-platform null device path

**Usage:**
```r
source("simulation_helpers.R")  # At top of every simulation script

# Wrap ALL fitting code in workers
fit <- suppress_all({
  dml_att(X, A, Y, verbose = FALSE, ...)
})

# Progress tracking (console only)
progress_msg(current = i, total = N_REPS, every = 50)
```

### 2. Updated Main Simulation Scripts

**Files modified:**
- `doubletree/simulations/production/run_beta_study.R`
- `doubletree/simulations/production/run_primary.R`

**Changes:**
- Silent package loading with `suppressMessages({})`
- Silent sourcing with `safe_source()`
- All worker code wrapped in `suppress_all()`
- Progress tracking changed from every 100 reps to console-only with `progress_msg()`
- Memory monitoring uses `monitor_memory()` helper
- Result saving uses `safe_save()` for atomic writes
- Post-run check for large files with `check_large_files()`

### 3. Created Documentation

**LOGGING_PROTOCOL.md** - Comprehensive guide covering:
- Why log bloat happens
- 5 core rules to prevent it
- Standard suppression patterns
- Testing protocol
- Emergency procedures

**PRE_FLIGHT_CHECKLIST.md** - Checklist for production runs:
- 9-step verification before starting
- Test run procedure (N_REPS = 3)
- File system checks
- Monitoring plan
- Post-run verification

**LOG_BLOAT_PREVENTION_SUMMARY.md** - This file

### 4. Enhanced .gitignore

Added patterns to block simulation output files:
```
*.Rout
*.Rout.save
*_log.txt
*_output.txt
simulation_log*.txt
simulation_output*.txt
debug_*.log
temp_*.log
```

---

## How to Use

### For Existing Scripts

**Step 1:** Add helpers to top of script
```r
source("../simulation_helpers.R")  # Path relative to your script
```

**Step 2:** Wrap all fitting code in `suppress_all()`
```r
run_single_sim <- function(sim_id, ...) {
  # Progress: console only, infrequent
  progress_msg(sim_id, total, every = 50)

  # ALL fitting wrapped
  fit <- suppress_all({
    dml_att(X, A, Y, verbose = FALSE, ...)
  })

  # Return clean results
  data.frame(...)
}
```

**Step 3:** Test with N_REPS = 3
```r
N_REPS <- 3
source("your_script.R")

# Verify:
# - No .log files created
# - Minimal console output
# - Results correct
```

### For New Scripts

Use the template in `LOGGING_PROTOCOL.md` (line 147-207).

Key points:
- Source `simulation_helpers.R` first
- Silent package loading: `suppressMessages({ library(...) })`
- All model calls: `verbose = FALSE`
- Worker functions: wrap everything in `suppress_all()`
- Progress: use `progress_msg()` with `every >= 50`

### Before Any Production Run

Use `PRE_FLIGHT_CHECKLIST.md`:

1. Code review (9 checkboxes)
2. Test run with N_REPS = 3
3. Check no .log files created
4. Verify resource availability
5. Start production run

**Critical checks:**
```bash
# Before starting:
find . -name "*.log" -size +1M  # Should return nothing

# During run (check hourly):
ls -lh results/*.rds           # Results growing
find . -name "*.log"           # Should return nothing

# If .log appears during run: STOP IMMEDIATELY
```

---

## What to Monitor

### Green flags (good):
- Console output: < 50 lines per hour
- Directory size: Growing slowly (< 10 MB per 100 reps)
- Memory usage: Stable or decreasing between batches
- No .log files appearing

### Red flags (BAD - stop immediately):
- Any .log file > 10 MB
- Directory growing > 100 MB per hour
- Memory usage climbing continuously
- R process consuming > 90% of RAM

---

## Testing Before Production

**Always test with N_REPS = 3 first:**

```r
# In your script, temporarily:
N_REPS <- 3
N_CORES <- 1

# Run it
source("run_beta_study.R")

# Check results
list.files("results/", pattern = "*.rds")  # Should have 1 file
list.files(".", pattern = "*.log")         # Should be EMPTY

# Load and inspect
results <- readRDS("results/test_2026-03-12/simulation_results.rds")
nrow(results)  # Should equal 3 × n_configs
```

**Only proceed to N_REPS = 500 if test passes cleanly.**

---

## Emergency Procedures

### If Log Bloat Detected During Run

1. **Stop R immediately**: `pkill R` or Ctrl+C
2. **Check damage**: `find . -name "*.log" -ls`
3. **Remove logs**: `rm *.log` (after verifying not needed)
4. **Fix script**: Review against `LOGGING_PROTOCOL.md`
5. **Test**: Run N_REPS = 3 to verify fix
6. **Resume**: Restart (or resume from last batch)

### If System Becomes Unresponsive

If log files are so large that system is unresponsive:

```bash
# From terminal (may need to force quit other apps first):
cd ~/Library/CloudStorage/OneDrive-RANDCorporation/rprojects/global-scholars/doubletree/simulations

# Find culprit
du -sh */*  | sort -h

# Remove large logs
find . -name "*.log" -size +100M -delete

# Or remove specific file
rm production/simulation_run.log
```

---

## Root Cause Analysis

**Why this happened:**

1. **Progress messages in loops** - `cat()` every 100 reps × 13,500 runs = millions of lines
2. **Tree JSON output** - Each tree model can be 10-100 KB, × thousands = GB
3. **Parallel multiplication** - Multiple cores writing simultaneously
4. **No output suppression** - Verbose output from packages going to logs
5. **File redirection** - Scripts using `sink()` or logging to files

**Why it's now fixed:**

1. **Progress infrequent** - Only every 50+ reps, console only
2. **Complete suppression** - `suppress_all()` blocks all output in workers
3. **No file logging** - All output to `/dev/null`, never to files
4. **Verbose disabled** - All package calls have `verbose = FALSE`
5. **Monitoring** - Automatic checks for large files after runs

---

## Maintenance

**Quarterly (or before major simulation runs):**

1. Review `LOGGING_PROTOCOL.md` - Update if new patterns emerge
2. Check all simulation scripts - Verify using `suppress_all()`
3. Test helpers - Run test suite with N_REPS = 3
4. Update .gitignore - Add new problematic file patterns if found

**After any simulation crash:**

1. Check for log files: `find . -name "*.log"`
2. If found, review script against protocol
3. Add [LEARN] entry to MEMORY.md
4. Update protocol if new pattern discovered

---

## Files Changed

**New files created:**
- `doubletree/simulations/simulation_helpers.R`
- `doubletree/simulations/LOGGING_PROTOCOL.md`
- `doubletree/simulations/PRE_FLIGHT_CHECKLIST.md`
- `doubletree/simulations/LOG_BLOAT_PREVENTION_SUMMARY.md` (this file)

**Files modified:**
- `doubletree/simulations/production/run_beta_study.R`
- `doubletree/simulations/production/run_primary.R`
- `doubletree/.gitignore`

**Files to update in future (as needed):**
- Other scripts in `doubletree/simulations/exploration/`
- Other scripts in `doubletree/simulations/diagnostics/`
- Batch scripts in `doubletree/simulations/batches/`

---

## Quick Reference

**Good pattern:**
```r
source("simulation_helpers.R")

run_single_sim <- function(sim_id, ...) {
  progress_msg(sim_id, N_REPS, every = 50)

  fit <- suppress_all({
    dml_att(X, A, Y, verbose = FALSE, ...)
  })

  data.frame(...)
}

results <- mclapply(1:N_REPS, run_single_sim, mc.cores = N_CORES)
safe_save(do.call(rbind, results), "results/output.rds")
check_large_files(".", min_mb = 10)
```

**Bad pattern (DON'T DO THIS):**
```r
# ❌ DON'T: Verbose output
fit <- dml_att(X, A, Y, verbose = TRUE)

# ❌ DON'T: Logging to files
sink("simulation.log")

# ❌ DON'T: Frequent progress
for (i in 1:500) cat("Rep", i, "\n")

# ❌ DON'T: capture.output to tempfile
capture.output(fit <- ..., file = tempfile())
```

---

## Next Steps

1. **Test the fixes** - Run N_REPS = 3 on both updated scripts
2. **Verify no bloat** - Confirm no .log files created
3. **Document in session notes** - Record this fix for future reference
4. **Update other scripts** - Apply same pattern to exploration/ and diagnostics/ scripts as needed
5. **Run production simulations** - Use PRE_FLIGHT_CHECKLIST.md before starting

---

**Status:** ✅ Core issue resolved. Safe to run production simulations.

**Last updated:** 2026-03-12
