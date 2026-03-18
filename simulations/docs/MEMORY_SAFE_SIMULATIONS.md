# Memory-Safe Simulation Guide

**Context:** After system restart from memory exhaustion during β smoothness simulations.

**Date:** 2026-03-11

---

## What Happened

Running large-scale tree-based DML simulations (50-500 reps × 5-fold CV × tree optimization) consumed all available memory (~16GB) and caused system restart.

**Memory hotspots:**
1. **Cross-validation:** 5 folds × model copies
2. **Tree optimization:** Candidate model generation (can hit 10k+ models)
3. **No garbage collection:** R doesn't clean up between reps
4. **Parallel processing:** Multiplies memory per core

---

## Memory-Safe Configuration

### 1. Start with Micro Tests (< 1 minute)

```r
N_REPS <- 3
N_VALUES <- c(400)  # Small n
K_FOLDS <- 2       # Fewer folds
use_rashomon <- FALSE
cv_regularization <- FALSE
```

**Expected memory:** ~500MB per rep × 3 reps = 1.5GB total

### 2. Small Test (2-5 minutes)

```r
N_REPS <- 10
N_VALUES <- c(800)
K_FOLDS <- 5
use_rashomon <- FALSE
cv_regularization <- FALSE
```

**Expected memory:** ~1GB per rep × 10 reps = 10GB total (sequential)

### 3. Medium Test (10-30 minutes)

```r
N_REPS <- 50
N_VALUES <- c(800)
K_FOLDS <- 5
use_rashomon <- FALSE
cv_regularization <- FALSE
```

**Expected memory:** ~1GB per rep × 50 reps = safe if sequential

### 4. Full Study (hours)

```r
N_REPS <- 500
N_VALUES <- c(400, 800, 1600)
K_FOLDS <- 5
```

**Expected memory:** Run in batches (see below)

---

## Memory Management Strategies

### Strategy 1: Explicit Garbage Collection

Add to simulation loop:

```r
for (rep in 1:N_REPS) {
  # ... run simulation ...

  # Force cleanup every 10 reps
  if (rep %% 10 == 0) {
    gc(verbose = FALSE, full = TRUE)
  }
}
```

### Strategy 2: Batch Processing

```r
# Instead of 500 reps at once, run 5 batches of 100
BATCH_SIZE <- 100
N_BATCHES <- 5

for (batch in 1:N_BATCHES) {
  batch_results <- run_batch(
    start_rep = (batch-1)*BATCH_SIZE + 1,
    end_rep = batch*BATCH_SIZE
  )

  # Save immediately
  saveRDS(batch_results,
          sprintf("results/batch_%d.rds", batch))

  # Clear memory before next batch
  rm(batch_results)
  gc(full = TRUE)
}

# Combine later
all_results <- lapply(1:N_BATCHES, function(b) {
  readRDS(sprintf("results/batch_%d.rds", b))
})
results <- do.call(rbind, all_results)
```

### Strategy 3: Monitor Memory During Runs

```r
monitor_memory <- function() {
  mem_info <- gc()
  used_mb <- sum(mem_info[, "used"]) * 0.001  # Convert to MB
  cat(sprintf("  [Memory: %.0f MB]\n", used_mb))

  # Warn if approaching limit
  if (used_mb > 12000) {  # 12GB
    warning("Approaching memory limit - forcing GC")
    gc(full = TRUE)
  }
}

# Add to loop
for (rep in 1:N_REPS) {
  cat(sprintf("Rep %d/%d...", rep, N_REPS))

  # ... run simulation ...

  if (rep %% 5 == 0) monitor_memory()
}
```

### Strategy 4: Disable Parallelization

**Problem:** Parallel processing multiplies memory usage

```r
# DON'T do this for large simulations:
mclapply(1:N_REPS, run_sim, mc.cores = 4)  # 4x memory!

# Instead, run sequentially:
lapply(1:N_REPS, run_sim)
```

### Strategy 5: Reduce Tree Model Complexity

If memory issues persist:

```r
dml_att(
  ...,
  regularization = 0.05,  # Higher = simpler trees
  model_limit = 10000,    # Limit candidate models
  max_depth = 10,         # Limit tree depth
  use_rashomon = FALSE,   # No Rashomon set
  cv_regularization = FALSE  # No CV for reg tuning
)
```

---

## Recommended Workflow

### Phase 1: Verify Setup (5 minutes)

Run `final_test.R` (5 reps):
- Confirms code works
- Estimates memory per rep
- Expected: ~100% convergence, ~95% coverage

### Phase 2: Pilot Study (30 minutes)

Run 50 reps × 1 sample size × 1 regime:
- Tests stability
- Estimates time per rep
- Expected memory: <10GB

### Phase 3: Full Study (batches)

Run in 5 batches of 100 reps each:
- Total: 500 reps × 3 regimes × 3 sample sizes
- Save after each batch
- Clear memory between batches

**Time estimate:** ~12-24 hours total

---

## Current Configuration Status

**Last successful run:** `final_test_20260311_1428.rds` (5 reps)
- Result: Should check if this completed successfully

**Blocked run:** Likely `test_beta_study.R` or full `run_beta_study.R`
- These were configured for 50-500 reps
- No batch processing
- No memory monitoring

---

## Quick Memory Check

```r
# Before starting large simulation:
cat("Available memory:\n")
system("vm_stat | head -10")

# Expected: "Pages free" should be >> 100k pages (~400MB)
# If "Pages free" < 50k, don't start large simulation
```

---

## Recovery After Crash

1. Check what completed:
```bash
ls -lh doubletree/simulations/production/results/*.rds
```

2. Load most recent result:
```r
last_result <- readRDS("results/final_test_20260311_1428.rds")
print(last_result)
```

3. Check if any batch files exist:
```bash
ls -lh doubletree/simulations/production/results/batch_*.rds
```

4. If partial results exist, load and combine
5. Resume from where it stopped (adjust SEED_OFFSET)

---

## Flags to Add to Simulation Scripts

**Essential flags:**
```r
MEMORY_SAFE_MODE <- TRUE

if (MEMORY_SAFE_MODE) {
  # Reduce batch size
  BATCH_SIZE <- 50

  # Force GC every N reps
  GC_INTERVAL <- 10

  # Monitor memory
  MEMORY_WARN_MB <- 12000

  # No parallelization
  USE_PARALLEL <- FALSE
}
```

---

## Next Steps

1. Check system memory: `vm_stat`
2. Verify last results completed: `final_test_20260311_1428.rds`
3. Choose configuration based on phases above
4. Start with Phase 1 (micro test) to ensure system stable
5. Progress to Phase 2 only if Phase 1 completes cleanly
6. Run Phase 3 in batches with explicit memory management

---

**Last updated:** 2026-03-11 after system restart
