# Runtime Estimates Based on Profiling

## Bottleneck Identified

**99.45% of time is spent in C++ tree optimization** (`treefarms_fit_with_config_cpp`)

Each DML-ATT estimate requires **15 tree optimizations**:
- 5 cross-validation folds
- 3 models per fold (propensity, outcome|treated, outcome|control)

## Empirical Timing

**Single tree optimization time:**
- n=400: ~0.04 seconds per tree
- n=800: ~0.60 seconds per tree (15x slower)
- n=1600: estimated ~2-3 seconds per tree (50-75x slower)

**Full DML-ATT estimate time:**
- n=400: 15 trees × 0.04s = ~0.6 seconds
- n=800: 15 trees × 0.60s = ~9 seconds
- n=1600: 15 trees × 2.5s = ~37.5 seconds

## Full Simulation Estimates

### Original Plan (24 configs × 500 reps = 12,000 runs)

**By sample size:**
- 8 configs at n=400: 8 × 500 × 0.6s = 2,400s (~40 min)
- 8 configs at n=800: 8 × 500 × 9s = 36,000s (~10 hours)
- 8 configs at n=1600: 8 × 500 × 37.5s = 150,000s (~42 hours)

**Total estimated time: ~52 hours** (2+ days running continuously)

### Problem

This is much slower than initial estimates because:
1. Tree optimization scales poorly with n (exponential in search space)
2. DML requires 15 tree fits per estimate (not just 1)
3. Initial timing tests used single trees, not full DML workflow

## Options

### Option 1: Reduce Sample Sizes (FASTEST)
**Use only n=400**
- 8 configs × 500 reps = 4,000 runs
- Estimated time: 40 minutes
- **Trade-off:** No evidence of scaling behavior

### Option 2: Reduce Replications
**Use 100 reps instead of 500**
- 24 configs × 100 reps = 2,400 runs
- Estimated time: ~10 hours (overnight run)
- n=400: 8 min
- n=800: 2 hours
- n=1600: 8.3 hours
- **Trade-off:** Less precise estimates (still adequate for 100 reps)

### Option 3: Mixed Strategy (RECOMMENDED)
**More reps for n=400, fewer for larger n:**
- n=400: 500 reps (8 configs × 500 = 4,000 runs, ~40 min)
- n=800: 200 reps (8 configs × 200 = 1,600 runs, ~4 hours)
- n=1600: 100 reps (8 configs × 100 = 800 runs, ~8.3 hours)
- **Total: 6,400 runs, ~12.6 hours** (overnight run)
- **Trade-off:** Balanced - high precision at n=400, adequate precision at larger n

### Option 4: Optimize Tree Algorithm
**Investigate why n=800 is 15x slower:**
- Is regularization too weak? (allows deeper trees)
- Can we increase regularization for larger n?
- Are there algorithmic optimizations available?

### Option 5: Drop n=1600, Focus on n=400 and n=800
**Keep two sample sizes:**
- 16 configs (2 DGPs × 2 n × 4 methods) × 500 reps = 8,000 runs
- Estimated time: ~11 hours
- **Trade-off:** Still shows scaling behavior, avoids slowest configs

## Recommendation

**Option 3 (Mixed Strategy)** or **Option 5 (Drop n=1600)**

Both provide:
- Evidence of finite-sample performance
- Evidence of scaling behavior
- Reasonable runtime (overnight run)
- Sufficient replications for statistical power

**Next step:** User decides which trade-off is acceptable.
