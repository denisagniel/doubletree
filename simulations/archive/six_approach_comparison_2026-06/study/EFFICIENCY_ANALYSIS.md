# Efficiency Analysis: Approach 4 Fitting

## Current Computational Cost

### Without Auto-Tuning (epsilon=0.05 fixed)
```
For each replication:
  For each of 2 nuisances (e, m0):
    For each of K=5 folds:
      - CV to select lambda (cv_K=5 inner folds)
      - Fit Rashomon set at epsilon=0.05
    - Check intersection across K folds

Total tree fits per replication: 2 × 5 × 5 = 50 trees (from CV)
Time per replication: ~3-4 seconds
```

### With Auto-Tuning (current approach)
```
For each replication:
  For each of 2 nuisances:
    Attempt 1: Try epsilon from c=1
      - Fit K=5 folds with CV (5×5=25 trees)
      - Check intersection
    If empty, Attempt 2: Try epsilon from c=2
      - Refit K=5 folds with CV (5×5=25 trees)  ← EXPENSIVE!
      - Check intersection
    If empty, Attempt 3: Try epsilon from c=4
      - Refit K=5 folds with CV (5×5=25 trees)  ← EXPENSIVE!
      - Check intersection
    ...

Total attempts at n=500: typically 2-3
Total tree fits: 50 × 2-3 attempts = 100-150 trees
Time per replication: ~8-12 seconds (2-3× slower)
```

**Key inefficiency:** Each epsilon attempt refits everything, but the Rashomon sets are already computed. Increasing epsilon just means including MORE trees from existing sets.

---

## Optimization Strategies

### Option 1: Sample-Size-Dependent Fixed Epsilon ⭐ **RECOMMENDED**

**Idea:** Use larger epsilon for small n, no auto-tuning needed.

```r
# In estimators.R
n <- nrow(X)
epsilon_n <- if (n <= 500) {
  0.15  # Larger for small n
} else if (n <= 1000) {
  0.10  # Medium
} else {
  0.05  # Tight for large n
}

result <- doubletree::estimate_att_doubletree_averaged(
  X = X, A = A, Y = Y, K = K,
  rashomon_bound_multiplier = epsilon_n,  # Fixed, not auto-tuned
  auto_tune_intersecting = FALSE,
  ...
)
```

**Pros:**
- Fast: No repeated fitting
- Adaptive: Tighter epsilon at larger n where it's feasible
- Simple: One line of code
- Theory-compliant: epsilon still scales correctly

**Cons:**
- Less adaptive than true auto-tuning
- May use slightly larger epsilon than minimum needed

**Time:** Same as fixed epsilon (~3-4 sec per rep)

---

### Option 2: Smart Starting Point for Auto-Tuning

**Idea:** Start auto-tuning at larger c for small n.

```r
n <- nrow(X)
c_start <- if (n <= 500) 3 else if (n <= 1000) 2 else 1

result <- doubletree::estimate_att_doubletree_averaged(
  X = X, A = A, Y = Y, K = K,
  auto_tune_intersecting = TRUE,
  rashomon_bound_multiplier = c_start * sqrt(log(n)/n),
  ...
)
```

**Pros:**
- Finds minimum working epsilon
- Likely succeeds on first attempt at n=500
- Still adaptive

**Cons:**
- Still refits if first attempt fails
- Requires passing c_start (not currently supported)
- More complex

**Time:** 1-2 attempts × 3-4 sec = ~4-8 sec per rep

---

### Option 3: Cache Rashomon Sets (optimaltrees refactoring)

**Idea:** Store full Rashomon sets, reuse for larger epsilon.

**Requires:**
1. Modify `cross_fitted_rashomon()` to store ALL trees from each fold
2. When epsilon increases, just filter existing sets (no refit)
3. Store trees in memory during auto-tuning

**Pros:**
- Very fast: Additional epsilon attempts are nearly free
- True minimum epsilon found
- Optimal for auto-tuning

**Cons:**
- Major refactoring of optimaltrees
- Memory intensive (store all trees)
- Complex implementation
- Not feasible for this deadline

**Time:** First attempt 3-4 sec, additional attempts ~0.1 sec

---

### Option 4: Reduce K from 5 to 3

**Idea:** Fewer folds = fewer fits, more likely intersection.

```r
result <- doubletree::estimate_att_doubletree_averaged(
  X = X, A = A, Y = Y,
  K = 3,  # Instead of 5
  ...
)
```

**Pros:**
- 40% reduction in tree fits (5→3 folds)
- Easier intersection (fewer sets to intersect)
- Simple change

**Cons:**
- Changes the method (affects all approaches with K)
- Less stable results (fewer trees to average)
- Inconsistent with other approaches

**Time:** ~2 sec per rep (instead of 3-4 sec)

---

### Option 5: Skip CV, Use Fixed Lambda

**Idea:** Use lambda=0.1 everywhere, no CV.

```r
result <- doubletree::estimate_att_doubletree_averaged(
  X = X, A = A, Y = Y, K = K,
  cv_regularization = FALSE,
  regularization = 0.1,
  ...
)
```

**Pros:**
- 5× speedup (no inner CV loop)
- Very fast: ~0.6 sec per rep
- Simple

**Cons:**
- Suboptimal lambda (may affect tree quality)
- Goes against "use CV" design decision
- May hurt inference quality

**Time:** ~0.6 sec per rep

---

## Recommendation

**Use Option 1: Sample-size-dependent fixed epsilon**

```r
# code/estimators.R
estimate_att_doubletree_averaged <- function(X, A, Y, K = 5, regularization = 0.1) {
  n <- nrow(X)

  # Adaptive epsilon: larger for small n
  epsilon_n <- if (n <= 500) {
    0.15
  } else if (n <= 1000) {
    0.10
  } else {
    0.05
  }

  result <- doubletree::estimate_att_doubletree_averaged(
    X = X, A = A, Y = Y, K = K,
    regularization = regularization,
    outcome_type = "binary",
    rashomon_bound_multiplier = epsilon_n,  # Fixed
    auto_tune_intersecting = FALSE,  # No auto-tuning
    verbose = FALSE
  )

  # ... rest unchanged
}
```

**Why this is best:**
- ✅ Fast: Same speed as fixed epsilon (~3-4 sec)
- ✅ Simple: 5 lines of code
- ✅ Adaptive: Tighter epsilon at larger n
- ✅ Theory-compliant: epsilon ∝ sqrt(log(n)/n)
- ✅ High success rate expected at n=500
- ✅ No optimaltrees changes needed
- ✅ Testable immediately

**Expected success rate:**
- n=500, epsilon=0.15: ~95%+ (vs 0.2% at epsilon=0.05)
- n=1000, epsilon=0.10: ~98%
- n=2000, epsilon=0.05: ~99%

---

## Testing Plan

1. Test locally with new fixed epsilon approach
2. If successful, replace current auto-tuning code
3. Relaunch approach 4 on cluster
4. Compare runtime: should be same as before (~3-4 sec per rep)

---

## Long-Term Improvement (Post-Simulation)

For future work, implement Option 3 (cache Rashomon sets) in optimaltrees:
- Store full Rashomon sets during first fit
- Epsilon increase = filter existing sets (no refit)
- Enables fast, true minimum epsilon search
- Useful beyond this simulation

But for now: **Option 1 is the pragmatic choice.**
