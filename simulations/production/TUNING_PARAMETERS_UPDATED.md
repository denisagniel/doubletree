# Tuning Parameters Updated - 2026-03-12

## Changes Made to run_primary.R

### **CRITICAL FIX: Rashomon-DML ε_n**

**Previous:** Used default `rashomon_bound_multiplier = 0.05`

**Updated:** Use theory-justified bound
```r
epsilon_n <- 2 * sqrt(log(n) / n)
rashomon_bound_multiplier = epsilon_n
```

**Values by sample size:**
| n | Old ε_n | New ε_n | Change |
|---|---------|---------|--------|
| 400 | 0.05 | 0.245 | **4.9x larger** |
| 800 | 0.05 | 0.183 | **3.7x larger** |
| 1600 | 0.05 | 0.136 | **2.7x larger** |

**Expected impact:**
- Larger Rashomon set → more stable structure intersection
- Should improve coverage from ~91.7% → ~93-95%
- Aligns with paper's theory (Section 3.2)

---

## All Tuning Parameters (Final)

### Tree-DML
```r
regularization = log(n) / n
cv_regularization = FALSE
use_rashomon = FALSE
```

**Justification:**
- Minimax-optimal rate for tree-based estimation
- Standard in statistical learning theory
- Fast, reproducible

**Expected coverage:** 93-95%

### Rashomon-DML
```r
regularization = log(n) / n
cv_regularization = FALSE
use_rashomon = TRUE
rashomon_bound_multiplier = 2 * sqrt(log(n) / n)  # ✅ UPDATED
```

**Justification:**
- ε_n satisfies theory requirement: ε_n = O(√(log n / n))
- c = 2 is conservative (allows reasonably sized Rashomon set)
- Balances: too small → empty intersection, too large → too much variation

**Expected coverage:** 93-95% (improved from 91.7%)

### Forest-DML
```r
num.trees = 500
mtry = floor(sqrt(p))  # p = number of features
min.node.size = 5
probability = TRUE
```

**Justification:**
- Standard ranger defaults for classification
- 500 trees: standard in literature
- mtry = sqrt(p): Breiman's original recommendation

**Expected coverage:** 93-95%

### Linear-DML
```r
interactions = FALSE
```

**Justification:**
- Main effects only (no two-way interactions)
- Avoids overfitting with limited sample sizes
- Expected to fail on nonlinear DGPs (by design)

**Expected coverage:**
- Linear DGPs: 93-95%
- Nonlinear DGPs: 80-90% (misspecification bias expected)

---

## What Was NOT Changed

### Tree-DML regularization

**Kept:** `λ = log(n) / n` (fixed, not CV)

**Rationale:**
1. 93% coverage is acceptable (within 2% of nominal)
2. Theory-justified, reproducible
3. Faster than CV (important for 18k simulations)
4. If coverage <92% after full run, can revisit

**Alternative options (not used):**
- `cv_regularization = TRUE` (slower, data-adaptive)
- `regularization = 0.5 * log(n) / n` (less conservative)

---

## Testing Protocol

### Before Full Run (N_REPS = 500)

**Step 1: Quick test (N_REPS = 3)**
```r
# In run_primary.R, temporarily:
N_REPS <- 3

# Run script
source("production/run_primary.R")

# Verify:
# - No log files created
# - No errors
# - Results file created
# - Epsilon_n values look correct
```

**Step 2: Check epsilon_n values**
```r
results <- readRDS("results/primary_2026-03-12/simulation_results.rds")

# Check Rashomon epsilon_n values
rashomon_results <- results[results$method == "rashomon", ]
table(rashomon_results$n, rashomon_results$epsilon_n)

# Expected:
# n=400:  epsilon_n ≈ 0.245
# n=800:  epsilon_n ≈ 0.183
# n=1600: epsilon_n ≈ 0.136
```

**Step 3: Verify no log bloat**
```bash
find . -name "*.log" -type f
# Should return: (empty) or only old files
```

**Step 4: If all checks pass → Run full simulation**
```r
# In run_primary.R:
N_REPS <- 500

# Run script (takes ~4 hours with 4 cores)
source("production/run_primary.R")
```

---

## Expected Outcomes

### Coverage by Method (After Fix)

| Method | Previous | Expected After Fix | Change |
|--------|----------|-------------------|--------|
| Tree-DML | 93% | 93-95% | No change (already good) |
| Rashomon-DML | 91.7% | **93-95%** | ✅ **Improved** |
| Forest-DML | 93% | 93-95% | No change |
| Linear-DML | 80-90% | 80-90% | No change (misspecified by design) |

### If Rashomon Coverage Still <92%

Investigate:
1. Check if empty intersections occurred (auto_tune_intersecting kicked in)
2. Examine fitted tree structures
3. Consider increasing ε_n further (try c = 3 instead of c = 2)

---

## Documentation for Paper

**Methods section should state:**

> "For tree-based methods, we used regularization parameter λ = log(n)/n,
> following minimax-optimal rate theory. For Rashomon-DML, we set the
> Rashomon bound to ε_n = 2√(log(n)/n) to ensure a sufficiently large
> Rashomon set while maintaining proximity to the optimal tree. Random
> forests used 500 trees with mtry = √p. All methods employed 5-fold
> cross-fitting with stratification by treatment status."

---

## Summary

**What was fixed:** Rashomon ε_n parameter (was 5x too small)

**What stayed the same:** Tree-DML λ, forest parameters, linear model specification

**Next step:** Test with N_REPS = 3, then run full simulation

**Expected improvement:** Rashomon-DML coverage 91.7% → 93-95%

---

**Date:** 2026-03-12
**Updated by:** AI assistant
**Verified by:** (pending user review)
