# Rashomon Intersection Failure Analysis

**Date:** 2026-05-04
**Status:** IN PROGRESS

---

## Problem Statement

Approach 4 (doubletree_singlefit) consistently fails with "Rashomon intersection failed - no common structure found" across multiple DGPs and `rashomon_bound_multiplier` values (0.01 to 0.20).

**Key observation:** Each fold finds only 1 tree in its Rashomon set, and these trees have different structures.

---

## Diagnostic Evidence

### Fold-Level Rashomon Sets

From diagnostic output:
```
Fold 1: Rashomon set size = 1 trees
Fold 2: Rashomon set size = 1 trees
Fold 3: Rashomon set size = 1 trees
Fold 4: Rashomon set size = 1 trees
Fold 5: Rashomon set size = 1 trees

Finding intersection across 5 Rashomon sets (by partition)...
Rashomon set sizes: 1, 1, 1, 1, 1
Note: Trees with same leaves but different split orders are considered equivalent
Starting with 1 unique partition(s) from fold 1
After intersecting with fold 2: 0 partitions remain
No common partitions appear in all folds
```

### Multiplier Sensitivity

Tested multipliers from 0.01 to 0.20:

| Multiplier | e_struct | m0_struct | theta | Notes |
|------------|----------|-----------|-------|-------|
| 0.01 | FALSE | FALSE | 0.1817 | Too tight |
| 0.05 | FALSE | FALSE | 0.1787 | Current default |
| 0.10 | FALSE | FALSE | 0.1799 | Moderate |
| 0.20 | FALSE | FALSE | 0.1706 | Loose |

**None succeed in finding common structure**, but all produce reasonable theta estimates via fallback to fold-specific trees.

---

## Root Causes

### 1. Small Per-Fold Sample Size

- Total n = 500
- Per fold: ~100 observations (K=5)
- For binary outcomes with log-loss, optimal tree is highly data-driven
- Small sample → high variance in optimal structure selection

### 2. DGP Complexity

**Simple DGP** (from code/dgps.R):
```r
# True propensity depends on X1, X2
e <- plogis(-1 + 2*X[,1] + X[,2])

# True outcome depends on X1, X3
m0 <- plogis(-0.5 + X[,1] + 1.5*X[,3])
```

With 3 binary features, there are multiple ways to approximate the same function:
- Split on X1 first vs X2 first
- Different split orders produce same partition
- Small sample differences push toward different structures

### 3. Optimization Landscape

For log-loss with binary outcomes:
- Many tree structures can achieve similar loss
- But loss is bumpy (discrete splits)
- Small sample differences → different local optima
- Even with 20% tolerance, Rashomon sets don't overlap

---

## Why Fallback Works

When intersection fails, doubletree falls back to **fold-specific trees**:
```
Rashomon intersection failed for: propensity, control outcome.
Using fold-specific trees for all nuisances.
```

This produces valid estimates because:
1. Each fold fits optimal tree on its training data
2. Predictions are made on held-out test data (cross-fitting)
3. EIF estimator is doubly robust
4. Only loses interpretability (multiple trees instead of one)

**Evidence:** Approach 3 (doubletree, uses fallback) produces theta=0.1745, very close to true ATT=0.15.

---

## Potential Solutions

### Solution 1: Increase Sample Size

**Idea:** Larger n → more stable optimal structures → higher chance of overlap

**Test:** Run with n = 1000, 2000, 5000

**Expected outcome:**
- Rashomon sets may grow (more trees near-optimal)
- Structures may stabilize (less cross-fold variance)
- But not guaranteed to help

**Tradeoff:** Simulation time increases linearly with n

### Solution 2: Increase Multiplier Aggressively

**Idea:** Use multiplier = 0.50 or 1.00 to include more sub-optimal trees

**Test:** See if very loose tolerance creates overlap

**Expected outcome:**
- Rashomon sets definitely grow
- More likely to find common structure
- But includes trees with 50-100% higher loss (undesirable)

**Tradeoff:** Loses optimality, may include bad trees

### Solution 3: Use Simpler DGPs

**Idea:** DGPs with single dominant split pattern → less ambiguity

**Example:**
```r
# Strong signal on X1 only
e <- plogis(-2 + 5*X[,1])  # Very strong X1 effect
m0 <- plogis(-1 + 4*X[,1])  # Same dominant feature
```

**Expected outcome:** All folds should split on X1 first → common structure

**Tradeoff:** Less realistic, doesn't stress-test method

### Solution 4: Post-hoc Structure Alignment

**Idea:** After finding fold-specific trees, identify most common structure

**Algorithm:**
1. Fit fold-specific trees (K trees)
2. Extract K structures
3. Find modal structure (most frequent partition)
4. Refit modal structure on full data

**Expected outcome:**
- Approach 4 becomes "modal structure + full-data refit"
- Similar to M-split approach (approach 6)
- Interpretability gain (one tree) without requiring intersection

**Tradeoff:** No longer uses Rashomon theory directly

### Solution 5: Adaptive epsilon_n Selection

**Idea:** Use `auto_tune_intersecting = TRUE` to automatically increase multiplier until intersection succeeds

**Current implementation:** Already available in `estimate_att()` parameter

**Test:**
```r
result <- estimate_att(...,
  use_rashomon = TRUE,
  auto_tune_intersecting = TRUE,
  rashomon_bound_multiplier = 0.05  # starting point
)
```

**Expected outcome:** Will find smallest epsilon that yields intersection

**Tradeoff:** May yield very large epsilon (e.g., 0.8) which includes poor trees

---

## Recommendations

### For This Simulation Study

**Option A: Accept the failure (recommended)**
- Approach 4 is an exploratory comparison
- Having 5/6 approaches work is sufficient
- Failure mode is informative (shows when Rashomon struggles)
- Fallback behavior is documented and valid

**Option B: Try Solution 1 (larger n)**
- Test with n = 2000 to see if structures stabilize
- Only for a subset of DGPs
- Document whether sample size helps

**Option C: Try Solution 5 (auto-tune)**
- Use existing parameter
- Document resulting epsilon values
- Compare to approaches without auto-tuning

### For Paper/Method

**Key points to discuss:**

1. **Rashomon intersection is not guaranteed**
   - Depends on: sample size, complexity, DGP structure
   - Empty intersection is a valid outcome (signals heterogeneity)

2. **Fallback is principled**
   - Fold-specific trees maintain cross-fit validity
   - Only loses interpretability, not statistical validity
   - This is by design (documented in code)

3. **Trade-offs are unavoidable**
   - Tight epsilon: better trees, fewer in Rashomon, less overlap
   - Loose epsilon: more overlap, worse trees, less appealing
   - No free lunch

4. **Alternative: M-split approach**
   - Approach 5/6 (M-split) finds modal structure across multiple splits
   - More robust to empty intersection
   - Different theoretical framework (multiple independent runs vs cross-fold intersection)

---

## Questions to Investigate

### Q1: Does n help?

**Hypothesis:** Larger sample → more stable structures → higher overlap probability

**Test:**
```r
test_sizes <- c(500, 1000, 2000, 5000)
for (n in test_sizes) {
  data <- generate_dgp_simple(n)
  result <- estimate_att_doubletree_singlefit(data$X, data$A, data$Y)
  # Check if structures found
}
```

**Prediction:** Might help at n=2000+, but not guaranteed

### Q2: Does simpler DGP help?

**Hypothesis:** Single dominant feature → all folds agree

**Test:**
```r
# Very strong X1 signal
generate_dgp_strong_signal <- function(n) {
  X <- data.frame(
    x1 = rbinom(n, 1, 0.5),
    x2 = rbinom(n, 1, 0.5),
    x3 = rbinom(n, 1, 0.5)
  )
  A <- rbinom(n, 1, plogis(-2 + 5*X$x1))  # Dominant X1
  Y0 <- rbinom(n, 1, plogis(-1 + 4*X$x1))  # Same dominant
  # ...
}
```

**Prediction:** Very likely to succeed (all split on X1)

### Q3: What epsilon is needed?

**Hypothesis:** auto_tune finds epsilon, but it's very large

**Test:**
```r
result <- estimate_att(..., auto_tune_intersecting = TRUE)
cat("Epsilon needed:", result$epsilon_n, "\n")
```

**Prediction:** epsilon ∈ [0.4, 1.0] if it succeeds

### Q4: How often does M-split succeed?

**Hypothesis:** M-split (approach 5/6) more robust than Rashomon

**Test:** Compare success rates across DGPs and sample sizes

**Prediction:** M-split succeeds more often (finds mode, not intersection)

---

## Action Items

**Immediate (to understand the issue):**
1. [x] Diagnose why intersection fails (completed - size=1 per fold)
2. [ ] Test Q1: Does larger n help?
3. [ ] Test Q3: What epsilon is needed with auto-tune?

**Short-term (for simulation):**
1. [ ] Decide: Accept failure, or try solution?
2. [ ] If trying: Implement and document
3. [ ] Update README with Rashomon limitations

**Long-term (for paper):**
1. [ ] Discuss empty intersection trade-offs
2. [ ] Compare Rashomon vs M-split approaches
3. [ ] Provide guidance on when each works

---

## Related Files

- `R/estimate_att.R`: Main doubletree implementation with Rashomon
- `R/estimate_att_msplit.R`: M-split alternative (modal structure)
- `code/estimators.R`: Simulation implementations
- `diagnose_approaches_3_4.R`: Diagnostic script
- This document: Analysis and recommendations

---

## Version History

- **2026-05-04**: Initial analysis and solution proposals
