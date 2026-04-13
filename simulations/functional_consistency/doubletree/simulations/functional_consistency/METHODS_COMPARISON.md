# Three Methods Comparison

**Goal:** Empirically test three approaches to achieving functional consistency

---

## Methods

### 1. Standard M-Split (Baseline)

**Algorithm:**
1. Run M splits with K-fold cross-fitting
2. Select modal tree structure
3. Refit on each fold's training data
4. Predict on test fold only (cross-fit)
5. Average predictions across M splits

**Properties:**
- ✓ No self-prediction bias (true cross-fitting)
- ✗ Imperfect functional consistency: max|η̂(X_i) - η̂(X_j)| = O(1/√n_x) when X_i = X_j
- ✓ Valid asymptotic inference

**Mechanism of FC gap:**
- Observations i, j with X_i = X_j may be in different folds
- Different training sets → different trees → different predictions
- Gap vanishes as n → ∞, but slowly

---

### 2. Averaged Tree

**Algorithm:**
1. Run M splits with K-fold cross-fitting
2. Select modal tree structure
3. Refit on each fold's training data
4. **Predict on FULL dataset** from each refit (breaks cross-fitting)
5. Average ALL M×K predictions

**Properties:**
- ✓ Perfect functional consistency: max|η̄(X_i) - η̄(X_j)| = 0
- ✗ Self-prediction bias: η̄(X_i) includes predictions trained with Y_i
- ? Asymptotic inference validity unclear

**Mechanism of bias:**
- For observation i, (K-1)/K of predictions used Y_i in training
- Direct self-prediction: E[η̂(X_i) | Y_i in training] - η(X_i) = O(1/n_ℓ)
- If n_ℓ ∝ n^β with β ≤ 1/2, bias may dominate variance

---

### 3. Pattern Aggregation

**Algorithm:**
1. Run M splits with K-fold cross-fitting (standard M-split)
2. Get cross-fit predictions η̂(X_i) for all i
3. **Group by covariate pattern** x
4. **Average within pattern**: η̃(x) = mean{η̂(X_i) : X_i = x}
5. Replace: η̃(X_i) = η̃(x) for all i with X_i = x

**Properties:**
- ✓ Perfect functional consistency: max|η̃(X_i) - η̃(X_j)| = 0 (by construction)
- ✗ Indirect contamination bias: η̃(X_i) includes predictions from Y_j with X_j = X_i
- ? Asymptotic inference validity unclear

**Mechanism of bias:**
- Each η̂(X_i) is individually unbiased (cross-fit)
- But η̃(X_i) = (1/n_x) Σ_{j: X_j=x} η̂(X_j)
- On average, (n_x-1) × (K-1)/K of these used Y_j in training
- Since X_j = X_i, the Y_j are informative about Y_i
- Bias ≈ O(1/n_x)

---

## Comparison

| Property | Standard M-Split | Averaged Tree | Pattern Aggregation |
|----------|-----------------|---------------|---------------------|
| **Functional Consistency** | O(1/√n_x) | 0 | 0 |
| **Bias mechanism** | None | Direct self-prediction | Indirect contamination |
| **Bias order** | 0 | O(1/n_ℓ) | O(1/n_x) |
| **Group size** | N/A | n_ℓ (tunable) | n_x (fixed by data) |
| **Asymptotic inference** | ✓ Valid | ? Unclear | ? Unclear |

---

## Key Questions

### Q1: Does group size matter?

**Hypothesis:** Pattern aggregation has worse coverage than averaged tree because n_x < n_ℓ

**Test:**
- Compare coverage at n=3200 between methods
- Check if difference correlates with n_x vs n_ℓ
- Sparse DGP (small n_x) should show larger gap

### Q2: Is there a practical difference?

**Two possibilities:**

**A. Methods are equivalent** (similar coverage):
- Both have bias O(1/group_size)
- Direct vs indirect contamination doesn't matter empirically
- → Use whichever is simpler (pattern aggregation)

**B. Averaged tree is better** (better coverage):
- n_ℓ > n_x, so O(1/n_ℓ) < O(1/n_x)
- Tunable group size helps
- → Use averaged tree with appropriate regularization

### Q3: Does either maintain valid inference?

**Three outcomes:**

**A. Both fail** (coverage degrades):
- Confirms bias = O(n^{-β}) with β ≤ 1/2 dominates
- √n-consistency fails for both
- → Trade-off: perfect FC vs valid inference

**B. Both succeed** (coverage stays 95%):
- Surprising! Maybe:
  - β > 1/2 (leaves grow fast)
  - Partial cross-fitting sufficient
  - Finite sample artifacts
- → Use either method, no cost for perfect FC!

**C. One succeeds, one fails**:
- Suggests mechanism matters (direct vs indirect)
- → Use the successful method

---

## DGP Impact

### Simple DGP (8 patterns, n_x ≈ n/8)

Expected:
- n_x relatively large → pattern aggregation bias small
- n_ℓ ≈ n_x in this setting (leaves often = patterns)
- **Methods should perform similarly**

### Complex DGP (16 patterns, n_x ≈ n/16)

Expected:
- n_x smaller → pattern aggregation bias larger
- n_ℓ can still span multiple patterns
- **Averaged tree should outperform pattern aggregation**

### Sparse DGP (32 patterns, n_x ≈ n/32)

Expected:
- n_x very small → pattern aggregation bias large
- May fail to maintain valid inference even at n=3200
- **Strong test of whether n_ℓ > n_x matters**

---

## Implications

### If pattern aggregation ≈ averaged tree:

**Conclusion:** Bias mechanism (direct vs indirect) doesn't matter much

**Recommendation:** Use pattern aggregation (simpler, no need to refit M×K times)

### If averaged tree > pattern aggregation:

**Conclusion:** Tunable group size n_ℓ matters

**Recommendation:** Use averaged tree with appropriate regularization to maximize n_ℓ

### If both maintain coverage:

**Conclusion:** Perfect FC comes "for free"!

**Impact:** Major practical advantage for applications requiring FC (fairness, policy targeting)

---

## Simulation Design

**Grid:**
- n: {200, 400, 800, 1600, 3200}
- K: {2, 3, 5}
- DGP: {simple, complex, sparse}
- method: {standard_msplit, averaged_tree, pattern_aggregation}

**Total:** 135 configurations, 500 reps each = 67,500 replications

**Critical comparisons:**
1. Coverage vs n by method (test degradation)
2. Coverage: averaged_tree vs pattern_aggregation (test if group size matters)
3. Functional consistency: all three methods (verify perfect FC)
4. Bias/SE ratio vs n (test if bias dominates)
