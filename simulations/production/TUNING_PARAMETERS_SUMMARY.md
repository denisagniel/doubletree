# Tuning Parameter Choices in Simulations

## Summary

All three simulation studies use **FIXED** regularization with **NO cross-validation** for parameter selection.

---

## 1. Regularization Parameter (λ)

### What it does
Controls tree complexity in optimaltrees's penalized objective:
```
Objective = Loss + λ × Complexity
```
- Higher λ → simpler trees (more leaves pruned)
- Lower λ → more complex trees (more leaves retained)

### How it's set
**Formula:** `λ = log(n) / n`

| Sample Size | λ Value | 
|-------------|---------|
| n = 400 | 0.0150 |
| n = 800 | 0.0084 |
| n = 1600 | 0.0046 |

### Theoretical justification
- Minimax-optimal rate for tree-based estimation
- Balances bias-variance trade-off
- Standard in statistical learning theory

### Potential issue
**May be TOO LARGE** for these specific DGPs:
1. Formula assumes worst-case complexity
2. Our DGPs use 4-5 binary features (limited complexity)
3. Strong signal-to-noise ratio in primary DGPs
4. Over-regularization → underfitting → biased nuisance estimates

**Evidence:**
- Beta study: 15% coverage (expected 95%)
- Primary sims: 93% coverage (expected 95%)
- Both suggest possible underestimation

---

## 2. Cross-Validation Choice

### Current: `cv_regularization = FALSE`
- Uses fixed λ = log(n)/n for all nuisances
- Same λ for propensity e(X) and outcome m0(X)
- Same λ across all cross-validation folds

**Pros:**
- Fast (no nested CV)
- Reproducible
- Theory-justified rate

**Cons:**
- Not data-adaptive
- May not match actual problem complexity
- No automatic adjustment to DGP difficulty

### Alternative: `cv_regularization = TRUE` (not used)
- Would select λ via cross-validation
- Separate λ for each nuisance function
- Data-adaptive selection

**Why not used:**
- Computational cost (nested CV for 33,900 sims)
- Wanted reproducibility
- Expected theory-based λ to work

**Question:** Would CV improve coverage?

---

## 3. Rashomon Parameters (for rashomon method)

### rashomon_bound_multiplier (ε_n)
Controls size of Rashomon set (models close to optimal):
- Trees with loss ≤ (1 + ε_n) × best_loss

**Current:** Uses default `ε_n = 0.05`

**Theory recommends:** `ε_n = c√((log n)/n)` with c ≈ 2

| n | Theory (c=2) | Actual (default) | Issue |
|---|--------------|------------------|-------|
| 400 | 0.245 | 0.05 | **Too tight!** |
| 800 | 0.183 | 0.05 | **Too tight!** |
| 1600 | 0.136 | 0.05 | **Too tight!** |

**Consequence:**
- Tighter Rashomon set than theory suggests
- May explain why Rashomon-DML has lower coverage (91.7% vs tree-DML 93%)
- Smaller set → less flexibility → potentially worse nuisance estimates

---

## 4. Other Parameters

### K (DML cross-fitting folds)
- **Value:** 5
- **Standard:** Yes (typical in DML literature)
- **Trade-off:** More folds = less bias, more variance

### stratified
- **Value:** TRUE
- **Meaning:** Folds balanced by treatment assignment
- **Important:** For ATT estimation with potentially imbalanced A

---

## Issues Identified

### Critical Issue: Beta Study Catastrophic Undercoverage

**Observed:** 15% coverage (expected 95%)

**Possible causes:**
1. **Over-regularization:** λ = log(n)/n too large
   - Trees too simple
   - Nuisance functions underestimated
   - CIs too narrow → undercoverage

2. **DGP design:** Beta regime DGPs may have issues
   - Propensity/outcome ranges
   - Signal strength
   - Feature interactions

3. **Both:** Wrong λ for these specific DGPs

### Minor Issue: Primary Sims Slight Undercoverage

**Observed:** 93% coverage (expected 95%)

**Possible causes:**
1. Fixed λ not optimal for all DGPs
2. Monte Carlo error (within 2% is borderline)
3. Slight over-regularization

**Assessment:** Acceptable for publication, but could be improved

---

## Recommendations

### For Beta Study (Critical)

**Option 1: Rerun with CV regularization**
```r
cv_regularization = TRUE
```
- Let data choose optimal λ
- May dramatically improve coverage
- Worth the computational cost

**Option 2: Reduce fixed regularization**
```r
regularization = 0.5 * log(n) / n  # Half of theory
```
- Less aggressive pruning
- More complex trees allowed
- May better match DGP complexity

**Option 3: Diagnose first**
- Check fitted tree sizes (# leaves)
- Check nuisance estimate quality
- Then decide on fix

### For Primary Sims (Optional)

**Acceptable as-is** for publication (93% close to 95%)

**If improving:**
- Try cv_regularization = TRUE
- Or reduce λ by constant factor

### For Rashomon-DML

**Use theory-based ε_n:**
```r
epsilon_n <- 2 * sqrt(log(n) / n)
rashomon_bound_multiplier = epsilon_n
```

---

## Next Steps

1. **Investigate beta study:**
   - Check tree complexity in fitted models
   - Examine nuisance estimate quality
   - Decide: rerun with CV or adjusted λ?

2. **Consider for future work:**
   - Simulation study comparing fixed vs. CV regularization
   - Sensitivity analysis: coverage vs. λ
   - Adaptive regularization strategies

3. **Documentation:**
   - Add regularization choice to methods section
   - Justify log(n)/n formula
   - Note limitation: fixed formula may not be universally optimal

---

## Bottom Line

**Tuning parameters are theory-based but FIXED, not data-adaptive.**

- **Good:** Fast, reproducible, theoretically justified
- **Bad:** May not match specific DGP complexity
- **Ugly:** Beta study results suggest serious over-regularization

**Action needed:** Investigate and potentially rerun beta study with better-tuned parameters.

