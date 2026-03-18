# β Smoothness Regime Study - Implementation Status

**Date:** 2026-03-11
**Status:** Implementation complete, initial testing reveals complexity mismatch

---

## Implementation Summary

### Files Created/Modified

1. **NEW: `dgps/dgps_beta_continuous.R`**
   - Continuous features (X ~ Uniform[0,1]^4) instead of manual discretization
   - Three DGPs: β=3 (cubic), β=2 (quadratic), β=1 (absolute value)
   - Theoretical s_n predictions: n^(d/(2β+d))

2. **MODIFIED: `run_beta_study.R`**
   - Sources new continuous DGP file
   - Adds `discretize_method="quantiles"`, `discretize_bins="adaptive"` to dml_att()
   - Captures tree complexity via `count_tree_leaves()` helper function
   - Returns `n_leaves_e`, `n_leaves_m0`, `theoretical_sn` in results

3. **MODIFIED: `analyze_beta_study.R`**
   - Adds Figure 0: Tree complexity verification (fitted vs theoretical s_n)
   - Includes s_n verification table at n=800
   - Extended summary statistics with tree complexity

4. **NEW: `verify_beta_dgps.R`**
   - Verification script for continuous DGPs
   - ✓ All checks passed (features continuous, ATT stable, overlap good)

5. **NEW: `test_beta_study.R`**
   - Small test (10 reps) to verify full pipeline

---

## Verification Results

### DGP Verification (n=800) ✓

| Regime | β | Theoretical s_n | X Continuous? | ATT  | Overlap |
|--------|---|----------------|---------------|------|---------|
| High   | 3 | 14.5 leaves    | ✓ (800 unique)| 0.10 | ✓ [0.63, 0.90] |
| Boundary | 2 | 28.3 leaves  | ✓ (800 unique)| 0.10 | ✓ [0.63, 0.85] |
| Low    | 1 | 86.2 leaves    | ✓ (800 unique)| 0.10 | ✓ [0.65, 0.89] |

All structural checks passed. Features are continuous, ATT is stable at 0.10, propensity overlap is good.

---

## Test Simulation Results (10 reps, n=800)

### Convergence
- **100%** (30/30 simulations converged)
- ✓ Tree complexity successfully captured (no longer NaN)

### Tree Complexity: Fitted vs Theoretical

| Regime | β | Theoretical s_n | Fitted e(X) | Fitted m0(X) | Ratio e | Ratio m0 |
|--------|---|----------------|-------------|--------------|---------|----------|
| High   | 3 | 14.5           | 61.0        | 13.7         | 4.21    | 0.948    |
| Boundary | 2 | 28.3         | 61.8        | 13.8         | 2.18    | 0.487    |
| Low    | 1 | 86.2           | 62.0        | 14.0         | 0.719   | 0.162    |

**Interpretation:**
- **Propensity trees (e):** Consistently ~61-62 leaves across ALL regimes (much larger than predicted)
- **Outcome trees (m0):** Consistently ~13-14 leaves across ALL regimes (smaller than predicted)
- **Fixed regularization λ = log(n)/n** does not produce theory-predicted s_n

**Why this happens:**
- Theory: s_n ~ n^(d/(2β+d)) is the OPTIMAL complexity given smoothness β
- Practice: λ = log(n)/n is a minimax-optimal choice but not tuned to specific β
- Result: Trees don't match theoretical predictions

**Is this a problem?**
- Not necessarily a bug - it's an empirical finding
- Shows gap between theory (optimal s_n) and practice (what fixed λ produces)
- Study can still test whether β < d/2 manifests, even if s_n doesn't match theory
- **BUT:** Propensity trees being ~4x too large for β=3 suggests possible overfitting

### Coverage (10 reps, n=800)

| Regime | β | Coverage | Mean CI Width |
|--------|---|----------|---------------|
| High   | 3 | 0%       | 0.0825        |
| Boundary | 2 | 30%    | 0.0855        |
| Low    | 1 | 10%      | 0.0837        |

⚠️ **CONCERNING:** Coverage is 0-30% (expected ~95%). Possible causes:
1. **Small sample:** Only 10 reps - could be noise
2. **Overfitting propensity:** e(X) trees with 61 leaves might be unstable
3. **Implementation issue:** Something fundamentally wrong with DGPs or estimation

**Next steps:** Run 50-rep test to see if coverage improves with more replications.

---

## Decision Point

### Option A: Proceed with Full Simulation (500 reps)
- **Pros:** Full results will clarify whether low coverage is noise or real
- **Cons:** 10-12 hours runtime; if there's a bug, we waste time

### Option B: Run 50-Rep Test First
- **Pros:** Better coverage estimates (1-2 hours), safer
- **Cons:** Delays full results

### Option C: Debug Coverage Issue First
- **Pros:** Ensures no fundamental bug before full run
- **Cons:** May be wild goose chase if it's just small-sample noise

---

## Recommended Next Step

**Run 50-rep test to diagnose coverage:**
1. Change `N_REPS <- 50` in `test_beta_study.R`
2. Run test (1-2 hours)
3. If coverage ≈ 90-95%: proceed with full 500-rep simulation
4. If coverage still < 80%: debug before full run

**If coverage is good with 50 reps:**
- Low coverage in 10-rep test was likely noise
- Complexity mismatch is real but doesn't prevent valid inference
- Document findings: "Fixed λ = log(n)/n produces trees larger/smaller than theory predicts"

**If coverage is still low with 50 reps:**
- Investigate: Are DGPs generating reasonable data? Is DML estimation correct?
- Check: Are propensity scores being used correctly?
- Consider: Reduce regularization for propensity trees (try λ/2)?

---

## Files Ready for Full Simulation

- ✓ `dgps/dgps_beta_continuous.R` - Continuous features, theory-aligned
- ✓ `run_beta_study.R` - Modified for continuous features + tree complexity
- ✓ `analyze_beta_study.R` - Includes s_n verification figure
- ✓ Tree complexity extraction working
- ⚠️ Coverage needs verification with larger test

---

## Expected Outcomes (Theory)

If implementation is correct and coverage stabilizes at ~95%:

| Regime | β | s_n Prediction | Expected Coverage | Expected Outcome |
|--------|---|----------------|-------------------|------------------|
| High   | 3 | n^0.4 ≈ 16     | ~95%              | Valid regime ✓   |
| Boundary | 2 | n^0.5 ≈ 28   | 92-96%            | Boundary case    |
| Low    | 1 | n^0.67 ≈ 70    | 85-92%            | Fails condition ✗|

If empirical s_n doesn't match predictions but coverage is good:
- Report: "Fixed regularization λ = log(n)/n produces consistent tree complexity across β regimes, not varying with smoothness as theory predicts"
- Study still valid: Tests whether performance degrades with β, even if s_n mechanism differs

---

## Summary

**Status:** Implementation technically complete, but coverage in small test is concerning.

**Next:** Run 50-rep test to determine if coverage issue is real or noise before committing to 500-rep simulation.
