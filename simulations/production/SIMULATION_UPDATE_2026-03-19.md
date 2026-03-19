# Simulation Update: CV-Based Tuning Parameter Selection

**Date:** 2026-03-19
**Issue:** Low coverage rates (84.7% fold-specific, 77.5% Rashomon vs 95% expected)
**Root Cause:** Fixed λ = log(n)/n may be overregularizing → underfit trees → biased estimates

---

## Changes Made

### Previous Approach (Theory-Only)
```r
estimate_att(
  X, A, Y, K = 5,
  regularization = log(n) / n,              # Fixed theory formula
  use_rashomon = FALSE
)

estimate_att(
  X, A, Y, K = 5,
  regularization = log(n) / n,              # Fixed theory formula
  use_rashomon = TRUE,
  rashomon_bound_multiplier = 2 * sqrt(log(n) / n),
  auto_tune_intersecting = FALSE            # No automatic tuning
)
```

**Problem:**
- λ = log(n)/n assumes smooth functions (minimax-optimal for Lipschitz functions)
- Our DGPs have **step functions** and **interactions**
- Fixed λ may be too large → trees too simple → bias → low coverage

---

### New Approach (CV + Auto-Tuning)
```r
estimate_att(
  X, A, Y, K = 5,
  cv_regularization = TRUE,     # AUTO-SELECT lambda via CV
  cv_K = 5,                     # 5-fold CV for lambda selection
  use_rashomon = FALSE
)

estimate_att(
  X, A, Y, K = 5,
  cv_regularization = TRUE,     # AUTO-SELECT lambda via CV
  cv_K = 5,                     # 5-fold CV for lambda selection
  use_rashomon = TRUE,
  rashomon_bound_multiplier = 2 * sqrt(log(n) / n),  # Initial guess
  auto_tune_intersecting = TRUE  # AUTO-TUNE to ensure intersections exist
)
```

**Benefits:**
1. **Data-adaptive λ:** CV finds the best penalty for each DGP/nuisance function
2. **Better tree fit:** Less bias in propensity/outcome models
3. **Improved coverage:** Confidence intervals should approach 95%
4. **Auto-tuning for Rashomon:** Automatically increases ε_n if intersections empty

---

## Files Updated

All 4 DGP batch files:
- `simulations/production/run_dgp1_batch.R` (binary features)
- `simulations/production/run_dgp2_batch.R` (continuous features)
- `simulations/production/run_dgp3_batch.R` (moderate complexity)
- `simulations/production/run_dgp4_batch.R` (weak overlap stress test)

---

## Expected Impact

### Coverage Improvement
- **Fold-specific:** 84.7% → ~93-95% (closer to nominal)
- **Rashomon:** 77.5% → ~90-93% (excluding DGP4)
- **DGP4:** Still challenging (weak overlap), but should improve

### RMSE
- May increase slightly (less regularization → more variance)
- But bias should decrease substantially
- Overall MSE should improve

### Computational Cost
- ~2-3× slower per simulation (nested CV: outer 5-fold DML, inner 5-fold lambda selection)
- From ~6 hours → ~12-18 hours for full batch
- Worth it for valid inference

---

## Running Updated Simulations

### Quick Test (n=400, 50 reps):
```bash
cd simulations/production
Rscript -e "
source('run_dgp2_batch.R')
# Edit N_REPS from 500 to 50 for quick test
"
```

### Full Batch (all DGPs, n=400/800/1600, 500 reps each):
```bash
cd simulations/production

# Run each DGP batch (can parallelize across DGPs)
nohup Rscript run_dgp1_batch.R > logs/dgp1_2026-03-19.log 2>&1 &
nohup Rscript run_dgp2_batch.R > logs/dgp2_2026-03-19.log 2>&1 &
nohup Rscript run_dgp3_batch.R > logs/dgp3_2026-03-19.log 2>&1 &
nohup Rscript run_dgp4_batch.R > logs/dgp4_2026-03-19.log 2>&1 &
```

---

## Validation Plan

After simulations complete:

1. **Check coverage:** Should be ~93-95% (fold-specific), ~90-93% (Rashomon)
2. **Check bias:** Should be minimal (< 0.01)
3. **Compare to oracle:** Gap should narrow (oracle has perfect models)
4. **Inspect selected λ values:** Should be smaller than log(n)/n for step functions

---

## Theoretical Justification

**Why CV is appropriate here:**

1. **DML framework protects from overfitting:**
   - Cross-fitting ensures nuisance functions trained on independent data
   - Neyman-orthogonality provides bias correction
   - CV for λ doesn't invalidate these guarantees

2. **Practical necessity:**
   - Theory assumes smooth functions (Lipschitz, Holder)
   - Our DGPs have **discontinuities** (step functions, thresholds)
   - Data-adaptive tuning is essential for non-smooth functions

3. **Evidence from oracle:**
   - Oracle coverage = 94.9% (validates DML implementation)
   - Gap to fold-specific = 10.2% (suggests nuisance model bias)
   - CV should close this gap

---

## Next Steps

1. **Monitor first batch:** Check convergence rate, coverage, bias
2. **If coverage still low:** Consider different tree implementations (BART, XGBoost)
3. **If coverage good:** Update manuscript with CV-based results
4. **Document selected λ:** Report median/IQR of CV-selected values per DGP

---

## Related Issue

This update addresses the concern that fixed λ = log(n)/n may be inappropriate for our DGPs. The code review fixes (Issue #1-15) ensure that when CV selects a better λ, the estimates will be **numerically stable and error-free** (no silent corruption, no catastrophic fallbacks).
