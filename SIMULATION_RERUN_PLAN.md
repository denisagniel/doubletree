# Simulation Re-Run Plan: CV-Based Parameter Selection

**Date:** 2026-03-19
**Motivation:** Low coverage (84.7% fold-specific, 77.5% Rashomon vs 95% expected)
**Solution:** Switch from fixed λ = log(n)/n to CV-based selection

---

## Problem Diagnosis

### Observed Issues
- **Oracle coverage:** 94.9% ✓ (validates DML implementation is correct)
- **Fold-specific coverage:** 84.7% ✗ (should be ~95%)
- **Rashomon coverage:** 77.5% ✗ (even worse, driven by DGP4)

### Root Cause
Fixed λ = log(n)/n is **theory-optimal for smooth functions** but our DGPs have:
- Step functions (threshold effects)
- Interactions between features
- Piecewise constant regions

Result: Trees are **underfit** → biased nuisance estimates → biased ATT → poor coverage

---

## Solution: CV + Auto-Tuning

### Changes Made

**Before (Theory-Only):**
```r
estimate_att(X, A, Y, K = 5,
  regularization = log(n) / n,              # Fixed
  use_rashomon = FALSE
)
```

**After (Data-Adaptive):**
```r
estimate_att(X, A, Y, K = 5,
  cv_regularization = TRUE,     # Let CV select best λ
  cv_K = 5,                     # 5-fold CV for λ
  use_rashomon = FALSE
)
```

**For Rashomon:**
```r
estimate_att(X, A, Y, K = 5,
  cv_regularization = TRUE,              # CV for λ
  cv_K = 5,
  use_rashomon = TRUE,
  rashomon_bound_multiplier = 2 * sqrt(log(n) / n),  # Initial
  auto_tune_intersecting = TRUE          # Auto-increase if needed
)
```

### Files Updated
✓ `simulations/production/run_dgp1_batch.R` (binary features)
✓ `simulations/production/run_dgp2_batch.R` (continuous features)
✓ `simulations/production/run_dgp3_batch.R` (moderate complexity)
✓ `simulations/production/run_dgp4_batch.R` (weak overlap stress test)

---

## Quick Test (Before Full Re-Run)

Run a quick diagnostic to confirm improvement:

```bash
cd /Users/dagniel/RAND/rprojects/global-scholars/doubletree

# Test CV improvement (100 sims, ~5-10 minutes)
Rscript simulations/test_cv_improvement.R
```

**Expected output:**
- Coverage improvement: +5% to +10%
- CV coverage: 90-95% (approaching nominal)
- Bias reduction: smaller |bias| with CV

**If test looks good:** Proceed to full re-run
**If test still poor:** May need different tree implementation (BART, XGBoost)

---

## Full Simulation Re-Run

### Setup
```bash
cd simulations/production

# Create logs directory
mkdir -p logs

# Verify code review fixes are in place
git log --oneline -1  # Should show: "Fix 15 critical code review issues..."
```

### Run Batches (Parallel)

Each batch: 500 reps × 3 sample sizes × 4 methods = 6,000 simulations

**Estimated time per batch:** 12-18 hours (was 6 hours, now 2-3× slower due to nested CV)

```bash
# DGP1: Binary features
nohup Rscript run_dgp1_batch.R > logs/dgp1_cv_2026-03-19.log 2>&1 &
PID1=$!

# DGP2: Continuous features
nohup Rscript run_dgp2_batch.R > logs/dgp2_cv_2026-03-19.log 2>&1 &
PID2=$!

# DGP3: Moderate complexity
nohup Rscript run_dgp3_batch.R > logs/dgp3_cv_2026-03-19.log 2>&1 &
PID3=$!

# DGP4: Weak overlap (stress test)
nohup Rscript run_dgp4_batch.R > logs/dgp4_cv_2026-03-19.log 2>&1 &
PID4=$!

echo "DGP1 PID: $PID1"
echo "DGP2 PID: $PID2"
echo "DGP3 PID: $PID3"
echo "DGP4 PID: $PID4"
```

### Monitor Progress

```bash
# Check logs
tail -f logs/dgp1_cv_2026-03-19.log

# Check convergence rates (should be high, >95%)
grep "Convergence rate" logs/*.log

# Check for errors
grep -i "error\|fail" logs/*.log | grep -v "squared_error"
```

---

## Expected Results

### Coverage Improvements
| Method | Old (Fixed λ) | New (CV λ) | Target |
|--------|---------------|------------|--------|
| Oracle | 94.9% | ~95% | 95% (benchmark) |
| Fold-specific | 84.7% | **92-95%** | 95% |
| Rashomon (DGPs 1-3) | ~86% | **90-93%** | 90-95% |
| Rashomon (DGP4) | 53.5% | **70-80%** | Difficult (weak overlap) |

### Other Metrics
- **Bias:** Should decrease (better fit nuisances)
- **RMSE:** May increase slightly (less regularization → more variance)
- **MSE:** Should improve overall (bias^2 reduction > variance increase)

---

## Validation Checklist

After simulations complete:

- [ ] **Coverage check:** Fold-specific ≥93%, Rashomon ≥90% (excluding DGP4)
- [ ] **Bias check:** Mean bias < 0.01 across DGPs
- [ ] **Oracle gap:** Gap between oracle and fold-specific narrows
- [ ] **Convergence:** >95% of simulations converge successfully
- [ ] **Lambda distribution:** Inspect CV-selected λ values (should be < log(n)/n)

---

## Post-Analysis

### Generate Updated Figures

```bash
cd simulations/production

# Combine results across DGPs
Rscript analyze_cv_results.R

# Expected outputs:
#   - simulations/figures_cv/figure1_coverage_comparison.pdf
#   - simulations/figures_cv/figure2_bias_rmse.pdf
#   - simulations/figures_cv/table1_summary_cv.csv
```

### Update Manuscript

1. Replace old simulation results (table1_summary.csv) with CV version
2. Update text: "Regularization parameters selected via 5-fold CV"
3. Add note: "CV selection ensures well-fit trees for step-function DGPs"
4. Report median CV-selected λ values per DGP in appendix

---

## Theoretical Justification for CV

**Why CV doesn't invalidate DML inference:**

1. **Cross-fitting protects against overfitting:**
   - Nuisance functions fitted on independent data (folds 1-4)
   - ATT estimated on held-out fold (fold 5)
   - Even if λ is "overfitted" via CV, cross-fitting prevents contamination

2. **Neyman-orthogonality provides bias correction:**
   - DML score is doubly robust
   - Small misspecification in nuisances → small bias in ATT
   - CV improves nuisance fit → reduces this misspecification

3. **Practical necessity:**
   - Theory assumes smooth functions (Hölder, Lipschitz)
   - Our DGPs violate this (step functions, interactions)
   - Data-adaptive tuning is essential for non-smooth settings

**Evidence:**
- Oracle (perfect nuisances) gets 95% coverage
- Gap to fold-specific (learned nuisances) is 10%
- CV should close this gap by improving nuisance fit

---

## If Coverage Still Poor After CV

If coverage doesn't improve to ≥93%, consider:

1. **More flexible models:**
   - BART (Bayesian Additive Regression Trees)
   - XGBoost (gradient boosting)
   - Random forests with probability calibration

2. **Sample splitting instead of cross-fitting:**
   - Fit nuisances on 50% of data
   - Estimate ATT on other 50%
   - May improve stability at cost of efficiency

3. **Standard error adjustment:**
   - Check if SE underestimated (Issue #15 validation helps)
   - May need robust sandwich estimator

---

## Contact for Issues

If simulations fail or results are unexpected:
- Check logs in `simulations/production/logs/`
- Verify code review fixes are in place (git log)
- Ensure optimaltrees package is up to date

---

## Timeline

- **Quick test:** 10 minutes (run now)
- **Full batch:** 12-18 hours per DGP (4 parallel streams)
- **Analysis:** 1-2 hours after completion
- **Total:** ~24 hours for full re-run + analysis

**Recommendation:** Run quick test first, then start overnight batch if results look promising.
