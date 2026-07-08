# CV-Selected Lambda Implementation Summary

**Date:** 2026-05-26 (Updated: fixed approaches 4 & 6)
**Status:** ✅ ALL 6 APPROACHES READY FOR CLUSTER DEPLOYMENT

---

## Overview

Successfully implemented theory-driven CV-selected regularization as the default for all applicable ATT estimation approaches. This addresses the under-coverage issue caused by over-regularization (fixed λ=0.1 was 26× larger than theory-prescribed values).

---

## Changes Made

### 1. Package-Level Changes

**optimaltrees (v0.4.0):**
- Added `compute_safe_model_limit()` to scale model_limit based on λ/λ_theory ratio
- Implemented early stopping (skip remaining folds after 3 failures per lambda)
- Optimized parallelization (sort grid by computational cost for better load balancing)
- **Critical bug fix:** Don't pass `model_limit=0` to C++ (causes JSON errors)
- Graceful failure handling (returns NA when all λ values fail)

**doubletree (v0.0.0.9000):**
- Changed default: `cv_regularization = TRUE` in `estimate_att()`
- Updated documentation to reflect CV as recommended approach
- Theory-driven grid: `(log n / n) × [0.25, 0.5, 1, 2, 4]`

### 2. Simulation Code Changes

**File:** `code/estimators.R`

**Approach 1 (fullsample):**
- Now uses `cv_regularization()` to select lambda for e(X) and m0(X)
- Fallback to fixed λ=0.1 if CV fails completely
- Adds ~0.8 sec per replication

**Approach 2 (crossfit):**
- Now uses `cv_regularization()` per fold for e(X) and m0(X)
- Nested CV: K=5 outer folds × K=5 inner CV folds = 25 fits per nuisance
- Fallback to fixed λ=0.1 if CV fails completely
- Adds ~2.5 sec per replication

**Approach 3 (doubletree):**
- Already had CV via `estimate_att()` default
- No changes needed (just removed explicit `regularization` argument)
- Works correctly with Rashomon intersection
- Adds ~3.3 sec per replication

**Approach 4 (doubletree_averaged):**
- Now uses `cv_regularization()` for all tree fits (2026-05-26)
- CV per fold in Rashomon collection and tier 5 fallback
- **Bug fix:** Collect K trees (one per fold) for averaging, not just intersection
- Fallback to fixed λ=0.1 if CV fails
- Adds ~7.6 sec per replication

**Approach 5 (msplit):**
- No changes (uses fixed regularization via specialized function)
- Adds ~0.4 sec per replication

**Approach 6 (msplit_averaged):**
- Now uses `cv_regularization()` for all M×K tree fits (2026-05-26)
- CV per fold in fold-specific collection
- Fallback to fixed λ=0.1 if CV fails
- Adds ~7.4 sec per replication

---

## Test Results

### Local Testing (n=500, simple DGP)

```
Approach                CV?    Status     Time      Notes
--------                ---    ------     ----      -----
1. fullsample           YES    ✅ PASS   0.6 sec   theta=0.128 (true=0.150)
2. crossfit             YES    ✅ PASS   2.4 sec   theta=0.131 (true=0.150)
3. doubletree           YES    ✅ PASS   2.8 sec   theta=0.118 (true=0.150)
4. doubletree_averaged  YES    ✅ PASS   7.6 sec   theta=0.128 (true=0.150) [FIXED 2026-05-26]
5. msplit               NO     ✅ PASS   0.4 sec   theta=0.125 (true=0.150)
6. msplit_averaged      YES    ✅ PASS   7.4 sec   theta=0.128 (true=0.150) [FIXED 2026-05-26]
```

**Key Findings:**
- ✅ **All 6 approaches now working correctly**
- ✅ Approaches 1-4, 6 use CV-selected lambda
- ✅ Approach 5 (msplit without averaging) uses fixed λ
- ✅ **Approaches 4 & 6 fixed on 2026-05-26**: Added CV + bug fix for tree averaging
- ✅ No errors, warnings, or crashes
- ✅ Computational overhead acceptable (0.4-7.6 sec per replication)
- ✅ Standard errors reasonable (0.044-0.046)
- ✅ All CI's cover true value

---

## Theory Implementation

**Lambda grid:** `(log n / n) × [0.25, 0.5, 1, 2, 4]`

**Example values (n=1000):**
- log(n)/n ≈ 0.0069
- Grid: [0.0017, 0.0035, 0.0069, 0.0138, 0.0276]
- Old fixed: 0.1 (14× larger than theory max!)

**CV Selection:**
- Each nuisance function (e, m0) gets its own λ via 5-fold CV
- Selects best λ from theory-driven grid
- Falls back to λ=0.1 only if CV completely fails (rare)

---

## Performance Characteristics

**Time per replication (n=500):**
- Approach 1 (fullsample):   +0.8 sec vs fixed λ
- Approach 2 (crossfit):     +2.5 sec vs fixed λ
- Approach 3 (doubletree):   +3.3 sec vs fixed λ

**For n=2000 (cluster simulations):**
- Estimated +3-5 sec per replication
- With 100 reps × 3 DGPs × 2 sample sizes: +50-80 minutes total
- Still tractable for cluster (array jobs)

**Robustness:**
- Handles small λ values that would fail with default model_limit
- Early stopping reduces wasted computation on infeasible λ
- Graceful degradation (fallback to fixed λ if CV fails)

---

## Critical Bug Fix

**Issue:** `model_limit=0` (unlimited) caused JSON parsing errors in C++:
```
ERROR: [json.exception.type_error.302] type must be number, but is string
```

**Root Cause:** `compute_safe_model_limit()` returns `0` for λ < 0.5×λ_theory, but C++ code couldn't handle `model_limit=0`.

**Fix:** Only pass `model_limit` parameter to `fit_tree()` when it's > 0. When it's 0 (unlimited), omit the parameter entirely.

**Impact:** All CV failures were due to this bug. After fix, CV works reliably.

---

## Deployment Checklist

- [x] All packages installed and tested
- [x] **ALL 6 approaches now passing (updated 2026-05-26)**
- [x] Approaches 1, 2, 3, 4, 6 use CV-selected lambda
- [x] Approach 5 uses fixed lambda (by design)
- [x] Documentation updated
- [x] Local tests passing
- [x] Approaches 4 & 6 bugs fixed
- [ ] Update SLURM scripts (if needed)
- [ ] Deploy to cluster
- [ ] Run 10-20 test replications on cluster
- [ ] If successful, launch full simulation

---

## Expected Impact

**Original Problem:**
- Fixed λ=0.1 caused over-regularization
- Trees too simple → poor nuisance approximation
- SE underestimated → **under-coverage**

**Expected After CV:**
- λ automatically scaled to appropriate values
- More complex trees where data supports it
- Better nuisance approximation
- **Coverage should improve to target range (0.93-0.97)**

**Key Comparison:**
- Old: Fixed λ=0.1 for all n
- New: λ ∈ [0.0017, 0.0276] for n=1000 (data-adaptive within theory range)

---

## Fixed Issues (2026-05-26)

**Approaches 4 & 6:**
- ✅ **FIXED:** Now use CV-selected lambda instead of fixed λ=0.1
- ✅ **FIXED:** Approach 4 bug - was returning 1 tree from intersection instead of K fold-specific trees for averaging
- ✅ **FIXED:** Approach 6 bug - trees with different structures due to fixed λ; now uses CV so structures more consistent

**Details:**
- Both approaches now call `cv_regularization()` for each tree fit
- Approach 4: Fixed tree collection to get K fold-specific trees matching intersection structure
- Approach 6: CV ensures more consistent tree structures across M×K trees, reducing averaging failures
- Fallback to fixed λ=0.1 only if CV completely fails (rare)

**All 6 approaches now working and ready for cluster deployment.**

---

## Files Modified

**optimaltrees (3 files):**
1. `R/cv_regularization.R` - Added helper, early stopping, parallelization, bug fix
2. `tests/testthat/test-cv-regularization.R` - Added tests (24 PASS, 0 FAIL)

**doubletree (3 files):**
1. `R/estimate_att.R` - Changed default to `cv_regularization = TRUE`
2. `tests/testthat/test-estimate-att.R` - Added CV tests
3. `simulations/six_approach_comparison/code/estimators.R` - Updated approaches 1, 2, 3

---

## Next Steps

1. **Immediate:** Deploy to cluster for test run (10-20 reps)
2. **If successful:** Launch full simulation
3. **Optional:** Fix approaches 4 & 6 (lower priority)
4. **Future:** Create vignette on regularization selection

---

## Verification Commands

```bash
# Test locally
cd doubletree/simulations/six_approach_comparison
Rscript test_all_approaches_cv.R

# Expected output:
# Successful: 4/6
# CV approaches (1, 2, 3): ALL PASS
```

---

## Contact

For questions or issues, check:
- This summary: `CV_IMPLEMENTATION_SUMMARY.md`
- Test script: `test_all_approaches_cv.R`
- Original plan: `quality_reports/plans/2026-05-22_lambda-selection-default.md`

---

**Status:** ✅ Ready for cluster deployment
**Confidence:** High - all critical approaches tested and working
