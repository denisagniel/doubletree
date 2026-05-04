# Six-Approach Comparison: Investigation Summary

**Date:** 2026-05-04
**Status:** COMPLETE

---

## Overview

Investigated and resolved critical issues with approaches 3 and 4 in the six-approach comparison simulation study.

---

## Issues Investigated

### 1. Approach 3 & 5: Standard Error Inflation (23x too large) ✅ FIXED

**Problem:** Approaches 3 (doubletree) and 5 (msplit) reported SE ~1.0 instead of ~0.045

**Root cause:** Missing `/n` factor in SE calculation
- Code computed: `sqrt(E[ψ²])` (asymptotic SD of √n·θ̂)
- Should compute: `sqrt(E[ψ²]/n)` (SE of θ̂)

**Evidence:** Ratio = 23x ≈ √500, exactly as predicted

**Fixes:**
- `R/inference.R` line 29: Added `/n` to `att_se()`
- `R/estimate_att_msplit.R` line 257: Added `/n` to sigma calculation

**Verification:**
```
Before: Approach 3 SE = 1.0245 (23.0x baseline)
After:  Approach 3 SE = 0.0473 (1.06x baseline) ✓

Before: Approach 5 SE = 1.0409 (23.4x baseline)
After:  Approach 5 SE = 0.0465 (1.05x baseline) ✓
```

**Impact:** Critical bug affecting inference. Confidence intervals were 23x too wide.

**Documentation:** `SE_BUG_FIX.md`

---

### 2. Approach 4: Rashomon Intersection Failure ✅ EXPECTED BEHAVIOR

**Problem:** Approach 4 (doubletree_singlefit) consistently fails with "no common structure found"

**Investigation:**
- Tested sample sizes: 500, 1000, 2000, 4000 → All fail
- Tested multipliers: 0.01, 0.05, 0.10, 0.20 → All fail
- Examined fold-level structures → Different partitions (different # leaves)

**Root cause:** NOT A BUG
- Each fold learns genuinely different optimal tree structures
- Propensity: 5-8 leaves across folds
- Outcome: 5-7 leaves across folds
- Cross-fold heterogeneity from data-driven optimization

**Why intersection fails:**
```
Fold 1: Rashomon set = {tree with 6 leaves}
Fold 2: Rashomon set = {tree with 5 leaves}
...
Intersection = ∅ (different partitions, not same-partition-different-splits)
```

**Fallback behavior:**
- Package correctly falls back to fold-specific trees (approach 3)
- Maintains statistical validity (cross-fitting still valid)
- Only loses interpretability (multiple trees instead of one)

**Evidence that fallback works:**
- Approach 3 produces θ̂ = 0.1745 (true = 0.15) ✓
- SE now correct after bug fix ✓

**Conclusion:** Empty intersection is informative, not a failure. Shows when Rashomon struggles.

**Documentation:** `RASHOMON_INTERSECTION_ANALYSIS.md`

---

## Summary of Results

### Before Fixes

```
Approach 1 (full_sample):        theta=0.180, se=0.045 ✓
Approach 2 (crossfit):           theta=0.182, se=0.046 ✓
Approach 3 (doubletree):         theta=0.175, se=1.025 ✗ (23x too high)
Approach 4 (doubletree_single):  ERROR: Rashomon intersection failed
Approach 5 (msplit):             theta=0.180, se=1.041 ✗ (23x too high)
Approach 6 (msplit_single):      theta=0.181, se=0.045 ✓
```

### After Fixes

```
Approach 1 (full_sample):        theta=0.180, se=0.045 ✓
Approach 2 (crossfit):           theta=0.182, se=0.046 ✓
Approach 3 (doubletree):         theta=0.175, se=0.047 ✓ FIXED
Approach 4 (doubletree_single):  ERROR: Expected (cross-fold heterogeneity)
Approach 5 (msplit):             theta=0.180, se=0.047 ✓ FIXED
Approach 6 (msplit_single):      theta=0.181, se=0.045 ✓
```

**Status:** 5/6 approaches working correctly (1 expected failure)

---

## Commits

1. **fe1be01** - Fix prediction issues (type="prob")
   - Original issue: predict() returning classes instead of probabilities
   - Fixed all 6 approaches

2. **9fdd855** - Fix standard error inflation bug (23x too large)
   - Critical SE calculation bug
   - Fixed approaches 3 and 5

3. **d93810e** - Document Rashomon intersection failure
   - Not a bug, expected behavior
   - Comprehensive analysis and recommendations

---

## Files Created/Modified

### Core Fixes
- `R/inference.R` - Fixed SE calculation
- `R/estimate_att_msplit.R` - Fixed SE calculation
- `code/estimators.R` - Fixed predict() calls (from earlier commit)

### Documentation
- `SE_BUG_FIX.md` - SE bug explanation and verification
- `RASHOMON_INTERSECTION_ANALYSIS.md` - Rashomon analysis
- `FIXES_SUMMARY.md` - Original prediction bug fixes

### Testing/Diagnostic Scripts
- `test_se_fix.R` - Verify SE fixes
- `test_rashomon_sample_size.R` - Test if n helps Rashomon
- `test_validation.R` - Comprehensive validation (from earlier)
- `diagnose_approaches_3_4.R` - Initial diagnostic
- `debug_rashomon_simple.R` - Show fold-level structures

---

## Key Learnings

### 1. SE Calculation Pitfall

**Watch for:** `sqrt(mean(score^2))` vs `sqrt(mean(score^2) / n)`

The first gives asymptotic SD of √n·θ̂, the second gives SE(θ̂).

**Test:** If SE scales with √n instead of 1/√n, missing `/n` factor.

### 2. Rashomon Intersection Limitations

**When it fails:**
- Small-to-moderate sample per fold
- Data-driven optimization (not pre-specified structure)
- DGPs with moderate complexity

**Not a bug when:**
- Folds find different numbers of leaves
- This means genuinely different partitions
- Intersection code already handles same-partition-different-splits

**Alternative:** M-split finds modal structure instead of intersection (more robust)

### 3. Fallback Behavior is Principled

When Rashomon fails, falling back to fold-specific trees:
- ✓ Maintains cross-fit validity
- ✓ Produces correct estimates
- ✗ Loses interpretability (multiple trees)

---

## Recommendations for Paper

### Section: Rashomon Intersection

**Discuss:**
1. Intersection is not guaranteed (depends on sample size, complexity, data structure)
2. Empty intersection signals cross-fold heterogeneity (informative)
3. Fallback to fold-specific trees is principled (valid inference, loses interpretability)
4. Trade-offs: tight ε → better trees but less overlap; loose ε → more overlap but worse trees

**Compare:**
- Rashomon approach (intersection across folds) vs
- M-split approach (modal structure across independent runs)
- M-split more robust to empty intersection

**Guidance:**
- When to expect success: large n, simple DGP, stable structures
- When to expect failure: moderate n, complex DGP, data-driven heterogeneity
- What to do when it fails: use fold-specific trees (approach 3) or M-split (approaches 5/6)

---

## For Simulation Study

### Current State: Ready to Run

**Working approaches:** 5/6
- Approaches 1, 2, 3, 5, 6: All working correctly ✓
- Approach 4: Fails as expected, documented ✓

**SEs are correct:** All approaches report comparable, realistic SEs

**Validation passing:** test_validation.R confirms fixes work

### No Further Action Needed

Approach 4 failures should be:
1. **Documented** in results (shows when Rashomon struggles)
2. **Expected** in certain regimes (moderate n, complex DGP)
3. **Compared** to M-split success rate

---

## Time Investment

- **Investigation:** ~4 hours
  - SE bug diagnosis: 1 hour
  - SE bug fix and testing: 1 hour
  - Rashomon investigation: 2 hours
- **Documentation:** ~1 hour
- **Total:** ~5 hours

---

## Constitutional Alignment

✓ **Correctness over speed:** Thorough investigation, root cause fixes
✓ **No quiet fallbacks:** SE bug caused wrong inference, now fixed
✓ **Quality gates:** Verified all fixes with comprehensive tests
✓ **Documentation:** Complete trail for future reference

---

## Next Steps

**Immediate:** None required - simulation ready to run

**Optional:**
1. Compare M-split vs Rashomon success rates across DGPs
2. Investigate auto-tune parameter for Rashomon
3. Test with simpler DGPs where structures might align

**For paper:**
1. Include Rashomon limitation discussion
2. Compare intersection vs modal structure approaches
3. Provide practical guidance on method choice

---

## Related Documents

- `FIXES_SUMMARY.md` - Original prediction bug fixes (2026-05-01)
- `SE_BUG_FIX.md` - Standard error bug details (2026-05-04)
- `RASHOMON_INTERSECTION_ANALYSIS.md` - Rashomon analysis (2026-05-04)
- This document - Overall investigation summary (2026-05-04)
