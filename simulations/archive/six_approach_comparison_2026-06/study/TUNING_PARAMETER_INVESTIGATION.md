# Tuning Parameter Selection Investigation

**Date:** 2026-05-22
**Investigation:** Theory vs Implementation for Regularization Parameter (λ)

---

## Executive Summary

**Major discrepancy found between theory, paper claims, and actual implementation:**

| Source | Lambda Selection Method | Value |
|--------|------------------------|-------|
| **Theory (manuscript)** | λ ≈ (log n)/n | n=500: 0.0124, n=1000: 0.0069, n=2000: 0.0038 |
| **Paper claims** | 5-fold cross-validation | Variable (data-driven) |
| **Actual implementation** | **Fixed** | **λ = 0.1 for all n** |

**Gap:** Simulations use λ = 0.1, which is **8-26× larger** than theory requires. This means **much stronger regularization** (simpler trees) than theory suggests.

**Potential impact on under-coverage:**
- Over-regularized trees → underestimated model uncertainty → underestimated SEs → narrow CIs → under-coverage

---

## Theory: What the Manuscript Says

**Source:** `doubletree/inst/paper/manuscript.tex`

### Theoretical Requirement

**Rate condition (appears in multiple theorems):**
- λ_n ≈ (log n)/n (asymptotic rate)
- Lines 115, 124, 156, 422, 492

**Rationale (line 439):**
- Matches Xu et al. (2026) equation (7): λ ≥ C(M+K)(log(nd)+u)/(δn)
- Ensures oracle inequalities and convergence rates hold
- Balances approximation error and estimation error

**Values implied by theory:**
```
n=500:  λ ≈ log(500)/500 ≈ 0.0124
n=1000: λ ≈ log(1000)/1000 ≈ 0.0069
n=2000: λ ≈ log(2000)/2000 ≈ 0.0038
```

### Paper's Implementation Claim

**Line 183:**
> "The regularization parameter λ_n is chosen by 5-fold cross-validation on each training fold."

**Interpretation:**
- Suggests data-driven selection via CV
- Should produce values consistent with (log n)/n rate
- Implies different λ for different n

### Acknowledged Limitations

**Line 300:**
> "Tuning choices—including the penalty λ, the Rashomon tolerance ε_n, and the maximum number of leaves—affect both interpretability and the validity of inference; sensitivity to these choices is a limitation in applications."

- Paper acknowledges tuning affects inference validity
- But does not quantify sensitivity or test robustness

---

## Implementation: What the Code Actually Does

**Source:** `doubletree/simulations/six_approach_comparison/code/estimators.R`

### Fixed Lambda Across All Approaches

**All six approaches use:**
```r
regularization = 0.1  # Fixed parameter
```

**Relevant code lines:**
- Approach 1 (full_sample): Line 179, `regularization = 0.1`
- Approach 2 (crossfit): Line 245, `regularization = 0.1`
- Approach 3 (doubletree): Line 315, `regularization = 0.1`
- Approach 4 (doubletree_averaged): Line 390, `regularization = 0.1`
- Approach 5 (msplit): Line 430, `regularization = 0.1`
- Approach 6 (msplit_averaged): Line 471, `regularization = 0.1`

**Characteristics:**
- λ = 0.1 **for all approaches**
- λ = 0.1 **for all DGPs** (simple, moderate, complex)
- λ = 0.1 **for all sample sizes** (n=500, 1000, 2000)
- **Not data-driven**
- **Not CV-selected**

### CV Option Exists But Is Not Used

**Source:** `doubletree/R/estimate_att.R`

**Parameters available:**
```r
estimate_att(...,
  regularization = 0.1,           # Default fixed value
  cv_regularization = FALSE,       # Default: don't use CV
  cv_K = 5                         # CV folds if CV enabled
)
```

**Documentation (line 33):**
> "Theory: Manuscript recommends λ ∝ (log n)/n for minimax-optimal trees."

**Implementation choice:**
- Simulations use `cv_regularization = FALSE` (default)
- Simulations use fixed `regularization = 0.1`
- CV option exists but is not invoked

---

## What CV Would Do If Enabled

**Source:** `optimaltrees/R/cv_regularization.R`

### Default CV Grid (line 94-96)

When `cv_regularization = TRUE` and no grid specified:
```r
lambda_grid <- (log(n) / n) * c(0.25, 0.5, 1, 2, 4)
```

**CV grid by sample size:**

| n | log(n)/n | CV Grid |
|---|----------|---------|
| 500 | 0.0124 | [0.0031, 0.0062, 0.0124, 0.0248, 0.0496] |
| 1000 | 0.0069 | [0.0017, 0.0035, 0.0069, 0.0138, 0.0276] |
| 2000 | 0.0038 | [0.00095, 0.0019, 0.0038, 0.0076, 0.0152] |

**Observations:**
- CV grid is **theory-driven** (centered around log(n)/n)
- All grid values are **much smaller** than λ = 0.1
- Grid scales appropriately with n
- CV would select λ ≈ 0.003-0.05 (depending on data)

### CV Procedure

**Method:**
1. For each λ in grid:
   - Fit tree on K-1 folds
   - Compute held-out loss on fold K
   - Average across K folds
2. Select λ that minimizes average held-out loss
3. Refit on full data with selected λ

**Loss function:**
- Misclassification: mean(1 - correct)
- Log-loss: mean cross-entropy (with clipping)

---

## Theory-Practice Gap Analysis

### Comparison Table

| Aspect | Theory | Paper Claims | Actual Implementation |
|--------|--------|--------------|----------------------|
| **Selection method** | λ ≈ (log n)/n | 5-fold CV | Fixed λ = 0.1 |
| **Scaling with n** | Yes (decreases) | Yes (via CV) | No (constant) |
| **Value at n=500** | 0.0124 | ~0.003-0.05 (CV) | 0.1 |
| **Value at n=1000** | 0.0069 | ~0.001-0.03 (CV) | 0.1 |
| **Value at n=2000** | 0.0038 | ~0.0009-0.02 (CV) | 0.1 |
| **Ratio (impl/theory)** | - | - | 8× to 26× |

### Gap Magnitude

**Simulations use λ = 0.1:**
- **8× larger** than theory at n=500 (0.1 / 0.0124 ≈ 8)
- **14× larger** than theory at n=1000 (0.1 / 0.0069 ≈ 14)
- **26× larger** than theory at n=2000 (0.1 / 0.0038 ≈ 26)

**Gap grows with n:**
- As n increases, theory requires smaller λ (more complexity allowed)
- Implementation keeps λ constant (same regularization)
- Gap widens: under-regularization relative to theory

**Implication for tree complexity:**
- Larger λ → stronger penalty on complexity → simpler trees
- Simulations fit **simpler trees** than theory requires
- Effect is **stronger at large n** (where gap is larger)

---

## Potential Impact on Under-Coverage

### Hypothesis: Over-Regularization Causes Under-Coverage

**Mechanism:**
1. λ = 0.1 is much larger than theory requires
2. Larger λ → stronger regularization → simpler trees (fewer leaves)
3. Simpler trees → less flexibility → higher bias, less variance captured
4. Influence function SE assumes "good enough" nuisance estimates
5. If trees are too simple (high misspecification), SE underestimates true uncertainty
6. Underestimated SE → narrow CIs → under-coverage

**Why it worsens with n:**
- Gap between λ=0.1 and theory grows with n (8× at n=500 → 26× at n=2000)
- At n=2000, trees are **most over-regularized relative to what theory allows**
- More misspecification → more SE underestimation → worse coverage

**Why specific to complex DGP:**
- Complex DGP requires more flexible trees to approximate well
- With λ=0.1, all DGPs get same regularization
- Simple DGP: even simple trees approximate well → SE ok
- Complex DGP: simple trees misspecified → SE underestimated

**Empirical pattern matches:**
- Under-coverage specific to complex DGP ✓
- Under-coverage worsens with n ✓
- All approaches affected (all use same λ) ✓
- Small bias but wrong SE (misspecification in variance, not mean) ✓

---

## Recommendations

### 1. Test Sensitivity to Lambda

**High priority:** Re-run simulations with theory-consistent λ

**Suggested values:**
```r
# Theory-driven
lambda_theory <- function(n) log(n) / n

# For comparison
lambda_values <- c(
  0.1,              # Current (baseline)
  0.05,             # Half current
  lambda_theory(n)  # Theory-consistent
)
```

**Compare:**
- Tree complexity (number of leaves)
- Prediction accuracy (held-out loss)
- Coverage rates by λ

**Expected:**
- λ = log(n)/n → more complex trees → better coverage (if hypothesis correct)
- Or: coverage still bad → issue is elsewhere (not λ)

### 2. Use CV-Selected Lambda

**Medium priority:** Enable CV in simulations

**Implementation:**
```r
estimate_att(...,
  cv_regularization = TRUE,  # Enable CV
  cv_K = 5                   # 5-fold CV for lambda
)
```

**Advantages:**
- Data-driven selection
- Matches paper's stated method
- Theory-consistent (CV grid centered on log(n)/n)

**Disadvantages:**
- Slower (adds CV overhead)
- Less reproducible (more variance across runs)

### 3. Update Paper Text

**Required:** Correct discrepancy between text and code

**Option A:** Change code to match paper
- Implement CV selection in simulations
- Update estimators.R to use `cv_regularization = TRUE`

**Option B:** Change paper to match code
- Update line 183: "λ = 0.1 (fixed)" instead of "chosen by CV"
- Add sensitivity analysis showing robustness to λ choice
- Or acknowledge limitation: "fixed λ may not be optimal"

### 4. Theoretical Investigation

**Low priority:** Prove sensitivity bounds

**Question:** How much does coverage degrade if λ is κ× too large?

**Approach:**
- Derive bound on SE bias as function of misspecification
- Express misspecification as function of λ deviation from optimal
- Predict coverage as function of λ

**Value:**
- Formal understanding of λ sensitivity
- Guidance on acceptable λ range

---

## Next Steps

### Immediate (Step 1)

**Re-run subset of simulations with theory-consistent λ:**
- Complex DGP only (where under-coverage occurs)
- n=2000 only (where gap is largest)
- Approaches 1-3 (fast)
- 100 reps (enough to detect coverage difference)

**Lambda values to test:**
```r
lambda_0.1   <- 0.1              # Current baseline
lambda_theory <- log(2000)/2000  # ≈ 0.0038
lambda_mid   <- 0.05             # Intermediate
```

**Compare:**
- Coverage rates
- Tree complexity
- Prediction accuracy

**Estimated time:** 1-2 hours (depending on cluster)

### If Hypothesis Confirmed

**Lambda is the culprit:**
1. Update simulation code to use theory-consistent λ
2. Re-run full simulations (all DGPs, all n)
3. Update paper with corrected results
4. Add section on λ sensitivity

### If Hypothesis Rejected

**Lambda is not the issue:**
1. Document that λ=0.1 vs theory λ doesn't explain under-coverage
2. Continue investigation with other hypotheses:
   - EIF assumption violations
   - Propensity score extremes
   - Tree complexity independent of λ
3. Update paper to correct λ selection statement regardless

---

## Files for Reference

### Theory
- `doubletree/inst/paper/manuscript.tex` (lines 94, 96, 115, 124, 156, 183, 300, 422, 492)

### Implementation
- `doubletree/simulations/six_approach_comparison/code/estimators.R` (all approaches)
- `doubletree/R/estimate_att.R` (cv_regularization parameter)
- `optimaltrees/R/cv_regularization.R` (CV grid and procedure)

### Results
- `results/combined/all_results.rds` (36,000 simulation results)
- `results/combined/summary_inference.csv` (coverage summary)

---

## Conclusion

**Definitive discrepancy identified:**
- Theory requires λ ≈ (log n)/n
- Paper claims CV selection
- Code uses fixed λ = 0.1 (8-26× too large)

**Plausible mechanism for under-coverage:**
- Over-regularization → tree misspecification → SE underestimation → narrow CIs → under-coverage
- Pattern matches: specific to complex DGP, worsens with n, affects all approaches

**Action required:**
1. Test λ sensitivity experimentally
2. Correct paper-code discrepancy
3. If confirmed, update simulations and results

**Status:** Investigation complete for tuning parameter selection. Hypothesis testable with targeted simulation runs.
