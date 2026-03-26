# Comprehensive Simulation Study Design Plan

**Goal:** Demonstrate tree-based DML-ATT applicability across diverse data characteristics beyond binary features.

---

## Current Limitations (DGP1-3)

**What we have:**
- ✅ Binary features only (4-5 features)
- ✅ Binary outcomes only
- ✅ Constant treatment effects
- ✅ Three sample sizes (400, 800, 1600)
- ✅ Four methods (tree, rashomon, forest, linear)

**What's missing:**
- ❌ Continuous features
- ❌ Continuous outcomes
- ❌ Mixed feature types
- ❌ Higher dimensions (10-20 features)
- ❌ Heterogeneous treatment effects
- ❌ Weak overlap scenarios
- ❌ Different noise levels

---

## Design Dimensions

### 1. Feature Types (Priority: HIGH)
- **Binary:** 4 features (current DGP1)
- **Continuous:** 4 continuous features
- **Mixed:** 2 binary + 2 continuous features
- **Higher-dim continuous:** 10 continuous features

**Why:** Trees handle continuous features differently (discretization). Need to show method works beyond binary case.

### 2. Outcome Types (Priority: HIGH)
- **Binary:** Current DGP1-3 (log-loss)
- **Continuous:** Gaussian outcomes (squared error loss)

**Why:** Continuous outcomes are common in practice (cost, test scores, etc.). Need to validate squared-error loss implementation.

### 3. Sample Sizes (Priority: MEDIUM)
- **Keep current:** n = 400, 800, 1600
- **Rationale:** Already covers small-to-medium practical range

### 4. Treatment Effect Heterogeneity (Priority: MEDIUM)
- **Constant ATT:** τ = 0.10 for all (current)
- **Heterogeneous ATT:** τ(X) varies by covariate subgroup

**Why:** Real effects often vary by covariates. Tests whether method can still estimate average effect correctly.

### 5. Overlap Quality (Priority: MEDIUM)
- **Good overlap:** e(X) ∈ [0.3, 0.7] (current)
- **Moderate overlap:** e(X) ∈ [0.1, 0.9]
- **Weak overlap:** Some e(X) near 0.05 or 0.95

**Why:** Tests propensity clipping robustness. Real data often has weak overlap.

### 6. Model Complexity (Priority: LOW)
- **Current:** Smooth functions, moderate interactions
- **Simple:** Linear/additive effects
- **Complex:** Strong interactions, non-monotonic

**Why:** Less critical - main focus is feature types and outcomes.

---

## Proposed DGP Suite (Prioritized)

### Phase 1: Feature Types & Outcomes (ESSENTIAL)

**DGP4: Continuous Features, Binary Outcome**
```r
generate_dgp_continuous_binary <- function(n, tau = 0.10, seed = NULL) {
  # 4 continuous features: X1, X2 ~ Unif(0,1); X3, X4 noise
  # Propensity: smooth function of X1, X2
  # Outcome: binary, depends on X1, X2
  # Trees will discretize continuous features automatically
}
```

**DGP5: Continuous Features, Continuous Outcome**
```r
generate_dgp_continuous_continuous <- function(n, tau = 0.10, seed = NULL) {
  # 4 continuous features
  # Propensity: smooth function
  # Outcome: continuous (Gaussian), depends on X1, X2
  # Tests squared-error loss for outcome model
}
```

**DGP6: Mixed Features, Binary Outcome**
```r
generate_dgp_mixed <- function(n, tau = 0.10, seed = NULL) {
  # 2 binary + 2 continuous features
  # Most realistic scenario
  # Tests handling of mixed feature types
}
```

### Phase 2: Higher Dimension (IMPORTANT)

**DGP7: High-Dimensional Continuous**
```r
generate_dgp_highdim_continuous <- function(n, tau = 0.10, seed = NULL) {
  # 10 continuous features: 4 signal, 6 noise
  # Tests whether trees can handle more features
  # Relevant for practical applications
}
```

### Phase 3: Robustness Checks (OPTIONAL)

**DGP8: Heterogeneous Treatment Effects**
```r
generate_dgp_heterogeneous_att <- function(n, tau_mean = 0.10, seed = NULL) {
  # Treatment effect varies: τ(X) = 0.05 + 0.10*X1
  # Average ATT still 0.10
  # Tests if method correctly estimates average despite heterogeneity
}
```

**DGP9: Weak Overlap**
```r
generate_dgp_weak_overlap <- function(n, tau = 0.10, seed = NULL) {
  # Propensity: e(X) ∈ [0.05, 0.95] with some near-extremes
  # Tests clipping robustness in harder settings
}
```

---

## Simulation Matrix

### Minimal (Phase 1 only): 108,000 replications
```
6 DGPs (DGP1-6) × 3 sample sizes × 4 methods × 1000 reps = 72,000
  Add DGP7 (high-dim): +12,000
  Add DGP8-9 (robustness): +24,000
```

### Recommended: Start with Phase 1 (DGP4-6)
```
3 new DGPs × 3 sample sizes × 4 methods × 1000 reps = 36,000 new
Total with current: 36,000 (current) + 36,000 (new) = 72,000
```

---

## Key Questions to Address

1. **Do trees work well with continuous features?**
   - Answer: Compare DGP4 (continuous) vs DGP1 (binary)
   - Metric: RMSE, coverage, convergence rates

2. **Does the method work for continuous outcomes?**
   - Answer: Analyze DGP5 (continuous Y)
   - Metric: Same as binary case

3. **What about mixed feature types (most realistic)?**
   - Answer: DGP6 results
   - Metric: Performance vs pure binary/continuous

4. **Can it handle higher dimensions?**
   - Answer: DGP7 (10 features) vs DGP4 (4 features)
   - Metric: Does performance degrade gracefully?

5. **How sensitive to overlap quality?**
   - Answer: DGP9 vs DGP1 (if including robustness phase)
   - Metric: Clipping frequency, bias, variance

---

## Implementation Strategy

### Step 1: Create New DGPs (DGP4-6)
- `dgps/dgps_continuous.R` - Continuous and mixed feature DGPs
- `dgps/dgps_continuous_outcome.R` - Continuous outcome support

### Step 2: Update Simulation Runner
- Extend `run_batch_replications.R` to handle:
  - `outcome_type = "continuous"` for DGP5
  - New DGP names (dgp4, dgp5, dgp6)

### Step 3: Test Locally
- Run 10 reps each of DGP4-6 to verify:
  - Continuous features discretize correctly
  - Continuous outcomes use squared-error loss
  - Mixed features work properly

### Step 4: Deploy to O2
- Add DGP4-6 to launch script
- Run 36,000 new replications (3 DGPs × 3 n × 4 methods × 1000)

### Step 5: Analysis & Decision
- Analyze Phase 1 results
- Decide whether Phase 2 (high-dim) needed
- Decide whether Phase 3 (robustness) needed

---

## Expected Findings

**Optimistic scenario:**
- Tree/rashomon work well for continuous features (auto-discretization works)
- Continuous outcomes work (squared-error loss validated)
- Mixed features combine best of both

**Realistic scenario:**
- Continuous features: small performance drop due to discretization approximation
- High-dim (DGP7): performance degrades, but forest/linear still work

**Worst case:**
- Continuous features fail badly → need better discretization strategy
- High-dim fails → dimensionality limit of tree methods

---

## Timeline Estimate

1. **DGP creation:** 2-3 hours (write + test locally)
2. **O2 deployment:** 1 hour (update scripts, launch)
3. **O2 execution:** ~12-18 hours (36k replications at 5-10 min/batch)
4. **Analysis:** 2-4 hours (combine results, generate tables/figures)

**Total: 1-2 days** for Phase 1

---

## Decision Points

**Before proceeding, decide:**

1. **Which phases?**
   - Phase 1 only (essential)?
   - Phase 1 + 2 (with high-dim)?
   - All phases (comprehensive)?

2. **How many replications per config?**
   - 1000 (current standard)?
   - 500 (faster, still adequate)?

3. **Same methods?**
   - Keep all 4 (tree, rashomon, forest, linear)?
   - Drop rashomon if it keeps failing?

4. **When to run?**
   - After current 36k finishes (wait for DGP1-3 results)?
   - Start now in parallel (separate configs)?

---

## Recommendation

**Start with Phase 1 (DGP4-6) immediately:**

**Why:**
1. **Essential validation** - Current simulations only show binary features work
2. **Manuscript needs it** - Can't claim "general method" without continuous features
3. **Manageable scope** - Only 36k new replications
4. **Fast turnaround** - 1-2 days to complete

**After Phase 1:**
- If continuous features work well → add Phase 2 (high-dim)
- If issues found → iterate on discretization/clipping
- Base final manuscript on Phases 1-2 combined

**Should we proceed with implementing Phase 1 (DGP4-6)?**
