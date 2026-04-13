# Functional Consistency Simulation Study

**Created:** 2026-04-13
**Purpose:** Test whether averaged tree approach achieves perfect functional consistency while maintaining valid inference

---

## Research Questions

### Primary Questions

**Q1:** Does averaged tree achieve perfect functional consistency?
- **Metric:** max|η̄(X_i) - η̄(X_j)| for X_i = X_j
- **Prediction:** ≈ 0 (machine precision) for all n, K

**Q2:** Does coverage degrade as n increases?
- **Theory:** If self-prediction bias = O(1/n_ℓ) and n_ℓ ∝ n^β with β ≤ 1/2, then √n × bias does not vanish
- **Implication:** Asymptotic CIs may be invalid
- **Test:** Coverage should degrade from 95% as n increases

**Q3:** Does K (cross-fitting proportion) matter?
- **Mechanism:** K=2 means 1/2 cross-fit, K=5 means 1/5 cross-fit
- **Theory:** More cross-fitting (larger K) should improve coverage
- **Test:** Coverage at n=3200 should be: K=5 > K=3 > K=2

**Q4:** Does standardized bias grow with n?
- **Metric:** √n × (θ̂ - θ_true)
- **Standard M-split:** Should stay constant (asymptotically unbiased)
- **Averaged tree:** Should grow (biased)

### Secondary Questions

**Q5:** Robustness across DGPs?
- Simple (3 binary covariates, linear)
- Complex (4 covariates, interactions)
- Sparse (5 covariates, many patterns)

**Q6:** How does leaf size scale?
- **Metric:** mean_leaf_size vs n
- **Theory:** Should scale as n^β for some β ∈ [0.3, 0.5]

---

## Design

### Parameter Grid

| Parameter | Values | Purpose |
|-----------|--------|---------|
| **n** | {200, 400, 800, 1600, 3200} | Test asymptotic behavior |
| **K** | {2, 3, 5} | Test cross-fitting proportion |
| **DGP** | {"simple", "complex", "sparse"} | Robustness |
| **method** | {"standard_msplit", "averaged_tree"} | Head-to-head |
| **M** | 10 (fixed) | Number of splits |

**Total configurations:** 5 × 3 × 3 × 2 = 90
**Replications per config:** 500
**Total replications:** 45,000

### DGP Specifications

**Simple:**
- 3 binary covariates (8 patterns, ~n/8 per pattern)
- Linear treatment and outcome models
- Strong HTE: τ = 0.3 + 0.4x₁ + 0.3x₂

**Complex:**
- 4 binary covariates (16 patterns, ~n/16 per pattern)
- Interactions in treatment and outcome
- τ = 0.3 + 0.4x₁ + 0.3x₂ + 0.2x₁x₃

**Sparse:**
- 5 binary covariates (32 patterns, ~n/32 per pattern)
- Only first 3 covariates matter (others are noise)
- Tests behavior with small n_x

### Estimation Methods

**Standard M-split:**
1. Run M splits with K-fold cross-fitting
2. Select modal tree structure
3. Refit structure on each fold
4. Average predictions within each split
5. Use cross-fit predictions only (observation i's prediction never uses Y_i)

**Averaged tree:**
1. Run M splits with K-fold cross-fitting
2. Select modal tree structure (same as standard M-split)
3. Refit structure on each fold
4. **Predict on FULL dataset from each refit** (breaks cross-fitting)
5. Average ALL M×K predictions
   - 1/K are truly cross-fit (i in test)
   - (K-1)/K have self-prediction bias (i in training)

### Key Difference

**Standard M-split:**
```
η̂(X_i) = average over {predictions where Y_i not in training}
```
**Averaged tree:**
```
η̄(X_i) = (1/MK) × Σ_{m,k} η̂_{m,k}(X_i)
        = (1/K) × [cross-fit] + [(K-1)/K] × [self-predictions]
```

---

## Expected Results

### Hypothesis 1: Perfect Functional Consistency

**Averaged tree:**
- max_diff ≈ 0 for all n, K, DGP
- Differences at machine precision (< 1e-10)

**Standard M-split:**
- max_diff ≈ O(1/√n_x) where n_x is pattern size
- At n=200: max_diff ≈ 0.05-0.20 (confirmed in pilot)
- At n=3200: max_diff ≈ 0.01-0.05 (should decrease)

### Hypothesis 2: Coverage Degradation

**If β ≤ 1/2 (typical for trees):**

| n | Standard M-split | Averaged tree K=2 | Averaged tree K=5 |
|---|------------------|-------------------|-------------------|
| 200 | 95% | 95% | 95% |
| 400 | 95% | 94% | 95% |
| 800 | 95% | 92% | 94% |
| 1600 | 95% | 88% | 91% |
| 3200 | 95% | 82% | 88% |

**If β > 1/2 (optimistic):**
- Coverage stays at 95% for all n (bias vanishes fast enough)

### Hypothesis 3: Standardized Bias

**Standard M-split:**
```
√n × bias ≈ constant (or → 0)
```

**Averaged tree:**
```
√n × bias = √n × O(n^{-β}) = O(n^{1/2-β})
```

If β = 0.4:
- n=200: √200 × 0.02 ≈ 0.28
- n=3200: √3200 × 0.05 ≈ 2.83 (10× larger!)

### Hypothesis 4: Bias/SE Ratio

**For valid inference, need: |bias|/SE << 1.96**

**Averaged tree:**
- bias = O(n^{-β})
- SE = O(n^{-1/2})
- ratio = O(n^{1/2-β})

If β ≤ 1/2: ratio grows with n → CIs become invalid

---

## Local Test Results

**Pilot test (n=200, K=2, M=10, 3 reps):**
```
Standard M-split:
  Mean ATT: 0.681
  Coverage: 100%
  Max FC diff: 0.21 (propensity), 0.05 (outcome)
```

**Interpretation:**
- Simulation works correctly
- Standard M-split shows expected imperfect FC
- Ready for full O2 deployment

---

## Deployment

**Quick test:**
```bash
cd doubletree/simulations/functional_consistency/slurm
bash quick_test.sh
```

**Full deployment (O2):**
```bash
# On O2
cd ~/global-scholars/doubletree/simulations/functional_consistency/slurm
bash launch_all_simulations.sh

# Monitor
bash check_progress.sh

# After completion (~6-8 hours)
Rscript combine_results.R
```

See `slurm/README_O2.md` for complete instructions.

---

## Analysis Plan

After results are collected, analyze:

### 1. Functional Consistency by Method
```r
# Should show: averaged tree ≈ 0, standard M-split > 0
plot(n, max_diff_e, color = method)
```

### 2. Coverage vs n by Method
```r
# Key test: does averaged tree coverage degrade?
coverage_table <- aggregate(coverage ~ method + n, mean)
```

### 3. Coverage vs K for Averaged Tree
```r
# Does more cross-fitting help?
coverage_by_K <- subset(coverage_table, method == "averaged_tree")
```

### 4. Standardized Bias vs n
```r
# Should grow for averaged tree, stay constant for standard
plot(n, sqrt(n) * bias, color = method)
```

### 5. Bias/SE Ratio vs n
```r
# Should grow for averaged tree if β ≤ 1/2
plot(n, abs(bias)/se, color = method)
abline(h = 1.96, lty = 2)  # Threshold for valid inference
```

### 6. Leaf Size Scaling
```r
# Verify n_ℓ ∝ n^β
plot(log(n), log(mean_leaf_size))
# Slope = β (should be 0.3-0.5)
```

---

## Implications

### If Coverage Stays Valid (H2 rejected):

**Good news!** Averaged tree provides:
- Perfect functional consistency
- Valid asymptotic inference
- No cost for breaking cross-fitting

**Possible explanations:**
1. β > 1/2 (leaves grow faster than expected)
2. Partial cross-fitting (1/K) sufficient for debiasing
3. Small sample artifacts masking asymptotic behavior

### If Coverage Degrades (H2 confirmed):

**Trade-off identified:**
- Averaged tree: perfect FC, invalid large-sample CIs
- Standard M-split: imperfect FC, valid CIs

**Resolution depends on application:**
- If FC critical (e.g., fairness): use averaged tree, acknowledge bias
- If inference critical: use standard M-split, acknowledge FC gap
- Possible hybrid: averaged tree for small n, standard for large n

### If K Matters (H3 confirmed):

**Practical recommendation:**
- Use K=5 or K=10 for averaged tree (more cross-fitting)
- Reduces self-prediction bias from (K-1)/K to ~0.8-0.9
- May preserve valid inference even if β ≤ 1/2

---

## Files

**Simulation code:**
- `run_fc_simulation.R` - Main simulation function with DGPs and methods
- `slurm/run_single_replication.R` - CLI wrapper for single rep
- `slurm/run_simulations.slurm` - SLURM batch script
- `slurm/launch_all_simulations.sh` - Submit all configs
- `slurm/quick_test.sh` - Local test (3 reps)
- `slurm/check_progress.sh` - Monitor O2 jobs
- `slurm/combine_results.R` - Aggregate results
- `slurm/README_O2.md` - Deployment documentation

**Output:**
- `results/*.rds` - Individual replications
- `results/combined_fc_simulations.rds` - Aggregated results
- `results/combined_fc_simulations.csv` - CSV format

---

## Version History

- **2026-04-13:** Initial simulation infrastructure created
- **2026-04-13:** Pilot test successful (n=200, 3 reps)
- **TBD:** Full O2 deployment
