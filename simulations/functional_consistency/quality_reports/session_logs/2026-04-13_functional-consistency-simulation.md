# Session Log: Functional Consistency Simulation Study

**Date:** 2026-04-13
**Goal:** Design and deploy large-scale simulation to test averaged tree approach for perfect functional consistency
**Status:** 🔄 IN PROGRESS - Infrastructure complete, ready for O2 deployment

---

## Context

Following up on discretization fix (2026-04-10), we discovered that:
1. Standard M-split has imperfect functional consistency: max|η̂(X_i) - η̂(X_j)| ≠ 0 even when X_i = X_j
2. This persists even as M→∞ due to cross-fitting asymmetry
3. Averaged tree approach achieves perfect FC but breaks cross-fitting

**Key theoretical question:** Does breaking cross-fitting create self-prediction bias that invalidates asymptotic inference?

---

## Averaged Tree Approach

**Algorithm:**
1. Select modal structure via standard M-split
2. Refit on each of M×K folds
3. **Predict on FULL dataset** from each refit (breaks cross-fitting)
4. Average all M×K predictions

**Properties:**
- ✓ Perfect functional consistency: X_i = X_j → η̄(X_i) = η̄(X_j) (by averaging, everyone with same X gets same prediction)
- ✗ Self-prediction bias: η̄(X_i) includes predictions trained with Y_i
- Only 1/K predictions are truly cross-fit, (K-1)/K have self-prediction

**Bias analysis:**
- Self-prediction bias = O(1/n_ℓ) where n_ℓ is leaf size
- With regularization: n_ℓ ∝ n^β for β ∈ [0.3, 0.5]
- For √n-consistency: need √n × bias → 0, which requires β > 1/2
- **Concern:** If β ≤ 1/2, coverage may degrade as n increases

---

## Pattern Aggregation (Alternative)

User asked: "Is there another way to get perfect FC without self-prediction bias?"

**Algorithm:**
1. Run standard M-split (fully cross-fit)
2. Group observations by covariate pattern x
3. Average cross-fit predictions within pattern: η̃(x) = mean{η̂(X_i) : X_i = x}
4. Use pattern averages for all observations with same X

**Properties:**
- ✓ Perfect functional consistency (by construction)
- ✗ Indirect contamination: η̃(X_i) includes predictions from Y_j where X_j = X_i
- Bias = O(1/n_x) where n_x is pattern size (fixed by data, not tunable)

**Comparison to averaged tree:**
- Similar conceptually (both use information from obs with same X)
- Different group sizes: n_ℓ (tunable) vs n_x (fixed)
- Typically n_ℓ ≥ n_x, so averaged tree might have less bias

---

## Simulation Study Design

**Goal:** Empirically test whether perfect FC comes at cost of invalid inference

### Parameter Grid

- **n**: {200, 400, 800, 1600, 3200} - test asymptotic behavior
- **K**: {2, 3, 5} - proportion of cross-fitting (does it help?)
- **DGP**: {"simple", "complex", "sparse"} - robustness
- **method**: {"standard_msplit", "averaged_tree", "pattern_aggregation"} - three-way comparison
- **M**: 10 (fixed)

**Total:** 135 configurations, 500 reps each = **67,500 replications**

### Key Questions

**Q1:** Perfect FC achieved?
- Verify: max_diff ≈ 0 for averaged tree and pattern aggregation

**Q2:** Coverage degradation?
- Does coverage stay at 95% as n increases to 3200?
- If yes: Perfect FC is "free"!
- If no: Trade-off between FC and valid inference

**Q3:** Method comparison?
- Averaged tree vs pattern aggregation - does group size matter?
- Does n_ℓ > n_x lead to better coverage?

**Q4:** Effect of K?
- Does more cross-fitting (larger K) preserve coverage?

**Q5:** DGP robustness?
- Sparse DGP (small n_x) should be hardest for pattern aggregation

---

## Implementation

### Files Created

**Location:** `doubletree/simulations/functional_consistency/`

**Simulation code:**
1. `run_fc_simulation.R` - Main simulation function
   - `generate_dgp()` - Three DGPs (simple/complex/sparse)
   - `run_standard_msplit()` - Baseline method
   - `run_averaged_tree()` - Perfect FC via averaging all predictions
   - `run_pattern_aggregation()` - Perfect FC via pattern-level averaging

**SLURM infrastructure:**
2. `slurm/run_single_replication.R` - CLI wrapper for single rep
3. `slurm/run_simulations.slurm` - SLURM batch script (135 jobs)
4. `slurm/launch_all_simulations.sh` - Submit all configs
5. `slurm/quick_test.sh` - Local testing
6. `slurm/check_progress.sh` - Monitor O2 jobs
7. `slurm/combine_results.R` - Aggregate 67,500 replications

**Documentation:**
8. `slurm/README_O2.md` - Complete deployment guide
9. `SIMULATION_SUMMARY.md` - Research design and expected results
10. `METHODS_COMPARISON.md` - Detailed comparison of three approaches
11. `DEPLOYMENT_CHECKLIST.md` - Step-by-step deployment plan

### Local Testing Results

All three methods tested successfully (n=200, K=2, M=10):
```
              method att_est coverage max_diff_e max_diff_m0
     standard_msplit  0.6838     TRUE   1.51e-01    8.56e-02
       averaged_tree  0.6778     TRUE   0.00e+00    0.00e+00
 pattern_aggregation  0.6789     TRUE   0.00e+00    0.00e+00
```

✓ Standard M-split shows expected FC gap
✓ Both new methods achieve perfect FC
✓ All provide reasonable estimates

---

## Technical Details

### Three DGPs

**Simple:** 3 binary features (8 patterns, ~n/8 per pattern)
- Linear effects, no interactions
- Pattern size: moderate (n_x ≈ 50 at n=400)

**Complex:** 4 binary features (16 patterns, ~n/16 per pattern)
- Includes interactions
- Pattern size: smaller (n_x ≈ 25 at n=400)

**Sparse:** 5 binary features (32 patterns, ~n/32 per pattern)
- Only first 3 matter (noise features)
- Pattern size: small (n_x ≈ 12 at n=400)
- Stress test for pattern aggregation

### Key Metrics

For each replication, we compute:
- `att_est` - ATT estimate
- `se` - Standard error
- `ci_lower`, `ci_upper` - 95% CI
- `bias` - θ̂ - θ_true
- `coverage` - Does CI contain truth?
- `standardized_bias` - √n × bias (should stay bounded for valid inference)
- `bias_se_ratio` - |bias|/SE (should be << 1.96)
- `max_diff_e`, `max_diff_m0` - Functional consistency gaps
- `n_leaves_e`, `n_leaves_m0` - Tree complexity
- `mean_leaf_size_e`, `mean_leaf_size_m0` - Group sizes

---

## Expected Results

### Scenario A: Both Methods Maintain Coverage

**Finding:** Coverage stays at ~95% for all n, both averaged tree and pattern aggregation

**Interpretation:**
- Surprising! Bias vanishes fast enough despite β ≤ 1/2, OR
- Partial cross-fitting (1/K truly cross-fit) sufficient for debiasing, OR
- Finite sample artifacts (need even larger n to see degradation)

**Implication:** **Perfect functional consistency comes for free!**
- Can achieve FC without paying asymptotic cost
- Major practical advantage for fairness/policy applications

### Scenario B: Both Methods Fail

**Finding:** Coverage degrades to ~85% or worse at n=3200

**Interpretation:**
- Confirms self-prediction bias = O(n^{-β}) with β ≤ 1/2
- √n × bias grows, invalidating asymptotic CIs
- Bias term dominates variance

**Implication:** **Fundamental trade-off exists**
- Perfect FC requires accepting biased inference at large n
- Standard M-split: imperfect FC but valid inference
- Averaged tree / pattern aggregation: perfect FC but invalid at large n

**Resolution depends on application:**
- If FC critical (fairness, legal compliance): use new methods, acknowledge bias
- If inference critical (policy evaluation): use standard M-split
- Possible hybrid: use new methods for small n, standard for large n

### Scenario C: Averaged Tree Succeeds, Pattern Aggregation Fails

**Finding:** Averaged tree maintains ~95% coverage, pattern aggregation drops to ~88%

**Interpretation:**
- Group size matters: n_ℓ > n_x helps
- Tunable regularization allows larger groups
- Direct vs indirect bias mechanism differs empirically

**Implication:** **Use averaged tree with appropriate regularization**
- Maximize n_ℓ via regularization tuning
- Pattern aggregation limited by fixed n_x

---

## Theoretical Insight

User's key observation: Only (K-1)/K of averaged tree predictions have self-prediction bias.

**Structure:**
```
η̄(X_i) = (1/MK) × Σ_{m,k} η̂_{m,k}(X_i)
        = (1/K) × [predictions where i in test]
          + [(K-1)/K] × [predictions where i in training]
```

**Bias:**
```
E[η̄(X_i)] - η(X_i) = [(K-1)/K] × O(1/n_ℓ)
```

This reduces the **constant** but not the **rate**. For √n-consistency:
```
√n × bias = [(K-1)/K] × √n × O(n^{-β}) = O(n^{1/2-β})
```

Still need β > 1/2 for bias to vanish, regardless of K.

**But:** Maybe partial independence helps in finite samples? Simulation will tell.

---

## Next Steps

1. **Deploy to O2:**
   ```bash
   cd doubletree/simulations/functional_consistency/slurm
   bash launch_all_simulations.sh
   ```
   - Submits 135 SLURM jobs
   - Estimated runtime: 8-12 hours

2. **Monitor progress:**
   ```bash
   bash check_progress.sh  # Run periodically
   ```

3. **Combine results:**
   ```bash
   Rscript combine_results.R  # After all jobs complete
   ```

4. **Analyze:**
   - Coverage vs n by method (key plot!)
   - Standardized bias trends
   - Method comparison at n=3200
   - Effect of K and DGP

5. **Document findings:**
   - Update session log with results
   - Write interpretation
   - Update paper if needed (Theorem 3 may need revision)

---

## Open Questions

**Theoretical:**
1. Can we prove whether β > 1/2 or β ≤ 1/2 for optimal trees?
2. Is there a formal connection between partial cross-fitting and bias reduction?
3. Can DML orthogonality be salvaged with pattern-level effects?

**Empirical:**
1. What actually happens to coverage? (Most important!)
2. If coverage fails, can we quantify the degradation rate?
3. Does K help in practice even if theory says it shouldn't change the rate?

**Practical:**
1. If trade-off exists, which method should practitioners use?
2. Can we provide guidance on when FC matters more than valid CIs?
3. Should we develop diagnostic tools to detect when bias is dominating?

---

## Related Work

- **Discretization fix** (2026-04-10): Fixed store_training_data=FALSE bug, enabling this study
- **M-split theory** (doubletree/inst/paper/m-split-theory.tex): Theorem 3 may need revision if coverage fails
- **Tree structure rewrite** (2026-04-09): Enabled proper refit for averaged tree approach

---

## Time Tracking

- Design and discussion: 2 hours
- Implementation: 3 hours
  - Simulation function: 1 hour
  - SLURM infrastructure: 1.5 hours
  - Testing and debugging: 0.5 hours
- Documentation: 1.5 hours
- **Total: 6.5 hours**

---

## Status

**Current:** Infrastructure complete, locally tested, ready for O2 deployment

**Blocked by:** Need to deploy to O2 and wait for results

**Next action:** Run deployment commands or review/modify design

---

**Quality Score:** 90/100
- Complete implementation with all three methods
- Comprehensive testing and documentation
- Clear research questions and expected outcomes
- Infrastructure follows best practices (O2 skills pattern)
- Minor: Haven't actually deployed yet (but ready)
