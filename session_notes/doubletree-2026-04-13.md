# Session Notes: Functional Consistency Simulation Study

**Date:** 2026-04-13
**Package:** doubletree
**Focus:** Design and deploy large-scale simulation to test averaged tree approach

---

## Summary

Designed and implemented complete simulation infrastructure (67,500 replications) to test whether "averaged tree" and "pattern aggregation" approaches can achieve perfect functional consistency without sacrificing valid asymptotic inference.

---

## Key Developments

### 1. Three-Method Comparison

Implemented and tested three approaches:

**Standard M-split (baseline):**
- True cross-fitting, no self-prediction bias
- Imperfect functional consistency: max|η̂(X_i) - η̂(X_j)| ≠ 0 when X_i = X_j
- Valid asymptotic inference

**Averaged tree:**
- Averages M×K predictions (some trained with Y_i)
- Perfect functional consistency
- Self-prediction bias = O(1/n_ℓ) where n_ℓ is leaf size
- Asymptotic validity unclear (testing this!)

**Pattern aggregation:**
- Averages cross-fit predictions within covariate patterns
- Perfect functional consistency
- Indirect contamination bias = O(1/n_x) where n_x is pattern size
- Asymptotic validity unclear (testing this!)

### 2. Theoretical Analysis

**Key insight:** Self-prediction bias is O(1/n_ℓ) where n_ℓ ∝ n^β with typical β ∈ [0.3, 0.5]

For √n-consistency: need √n × bias → 0
- This requires β > 1/2
- Typical trees have β ≤ 1/2

**Implication:** Coverage may degrade as n increases if β ≤ 1/2

**User's observation:** Only (K-1)/K of averaged tree predictions have self-prediction bias
- Reduces constant but not rate
- Still need β > 1/2 for asymptotic validity

### 3. Simulation Design

**Parameter grid:**
- n: {200, 400, 800, 1600, 3200} - test asymptotic behavior
- K: {2, 3, 5} - proportion of cross-fitting
- DGP: {simple, complex, sparse} - robustness
- method: {standard_msplit, averaged_tree, pattern_aggregation}

**Total:** 135 configs × 500 reps = 67,500 replications

**Key questions:**
1. Perfect FC achieved? (should be 0 for both new methods)
2. Coverage degrades? (critical test!)
3. Averaged tree vs pattern aggregation - does group size matter?
4. Effect of K?
5. DGP robustness?

### 4. Implementation

**Location:** `doubletree/simulations/functional_consistency/`

**Created:**
- `run_fc_simulation.R` - Main simulation function with all 3 methods and 3 DGPs
- Complete SLURM infrastructure (8 files) for O2 deployment
- Comprehensive documentation (4 markdown files)

**Local testing:** All three methods work correctly
```
              method att_est coverage max_diff_e max_diff_m0
     standard_msplit  0.6838     TRUE   1.51e-01    8.56e-02
       averaged_tree  0.6778     TRUE   0.00e+00    0.00e+00
 pattern_aggregation  0.6789     TRUE   0.00e+00    0.00e+00
```

✓ Both new methods achieve perfect FC as expected

---

## Next Steps

### Deployment to O2

Code pushed to GitHub (commit c6659e4). To deploy:

```bash
# On O2
ssh username@o2.hms.harvard.edu
cd ~/global-scholars
git pull  # In optimaltrees and doubletree repos

# Install packages if needed
module load gcc/9.2.0 R/4.2.1
R CMD INSTALL optimaltrees
R CMD INSTALL doubletree

# Test
cd doubletree/simulations/functional_consistency/slurm
bash quick_test.sh

# Deploy (submit 135 jobs)
bash launch_all_simulations.sh
```

Expected runtime: 8-12 hours

### Analysis

After results complete:
1. Combine: `Rscript combine_results.R`
2. Key analyses:
   - Coverage vs n by method (CRITICAL)
   - Functional consistency verification
   - Standardized bias trends
   - Method comparison at n=3200

### Outcomes

**If both new methods maintain coverage:**
- Perfect FC comes "for free"!
- Major practical advantage

**If both fail:**
- Confirms fundamental trade-off
- Need to choose: perfect FC vs valid inference

**If one succeeds:**
- Use the successful method
- Group size or mechanism matters

---

## Files Modified/Created

**doubletree:**
- `simulations/functional_consistency/` (entire directory, 17 files)
- Committed and pushed to GitHub

**Root level:**
- `quality_reports/session_logs/2026-04-13_functional-consistency-simulation.md`
- `session_notes/doubletree-2026-04-13.md` (this file)

---

## Time

- Design and discussion: 2 hours
- Implementation: 3 hours
- Testing and documentation: 1.5 hours
- Session logging: 0.5 hours
- **Total: 7 hours**

---

## Links

- **Full session log:** `quality_reports/session_logs/2026-04-13_functional-consistency-simulation.md`
- **Simulation summary:** `doubletree/simulations/functional_consistency/SIMULATION_SUMMARY.md`
- **Methods comparison:** `doubletree/simulations/functional_consistency/METHODS_COMPARISON.md`
- **Deployment guide:** `doubletree/simulations/functional_consistency/slurm/README_O2.md`

---

## Context for Future

This simulation will definitively answer whether we can achieve perfect functional consistency without paying an asymptotic cost. Results will inform:
- Method recommendations in doubletree package
- Potential revision of Theorem 3 in paper
- Practical guidance for when FC matters most

The averaged tree approach was motivated by discovering that standard M-split has persistent FC gaps due to cross-fitting asymmetry. If it works (maintains valid inference), it's a significant methodological contribution.
