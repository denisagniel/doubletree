# Diagnostic Suite for Complex DGP Calibration Failure

**Created:** 2026-05-27
**Purpose:** Systematic diagnosis of tree-based nuisance function estimation issues

## Problem

The six-approach comparison simulation shows severe inference problems in the complex DGP (dgp=3):

- **Bias:** All approaches underestimate ATT by 0.017-0.026 (true = 0.15)
- **Undercoverage:** 83-90% coverage (should be 95%)
- **Consistent across approaches:** Full-sample, crossfit, doubletree, msplit all fail similarly
- **Persists with n:** Even at n=2000, coverage only 88%

**Hypothesis:** The problem is in the **tree-based nuisance estimation**, not the inference approach.

---

## Diagnostic Scripts

### Phase 1: Nuisance Function Quality

#### `01_propensity_diagnostics.R`

**What it does:**
- Fits propensity score trees to simulated data
- Compares estimated e(X) to true propensity function
- Analyzes:
  - Bias, RMSE, calibration
  - Overlap quality and extreme weights
  - Tree complexity (number of leaves, depth)

**Key metrics:**
- RMSE: How well trees estimate e(X)?
- Calibration slope: Are predicted probabilities well-calibrated?
- Extreme weights: Overlap violations?

**Run time:** ~15-30 min (100 reps × 3 sample sizes × 4 DGPs)

**Output:**
- `diagnostics/results/propensity/propensity_diagnostics.rds`
- Plots: distribution, calibration, error, RMSE by DGP

---

#### `02_outcome_diagnostics.R`

**What it does:**
- Fits outcome model trees (on control units)
- Evaluates performance on:
  - Control units (in-sample)
  - Treated units (extrapolation)
- Compares data-fit trees to oracle trees (fit to true function)

**Key metrics:**
- RMSE (control vs treated): Extrapolation error?
- Oracle RMSE: Can trees represent function at all?
- Expressiveness gap: Oracle vs data performance

**Critical diagnostic:**
- If oracle RMSE > 0.10: Trees cannot represent function even with infinite data
- If data RMSE >> oracle RMSE: Estimation noise is the problem

**Run time:** ~15-30 min

**Output:**
- `diagnostics/results/outcome/outcome_diagnostics.rds`
- Plots: error by group, oracle vs data, RMSE by DGP

---

### Phase 5: EIF Decomposition

#### `05_eif_decomposition.R`

**What it does:**
- Fits full ATT estimation pipeline
- Decomposes bias into two components:
  - **Component 1:** Outcome model error on treated units
  - **Component 2:** Propensity-weighted residuals on control units
- Identifies which nuisance function causes most bias

**Key insight:**
- If Component 1 dominates: Outcome model is the problem
- If Component 2 dominates: Propensity score is the problem
- If both contribute: Both nuisance functions inadequate

**Run time:** ~20-40 min (fits both PS and outcome trees per replication)

**Output:**
- `diagnostics/results/eif_decomposition/eif_decomposition.rds`
- Plots: bias decomposition, component comparison

---

## Usage

### Quick Start (run all diagnostics)

```bash
cd /path/to/doubletree/simulations/six_approach_comparison
Rscript diagnostics/run_all_diagnostics.R
```

This will:
1. Run propensity diagnostics
2. Run outcome diagnostics
3. Run EIF decomposition
4. Generate summary report

**Total run time:** ~1-2 hours

---

### Run Individual Diagnostics

```bash
# Propensity score analysis
Rscript diagnostics/01_propensity_diagnostics.R

# Outcome model analysis
Rscript diagnostics/02_outcome_diagnostics.R

# EIF decomposition
Rscript diagnostics/05_eif_decomposition.R
```

---

### Interpret Results

After running diagnostics, check:

1. **Console output:** Summary statistics printed at end of each script
2. **RDS files:** Full results saved in `diagnostics/results/*/`
3. **Plots:** Visualizations in same directories

**Key questions to answer:**

| Question | Where to look | What to check |
|----------|---------------|---------------|
| Are propensity trees adequate? | `01_propensity_diagnostics.R` | RMSE, calibration slope |
| Are outcome trees adequate? | `02_outcome_diagnostics.R` | RMSE, oracle performance |
| Can trees represent functions? | `02_outcome_diagnostics.R` | Oracle RMSE (should be < 0.05) |
| Which nuisance function causes bias? | `05_eif_decomposition.R` | Component 1 vs Component 2 |
| Is extrapolation the problem? | `02_outcome_diagnostics.R` | RMSE (treated) vs RMSE (control) |

---

## Expected Findings

Based on the diagnostic plan, we expect to find:

**If trees are too simple:**
- High RMSE for both propensity and outcome models
- Large calibration slope deviation from 1.0
- Few leaves (< 6-8 for complex DGP)
- Oracle RMSE may be moderate

**If trees can't represent functions:**
- High oracle RMSE (> 0.10)
- Even with infinite data, trees struggle
- Suggests need for deeper trees or ensemble methods

**If regularization is too strong:**
- Trees consistently select max lambda
- Very small trees (2-3 leaves for complex DGP)
- Oracle trees perform much better than data trees

**If sample size is too small:**
- Large variance in tree structures across reps
- Performance improves substantially with n
- CV is unstable (high fold-to-fold variance)

---

## Next Steps After Diagnostics

Based on diagnostic findings, proceed to calibration experiments:

1. **If propensity trees inadequate:** Test weaker regularization for PS trees
2. **If outcome trees inadequate:** Test weaker regularization for outcome trees
3. **If both inadequate:** Test weaker regularization globally
4. **If oracle trees fail:** Consider ensemble methods (random forests, gradient boosting)
5. **If extrapolation fails:** Investigate overlap, consider different loss functions

See: `diagnostics/06_calibration_experiments.R` (to be created)

---

## File Structure

```
diagnostics/
├── README.md                      # This file
├── run_all_diagnostics.R          # Master script
├── 01_propensity_diagnostics.R    # Phase 1.1
├── 02_outcome_diagnostics.R       # Phase 1.2
├── 05_eif_decomposition.R         # Phase 5
├── utils/
│   ├── tree_diagnostics.R         # Tree analysis helpers
│   ├── eif_components.R           # EIF decomposition helpers
│   └── plotting.R                 # Visualization utilities
└── results/
    ├── propensity/                # PS diagnostic results
    ├── outcome/                   # Outcome diagnostic results
    └── eif_decomposition/         # EIF decomposition results
```

---

## Configuration

All scripts use the same configuration:

- **Replications:** 100 (default)
- **Sample sizes:** 500, 1000, 2000
- **DGPs:** 1-4 (simple, moderate, complex, continuous)
- **Regularization:** Fixed λ = log(n)/n (no CV by default)
- **Random seed:** 20260527 (for reproducibility)

To modify:
- Edit configuration sections at top of each script
- Or use command-line arguments (if implemented)

---

## Computational Requirements

- **Memory:** ~4-8 GB per script
- **Time:** 1-2 hours total for all diagnostics
- **Storage:** ~100-200 MB for all results

Can be parallelized by:
- Running each diagnostic script separately
- Running different DGPs in parallel (modify scripts)
- Using cluster (see `DEPLOYMENT_INSTRUCTIONS.md` for SLURM setup)

---

## Troubleshooting

**Script fails with package error:**
```
Error: Package 'optimaltrees' not found
```

**Solution:**
```bash
cd /path/to/optimaltrees
R CMD INSTALL .

cd /path/to/doubletree
R CMD INSTALL .
```

**Script fails with path error:**
```
Error: cannot open file 'diagnostics/utils/...'
```

**Solution:** Ensure you're running from the correct directory:
```bash
cd doubletree/simulations/six_approach_comparison
```

**Results look strange:**

1. Check random seed is set correctly
2. Verify DGP functions haven't changed
3. Compare to baseline results in `../results/combined/`

---

## References

- **Diagnostic plan:** See root-level plan document
- **DGP definitions:** `code/dgps.R`
- **Original simulation:** `code/run_single_replication.R`
- **Existing results:** `results/combined/all_results.rds`

---

## Contact

For questions or issues, document in session notes or quality reports.
