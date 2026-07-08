# Diagnostic Infrastructure Implementation Summary

**Date:** 2026-05-27
**Status:** Phase 1 & 5 Complete (Core diagnostics implemented)

---

## What Was Implemented

### 1. Utility Functions (`utils/`)

#### `tree_diagnostics.R`
**Purpose:** Helper functions for analyzing tree structure and quality

**Functions:**
- `count_leaves()`: Count terminal nodes in tree
- `max_depth()`: Compute maximum tree depth
- `extract_tree_structure()`: Get split information as data frame
- `compute_prediction_metrics()`: Accuracy, calibration, log loss
- `analyze_overlap()`: Propensity score overlap diagnostics
- `compare_tree_structures()`: Cross-fold structure comparison
- `compute_oracle_performance()`: Fit tree to true function (no noise)

**Use case:** All diagnostic scripts use these for consistent tree analysis

---

#### `eif_components.R`
**Purpose:** EIF-based bias decomposition

**Functions:**
- `compute_true_propensity()`: True e(X) for each DGP
- `compute_true_outcome()`: True E[Y(0)|X] for each DGP
- `decompose_eif_components()`: Break ATT bias into outcome vs propensity contributions
- `compute_eif_values()`: Calculate influence function values
- `compute_eif_se()`: Standard error from EIF

**Use case:** Identifies which nuisance function (propensity or outcome) causes bias

---

#### `plotting.R`
**Purpose:** Visualization utilities

**Functions:**
- `plot_propensity_distribution()`: Histogram of e(X)
- `plot_calibration()`: Calibration curve
- `plot_prediction_error()`: Error vs true values
- `plot_tree_sizes()`: Tree complexity across folds
- `plot_bias_decomposition()`: Component contributions
- `plot_lambda_selection()`: CV curve
- `plot_diagnostics_grid()`: Multi-panel diagnostic plot

**Use case:** Consistent visualization across all diagnostics

---

### 2. Diagnostic Scripts

#### `01_propensity_diagnostics.R`
**Phase:** 1.1 - Propensity Score Quality

**What it analyzes:**
1. **Error metrics:** Bias, RMSE, MAE, max error vs true e(X)
2. **Calibration:** Does P(A=1|e(X)=p) ≈ p? (calibration slope should be ~1.0)
3. **Overlap:** Extreme weights, positivity violations, effective sample size
4. **Tree complexity:** Number of leaves, depth

**Simulation settings:**
- 100 reps × 3 sample sizes × 4 DGPs
- Fixed lambda = log(n)/n (no CV by default)
- Focus on complex DGP (dgp=3)

**Output:**
- `results/propensity/propensity_diagnostics.rds`
- Plots: distribution, calibration, RMSE by DGP, tree sizes

**Key diagnostic:**
- If RMSE > 0.10: Propensity trees inadequate
- If calibration slope ≠ 1.0: Poor probability estimates
- If extreme_weights > 10%: Overlap problems

---

#### `02_outcome_diagnostics.R`
**Phase:** 1.2 - Outcome Model Quality

**What it analyzes:**
1. **In-sample (control units):** RMSE, calibration, correlation
2. **Extrapolation (treated units):** RMSE, error increase vs controls
3. **Oracle performance:** Fit tree to TRUE outcome function (no noise)
4. **Expressiveness gap:** Oracle RMSE vs data RMSE

**Simulation settings:**
- Same as propensity diagnostics
- Fits outcome trees on control units only
- Evaluates on both control and treated units

**Output:**
- `results/outcome/outcome_diagnostics.rds`
- Plots: error by treatment group, oracle vs data, RMSE comparisons

**Key diagnostic:**
- **CRITICAL:** If oracle RMSE > 0.10, trees cannot represent function even with infinite data
  - This indicates fundamental expressiveness problem
  - Suggests need for deeper trees or ensemble methods
- If extrapolation_ratio > 1.5: Poor generalization to treated units
- If expressiveness_gap is large: Estimation noise (not expressiveness) is problem

---

#### `05_eif_decomposition.R`
**Phase:** 5 - Bias Decomposition

**What it analyzes:**

The EIF for ATT has two components:
1. **Component 1:** Outcome model error on treated units
   - How well does m̂0(X) predict E[Y(0)|X] for treated?
   - Directly affects ATT estimate

2. **Component 2:** Propensity-weighted residuals on control units
   - Correction term: weights control residuals by ê(X)/(1-ê(X))
   - Depends on both propensity and outcome model

**Identifies primary bias source:**
- If |Component 1| >> |Component 2|: Outcome model is problem
- If |Component 2| >> |Component 1|: Propensity score is problem
- If both large: Both nuisance functions inadequate

**Simulation settings:**
- Same 100 reps × 3 n × 4 DGPs
- Fits FULL pipeline: propensity + outcome + ATT estimation
- Compares to true nuisance functions

**Output:**
- `results/eif_decomposition/eif_decomposition.rds`
- Plots: bias by component, component comparison, scatterplots

**Key diagnostic:**
- Tells you WHERE to focus calibration efforts
- If Component 1 dominates: Fix outcome model first
- If Component 2 dominates: Fix propensity score first

---

### 3. Master Script

#### `run_all_diagnostics.R`
**Purpose:** Execute all diagnostics in sequence, generate summary

**What it does:**
1. Runs propensity diagnostics
2. Runs outcome diagnostics
3. Runs EIF decomposition
4. Combines results
5. Generates summary report with findings and recommendations

**Output:**
- `results/summary/diagnostic_summary.txt`
  - Status of each diagnostic
  - Key findings (RMSE, bias, tree sizes)
  - Automated interpretation
  - Recommended next steps

**Run time:** ~1-2 hours for all diagnostics

**Usage:**
```bash
cd doubletree/simulations/six_approach_comparison
Rscript diagnostics/run_all_diagnostics.R
```

---

## Implementation Completeness

### ✅ Implemented (Phases 1 & 5)

- [x] Utility functions (tree analysis, EIF decomposition, plotting)
- [x] Propensity score diagnostics
- [x] Outcome model diagnostics
- [x] Oracle tree analysis (expressiveness)
- [x] EIF bias decomposition
- [x] Master script with summary report
- [x] Documentation (README, this summary)

### ⏭️ Not Yet Implemented (Deferred)

**Phase 2: Cross-Validation Behavior**
- `02_cv_lambda_diagnostics.R`: Analyze CV selection of lambda
- Rationale: Using fixed lambda initially; add CV analysis if needed

**Phase 3: Rashomon Set Analysis**
- `03_rashomon_diagnostics.R`: Characterize Rashomon sets, intersection behavior
- Rationale: Lower priority; focus on core nuisance quality first

**Phase 4: Tree Structure Adequacy (Extended)**
- Oracle analysis is included in `02_outcome_diagnostics.R`
- Extended lambda grid experiments can be added later

**Phase 6: Calibration Experiments**
- `06_calibration_experiments.R`: Test interventions based on findings
- Rationale: Run after diagnostics to see what needs fixing

---

## File Structure

```
diagnostics/
├── README.md                       # User guide
├── IMPLEMENTATION_SUMMARY.md       # This file
├── run_all_diagnostics.R           # Master script
│
├── 01_propensity_diagnostics.R     # Phase 1.1 ✓
├── 02_outcome_diagnostics.R        # Phase 1.2 ✓
├── 05_eif_decomposition.R          # Phase 5 ✓
│
├── utils/
│   ├── tree_diagnostics.R          # ✓
│   ├── eif_components.R            # ✓
│   └── plotting.R                  # ✓
│
└── results/                        # Created on first run
    ├── propensity/
    ├── outcome/
    ├── eif_decomposition/
    └── summary/
```

---

## How to Use

### Quick Start

```bash
# From simulation root
cd doubletree/simulations/six_approach_comparison

# Run all diagnostics
Rscript diagnostics/run_all_diagnostics.R

# Check summary
cat diagnostics/results/summary/diagnostic_summary.txt

# View plots
open diagnostics/results/*/
```

### Run Individual Diagnostics

```bash
# Just propensity analysis
Rscript diagnostics/01_propensity_diagnostics.R

# Just outcome analysis
Rscript diagnostics/02_outcome_diagnostics.R

# Just EIF decomposition
Rscript diagnostics/05_eif_decomposition.R
```

### Customize Settings

Edit configuration section in each script:
```r
# Simulation parameters
n_reps <- 100              # Number of replications
sample_sizes <- c(500, 1000, 2000)  # Sample sizes to test
dgps <- 1:4                # Which DGPs (1=simple, 2=moderate, 3=complex, 4=continuous)

# Regularization
use_cv <- FALSE            # Use CV to select lambda?
fixed_lambda_multiplier <- 1.0  # Multiply theory lambda by this
```

---

## Expected Run Time

| Diagnostic | Reps | Time | Output Size |
|------------|------|------|-------------|
| Propensity | 100 × 3 × 4 | 15-30 min | ~20 MB |
| Outcome | 100 × 3 × 4 | 15-30 min | ~25 MB |
| EIF | 100 × 3 × 4 | 20-40 min | ~30 MB |
| **Total** | 1,200 sims | **1-2 hours** | **~75 MB** |

**Note:** Can be parallelized:
- Run each diagnostic on separate cores
- Or run different DGPs in parallel (modify scripts)
- Or run on cluster (see `DEPLOYMENT_INSTRUCTIONS.md`)

---

## Interpreting Results

### Step 1: Check Summary Report

```bash
cat diagnostics/results/summary/diagnostic_summary.txt
```

Look for:
- **Propensity RMSE:** Should be < 0.10 for adequate estimation
- **Outcome oracle RMSE:** Should be < 0.05 if trees can represent function
- **EIF component dominance:** Which component (1 or 2) is larger?

### Step 2: Review Plots

**Propensity plots:**
- `ps_distribution_example.png`: Are estimates close to truth?
- `ps_calibration_example.png`: Is calibration slope ≈ 1.0?
- `ps_rmse_by_dgp.png`: Which DGPs have worst estimation?

**Outcome plots:**
- `outcome_error_example.png`: Error pattern?
- `outcome_by_treatment_example.png`: Extrapolation to treated?
- `oracle_vs_data_rmse.png`: Is gap large? (estimation noise)

**EIF plots:**
- `bias_decomposition_by_dgp.png`: Which component dominates?
- `complex_bias_by_n.png`: Does bias decrease with n?

### Step 3: Form Hypothesis

Based on findings, determine:
1. **What's the root cause?**
   - Trees too simple? (High oracle RMSE)
   - Regularization too strong? (Very small trees)
   - Sample size too small? (High variance)

2. **Which nuisance function is worse?**
   - Propensity? (Component 2 dominates)
   - Outcome? (Component 1 dominates)
   - Both? (Both components large)

3. **What intervention to test?**
   - Weaker regularization? (λ × 0.5 or λ × 0.25)
   - Deeper trees? (reduce lambda substantially)
   - Ensemble methods? (if oracle RMSE high)
   - Separate lambda for PS vs outcome?

### Step 4: Run Calibration Experiments

Based on hypothesis, test interventions in Phase 6 calibration experiments.

---

## Troubleshooting

### Issue: Script fails with package error

```
Error: there is no package called 'optimaltrees'
```

**Solution:**
```bash
cd optimaltrees
R CMD INSTALL .

cd doubletree
R CMD INSTALL .
```

### Issue: Script fails with path error

```
Error: cannot open file 'diagnostics/utils/...'
```

**Solution:** Ensure you're in the correct directory:
```bash
cd doubletree/simulations/six_approach_comparison
pwd  # Should end in six_approach_comparison
```

### Issue: Out of memory

**Solution:** Reduce n_reps:
```r
n_reps <- 50  # Instead of 100
```

Or run one DGP at a time:
```r
dgps <- 3  # Just complex DGP
```

### Issue: Results don't match expectations

1. Check random seed is set correctly
2. Verify DGP functions haven't changed (compare to `code/dgps.R`)
3. Check regularization settings
4. Compare to baseline results in `results/combined/all_results.rds`

---

## Next Steps

### Immediate (After Running Diagnostics)

1. Review summary report
2. Examine plots to confirm findings
3. Document key findings in session notes

### Short-term (Based on Findings)

If diagnostics reveal:
- **High oracle RMSE:** Need deeper trees or ensemble methods
- **Large expressiveness gap:** Estimation noise, calibrate regularization
- **Component 1 dominates:** Focus on outcome model first
- **Component 2 dominates:** Focus on propensity score first

### Long-term (Calibration Phase)

1. Create `06_calibration_experiments.R` based on findings
2. Test interventions:
   - Weaker regularization grid: λ × [0.1, 0.25, 0.5, 1, 2]
   - CV for lambda selection
   - Separate lambda for PS vs outcome
   - Ensemble methods if needed
3. Validate solution on full simulation (500 reps)

---

## Contact / Questions

For questions about this implementation:
- Review code comments in each script
- Check README.md for usage guidance
- Document issues in session notes
- Refer to original diagnostic plan document

---

## Version History

- **2026-05-27:** Initial implementation (Phases 1 & 5)
  - Core diagnostic infrastructure
  - Propensity, outcome, and EIF decomposition
  - Master script and documentation
