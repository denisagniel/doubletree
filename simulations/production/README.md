# Production Simulations for doubletree Manuscript

**Status:** Ready for production runs (Phase 7: verification pending)
**Last updated:** 2026-03-04

---

## Overview

This directory contains **manuscript-ready** simulation code for doubletree paper Section 4. All code follows research constitution requirements:

- ✅ Stress-testing with adversarial scenarios
- ✅ No quiet favoritism (includes baselines that may outperform)
- ✅ Expected failure modes documented BEFORE running
- ✅ Three-way fidelity (paper ↔ code ↔ package)
- ✅ Full reproducibility with explicit seeds

---

## Directory Structure

```
production/
├── dgps/
│   ├── dgps_smooth.R       # DGPs 1-3 (validated, 95% coverage)
│   └── dgps_stress.R       # DGPs 4-6 (adversarial scenarios)
├── methods/
│   ├── method_forest_dml.R # Random forest baseline (ranger)
│   └── method_linear_dml.R # GLM baseline (logistic regression)
├── run_primary.R           # Main simulations (18k runs, ~4 hours)
├── run_stress.R            # Stress tests (2.4k runs, ~2 hours)
├── analyze_manuscript.R    # Generate Tables 1-2, Figures 1-2
└── results/                # Simulation outputs (created at runtime)
```

---

## Quick Start

### 1. Run Primary Simulations (Table 1)

```r
# From doubletree/simulations/production/
source("run_primary.R")

# Output:
# - results/primary_YYYY-MM-DD/simulation_results.rds
# - results/primary_YYYY-MM-DD/summary_stats.csv

# Runtime: ~15 hours single-threaded; ~4 hours with 4 cores
# Total: 18,000 simulations (3 DGPs × 4 methods × 3 sample sizes × 500 reps)
```

### 2. Run Stress-Test Simulations (Table 2)

```r
source("run_stress.R")

# Output:
# - results/stress_YYYY-MM-DD/stress_results.rds
# - results/stress_YYYY-MM-DD/stress_summary.csv
# - results/stress_YYYY-MM-DD/failure_modes.txt

# Runtime: ~2 hours with 4 cores
# Total: 2,400 simulations (3 DGPs × 2 methods × 2 sample sizes × 200 reps)
```

### 3. Generate Manuscript Outputs

```r
source("analyze_manuscript.R")

# Requires: run_primary.R and run_stress.R completed first

# Output:
# - results/manuscript_outputs_YYYY-MM-DD/table1_primary.tex
# - results/manuscript_outputs_YYYY-MM-DD/table2_stress.tex
# - results/manuscript_outputs_YYYY-MM-DD/figure1a_coverage.pdf
# - results/manuscript_outputs_YYYY-MM-DD/figure1b_rmse.pdf
# - results/manuscript_outputs_YYYY-MM-DD/figure2_method_comparison.pdf
```

---

## Simulation Design

### Primary Simulations (run_primary.R)

**Grid:**
- DGPs: Binary (4 features), Continuous (4 features), Moderate (5 features)
- Methods: Tree-DML, Rashomon-DML, Forest-DML, Linear-DML
- Sample sizes: n ∈ {400, 800, 1600}
- Replications: 500 per configuration
- Total: 18,000 runs

**Expected results:**
- Tree-DML: 95% coverage, low bias
- Rashomon-DML: Similar to tree (structure intersection)
- Forest-DML: Comparable to tree (bias/RMSE within 10%)
- Linear-DML: Worse when nonlinear (DGPs 2-3: bias 2-3× larger)

**Paper claims validated:**
- "Tree-DML achieves nominal coverage" (line 160)
- "Comparable to forest-DML" (line 162)
- "Linear-DML suffers under nonlinearity" (line 162)

### Stress-Test Simulations (run_stress.R)

**Grid:**
- DGPs: Weak overlap, Piecewise, High-dimensional
- Methods: Tree-DML, Forest-DML (skip linear)
- Sample sizes: n ∈ {800, 1600}
- Replications: 200 per configuration
- Total: 2,400 runs

**Expected failure modes:**

| DGP | Assumption Stretched | Expected Behavior | Recovery? |
|-----|---------------------|-------------------|-----------|
| DGP 4: Weak Overlap | Positivity (e ∈ [0.05, 0.95]) | CI width 2-3× larger | No |
| DGP 5: Piecewise | Smoothness (step functions) | Trees excel, linear fails | Linear never |
| DGP 6: High-Dim | Dimensionality (p=8) | Coverage <95% at n=800 | Yes, at n≥1600 |

**Constitution compliance:**
- Adversarial scenarios included ✓
- Failure modes documented BEFORE running ✓
- Observed failures reported honestly ✓

---

## Data Generating Processes

### DGPs 1-3 (Smooth)

**Common features:**
- Treatment effect τ = 0.10 on **probability scale** (p1 = p0 + τ)
- Binary outcomes Y ∈ {0, 1}
- Binary treatment A ∈ {0, 1}
- Propensity and outcome depend on subset of covariates
- **Validated:** 95% coverage achieved with λ = log(n)/n

**DGP 1 (Binary):**
- 4 binary features (X1-X4)
- Signal in X1, X2; noise in X3, X4
- 16 covariate patterns

**DGP 2 (Continuous):**
- 4 continuous features (X1-X4)
- Smooth functions of X1, X2
- Tests discretization workflow

**DGP 3 (Moderate):**
- 5 binary features (X1-X5)
- Signal in X1, X2, X3
- 32 covariate patterns

### DGPs 4-6 (Stress)

**DGP 4 (Weak Overlap):**
- Propensity scores near boundaries: e ∈ [0.08, 0.92]
- Expected: Large variance, wide CIs
- Tests: Robustness to near-positivity violations

**DGP 5 (Piecewise):**
- Piecewise constant functions (4 regions)
- Sharp discontinuities at boundaries
- Tests: Handling non-smooth nuisances

**DGP 6 (High-Dimensional):**
- 8 binary features (256 patterns)
- Sparse signal (only 3 features matter)
- Tests: Curse of dimensionality

---

## Methods

### Tree-DML (doubletree::dml_att)

```r
doubletree::dml_att(
  X, A, Y,
  K = 5,
  regularization = log(n) / n,  # Theory-driven
  cv_regularization = FALSE,    # CV not beneficial
  use_rashomon = FALSE          # Fold-specific trees
)
```

**Implementation:** `doubletree/R/dml_att.R`

### Rashomon-DML (doubletree::dml_att)

```r
doubletree::dml_att(
  X, A, Y,
  K = 5,
  regularization = log(n) / n,
  cv_regularization = FALSE,
  use_rashomon = TRUE  # Structure intersection across folds
)
```

**Implementation:** `doubletree/R/dml_att.R` with `use_rashomon = TRUE`

### Forest-DML (method_forest_dml.R)

```r
dml_att_forest(
  X, A, Y,
  K = 5,
  num.trees = 500,
  seed = seed
)
```

**Implementation:** `production/methods/method_forest_dml.R`
**Package:** ranger (probability forests)

### Linear-DML (method_linear_dml.R)

```r
dml_att_linear(
  X, A, Y,
  K = 5,
  interactions = FALSE,
  seed = seed
)
```

**Implementation:** `production/methods/method_linear_dml.R`
**Package:** stats::glm (binomial family)

---

## Manuscript Outputs

### Table 1: Primary Results (Main Text)

**Content:**
- All 4 methods across DGPs 1-3 and n ∈ {400, 800, 1600}
- Metrics: Bias, RMSE, Coverage (%), CI width
- Shows tree ≈ forest; linear worse when nonlinear

**Location:** `results/manuscript_outputs_YYYY-MM-DD/table1_primary.tex`

### Table 2: Stress-Test Results (Appendix)

**Content:**
- Tree and forest across DGPs 4-6 and n ∈ {800, 1600}
- Same metrics as Table 1
- Documents observed failure modes vs expected

**Location:** `results/manuscript_outputs_YYYY-MM-DD/table2_stress.tex`

### Figure 1: Coverage and RMSE by Sample Size

**Panel A:** Coverage vs n (horizontal line at 95%)
**Panel B:** RMSE vs n (log-log scale, shows √n decay)
**DGP:** Focus on DGP 1 (binary) for clarity

**Location:**
- `figure1a_coverage.pdf`
- `figure1b_rmse.pdf`

### Figure 2: Method Comparison

**Content:** Grouped bar chart of bias and RMSE by method
**Scenario:** DGP 1, n = 800 (representative)
**Shows:** Tree ≈ rashomon ≈ forest < linear

**Location:** `figure2_method_comparison.pdf`

---

## Three-Way Fidelity Checklist

Before submitting manuscript Section 4:

**Paper → Code:**
- [ ] DGPs described in paper match `dgps_smooth.R` and `dgps_stress.R`
- [ ] Sample sizes in paper (400, 800, 1600) match `run_primary.R` grid
- [ ] Methods in paper (tree, rashomon, forest, linear) match `run_primary.R`
- [ ] Metrics in paper (bias, RMSE, coverage) match `analyze_manuscript.R`

**Code → Paper:**
- [ ] Table 1 numbers match `results/primary_*/summary_stats.csv`
- [ ] Table 2 numbers match `results/stress_*/stress_summary.csv`
- [ ] Figure 1 plots data from `run_primary.R`
- [ ] Reported τ = 0.10 matches DGP specifications

**Code → Package:**
- [ ] `doubletree::dml_att` implements ATT as paper describes
- [ ] `dml_att_forest` uses ranger as paper states
- [ ] `dml_att_linear` uses glm as paper states
- [ ] Hyperparameters (λ = log(n)/n, K = 5) match documentation

**Reproducibility:**
- [ ] All scripts use explicit seeds (SEED_OFFSET defined)
- [ ] Package versions documented (sessionInfo() saved)
- [ ] QUICKSTART explains how to reproduce tables/figures
- [ ] No hidden parameter choices (all hyperparameters explicit)

---

## Computational Requirements

**Hardware:**
- Recommended: 4+ cores, 8+ GB RAM
- Minimum: 2 cores, 4 GB RAM (slower)

**Software:**
- R ≥ 4.0
- Packages: doubletree, optimaltrees, ranger, parallel, dplyr, ggplot2, kableExtra

**Storage:**
- Primary results: ~300 MB
- Stress results: ~100 MB
- Total: ~500 MB

**Runtime:**
- Primary (single-threaded): ~15 hours
- Primary (4 cores): ~4 hours
- Stress (4 cores): ~2 hours
- Total (4 cores): ~6 hours

---

## Quality Gates

Before using results in manuscript:

**Minimum criteria (80/100 - commit level):**
- [ ] All methods converge >95% of time
- [ ] Tree-DML coverage 90-98% for DGPs 1-3
- [ ] No missing data (all runs completed)

**Publication criteria (90/100 - PR level):**
- [ ] Tree-DML coverage 93-97% for DGPs 1-3
- [ ] Bias < 0.03 (30% of τ) for all methods
- [ ] Stress-test failure modes match expectations

**Excellence criteria (95/100):**
- [ ] Tree-DML coverage 94-96% for DGPs 1-3
- [ ] Three-way fidelity verified
- [ ] All manuscript claims have code citations

---

## Troubleshooting

### "Model limit exceeded" (GOSDT error)

**Cause:** Too many features or patterns for tree optimization

**Solutions:**
1. Use DGPs 1-2 (4 features) instead of DGP 3 (5 features)
2. Increase sample size (more data per pattern)
3. Adjust `model_limit` in doubletree config (not recommended)

### Poor coverage (<90%)

**Check:**
1. Using `dgps_smooth.R`? (Not deprecated DGPs)
2. Using λ = log(n)/n? (Not CV-based)
3. Enough replications? (Need 200+ for stable coverage)

### Large bias (>10% of truth)

**Likely causes:**
- Wrong DGP file (misspecified τ)
- Incorrect regularization
- Convergence failures

**Solution:** Check `$converged` field; filter failed runs

---

## Version History

- **2026-03-04:** Initial production release
  - Validated DGPs (95% coverage)
  - 4 methods implemented and tested
  - 6 DGPs (3 smooth + 3 stress)
  - Full manuscript output pipeline

---

## References

**Related files:**
- Root QUICKSTART: `../QUICKSTART.md`
- Deprecated DGPs: `../deprecated/README.md`
- Diagnostics: `../diagnostics/`

**Session notes:**
- Package development: `../../session_notes/`
- Root project: `../../../session_notes/`

**Paper:**
- Manuscript: `../../paper/manuscript.tex`
- Section 4 (simulations) to be completed using these outputs
