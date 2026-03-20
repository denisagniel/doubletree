# doubletree: Causal Estimation with Interpretable Trees

Implements causal inference for the **Average Treatment Effect on the Treated (ATT)** using efficient influence function-based estimation with cross-fitting and interpretable optimal decision trees. Doubly robust, semiparametric estimator with optional Rashomon set integration for interpretable, stable tree selection.

**Current version:** 0.0.0.9000 (development)
**Repository:** [github.com/denisagniel/doubletree](https://github.com/denisagniel/doubletree)
**Depends on:** [optimaltrees](https://github.com/denisagniel/treefarmr) v0.4.0+

## Features

- **DML-ATT estimation** with tree-based nuisance functions (propensity and outcome models)
- **Rashomon-DML integration** for interpretable tree selection via cross-validated structure intersection
- **Binary and continuous outcomes** (via `outcome_type`)
- **Parallel execution** on O2 cluster for large-scale simulations
- **Theory-aligned implementation** with comprehensive simulation infrastructure

## Recent Updates (March 2026)

- **O2/SLURM Infrastructure:** Complete setup for distributed simulations (18,000+ replications in 30-60 minutes)
- **S7 Integration:** Full compatibility with optimaltrees S7 class system
- **Simulation Grid:** 3 DGPs × 4 methods × 3 sample sizes × 500 replications
- **Methods:** tree-DML, Rashomon-DML, forest-DML (ranger), linear-DML

## Project Structure

```
doubletree/
├── .Rprofile              # Auto-loads optimaltrees from local dev
├── .gitignore             # Standard R/research project ignores
├── README.md              # This file
├── DESCRIPTION            # R package metadata
├── NAMESPACE              # R package namespace
├── R/                     # R source code
├── man/                   # Package documentation
├── tests/                 # Test files
├── paper/                 # Manuscript and related files
│   ├── manuscript.tex
│   ├── figures/
│   └── references.bib
├── simulations/           # Simulation scripts and results
│   ├── run_simulations.R
│   └── results/
└── data/                  # Data files
```

## Installation

### From GitHub

```r
# Install optimaltrees first
devtools::install_github("denisagniel/treefarmr")

# Install doubletree
devtools::install_github("denisagniel/doubletree")
```

### For Development

```bash
# Clone repositories
git clone git@github.com:denisagniel/treefarmr.git optimaltrees
git clone git@github.com:denisagniel/doubletree.git

# Install optimaltrees
cd optimaltrees
R CMD INSTALL .

# Load doubletree in development mode
cd ../doubletree
R
```

In R:
```r
devtools::load_all()  # Development mode
# OR
devtools::install()   # Install locally
library(doubletree)
```

### Dependencies

Required R packages:
- `optimaltrees` (>= 0.4.0)
- `dplyr`
- `ranger` (for forest-DML baseline)

For O2 cluster simulations, also install:
- `optparse` (command-line arguments)

## Methods

### Tree-DML (Standard)

Fits one optimal tree per fold for each nuisance function (propensity and outcome models):

```r
fit <- estimate_att(X, A, Y, K = 5, use_rashomon = FALSE)
```

### Rashomon-DML (Interpretable)

Selects a **single interpretable tree per nuisance** via the intersection of Rashomon sets across cross-fitting folds, then refits that structure per fold for valid cross-fitted estimation:

```r
fit <- estimate_att(
  X, A, Y,
  K = 5,
  use_rashomon = TRUE,
  rashomon_bound_multiplier = 0.05
)
```

**How it works:**
1. Fit Rashomon sets (near-optimal trees) for each nuisance in each fold
2. Find structural intersection across folds (stable tree structures)
3. Refit the intersecting structure per fold for valid cross-fitting
4. Use fold-specific predictions for DML estimation

**Benefits:** Interpretable trees (single structure for each nuisance) with valid statistical inference.

### Baseline Methods

For comparison:
- **forest-DML:** Random forests via `ranger` (see `simulations/production/methods/method_forest.R`)
- **linear-DML:** Logistic regression (see `simulations/production/methods/method_linear.R`)

## Minimal example

**Binary outcome (default):**

```r
devtools::load_all()
# X: data.frame of binary (0/1) covariates; A, Y: binary (0/1) treatment and outcome
set.seed(42)
n <- 300
X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5))
A <- rbinom(n, 1, plogis(0.5 * X$X1 - 0.2))
Y <- rbinom(n, 1, 0.3 + 0.2 * X$X1 + 0.15 * A)
fit <- estimate_att(X, A, Y, K = 5)
fit$theta   # point estimate
fit$ci_95   # 95% Wald CI
```

**Continuous outcome:** Use `outcome_type = "continuous"` and numeric Y. Requires optimaltrees to support `squared_error` loss for the outcome trees.

## Running Simulations

### Local Simulations

```r
# See simulations/run_simulations.R for examples
source("simulations/production/dgps/dgps_smooth.R")

# Generate data
d <- generate_dgp_binary_att(n = 400, tau = 0.10, seed = 123)

# Estimate ATT
fit <- estimate_att(
  X = d$X, A = d$A, Y = d$Y,
  K = 5,
  regularization = log(400) / 400,
  use_rashomon = FALSE
)

print(fit$theta)  # Point estimate
print(fit$ci)     # 95% CI
```

### Large-Scale Simulations on O2 Cluster

Complete SLURM infrastructure for distributed simulations:

```bash
# On O2
cd doubletree/simulations/production

# Quick test (30 seconds)
bash slurm/quick_test.sh

# Launch all 18,000 simulations
bash slurm/launch_all_simulations.sh

# Monitor progress
bash slurm/check_progress.sh

# Combine results (after completion)
Rscript slurm/combine_results.R
```

**Simulation grid:**
- 3 DGPs × 4 methods × 3 sample sizes × 500 replications = 18,000 runs
- Estimated time: 30-60 minutes (vs. 15 hours sequential)
- Output: Individual .rds files → combined dataset + summary statistics

See `simulations/production/slurm/README_O2.md` for complete documentation.

## Package Development

This project is structured as an R package at the root level, making it easy to extract and publish to GitHub as a standalone package in the future.

### Development Workflow

1. Add functions to `R/`
2. Document functions using roxygen2 comments
3. Generate documentation: `devtools::document()`
4. Write tests in `tests/testthat/`
5. Run tests: `devtools::test()`

## Links

- **GitHub Repository:** [github.com/denisagniel/doubletree](https://github.com/denisagniel/doubletree)
- **optimaltrees (treefarmr):** [github.com/denisagniel/treefarmr](https://github.com/denisagniel/treefarmr)
- **Manuscript:** `paper/manuscript.tex` (theory complete, restructuring in progress)
- **O2 Simulation Docs:** `simulations/production/slurm/README_O2.md`

## Citation

If you use doubletree in your research, please cite:

```bibtex
@software{doubletree2026,
  title = {doubletree: Causal Estimation with Interpretable Trees},
  author = {Denis Agniel},
  year = {2026},
  url = {https://github.com/denisagniel/doubletree}
}
```

## License

MIT License - see LICENSE file for details.
