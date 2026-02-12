# dmltree: DML Causal Estimation with Interpretable Trees

This package implements causal inference for the **Average Treatment effect on the Treated (ATT)** using Double Machine Learning (DML) with interpretable optimal decision trees. It depends on [treefarmr](https://github.com/) for fitting the nuisance functions (propensity and outcome trees). Theory-aligned API expectations for the tree side are described in `paper/Implementation-requirements-Rashomon-DML.md`.

**Note:** Only **binary outcome Y** is currently supported. Continuous outcome is not implemented yet (no fallback).

## Project Structure

```
dmltree/
├── .Rprofile              # Auto-loads treefarmr from local dev
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

## Setup

### Prerequisites

1. **R** (version 4.0 or higher recommended)
2. **treefarmr package**: This project depends on `treefarmr` for tree-based modeling. The `.Rprofile` file will automatically attempt to load it from a sibling directory (`../treefarmr`) if available.

### Installing Dependencies

```r
# Install required packages
install.packages(c("devtools", "testthat"))

# If treefarmr is not in ../treefarmr, install it from its source
# (Update path as needed)
devtools::install("../treefarmr")
```

### Loading the Package

Since this is both a research project and an R package, you can load the package in development mode:

```r
devtools::load_all()
```

Or install it locally:

```r
devtools::install()
library(dmltree)
```

## Using treefarmr

This project builds on the `treefarmr` package. If `treefarmr` is located in a sibling directory (`../treefarmr`), it will be automatically loaded when you start R in this directory (via `.Rprofile`).

If `treefarmr` is located elsewhere:
- Update the path in `.Rprofile`, or
- Manually load it: `devtools::load_all("/path/to/treefarmr")`

### Rashomon DML (future)

The manuscript selects a **single interpretable tree per nuisance** via the intersection of Rashomon sets across cross-fitting folds, then refits that structure per fold for valid DML. The current package fits one optimal tree per fold (no Rashomon or intersection). To support the full Rashomon-DML workflow, **treefarmr** would need to implement the API described in `paper/Implementation-requirements-Rashomon-DML.md`: Rashomon set per fold (tolerance, optional max_leaves), structure comparison and intersection across folds, and refit-by-structure for fold-specific nuisances. dmltree would then call those treefarmr APIs; no dmltree code changes are required for estimation validity (fold-specific nuisances are already used).

## Minimal example

```r
devtools::load_all()
# X: data.frame of binary (0/1) covariates; A, Y: binary (0/1) treatment and outcome
set.seed(42)
n <- 300
X <- data.frame(X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5))
A <- rbinom(n, 1, plogis(0.5 * X$X1 - 0.2))
Y <- rbinom(n, 1, 0.3 + 0.2 * X$X1 + 0.15 * A)
fit <- dml_att(X, A, Y, K = 5)
fit$theta   # point estimate
fit$ci_95   # 95% Wald CI
```

## Running Simulations

See `simulations/run_simulations.R` for a full example (DGP and replications).

## Package Development

This project is structured as an R package at the root level, making it easy to extract and publish to GitHub as a standalone package in the future.

### Development Workflow

1. Add functions to `R/`
2. Document functions using roxygen2 comments
3. Generate documentation: `devtools::document()`
4. Write tests in `tests/testthat/`
5. Run tests: `devtools::test()`

## Links

- **treefarmr**: [Link to treefarmr repository/package]

## License

[To be determined]
