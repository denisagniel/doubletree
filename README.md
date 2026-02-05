# dmltree: DML Causal Estimation with Interpretable Trees

This research project implements and evaluates methods for causal inference using Double Machine Learning (DML) with interpretable tree-based models.

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

## Running Simulations

See `simulations/run_simulations.R` for the main simulation script template.

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
