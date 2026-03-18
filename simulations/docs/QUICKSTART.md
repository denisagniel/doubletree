# Quick Start: DML-ATT Simulations

**Last updated:** 2026-03-04 (restructured for manuscript production)

---

## Directory Structure

```
simulations/
├── production/          # Manuscript simulations (use these)
│   ├── dgps/           # Validated DGPs
│   ├── methods/        # Baseline comparison methods
│   ├── run_primary.R   # Main simulation
│   ├── run_stress.R    # Stress tests
│   └── analyze_manuscript.R
├── diagnostics/        # Testing and verification
├── deprecated/         # Misspecified DGPs (do not use)
└── exploration/        # Archived test scripts
```

---

## Minimal Working Example

```r
# Load packages
devtools::load_all("/path/to/optimaltrees")
devtools::load_all("/path/to/doubletree")
source("production/dgps/dgps_smooth.R")

# Generate data
n <- 400
tau <- 0.10
d <- generate_dgp_binary_att(n, tau = tau, seed = 123)

# Fit DML-ATT with fixed regularization
result <- doubletree::dml_att(
  X = d$X,
  A = d$A,
  Y = d$Y,
  K = 5,
  regularization = log(n) / n,  # Theory-driven choice
  cv_regularization = FALSE,     # CV not beneficial (see SUMMARY)
  verbose = TRUE
)

# Extract results
cat(sprintf("Estimated ATT: %.4f\n", result$theta))
cat(sprintf("Standard error: %.4f\n", result$sigma))
cat(sprintf("95%% CI: [%.4f, %.4f]\n", result$ci[1], result$ci[2]))
cat(sprintf("True ATT: %.4f\n", tau))
```

**Expected output:**
- Estimated ATT ≈ 0.10 (close to true value)
- 95% CI covers true ATT in ~95% of replications

---

## Available DGPs

### 1. Binary Features (Recommended)

```r
generate_dgp_binary_att(n, tau = 0.10, seed = NULL)
```

- **Features:** 4 binary (X1, X2, X3, X4)
- **Signal:** X1, X2 in both propensity and outcome
- **Noise:** X3, X4
- **Patterns:** 2^4 = 16 covariate patterns
- **Use case:** Fast, reliable, good for testing

### 2. Continuous Features

```r
generate_dgp_continuous_att(n, tau = 0.10, seed = NULL)
```

- **Features:** 4 continuous (X1, X2, X3, X4)
- **Signal:** X1, X2 in both functions
- **Discretization:** Automatic via `fit_tree()` (adaptive bins)
- **Use case:** More realistic, tests discretization workflow

### 3. Moderate Complexity

```r
generate_dgp_moderate_att(n, tau = 0.10, seed = NULL)
```

- **Features:** 5 binary (X1, X2, X3, X4, X5)
- **Signal:** X1, X2, X3
- **Noise:** X4, X5
- **Patterns:** 2^5 = 32 covariate patterns
- **Use case:** Higher complexity, more challenging

---

## Recommended Simulation Setup

### Sample Sizes

```r
n_values <- c(400, 800, 1600)  # Standard choices
tau <- 0.10                    # Effect size (10pp)
K <- 5                         # DML folds
n_reps <- 500                  # For publication-quality results
```

### Regularization

**Use fixed λ = log(n)/n** (not CV-based):
- Achieves 95% coverage
- Theory-aligned
- Computationally efficient

### Metrics to Report

```r
# For each replication:
theta_hat <- result$theta
sigma_hat <- result$sigma
ci_lower <- result$ci[1]
ci_upper <- result$ci[2]

# Summary statistics:
bias <- mean(theta_hats) - tau
coverage <- mean(ci_lower <= tau & ci_upper >= tau)
rmse <- sqrt(mean((theta_hats - tau)^2))
mean_ci_width <- mean(ci_upper - ci_lower)
```

---

## Full Simulation Script Template

```r
# setup
devtools::load_all("/path/to/optimaltrees")
devtools::load_all("/path/to/doubletree")
source("simulations/dgps_att_correct.R")

# parameters
n_values <- c(400, 800, 1600)
tau <- 0.10
K <- 5
n_reps <- 500

# storage
results <- expand.grid(
  n = n_values,
  rep = 1:n_reps,
  theta = NA,
  sigma = NA,
  ci_lower = NA,
  ci_upper = NA
)

# simulation loop
for (i in 1:nrow(results)) {
  n <- results$n[i]

  # Generate data
  d <- generate_dgp_binary_att(n, tau = tau, seed = 10000 + i)

  # Fit DML
  capture.output({  # Suppress GOSDT output
    fit <- tryCatch({
      doubletree::dml_att(
        d$X, d$A, d$Y,
        K = K,
        regularization = log(n) / n,
        cv_regularization = FALSE,
        verbose = FALSE
      )
    }, error = function(e) NULL)
  }, file = tempfile())

  # Store results
  if (!is.null(fit)) {
    results$theta[i] <- fit$theta
    results$sigma[i] <- fit$sigma
    results$ci_lower[i] <- fit$ci[1]
    results$ci_upper[i] <- fit$ci[2]
  }

  # Progress
  if (i %% 50 == 0) cat(sprintf("Completed %d/%d\n", i, nrow(results)))
}

# Compute summary statistics
library(dplyr)
summary_stats <- results %>%
  filter(!is.na(theta)) %>%
  group_by(n) %>%
  summarize(
    n_valid = n(),
    bias = mean(theta) - tau,
    rmse = sqrt(mean((theta - tau)^2)),
    coverage = mean(ci_lower <= tau & ci_upper >= tau),
    mean_ci_width = mean(ci_upper - ci_lower),
    .groups = "drop"
  )

print(summary_stats)
```

**Expected results:**
- Bias ≈ 0 (within ±0.02)
- Coverage ≈ 95% (within 93-97%)
- RMSE decreases with n (~ 1/√n)

---

## Troubleshooting

### Issue: "Model limit exceeded"

**Cause:** Too many features or patterns for GOSDT

**Solution:**
- Use binary_att (4 features) instead of moderate_att (5 features)
- Or increase `model_limit` (not recommended)

### Issue: Poor coverage (<90%)

**Check:**
1. Using corrected DGPs from `dgps_att_correct.R`? (Not dgps_realistic.R)
2. True ATT = specified tau? Run `check_dgp_att.R` to verify
3. Regularization = log(n)/n?

### Issue: Large bias (>10% of truth)

**Likely causes:**
- Using wrong DGP (misspecified tau)
- Too-strong regularization (underfitting)
- Too-weak regularization (overfitting)

**Solution:** Use fixed λ = log(n)/n, verify with `verify_att_dgp.R`

---

## What NOT to Use

❌ **deprecated/** directory - Contains misspecified DGPs with tau on logit scale
❌ **exploration/** directory - 25+ test scripts for development only

These files are archived for reference but should not be used for manuscript simulations.

---

## Production Scripts (Coming Soon)

- `production/run_primary.R` - Main simulations for manuscript Table 1
- `production/run_stress.R` - Stress tests for manuscript Table 2
- `production/analyze_manuscript.R` - Generate tables and figures

## References

- **Verified DGPs:** `production/dgps/dgps_smooth.R` (95% coverage validated)
- **Deprecated DGPs:** `deprecated/README.md` (explanation of tau bug)
- **Verification script:** `diagnostics/verify_dgp_att.R`
