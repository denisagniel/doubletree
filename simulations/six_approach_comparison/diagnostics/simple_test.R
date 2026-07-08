#!/usr/bin/env Rscript
# Simple test of diagnostic utilities
library(optimaltrees)
library(ggplot2)

source('diagnostics/utils/tree_diagnostics.R')
source('diagnostics/utils/eif_components.R')
source('code/dgps.R')

cat('Testing basic functionality...\n\n')

# Test 1: Generate data
data <- generate_dgp_complex(n = 100)
cat('✓ DGP generation works\n')

# Test 2: Fit tree
tree <- optimaltrees::fit_tree(
  X = data$X,
  y = data$A,
  loss_function = 'log_loss',
  regularization = 0.01
)
cat('✓ Tree fitting works\n')

# Test 3: Predict
e_hat <- predict(tree, data$X)
cat('✓ Prediction works\n')

# Test 4: Count leaves
n_leaves <- count_leaves(tree)
cat(sprintf('✓ Count leaves works: %d leaves\n', n_leaves))

# Test 5: Metrics
rmse <- sqrt(mean((e_hat - data$e_true)^2))
cat(sprintf('✓ Metrics work: RMSE = %.4f\n', rmse))

# Test 6: EIF components
e_true <- compute_true_propensity(data$X, dgp = 3)
mu0_true <- compute_true_outcome(data$X, dgp = 3)
cat('✓ True nuisance functions work\n')

cat('\n✓✓✓ All basic tests passed! ✓✓✓\n')
cat('\nDiagnostic infrastructure is functional.\n')
cat('Ready to run full diagnostics.\n')
