# Diagnostic script for beta study DGPs
# Analyzes tree fit quality and fundamental information loss

library(dplyr)
suppressPackageStartupMessages({
  devtools::load_all("../../../optimaltrees")
  source("dgps/dgps_beta_regimes.R")
})

cat('\n')
cat(strrep('=', 75), '\n')
cat('BETA STUDY DGP DIAGNOSTIC\n')
cat(strrep('=', 75), '\n\n')

# Generate example data
set.seed(12345)
n <- 800
lambda <- log(n) / n

cat('1. FUNDAMENTAL PROBLEM: INFORMATION LOSS\n')
cat(strrep('-', 75), '\n\n')

d <- generate_dgp_beta_high(n = n, tau = 0.10)

# Analyze variation within binary groups
variation <- data.frame(
  X1 = d$X$X1,
  X2 = d$X$X2,
  true_e = d$true_e,
  true_m0 = d$true_m0
) %>%
  group_by(X1, X2) %>%
  summarize(
    n_obs = n(),
    e_mean = mean(true_e),
    e_range = max(true_e) - min(true_e),
    m0_mean = mean(true_m0),
    m0_range = max(true_m0) - min(true_m0),
    .groups = 'drop'
  )

cat('Within-group variation (beta-high, n=800):\n\n')
print(variation)

cat('\nSummary:\n')
cat(sprintf('  • Propensity e varies by %.3f within groups (mean)\n',
            mean(variation$e_range)))
cat(sprintf('  • Outcome m0 varies by %.3f within groups (mean)\n',
            mean(variation$m0_range)))
cat('  • Trees can only predict group averages\n')
cat('  • This creates UNAVOIDABLE approximation error\n\n')

# Fit trees quietly
cat('2. TREE FIT QUALITY (regularization = log(n)/n = ',
    sprintf('%.6f', lambda), ')\n')
cat(strrep('-', 75), '\n\n')

cat('Fitting propensity tree...\n')
e_fit <- optimaltrees::fit_tree(
  X = d$X, y = d$A, loss = 'log_loss',
  regularization = lambda, verbose = FALSE
)

cat('Fitting outcome tree (controls only)...\n')
m0_fit <- optimaltrees::fit_tree(
  X = d$X[d$A == 0, ],
  y = d$Y[d$A == 0],
  loss = 'log_loss',
  regularization = lambda,
  verbose = FALSE
)

# Get predictions
e_pred_obj <- predict(e_fit, d$X)
if (is.list(e_pred_obj)) {
  e_pred <- e_pred_obj$probabilities[, 2]
} else {
  e_pred <- e_pred_obj[, 2]
}

m0_pred_obj <- predict(m0_fit, d$X)
if (is.list(m0_pred_obj)) {
  m0_pred <- m0_pred_obj$probabilities[, 2]
} else {
  m0_pred <- m0_pred_obj[, 2]
}

cat('\n')
cat('PROPENSITY TREE (e):\n')
cat(sprintf('  Leaves: %d\n', e_fit$n_leaves))
cat(sprintf('  Depth: %d\n', e_fit$depth))
cat(sprintf('  MAE: %.4f\n', mean(abs(d$true_e - e_pred))))
cat(sprintf('  RMSE: %.4f\n', sqrt(mean((d$true_e - e_pred)^2))))
cat(sprintf('  Correlation: %.4f\n', cor(d$true_e, e_pred)))

cat('\nOUTCOME TREE (m0):\n')
cat(sprintf('  Leaves: %d\n', m0_fit$n_leaves))
cat(sprintf('  Depth: %d\n', m0_fit$depth))
cat(sprintf('  MAE: %.4f\n', mean(abs(d$true_m0 - m0_pred))))
cat(sprintf('  RMSE: %.4f\n', sqrt(mean((d$true_m0 - m0_pred)^2))))
cat(sprintf('  Correlation: %.4f\n', cor(d$true_m0, m0_pred)))

# Analyze predictions vs. truth by group
cat('\n3. PREDICTIONS BY (X1, X2) GROUP\n')
cat(strrep('-', 75), '\n\n')

pred_analysis <- data.frame(
  X1 = d$X$X1,
  X2 = d$X$X2,
  true_e = d$true_e,
  pred_e = e_pred,
  true_m0 = d$true_m0,
  pred_m0 = m0_pred
) %>%
  group_by(X1, X2) %>%
  summarize(
    n = n(),
    true_e_mean = mean(true_e),
    pred_e_mean = mean(pred_e),
    e_error = mean(abs(true_e - pred_e)),
    true_m0_mean = mean(true_m0),
    pred_m0_mean = mean(pred_m0),
    m0_error = mean(abs(true_m0 - pred_m0)),
    .groups = 'drop'
  )

print(pred_analysis, n = Inf)

cat('\n4. TESTING ALL THREE BETA REGIMES\n')
cat(strrep('-', 75), '\n\n')

test_regime <- function(dgp_func, name, beta) {
  set.seed(12345)
  d <- dgp_func(n = 800, tau = 0.10)

  # Fit trees
  e_fit <- optimaltrees::fit_tree(
    X = d$X, y = d$A, loss = 'log_loss',
    regularization = log(800)/800, verbose = FALSE
  )
  m0_fit <- optimaltrees::fit_tree(
    X = d$X[d$A == 0, ], y = d$Y[d$A == 0],
    loss = 'log_loss', regularization = log(800)/800, verbose = FALSE
  )

  # Predictions
  e_pred <- predict(e_fit, d$X)$probabilities[, 2]
  m0_pred <- predict(m0_fit, d$X)$probabilities[, 2]

  # Errors
  e_rmse <- sqrt(mean((d$true_e - e_pred)^2))
  m0_rmse <- sqrt(mean((d$true_m0 - m0_pred)^2))
  e_cor <- cor(d$true_e, e_pred)
  m0_cor <- cor(d$true_m0, m0_pred)

  cat(sprintf('%-20s (β=%d):\n', name, beta))
  cat(sprintf('  Propensity:  RMSE=%.4f, Cor=%.3f, %d leaves\n',
              e_rmse, e_cor, e_fit$n_leaves))
  cat(sprintf('  Outcome:     RMSE=%.4f, Cor=%.3f, %d leaves\n',
              m0_rmse, m0_cor, m0_fit$n_leaves))
  cat('\n')
}

test_regime(generate_dgp_beta_high, "Beta-high", 3)
test_regime(generate_dgp_beta_boundary, "Beta-boundary", 2)
test_regime(generate_dgp_beta_low, "Beta-low", 1)

cat(strrep('=', 75), '\n')
cat('DIAGNOSIS COMPLETE\n')
cat(strrep('=', 75), '\n')
