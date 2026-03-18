# debug_cv_lambdas.R
# Investigate what lambda values CV is selecting and why coverage is so bad

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools required")
devtools::load_all()
if (!requireNamespace("optimaltrees", quietly = TRUE)) stop("optimaltrees required")

source("simulations/run_simulations_extended.R", local = TRUE)

n <- 400
tau <- 0.15
K <- 5

dgp_fn <- generate_data_dgp1
d <- dgp_fn(n, tau, seed = 123)

message("=== Debugging CV Lambda Selection ===\n")

# 1. Check what CV selects for each nuisance on one dataset
message("1. Lambda selection via CV for one dataset:")

# Propensity
message("\n  Propensity e(X):")
cv_e <- optimaltrees::cv_regularization(d$X, d$A, loss_function = "log_loss",
                                     K = 5, refit = FALSE, verbose = FALSE)
message(sprintf("    Lambda grid: %s", paste(round(cv_e$lambda_grid, 5), collapse=", ")))
message(sprintf("    CV loss: %s", paste(round(cv_e$cv_loss, 4), collapse=", ")))
message(sprintf("    Selected lambda: %.5f", cv_e$best_lambda))
message(sprintf("    Fixed lambda (log(n)/n): %.5f", log(n)/n))

# m0
idx0 <- which(d$A == 0)
message("\n  Outcome m0(X) (A=0, n0=", length(idx0), "):")
cv_m0 <- optimaltrees::cv_regularization(d$X[idx0, , drop=FALSE], d$Y[idx0],
                                      loss_function = "log_loss",
                                      K = 5, refit = FALSE, verbose = FALSE)
message(sprintf("    Lambda grid: %s", paste(round(cv_m0$lambda_grid, 5), collapse=", ")))
message(sprintf("    CV loss: %s", paste(round(cv_m0$cv_loss, 4), collapse=", ")))
message(sprintf("    Selected lambda: %.5f", cv_m0$best_lambda))

# m1
idx1 <- which(d$A == 1)
message("\n  Outcome m1(X) (A=1, n1=", length(idx1), "):")
cv_m1 <- optimaltrees::cv_regularization(d$X[idx1, , drop=FALSE], d$Y[idx1],
                                      loss_function = "log_loss",
                                      K = 5, refit = FALSE, verbose = FALSE)
message(sprintf("    Lambda grid: %s", paste(round(cv_m1$lambda_grid, 5), collapse=", ")))
message(sprintf("    CV loss: %s", paste(round(cv_m1$cv_loss, 4), collapse=", ")))
message(sprintf("    Selected lambda: %.5f", cv_m1$best_lambda))

# 2. Compare estimates with fixed vs CV lambda
message("\n\n2. Comparing estimates (1 rep only):")

result_fixed <- estimate_att(
  d$X, d$A, d$Y, K = K,
  use_rashomon = FALSE,
  regularization = log(n)/n,
  verbose = FALSE,
  seed = 123
)

result_cv <- estimate_att(
  d$X, d$A, d$Y, K = K,
  use_rashomon = FALSE,
  cv_regularization = TRUE,
  cv_K = 5,
  verbose = FALSE,
  seed = 123
)

message(sprintf("\n  Fixed (lambda=%.5f):", log(n)/n))
message(sprintf("    theta: %.4f (true: %.2f, bias: %.4f)", result_fixed$theta, tau, result_fixed$theta - tau))
message(sprintf("    sigma: %.4f (sqrt(n) scale)", result_fixed$sigma))
message(sprintf("    CI: [%.4f, %.4f]", result_fixed$ci_95[1], result_fixed$ci_95[2]))
message(sprintf("    Contains truth? %s", ifelse(result_fixed$ci_95[1] <= tau & tau <= result_fixed$ci_95[2], "YES", "NO")))

message("\n  CV (lambda selected per nuisance):")
message(sprintf("    theta: %.4f (true: %.2f, bias: %.4f)", result_cv$theta, tau, result_cv$theta - tau))
message(sprintf("    sigma: %.4f (sqrt(n) scale)", result_cv$sigma))
message(sprintf("    CI: [%.4f, %.4f]", result_cv$ci_95[1], result_cv$ci_95[2]))
message(sprintf("    Contains truth? %s", ifelse(result_cv$ci_95[1] <= tau & tau <= result_cv$ci_95[2], "YES", "NO")))

message("\n=== Diagnosis ===")
if (abs(result_cv$theta - tau) > abs(result_fixed$theta - tau)) {
  message("CV has HIGHER bias than fixed lambda")
  message("Likely cause: CV selecting too large lambda (underfitting)")
} else {
  message("CV has same or lower bias than fixed lambda")
}

message("\nPossible issues:")
message("1. CV grid may be poorly chosen for DML (default grid assumes n for full data, but nuisances use n/K)")
message("2. For m0/m1, sample sizes are smaller (n0, n1) and CV should account for this")
message("3. Default grid (log(n)/n) * c(0.25, 0.5, 1, 2, 4) may be too coarse")
