# Local coverage check: 20 reps × 6 approaches × 4 DGPs
#
# Purpose: Verify coverage rates after CI formula fixes (approaches 4, 5, 6).
# Uses INSTALLED packages (library()) to match cluster environment exactly.
#
# Expected (n=500, 20 reps):
#   theta:    near 0.15 (true ATT)
#   bias:     |mean(theta) - 0.15| < 0.05
#   coverage: roughly 15-20/20 (75-100%) — with 20 reps, expect ~18/20 at 95% CI
#   width:    ~0.17-0.22, NOT ~0.007
#
# Runtime: ~30-60 min (approaches 5/6 slow; M=3, K=3 to keep manageable)
#
# Run from doubletree package root:
#   Rscript simulations/six_approach_comparison/local_coverage_check.R

suppressPackageStartupMessages({
  library(optimaltrees)
  library(doubletree)
})

cat("optimaltrees:", as.character(packageVersion("optimaltrees")), "\n")
cat("doubletree:  ", as.character(packageVersion("doubletree")), "\n\n")

SIM_DIR <- tryCatch(
  dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
    grep("--file=", commandArgs(trailingOnly = FALSE))
  ][1], mustWork = FALSE)),
  error = function(e) "."
)
if (!file.exists(file.path(SIM_DIR, "code/dgps.R"))) {
  SIM_DIR <- "/Users/dagniel/RAND/rprojects/global-scholars/doubletree/simulations/six_approach_comparison"
}
source(file.path(SIM_DIR, "code/dgps.R"))

# ---- parameters ---------------------------------------------------------------
# Production vs local differences:
#   Production (cluster): M=10, K_MSPLIT=5  (see code/estimators.R SIM_M, SIM_K_MSPLIT)
#   Local validation:     M=3,  K_MSPLIT=3  (reduced for manageable runtime ~30-60 min)
# All other parameters match production: K=5, eps_n=2*sqrt(log(n)/n) for appr 3/4,
# cv_regularization_adaptive with max_lambda=20*log(n)/n for all 6 approaches

N        <- 500
TRUE_ATT <- 0.15
N_REPS   <- 20
K        <- 5    # folds for approaches 1-4 (matches production)
M        <- 3    # splits for approaches 5/6 (production: M=10)
K_MSPLIT <- 3    # folds per split for approaches 5/6 (production: K_MSPLIT=5)

set.seed(42)

dgp_fns <- list(
  simple     = generate_dgp_simple,
  moderate   = generate_dgp_moderate,
  complex    = generate_dgp_complex,
  continuous = generate_dgp_continuous
)

approach_names <- c(
  "1_fullsample",
  "2_crossfit",
  "3_doubletree",
  "4_dt_averaged",
  "5_msplit",
  "6_msplit_avg"
)

# ---- run one rep for one approach on one dataset ------------------------------

run_one <- function(approach_idx, X, A, Y) {
  tryCatch({
    switch(approach_idx,
      # 1: fullsample
      {
        n_all <- length(Y)
        cv_e  <- cv_regularization_adaptive(X, A, loss_function = "log_loss",
                                            K = 5,
                                            max_lambda = 20 * log(n_all) / n_all,
                                            refit = TRUE, verbose = FALSE)
        e_hat <- predict(cv_e$model, X, type = "prob")[, 2L]
        X0    <- X[A == 0, , drop = FALSE]; Y0 <- Y[A == 0]
        n0    <- length(Y0)
        cv_m0  <- cv_regularization_adaptive(X0, Y0, loss_function = "log_loss",
                                             K = 5,
                                             max_lambda = 20 * log(n0) / n0,
                                             refit = TRUE, verbose = FALSE)
        m0_hat <- predict(cv_m0$model, X, type = "prob")[, 2L]
        psi    <- doubletree:::psi_att(Y, A, 0,
                                       list(e = e_hat, m0 = m0_hat, m1 = NULL),
                                       mean(A))
        theta  <- sum(psi) / sum(A / mean(A))
        sigma  <- att_se(doubletree:::psi_att(Y, A, theta,
                                              list(e = e_hat, m0 = m0_hat, m1 = NULL),
                                              mean(A)))
        list(theta = theta, sigma = sigma, ci = att_ci(theta, sigma))
      },
      # 2: crossfit
      {
        r <- estimate_att(X, A, Y, K = K, outcome_type = "binary",
                          use_rashomon = FALSE, cv_regularization = TRUE,
                          verbose = FALSE)
        list(theta = r$theta, sigma = r$sigma, ci = att_ci(r$theta, r$sigma))
      },
      # 3: doubletree
      {
        eps_n <- 2 * sqrt(log(length(Y)) / length(Y))
        r <- estimate_att(X, A, Y, K = K, outcome_type = "binary",
                          use_rashomon = TRUE, cv_regularization = TRUE,
                          rashomon_bound_multiplier = eps_n,
                          auto_tune_intersecting = TRUE,
                          verbose = FALSE)
        if (!isTRUE(r$converged)) {
          stop("Approach 3: Rashomon intersection failed after auto-tuning. ",
               "Hard failure -- do not silently proceed. converged=", r$converged)
        }
        list(theta = r$theta, sigma = r$sigma, ci = att_ci(r$theta, r$sigma))
      },
      # 4: doubletree averaged
      {
        eps_n <- 2 * sqrt(log(length(Y)) / length(Y))
        r <- estimate_att_doubletree_averaged(X, A, Y, K = K,
                                              outcome_type = "binary",
                                              rashomon_bound_multiplier = eps_n,
                                              auto_tune_intersecting = TRUE,
                                              verbose = FALSE)
        list(theta = r$theta, sigma = r$sigma, ci = r$ci_95)
      },
      # 5: msplit
      {
        r <- estimate_att_msplit(X, A, Y, M = M, K = K_MSPLIT,
                                 outcome_type = "binary", verbose = FALSE)
        list(theta = r$theta, sigma = r$sigma, ci = r$ci_95)
      },
      # 6: msplit averaged
      {
        r <- estimate_att_msplit_averaged(X, A, Y, M = M, K = K_MSPLIT,
                                          outcome_type = "binary",
                                          seed_base = 100, verbose = FALSE)
        list(theta = r$theta, sigma = r$sigma, ci = r$ci_95)
      }
    )
  }, error = function(e) list(theta = NA_real_, sigma = NA_real_, ci = c(NA, NA),
                               err = conditionMessage(e)))
}

# ---- main loop ----------------------------------------------------------------

all_results <- list()

for (dgp_name in names(dgp_fns)) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("DGP: ", toupper(dgp_name), " (n=", N, ", true ATT=", TRUE_ATT, ")\n", sep = "")
  cat(strrep("=", 70), "\n", sep = "")

  key <- dgp_name
  all_results[[key]] <- vector("list", 6)
  names(all_results[[key]]) <- approach_names

  for (appr in seq_along(approach_names)) {
    aname <- approach_names[appr]
    t0 <- proc.time()["elapsed"]

    thetas  <- numeric(N_REPS)
    sigmas  <- numeric(N_REPS)
    covers  <- logical(N_REPS)
    widths  <- numeric(N_REPS)
    errors  <- character(N_REPS)

    for (r in seq_len(N_REPS)) {
      set.seed(1000 + appr * 10000 + match(dgp_name, names(dgp_fns)) * 1000 + r)
      d <- dgp_fns[[dgp_name]](N)
      res <- run_one(appr, d$X, d$A, d$Y)

      thetas[r] <- res$theta
      sigmas[r] <- res$sigma
      ci <- res$ci
      if (!is.null(ci) && length(ci) == 2 && all(is.finite(ci))) {
        covers[r] <- (ci[1] <= TRUE_ATT && TRUE_ATT <= ci[2])
        widths[r] <- diff(ci)
      } else {
        covers[r] <- NA
        widths[r] <- NA
      }
      if (!is.null(res$err)) errors[r] <- res$err
    }

    elapsed <- as.integer(proc.time()["elapsed"] - t0)
    n_ok     <- sum(!is.na(thetas))
    n_err    <- sum(nchar(errors) > 0)
    coverage <- mean(covers, na.rm = TRUE)
    bias     <- mean(thetas, na.rm = TRUE) - TRUE_ATT
    mean_w   <- mean(widths, na.rm = TRUE)

    cat(sprintf("  %-16s  n_ok=%2d  bias=%+.3f  coverage=%.0f%%  width=%.3f  (%ds)\n",
                aname, n_ok, bias, 100 * coverage, mean_w, elapsed))
    if (n_err > 0) {
      cat(sprintf("    ERRORS (%d): %s\n", n_err,
                  substr(paste(unique(errors[nchar(errors)>0]), collapse=" | "), 1, 80)))
    }

    all_results[[key]][[aname]] <- list(
      thetas = thetas, sigmas = sigmas, covers = covers, widths = widths,
      errors = errors, n_ok = n_ok, n_err = n_err,
      coverage = coverage, bias = bias, mean_width = mean_w
    )
  }
}

# ---- summary ------------------------------------------------------------------

cat("\n\n", strrep("=", 80), "\n", sep = "")
cat("SUMMARY TABLE  (n=", N, ", ", N_REPS, " reps per cell)\n", sep = "")
cat("Expected: coverage ~90-95%, width ~0.17-0.22, |bias| < 0.05\n")
cat(strrep("=", 80), "\n", sep = "")
cat(sprintf("%-12s  %-16s  %6s  %8s  %6s  %5s\n",
            "DGP", "Approach", "n_ok", "coverage", "width", "bias"))
cat(strrep("-", 70), "\n", sep = "")

flag_width_ok <- function(w) if (!is.na(w) && w < 0.05) " *** TOO NARROW" else ""
flag_cov      <- function(c) if (!is.na(c) && c < 0.70) " *** LOW" else ""

for (dgp_name in names(all_results)) {
  for (aname in approach_names) {
    r <- all_results[[dgp_name]][[aname]]
    cat(sprintf("%-12s  %-16s  %6d  %7.1f%%  %6.3f  %+5.3f%s%s\n",
                dgp_name, aname, r$n_ok,
                100 * r$coverage, r$mean_width, r$bias,
                flag_width_ok(r$mean_width),
                flag_cov(r$coverage)))
  }
}

cat("\nDone.\n")
