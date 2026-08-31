# Local smoke test: all six approaches across all four DGPs
#
# Purpose: Verify correctness (not coverage) after CI-formula and
# silent-fallback fixes. One rep per DGP; M=3, K=3 for approaches 5/6
# to keep runtime manageable (~10-15 min total).
#
# Checks:
#   - All approaches run without error
#   - theta is finite and in plausible range
#   - sigma is positive and finite
#   - CI width reasonable (approaches 4/6: ~0.15-0.20, NOT ~0.007)
#   - CI contains true ATT (one rep, so this may fail by chance)
#
# Run from the doubletree package root:
#   Rscript simulations/six_approach_comparison/test_all_six_local.R

PKGROOT <- "/Users/dagniel/RAND/rprojects/global-scholars"
devtools::load_all(file.path(PKGROOT, "optimaltrees"), quiet = TRUE)
devtools::load_all(file.path(PKGROOT, "doubletree"),    quiet = TRUE)

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
    grep("--file=", commandArgs(trailingOnly = FALSE))
  ][1], mustWork = FALSE)),
  error = function(e) "."
)
# Fallback: known absolute path
if (!file.exists(file.path(SCRIPT_DIR, "code/dgps.R"))) {
  SCRIPT_DIR <- "/Users/dagniel/RAND/rprojects/global-scholars/doubletree/simulations/six_approach_comparison"
}
source(file.path(SCRIPT_DIR, "code/dgps.R"))

set.seed(2026)

N        <- 500
TRUE_ATT <- 0.15
K        <- 5   # folds (approaches 1-4)
M_SMALL  <- 3   # splits for approaches 5/6 (reduced for local speed)
K_SMALL  <- 3   # folds per split for approaches 5/6

dgp_fns <- list(
  simple     = generate_dgp_simple,
  moderate   = generate_dgp_moderate,
  complex    = generate_dgp_complex,
  continuous = generate_dgp_continuous
)

# ---- helper: run one approach, print results ---------------------------------

run_approach <- function(label, fn) {
  t0 <- proc.time()["elapsed"]
  result <- tryCatch(fn(), error = function(e) list(err = conditionMessage(e)))
  elapsed <- as.integer(proc.time()["elapsed"] - t0)

  if (!is.null(result$err)) {
    cat(sprintf("  %-22s  ERROR: %s\n", label, substr(result$err, 1, 70)))
    return(invisible(NULL))
  }

  theta <- result$theta
  sigma <- result$sigma
  ci    <- result$ci_95

  if (is.null(ci) || length(ci) != 2 || any(!is.finite(ci))) {
    width  <- NA_real_
    covers <- NA
  } else {
    width  <- diff(ci)
    covers <- ci[1] <= TRUE_ATT && TRUE_ATT <= ci[2]
  }

  cat(sprintf(
    "  %-22s  theta=%+.3f  SE=%.3f  CI=[%+.3f, %+.3f]  width=%.3f  covers=%s  (%ds)\n",
    label,
    if (is.finite(theta)) theta else NaN,
    if (is.finite(sigma)) sigma else NaN,
    if (!is.na(width)) ci[1] else NaN,
    if (!is.na(width)) ci[2] else NaN,
    if (!is.na(width)) width else NaN,
    if (isTRUE(covers)) "YES" else if (isFALSE(covers)) "NO " else "???",
    elapsed
  ))

  invisible(list(theta = theta, sigma = sigma, ci = ci, width = width, covers = covers))
}

# ---- main loop --------------------------------------------------------------

all_results <- list()

for (dgp_name in names(dgp_fns)) {
  cat("\n", strrep("=", 75), "\n", sep = "")
  cat("DGP: ", toupper(dgp_name), "  (n=", N, ", true ATT=", TRUE_ATT, ")\n", sep = "")
  cat(strrep("=", 75), "\n", sep = "")

  d <- dgp_fns[[dgp_name]](N)
  X <- d$X; A <- d$A; Y <- d$Y

  res_dgp <- list()

  # ---- Approach 1: Fullsample (no cross-fitting) ----------------------------
  # Implemented via estimators.R logic: CV-select lambda on all data, predict all.
  res_dgp[["1_fullsample"]] <- run_approach("1_fullsample", function() {
    cv_e <- optimaltrees::cv_regularization(X, A, loss_function = "log_loss",
                                            K = 5, refit = TRUE, verbose = FALSE)
    e_hat <- predict(cv_e$model, X, type = "prob")[, 2L]

    X0 <- X[A == 0, , drop = FALSE]; Y0 <- Y[A == 0]
    cv_m0 <- optimaltrees::cv_regularization(X0, Y0, loss_function = "log_loss",
                                             K = 5, refit = TRUE, verbose = FALSE)
    m0_hat <- predict(cv_m0$model, X, type = "prob")[, 2L]

    psi <- doubletree:::psi_att(Y, A, theta = 0,
                                list(e = e_hat, m0 = m0_hat, m1 = NULL), mean(A))
    theta <- sum(psi) / sum(A / mean(A))
    sigma <- att_se(doubletree:::psi_att(Y, A, theta, list(e=e_hat, m0=m0_hat, m1=NULL), mean(A)))
    list(theta = theta, sigma = sigma, ci_95 = att_ci(theta, sigma))
  })

  # ---- Approach 2: Crossfit (standard) -------------------------------------
  res_dgp[["2_crossfit"]] <- run_approach("2_crossfit", function() {
    r <- estimate_att(X, A, Y, K = K, outcome_type = "binary",
                      use_rashomon = FALSE, cv_regularization = TRUE, verbose = FALSE)
    list(theta = r$theta, sigma = r$sigma, ci_95 = att_ci(r$theta, r$sigma))
  })

  # ---- Approach 3: Doubletree (Rashomon intersection + cross-fit) -----------
  res_dgp[["3_doubletree"]] <- run_approach("3_doubletree", function() {
    r <- estimate_att(X, A, Y, K = K, outcome_type = "binary",
                      use_rashomon = TRUE, cv_regularization = TRUE,
                      auto_tune_intersecting = TRUE, verbose = FALSE)
    list(theta = r$theta, sigma = r$sigma, ci_95 = att_ci(r$theta, r$sigma))
  })

  # ---- Approach 4: Doubletree averaged -------------------------------------
  res_dgp[["4_dt_averaged"]] <- run_approach("4_dt_averaged", function() {
    r <- estimate_att_doubletree_averaged(X, A, Y, K = K,
                                          outcome_type = "binary",
                                          auto_tune_intersecting = TRUE,
                                          verbose = FALSE)
    list(theta = r$theta, sigma = r$sigma, ci_95 = r$ci_95)
  })

  # ---- Approach 5: M-split -------------------------------------------------
  res_dgp[["5_msplit"]] <- run_approach("5_msplit", function() {
    r <- estimate_att_msplit(X, A, Y, M = M_SMALL, K = K_SMALL,
                             outcome_type = "binary", verbose = FALSE)
    list(theta = r$theta, sigma = r$sigma, ci_95 = r$ci_95)
  })

  # ---- Approach 6: M-split averaged ----------------------------------------
  res_dgp[["6_msplit_averaged"]] <- run_approach("6_msplit_averaged", function() {
    r <- estimate_att_msplit_averaged(X, A, Y, M = M_SMALL, K = K_SMALL,
                                      outcome_type = "binary",
                                      seed_base = 42, verbose = FALSE)
    list(theta = r$theta, sigma = r$sigma, ci_95 = r$ci_95)
  })

  all_results[[dgp_name]] <- res_dgp
}

# ---- summary ----------------------------------------------------------------

cat("\n", strrep("=", 75), "\n", sep = "")
cat("SUMMARY  (CI width: approaches 4 & 6 should be ~0.15-0.22, NOT ~0.007)\n")
cat(strrep("=", 75), "\n", sep = "")
cat(sprintf("%-12s  %-22s  %8s  %8s  %s\n",
            "DGP", "Approach", "theta", "width", "covers"))
cat(strrep("-", 65), "\n", sep = "")

for (dgp_name in names(all_results)) {
  for (appr in names(all_results[[dgp_name]])) {
    r <- all_results[[dgp_name]][[appr]]
    if (!is.null(r)) {
      cat(sprintf("%-12s  %-22s  %+8.3f  %8.4f  %s\n",
                  dgp_name, appr,
                  if (is.finite(r$theta)) r$theta else NaN,
                  if (!is.null(r$width) && !is.na(r$width)) r$width else NaN,
                  if (isTRUE(r$covers)) "YES" else if (isFALSE(r$covers)) "NO" else "???"))
    } else {
      cat(sprintf("%-12s  %-22s  %8s  %8s  %s\n", dgp_name, appr, "ERROR", "ERROR", "---"))
    }
  }
}
cat("\nDone.\n")
