# ============================================================================
# Validate the full-sample two-tree ATT estimator, estimate_att()
#   (doubletree/inst/paper/theory.tex Part II, sec:main / thm:main)
#
# Purpose: check empirically that estimate_att() behaves as Part II predicts,
#   and that it fails as Part II predicts when its identifying assumption
#   (ass:sparsity) is violated. Compared against estimate_att_crossfit()
#   (Part I fallback) on the SAME simulated datasets.
#
# Two studies:
#   A. sparsity-SATISFYING DGP  -- ass:sparsity holds (delta_e = delta_mu = 0)
#      => thm:main applies; expect agreement with cross-fit, full-sample mean SE
#         not exceeding cross-fit's, coverage near nominal 0.95.
#   B. sparsity-VIOLATING DGP   -- delta_e, delta_mu > 0 at any small budget
#      => thm:main does NOT apply; degraded coverage is the PREDICTED outcome,
#         not a bug.
#
# Estimator arms (all on identical data within a replication):
#   fullsample   estimate_att(leaf_budget = 4, lambda_n = log(n)/n)
#   crossfit_cv  estimate_att_crossfit(K = 5)                 [package defaults]
#   crossfit_fix estimate_att_crossfit(K = 5, cv_regularization = FALSE,
#                                      regularization = log(n)/n)
#     -- the lambda-MATCHED cross-fit arm. Without it, fullsample-vs-crossfit_cv
#        confounds two differences at once (sample splitting AND fixed-vs-CV
#        penalty); crossfit_fix isolates the splitting effect alone.
#   DGP B only, fullsample_L16: estimate_att(leaf_budget = 16) -- probes whether
#     B's failure is the BUDGET or the additive structure itself.
#   oracle       the TRUE e(X), m0(X) pushed through the SAME eif_att_solve().
#     -- the implementation-correctness gate. If `oracle` is unbiased and covers
#        at ~0.95 on a DGP while `fullsample` does not, the shared EIF/variance
#        plumbing is sound and the gap is nuisance approximation error (which is
#        what theory predicts under a sparsity violation), NOT a coding bug.
#        Costs nothing: no tree is fitted.
#
# Inputs: none.
# Outputs: dev-scripts/output/2026-08-21_validate_fullsample/
#            results_<DGP>.rds  (per-replication rows)
#            summary_<DGP>.csv  (the reported table)
#
# Usage (defaults run BOTH studies serially, ~40 min):
#   Rscript dev-scripts/2026-08-21_validate_fullsample_estimate_att.R
#   Rscript dev-scripts/2026-08-21_validate_fullsample_estimate_att.R --dgp A
#   Rscript dev-scripts/2026-08-21_validate_fullsample_estimate_att.R --dgp B --reps 50
# The two studies are independent, so running --dgp A and --dgp B as two
# concurrent processes halves wall time (~20 min); that is how the reported
# n = 400 numbers were produced.
#
# Supplementary large-n run (disentangles "sparsity violated" from "signal below
# lambda_n at n = 400"; at n = 2000 the trees do split, so the residual bias is
# genuine approximation error delta > 0 rather than a stump):
#   Rscript dev-scripts/2026-08-21_validate_fullsample_estimate_att.R \
#     --dgp B --n 2000 --reps 200 --skip-cv
# ============================================================================

# ---- 0. Setup --------------------------------------------------------------

suppressPackageStartupMessages({
  library(optimaltrees)
  devtools::load_all(".", quiet = TRUE)
})

args <- commandArgs(trailingOnly = TRUE)
opt_arg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}

WHICH_DGP <- opt_arg("--dgp", "both")          # "A", "B", or "both"
REPS      <- as.integer(opt_arg("--reps", "200"))
N_OBS     <- as.integer(opt_arg("--n", "400"))
# --skip-cv drops the crossfit_cv arm. Needed only for the large-n supplementary
# run: estimate_att_crossfit() at its CV defaults costs ~10 s/replication at
# n = 2000 (vs 1.6 s for estimate_att and 0.15 s for the fixed-lambda arm), i.e.
# ~34 min for that arm alone. Never used for the headline n = 400 studies.
SKIP_CV   <- "--skip-cv" %in% args
LEAF_BUDGET <- 4L        # Lbar: FIXED analyst choice (theory.tex ass:budget)
K_FOLDS   <- 5L
TRUE_ATT  <- 0.15        # constant additive effect => ATT known exactly
MASTER_SEED <- 20260821L

OUT_DIR <- file.path("dev-scripts", "output", "2026-08-21_validate_fullsample")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

stopifnot(WHICH_DGP %in% c("A", "B", "both"), REPS >= 1L, N_OBS >= 50L)

# One seed governs the whole study. Per-replication seeds are DRAWN from it and
# re-set before each replication's data generation, because estimate_att_crossfit
# (seed =) calls set.seed() internally and would otherwise perturb the stream
# that generates subsequent replications' data.
set.seed(MASTER_SEED)
REP_SEEDS <- sample.int(.Machine$integer.max, REPS)

# ---- 1. Data / DGP ---------------------------------------------------------

#' Sparsity-SATISFYING DGP: hierarchical structure, e on X1 only, m0 on X2 only.
#' Each nuisance takes 2 distinct values, so a 2-leaf tree represents it
#' EXACTLY: delta_e = delta_mu = 0 at any leaf_budget >= 2 (ass:sparsity holds).
#' theory.tex ass:sparsity discussion: sparsity is CHEAP for region-varying
#' active coordinate sets.
generate_dgp_sparse <- function(n) {
  X <- data.frame(
    X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5),
    X3 = rbinom(n, 1, 0.5), X4 = rbinom(n, 1, 0.5)   # inert noise covariates
  )
  e_true  <- 0.35 + 0.30 * X$X1
  m0_true <- 0.25 + 0.30 * X$X2
  A <- rbinom(n, 1, e_true)
  Y <- rbinom(n, 1, m0_true + TRUE_ATT * A)
  list(X = X, A = A, Y = Y, e_true = e_true, m0_true = m0_true)
}

#' Sparsity-VIOLATING DGP: both nuisances ADDITIVE in 6 binary covariates.
#' theory.tex ass:sparsity discussion: "The assumption is expensive precisely
#' for additive structure, where a leaf is required for each configuration of
#' the active coordinates." Level sets of a sum are not axis-aligned rectangles,
#' so exact representation needs 2^6 = 64 leaves -- and at lambda_n = log(n)/n
#' the penalised selector stops well short of that, leaving delta > 0 at ANY
#' budget (verified: RMSE(e_hat, e) plateaus at ~0.098 from budget 8 to 64 at
#' n = 4000, and does not shrink).
generate_dgp_additive <- function(n) {
  X <- data.frame(
    X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5), X3 = rbinom(n, 1, 0.5),
    X4 = rbinom(n, 1, 0.5), X5 = rbinom(n, 1, 0.5), X6 = rbinom(n, 1, 0.5)
  )
  S <- rowSums(X)
  e_true  <- 0.20 + 0.10 * S          # 0.20 .. 0.80, interior (ass:causal ok)
  m0_true <- 0.15 + 0.10 * S          # 0.15 .. 0.75
  A <- rbinom(n, 1, e_true)
  Y <- rbinom(n, 1, m0_true + TRUE_ATT * A)
  list(X = X, A = A, Y = Y, e_true = e_true, m0_true = m0_true)
}

# ---- 2. Sparsity verification (run BEFORE the Monte Carlo) -----------------

#' Verify a DGP's sparsity status directly, by oracle approximation error.
#'
#' certified_*/n_leaves_* alone cannot settle this: `certified` means the fit
#' provably solves eq:select, NOT that ass:sparsity holds. The decisive evidence
#' is whether the fitted nuisance's error against the KNOWN truth shrinks toward
#' 0 as the leaf budget grows (delta = 0) or plateaus at a positive value
#' (delta > 0). Run at large n so estimation noise is small relative to delta.
verify_sparsity <- function(gen, label, budgets, n_big = 4000L, seed = 20260821L) {
  set.seed(seed)
  d <- gen(n_big)
  rows <- lapply(budgets, function(L) {
    fit <- estimate_att(d$X, d$A, d$Y, leaf_budget = L, outcome_type = "binary")
    data.frame(
      budget      = L,
      n_leaves_e  = fit$n_leaves_e,
      n_leaves_m0 = fit$n_leaves_m0,
      certified_e = fit$certified_e,
      certified_m0 = fit$certified_m0,
      used_search_e = fit$used_search_e,
      rmse_e      = sqrt(mean((fit$nuisance_fits$propensity - d$e_true)^2)),
      rmse_m0     = sqrt(mean((fit$nuisance_fits$outcome_control - d$m0_true)^2))
    )
  })
  tab <- do.call(rbind, rows)
  cat(sprintf("\n--- sparsity verification: %s (n = %d) ---\n", label, n_big))
  print(tab, row.names = FALSE, digits = 4)
  cat(sprintf("  oracle RMSE floor  e: %.4f   m0: %.4f   (sampling noise ~ %.4f)\n",
              min(tab$rmse_e), min(tab$rmse_m0), 1 / sqrt(n_big)))
  flush.console()
  tab
}

# ---- 3. Estimation arms ----------------------------------------------------

#' Run every estimator arm on one dataset. Returns a one-row data.frame.
#' A failing arm records its CONDITION MESSAGE (never a silent NA substitution)
#' so the run can continue across 200 replications without discarding the reason.
run_arms <- function(d, rep_id, fold_seed, extra_budget = NULL) {
  n <- nrow(d$X)
  lam <- log(n) / n
  errs <- character(0)

  grab <- function(expr, arm) {
    withCallingHandlers(
      tryCatch(expr, error = function(e) {
        errs[[arm]] <<- conditionMessage(e); NULL
      }),
      warning = function(w) invokeRestart("muffleWarning")
    )
  }

  fs <- grab(estimate_att(d$X, d$A, d$Y, leaf_budget = LEAF_BUDGET,
                          outcome_type = "binary", lambda_n = lam), "fullsample")
  cf_cv <- if (SKIP_CV) NULL else
    grab(estimate_att_crossfit(d$X, d$A, d$Y, K = K_FOLDS,
                               seed = fold_seed, outcome_type = "binary"),
         "crossfit_cv")
  cf_fix <- grab(estimate_att_crossfit(d$X, d$A, d$Y, K = K_FOLDS,
                                       seed = fold_seed, outcome_type = "binary",
                                       cv_regularization = FALSE,
                                       regularization = lam), "crossfit_fix")
  fs_big <- if (is.null(extra_budget)) NULL else
    grab(estimate_att(d$X, d$A, d$Y, leaf_budget = extra_budget,
                      outcome_type = "binary", lambda_n = lam), "fullsample_big")

  # Implementation-correctness gate: the TRUE nuisances through the SAME solver
  # estimate_att() uses. Clamped exactly as estimate_att() clamps, so the only
  # difference from the `fullsample` arm is the nuisance values themselves.
  orc <- grab(doubletree:::eif_att_solve(
    d$Y, d$A,
    pmax(doubletree:::.PROPENSITY_LOWER_BOUND,
         pmin(doubletree:::.PROPENSITY_UPPER_BOUND, d$e_true)),
    d$m0_true, n), "oracle")

  covered <- function(f) if (is.null(f)) NA else
    as.numeric(TRUE_ATT >= f$ci_95[1] && TRUE_ATT <= f$ci_95[2])
  pull <- function(f, fld) if (is.null(f)) NA_real_ else as.numeric(f[[fld]])

  out <- data.frame(
    rep = rep_id,
    theta_fullsample   = pull(fs, "theta"),
    sigma_fullsample   = pull(fs, "sigma"),
    cover_fullsample   = covered(fs),
    theta_crossfit_cv  = pull(cf_cv, "theta"),
    sigma_crossfit_cv  = pull(cf_cv, "sigma"),
    cover_crossfit_cv  = covered(cf_cv),
    theta_crossfit_fix = pull(cf_fix, "theta"),
    sigma_crossfit_fix = pull(cf_fix, "sigma"),
    cover_crossfit_fix = covered(cf_fix),
    theta_oracle = pull(orc, "theta"),
    sigma_oracle = pull(orc, "sigma"),
    cover_oracle = covered(orc),
    n_leaves_e  = pull(fs, "n_leaves_e"),
    n_leaves_m0 = pull(fs, "n_leaves_m0"),
    # Spread of the FITTED nuisances. The decisive diagnostic for DGP B: when the
    # per-split risk gain falls below lambda_n, eq:select correctly returns a
    # STUMP, sd(e_hat) collapses to 0, and the doubly-robust correction has
    # nothing to correct with -- the estimate slides toward the unadjusted
    # difference. Compare against sd(e_true) for the DGP.
    sd_e_hat  = if (is.null(fs)) NA_real_ else
      stats::sd(fs$nuisance_fits$propensity),
    sd_m0_hat = if (is.null(fs)) NA_real_ else
      stats::sd(fs$nuisance_fits$outcome_control),
    sd_e_true = stats::sd(d$e_true),
    certified_e   = if (is.null(fs)) NA else fs$certified_e,
    certified_m0  = if (is.null(fs)) NA else fs$certified_m0,
    used_search_e = if (is.null(fs)) NA else fs$used_search_e,
    used_search_m0 = if (is.null(fs)) NA else fs$used_search_m0,
    n_errors = length(errs),
    err_msg = if (length(errs)) paste(names(errs), errs, sep = ": ",
                                      collapse = " | ") else NA_character_,
    stringsAsFactors = FALSE
  )
  if (!is.null(extra_budget)) {
    out$theta_fullsample_big <- pull(fs_big, "theta")
    out$sigma_fullsample_big <- pull(fs_big, "sigma")
    out$cover_fullsample_big <- covered(fs_big)
    out$n_leaves_e_big  <- pull(fs_big, "n_leaves_e")
    out$n_leaves_m0_big <- pull(fs_big, "n_leaves_m0")
  }
  out
}

# ---- 4. Monte Carlo driver -------------------------------------------------

summarise_arm <- function(res, arm) {
  th <- res[[paste0("theta_", arm)]]
  sg <- res[[paste0("sigma_", arm)]]
  cv <- res[[paste0("cover_", arm)]]
  ok <- !is.na(th)
  data.frame(
    arm         = arm,
    n_ok        = sum(ok),
    mean_theta  = mean(th, na.rm = TRUE),
    bias        = mean(th, na.rm = TRUE) - TRUE_ATT,
    sd_theta    = stats::sd(th, na.rm = TRUE),
    mean_sigma  = mean(sg, na.rm = TRUE),
    se_ratio    = mean(sg, na.rm = TRUE) / stats::sd(th, na.rm = TRUE),
    coverage    = mean(cv, na.rm = TRUE),
    cov_mcse    = sqrt(mean(cv, na.rm = TRUE) * (1 - mean(cv, na.rm = TRUE)) /
                         sum(!is.na(cv))),
    stringsAsFactors = FALSE
  )
}

run_study <- function(gen, label, tag, extra_budget = NULL) {
  if (N_OBS != 400L) tag <- paste0(tag, "_n", N_OBS)
  cat(sprintf("\n================ STUDY %s: %s ================\n", tag, label))
  cat(sprintf("n = %d, reps = %d, leaf_budget = %d, K = %d, true ATT = %.3f%s\n",
              N_OBS, REPS, LEAF_BUDGET, K_FOLDS, TRUE_ATT,
              if (SKIP_CV) "  [crossfit_cv arm SKIPPED, see --skip-cv]" else ""))
  arms <- c("fullsample", "crossfit_cv", "crossfit_fix")
  if (SKIP_CV) arms <- setdiff(arms, "crossfit_cv")
  if (!is.null(extra_budget)) arms <- c(arms, "fullsample_big")
  arms <- c(arms, "oracle")

  t_start <- Sys.time()
  rows <- vector("list", REPS)
  for (r in seq_len(REPS)) {
    set.seed(REP_SEEDS[r])
    d <- gen(N_OBS)
    rows[[r]] <- run_arms(d, rep_id = r, fold_seed = REP_SEEDS[r] %% 100000L,
                          extra_budget = extra_budget)
    if (r %% 25L == 0L) {
      cat(sprintf("  ... rep %d/%d  (%.1f min elapsed)\n", r, REPS,
                  as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
      flush.console()
    }
  }
  res <- do.call(rbind, rows)

  summ <- do.call(rbind, lapply(arms, function(a) summarise_arm(res, a)))
  cat("\n-- results (true ATT = ", TRUE_ATT, ") --\n", sep = "")
  print(summ, row.names = FALSE, digits = 4)

  cat(sprintf("\n-- full-sample nuisance diagnostics (leaf_budget = %d) --\n",
              LEAF_BUDGET))
  cat(sprintf("  mean n_leaves: e = %.2f  m0 = %.2f   (modal: e = %s, m0 = %s)\n",
              mean(res$n_leaves_e, na.rm = TRUE), mean(res$n_leaves_m0, na.rm = TRUE),
              names(sort(table(res$n_leaves_e), decreasing = TRUE))[1],
              names(sort(table(res$n_leaves_m0), decreasing = TRUE))[1]))
  cat(sprintf("  certified: e = %.3f  m0 = %.3f | used_search: e = %.3f  m0 = %.3f\n",
              mean(res$certified_e, na.rm = TRUE), mean(res$certified_m0, na.rm = TRUE),
              mean(res$used_search_e, na.rm = TRUE), mean(res$used_search_m0, na.rm = TRUE)))
  cat(sprintf("  fitted spread: mean sd(e_hat) = %.4f vs sd(e_true) = %.4f | mean sd(m0_hat) = %.4f\n",
              mean(res$sd_e_hat, na.rm = TRUE), mean(res$sd_e_true, na.rm = TRUE),
              mean(res$sd_m0_hat, na.rm = TRUE)))
  cat(sprintf("  reps with CONSTANT e_hat (stump): %.3f | constant m0_hat: %.3f\n",
              mean(res$sd_e_hat < 1e-10, na.rm = TRUE),
              mean(res$sd_m0_hat < 1e-10, na.rm = TRUE)))
  if (!is.null(extra_budget)) {
    cat(sprintf("  budget %d arm mean n_leaves: e = %.2f  m0 = %.2f\n", extra_budget,
                mean(res$n_leaves_e_big, na.rm = TRUE),
                mean(res$n_leaves_m0_big, na.rm = TRUE)))
  }
  n_fail <- sum(res$n_errors > 0)
  cat(sprintf("  replications with >=1 arm failure: %d/%d\n", n_fail, REPS))
  if (n_fail > 0) print(utils::head(stats::na.omit(res$err_msg), 5))
  cat(sprintf("  wall time: %.1f min\n",
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

  # ---- 5. Export ----
  saveRDS(list(per_rep = res, summary = summ, label = label,
               n = N_OBS, reps = REPS, leaf_budget = LEAF_BUDGET,
               true_att = TRUE_ATT, master_seed = MASTER_SEED),
          file.path(OUT_DIR, paste0("results_", tag, ".rds")))
  utils::write.csv(summ, file.path(OUT_DIR, paste0("summary_", tag, ".csv")),
                   row.names = FALSE)
  invisible(list(res = res, summ = summ))
}

# ---- Run -------------------------------------------------------------------

if (WHICH_DGP %in% c("A", "both")) {
  verify_sparsity(generate_dgp_sparse, "DGP A (sparsity-SATISFYING)",
                  budgets = c(2L, 4L, 8L, 16L, 32L))
  run_study(generate_dgp_sparse, "DGP A -- hierarchical, ass:sparsity HOLDS", "A")
}

if (WHICH_DGP %in% c("B", "both")) {
  verify_sparsity(generate_dgp_additive, "DGP B (sparsity-VIOLATING)",
                  budgets = c(2L, 4L, 8L, 16L, 32L, 64L))
  run_study(generate_dgp_additive, "DGP B -- additive, ass:sparsity FAILS", "B",
            extra_budget = 16L)
}

cat("\nOutputs written to: ", normalizePath(OUT_DIR), "\n", sep = "")
