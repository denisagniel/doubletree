## ============================================================================
## 07b_compare_ols.R
##
## Compare the tree estimators with OLS/logit on the SAME winsorized outcome and
## the SAME folds as 07's cross-fit anchor. Reads 07's outputs, so run 07 first.
##
## Written to output/export/ (shareable aggregates; reads from output/private/):
##   ols_att_estimates.parquet     ATT for each (e learner x m0 learner) pair
##       full_sample  e: flagship tree | logit     m0: flagship tree | OLS main | OLS +pairwise
##       cross_fit    e: tree_cv | anchor | logit   m0: tree_cv | anchor | OLS main | OLS +pairwise
##   ols_nuisance_metrics.parquet  out-of-fold m0 MSE/MAE (controls) and
##                                 e Brier/log-loss, with fold-wise SEs
##   ols_m0_paired_cv.parquet      held-out outcome model, tree vs OLS: paired
##                                 MSE differences and skill vs a constant
##   ols_m0_fold_mse.parquet       the per-fold MSEs behind it
##   ols_complexity.csv            parameter counts (the interpretability side)
##
## `tree_cv` is the flagship pipeline (estimate_att() with the flagship's
## leaf_budget, depth, m_n and lambda_n) refit on each fold's TRAINING rows and
## evaluated on the held-out rows, so its out-of-fold numbers are comparable
## with the other learners'. Repetition 1 of the CV reuses the anchor's folds;
## later repetitions draw fresh stratified folds (config_oof_repeats).
##
## How to read the comparison. Real data, no ground truth: a narrower CI is not
## "beating" OLS. The defensible statements are about out-of-fold predictive
## risk and parameter count. Claims about ATT bias or coverage belong to
## simulations with a known truth. The `full_sample` ATT rows use in-sample
## nuisances (SEs downward-biased, O(p/n)); the `cross_fit` rows do not.
##
## Run time: dominated by the tree refits, config_oof_repeats * K flagship fits
## (each ~7 minutes at full size, so roughly 1.5-2 hours at the defaults). The
## fast parametric fits run FIRST so a separation error surfaces immediately.
##
## SERVER-ONLY: reads patient-level files. The outputs are aggregates.
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))
source(file.path(app_dir, "helpers", "covariate_blocks.R"))
source(file.path(app_dir, "helpers", "design_matrix.R"))
source(file.path(app_dir, "helpers", "ols_comparison.R"))
if (config_echo_decisions) echo_open_decisions()

cli::cli_h1("07b -- tree vs OLS/logit comparison")

## internals this script depends on (no exported equivalents yet)
stopifnot(
  exists("eif_att_solve", envir = asNamespace("doubletree")),
  exists("predict_nuisances_fold", envir = asNamespace("doubletree"))
)

analytic_path <- output_file("private", "analytic_cohort.parquet", app_dir = app_dir)
fits_path <- output_file("private", "att_fits.rds", app_dir = app_dir)
manifest_path <- output_file("private", "analytic_cohort_manifest.csv", app_dir = app_dir)
prov_path <- output_file("export", "att_provenance.csv", app_dir = app_dir)
out_dir <- dirname(output_file("export", "x", app_dir = app_dir))
for (p in c(analytic_path, fits_path, manifest_path, prov_path)) {
  if (!file.exists(p)) cli::cli_abort("Missing {.file {p}}; run 06 and 07 first.")
}
manifest <- utils::read.csv(manifest_path)
built_under <- manifest$value[manifest$field == "config_fingerprint"]
if (length(built_under) != 1L || !identical(built_under, config_fingerprint())) {
  cli::cli_abort("{.file {analytic_path}} was built under a different config; re-run 06 and 07.")
}
## att_fits.rds must come from THIS parquet and THIS config, or the nuisances
## and folds below would be mixed with a different Y.
prov <- readr::read_csv(prov_path, show_col_types = FALSE)
pv <- stats::setNames(as.character(prov$value), prov$field)
if (!identical(unname(pv[["analytic_md5"]]), unname(tools::md5sum(analytic_path))) ||
    !identical(unname(pv[["crossfit_learner"]]), config_crossfit_learner) ||
    !identical(unname(pv[["seed"]]), as.character(config_seed))) {
  cli::cli_abort(c(
    "{.file {fits_path}} is from a different 07 run than the current data/config.",
    "i" = "Re-run {.file 07_estimate_att.R} (data hash, crossfit learner or seed differ)."
  ))
}

analytic_data <- arrow::read_parquet(analytic_path)
X <- as.data.frame(analytic_data[, design_matrix_columns(), drop = FALSE])
A <- as.integer(analytic_data$A)
Y <- as.numeric(analytic_data$Y)
assert_binary_design_matrix(X)
n <- nrow(X)

fits <- readr::read_rds(fits_path)
fit_att <- fits$fit_att
fit_cf <- fits$fit_cf
stopifnot(identical(fit_cf$n, n), identical(fit_att$n, n))
fold_indices <- fit_cf$fold_indices
K <- length(unique(fold_indices))
clip_c <- fit_att$clip_c

## ---- folds for each CV repetition -------------------------------------------
cv_folds <- lapply(seq_len(config_oof_repeats), function(r) {
  if (r == 1L) fold_indices else doubletree::create_folds(n, K, strata = A, seed = config_seed + r)
})

## ---- parametric nuisances FIRST (fast; separation errors surface now) -------
par_reps <- lapply(cv_folds, function(f) crossfit_parametric(X, A, Y, f))
par_cf <- par_reps[[1]]
par_full <- full_sample_parametric(X, A, Y)

## ---- tree refit per fold (the slow part) ------------------------------------
tree_oof <- function(folds) {
  e <- m0 <- numeric(n)
  for (k in sort(unique(folds))) {
    te <- which(folds == k); tr <- which(folds != k)
    cli::cli_alert_info("tree_cv fold {k}: flagship refit on {length(tr)} training rows")
    f <- doubletree::estimate_att(
      X = X[tr, , drop = FALSE], A = A[tr], Y = Y[tr],
      leaf_budget = config_leaf_budget, outcome_type = "continuous",
      max_depth = config_max_depth, depth_restricted = TRUE, model_limit = 0L,
      m_n = fit_att$m_n, lambda_n = fit_att$lambda_n   # match the flagship exactly
    )
    models <- list(e_model = f$nuisance_fits$e_model, m0_model = f$nuisance_fits$m0_model,
                   outcome_type = "continuous")
    pr <- doubletree:::predict_nuisances_fold(models, X, te)
    e[te] <- pr$e
    m0[te] <- pr$m0
  }
  list(e = e, m0 = m0)
}
cv_reps <- vector("list", config_oof_repeats)
for (r in seq_len(config_oof_repeats)) {
  cli::cli_alert_info("outcome-model CV repetition {r}/{config_oof_repeats}")
  cv_reps[[r]] <- list(folds = cv_folds[[r]], tree = tree_oof(cv_folds[[r]]), par = par_reps[[r]])
}
tree_cv <- cv_reps[[1]]$tree

## ---- paired held-out MSE: tree vs OLS ---------------------------------------
const_oof <- function(folds) {
  out <- numeric(n)
  for (k in sort(unique(folds))) {
    out[folds == k] <- mean(Y[folds != k & A == 0L])
  }
  out
}
fm <- do.call(rbind, lapply(seq_len(config_oof_repeats), function(r) {
  cr <- cv_reps[[r]]
  fold_mse_m0(A, Y, cr$folds, list(tree_cv = cr$tree$m0, ols_main = cr$par$m0_main,
                                   ols_pairwise = cr$par$m0_inter,
                                   constant = const_oof(cr$folds)), rep_id = r)
}))
paired <- paired_mse_summary(fm, reference = "tree_cv", baseline = "constant")
cli::cli_h2("Held-out outcome model: flagship tree vs OLS (controls, winsorized Y)")
print(paired, digits = 4, row.names = FALSE)

## ---- ATT grids ----------------------------------------------------------------
e_cf <- list(tree_cv = tree_cv$e, anchor = fit_cf$nuisance_fits$propensity, logit = par_cf$e)
m0_cf <- list(tree_cv = tree_cv$m0, anchor = fit_cf$nuisance_fits$outcome_control,
              ols_main = par_cf$m0_main, ols_pairwise = par_cf$m0_inter)
e_full <- list(tree_flagship = fit_att$nuisance_fits$propensity, logit = par_full$e)
m0_full <- list(tree_flagship = fit_att$nuisance_fits$outcome_control,
                ols_main = par_full$m0_main, ols_pairwise = par_full$m0_inter)
n_par_m0 <- c(tree_flagship = fit_att$n_leaves_m0, tree_cv = NA_integer_, anchor = NA_integer_,
              ols_main = ncol(X) + 1L, ols_pairwise = ncol(X) + choose(ncol(X), 2L) + 1L)

grid_rows <- function(e_l, m0_l, scheme, se_method) {
  do.call(rbind, lapply(names(e_l), function(en) {
    do.call(rbind, lapply(names(m0_l), function(mn) {
      cbind(scheme = scheme, e_learner = en, m0_learner = mn,
            att_from_nuisances(Y, A, e_l[[en]], m0_l[[mn]], clip_c),
            se_method = se_method, n_parameters_m0 = unname(n_par_m0[mn]))
    }))
  }))
}
att_grid <- rbind(
  grid_rows(e_full, m0_full, "full_sample",
            "in-sample EIF plug-in; downward-biased, O(p/n)"),
  grid_rows(e_cf, m0_cf, "cross_fit", "cross-fit EIF")
)
att_grid$se_method <- paste0(att_grid$se_method,
                             "; treats the estimated winsorization cap as fixed")
att_grid$n <- n
att_grid$n_treated <- sum(A == 1L)
att_grid$n_control <- sum(A == 0L)
print(att_grid, digits = 4, row.names = FALSE)

## ---- out-of-fold metrics (repetition 1) ----------------------------------------
metrics <- nuisance_metrics(A, Y, fold_indices, e_cf, m0_cf, clip_c)
metrics$scheme <- "cross_fit_out_of_fold"
print(metrics, digits = 4, row.names = FALSE)

## ---- complexity (the interpretability side of the comparison) ----------------
complexity <- data.frame(
  learner = c("tree_flagship_e", "tree_flagship_m0", "ols_main", "ols_pairwise", "logit_main"),
  n_parameters = c(fit_att$n_leaves_e, fit_att$n_leaves_m0,
                   ncol(X) + 1L, ncol(X) + choose(ncol(X), 2L) + 1L, ncol(X) + 1L),
  n_aliased_full_sample = c(NA, NA, par_full$n_aliased[["main"]],
                            par_full$n_aliased[["inter"]], NA)
)
print(complexity, row.names = FALSE)

arrow::write_parquet(att_grid, file.path(out_dir, "ols_att_estimates.parquet"))
arrow::write_parquet(metrics, file.path(out_dir, "ols_nuisance_metrics.parquet"))
arrow::write_parquet(paired, file.path(out_dir, "ols_m0_paired_cv.parquet"))
arrow::write_parquet(fm, file.path(out_dir, "ols_m0_fold_mse.parquet"))
utils::write.csv(complexity, file.path(out_dir, "ols_complexity.csv"), row.names = FALSE)
cli::cli_alert_success("Wrote the 07b result files to {.file {out_dir}}.")

if (!interactive() && sys.nframe() == 0L) {
  cli::cli_alert_info("07b_compare_ols.R complete.")
}
