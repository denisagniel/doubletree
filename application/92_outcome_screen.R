## ============================================================================
## 92_outcome_screen.R
##
## SCREEN candidate alternative OUTCOMES for the unchanged exposure A (SAA
## achievement) on their NUISANCE behaviour only. It estimates no ATT and reads
## no treatment effect: nothing here is a result, it is a feasibility screen,
## in the same spirit as 91_leaf_budget_feasibility.R. STANDALONE: not part of
## the 01-07 cascade, not gated by assert_no_blocking_open().
##
## WHY. The cost outcome is heavy-tailed (MSE ~6e9), so out-of-fold comparisons
## of the tree against OLS are dominated by a few very large claims, and on the
## 15 binary covariates even xgboost beats additive OLS by ~0.2% of the MSE: there
## is almost no interaction signal for a tree to find. Before redesigning the
## outcome, this asks of each candidate: is it defined for (nearly) everyone, how
## common is it, and is there any gap between a flexible learner (xgboost) and an
## additive logit on the same covariates? If the gap is ~0, no tree can win on it.
##
## WINDOW. Every candidate is read from the measure period that starts EXACTLY at
## s + 12 months (T0 + 12mo), i.e. the 12 months [s+12mo, s+24mo) that the cost
## outcome uses. Rows are selected with smidata::smi_select_measure_row(anchor =
## "ceiling_month"), as 02 does for the exposure: no row at that start = no match.
##
## CANDIDATES (see the `estimand_status` column; THE SCREEN DOES NOT MAKE A
## CANDIDATE A VALID ATT TARGET):
##   hosp_any_fuh30  Y = 1{matched row with MSR_DEN >= 1} in msr_fuh30_v2, i.e. any
##                   MH discharge in the window (MSR_DEN counts discharges).
##                   Defined for the WHOLE analytic cohort ONLY IF a missing row
##                   means "no discharge" -- smidata's docs say treat as DEN = 0,
##                   but also that some CDM months are missing. The coverage table
##                   reports the share of the cohort with a matched row so that
##                   assumption can be judged, not assumed.
##   fuh30_followup, fuh07_followup
##                   Y = MSR_NUM / MSR_DEN >= threshold among patients WITH a
##                   discharge. This conditions on a post-exposure event (SAA
##                   adherence can change hospitalisation), so it is a
##                   selected-outcome screen for prediction only, NOT an ATT.
##   aap_access      Y = AAP achieved among rows with MSR_DEN >= 1. Unconditional
##                   only if that covers ~all of the cohort; the status column
##                   says which, from the coverage actually observed.
##
## The exposure model e does not depend on the outcome, so it is NOT refit here
## (07b's ols_nuisance_metrics.parquet already has it).
##
## LEARNERS for m0 = P(Y = 1 | X, A = 0), fit on training controls, scored on
## held-out controls (Brier and log-loss): constant, additive logit, logit with
## all pairwise terms, xgboost (doubletree's own early-stopped fit, same default
## parameters as estimate_att_crossfit_xgb()). No optimal tree is fit: one
## flagship fit takes ~7 minutes at this n, so it is left out of the screen on
## purpose. The screening statistic is
##   gap = skill(xgboost) - skill(additive logit),   skill = 1 - Brier / Brier(constant)
## averaged over config_oof_repeats repetitions of K-fold CV, with the SE of the
## paired per-fold Brier difference from repetition 1 (folds are not independent;
## it is a rough scale, not a test).
##
## Written to output/export/ (aggregates). Any row with a count below
## config_report_min_cell is suppressed WHOLE (every number in it is NA), which
## leaves nothing to difference. This is simpler than, and not audited by,
## smidata::smi_report(): look at the files before pulling them back.
##   outcome_screen_coverage.parquet  per table: rows read, distinct MSR values,
##                                    cohort share with a matched row / DEN >= 1
##   outcome_screen_summary.csv       per candidate: counts, prevalence by arm,
##                                    skills, gap and its SE, estimand status
##   outcome_screen_metrics.parquet   per candidate x learner: Brier, log-loss, skill
##
## SERVER-ONLY: reads patient-level files. Run time: seconds to minutes (xgboost
## and logits only).
##   Rscript application/92_outcome_screen.R
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))
source(file.path(app_dir, "helpers", "covariate_blocks.R"))
source(file.path(app_dir, "helpers", "design_matrix.R"))
source(file.path(app_dir, "helpers", "ols_comparison.R"))
source(file.path(app_dir, "helpers", "disclosure.R"))   # config_report_min_cell
if (config_echo_decisions) echo_open_decisions()

cli::cli_h1("92 -- outcome screen (NOT a reportable estimate)")

if (!config_has_smidata || identical(smidata::smi_env(), "local")) {
  cli::cli_abort(c(
    "This screen reads the real measure tables and must run on the server.",
    "i" = "Locally {.fn smidata::smi_read} would return a synthetic fixture and the screen would be meaningless."
  ))
}
stopifnot(requireNamespace("xgboost", quietly = TRUE))
## internals this script depends on (no exported equivalents), as in 07b
stopifnot(
  exists(".xgb_fit_one", envir = asNamespace("doubletree")),
  exists(".xgb_predict", envir = asNamespace("doubletree"))
)

## ---- the analytic cohort, from the same build 07 used --------------------------
analytic_path <- output_file("private", "analytic_cohort.parquet", app_dir = app_dir)
manifest_path <- output_file("private", "analytic_cohort_manifest.csv", app_dir = app_dir)
out_dir <- dirname(output_file("export", "x", app_dir = app_dir))
for (p in c(analytic_path, manifest_path)) {
  if (!file.exists(p)) cli::cli_abort("Missing {.file {p}}; run 06 first.")
}
manifest <- utils::read.csv(manifest_path)
built_under <- manifest$value[manifest$field == "config_fingerprint"]
if (length(built_under) != 1L || !identical(built_under, config_fingerprint())) {
  cli::cli_abort("{.file {analytic_path}} was built under a different config; re-run 06.")
}
analytic_data <- arrow::read_parquet(analytic_path)
X_all <- as.data.frame(analytic_data[, design_matrix_columns(), drop = FALSE])
assert_binary_design_matrix(X_all)
A_all <- as.integer(analytic_data$A)
ids <- analytic_data$ID
n_cohort <- nrow(X_all)
index_dt <- data.frame(ID = ids, INDEX_DT = as.Date(analytic_data$T0))
cli::cli_alert_info("Analytic cohort: n = {n_cohort} ({sum(A_all == 1L)} treated, {sum(A_all == 0L)} control).")

## ---- candidate measure tables ---------------------------------------------------
## Dataset keys are not in SMI_KEYS (this screen is outside the pipeline), so they
## are named here. Reads still go through smi_read_pinned(), never smi_read().
screen_tables <- c(fuh30 = "msr_fuh30_v2", fuh07 = "msr_fuh07_v2", aap = "msr_aap")

read_window_row <- function(key) {
  raw <- smi_read_pinned(key, columns = SMI_COLS$saa)
  raw <- raw[raw$ID %in% ids, , drop = FALSE]
  msr_values <- sort(unique(as.character(raw$MSR)))
  ## One measure per file is a premise, not a fact (see OPEN_DECISIONS$msr_code):
  ## smi_select_measure_row() groups on ID, so a second MSR code could win silently.
  if (length(msr_values) > 1L) {
    cli::cli_abort(c(
      "{.val {key}} has {length(msr_values)} distinct {.field MSR} values within the cohort: {.val {msr_values}}.",
      "i" = "Filter to the intended code before screening it (see helpers/msr_filter.R)."
    ))
  }
  sel <- smidata::smi_select_measure_row(
    raw, index_dt = index_dt, target_offset_months = 12L,
    achievement_threshold = config_saa_achievement_threshold,
    anchor = "ceiling_month"
  )
  stopifnot(identical(sel$ID, ids))
  attr(sel, "msr_values") <- paste(msr_values, collapse = ";")
  attr(sel, "rows_read") <- nrow(raw)
  sel
}

windows <- lapply(screen_tables, function(k) {
  cli::cli_alert_info("Reading {.val {k}}.")
  read_window_row(k)
})

coverage <- do.call(rbind, lapply(names(windows), function(nm) {
  w <- windows[[nm]]
  data.frame(
    table = screen_tables[[nm]],
    rows_in_cohort = attr(w, "rows_read"),
    msr_values = attr(w, "msr_values"),
    n_cohort = n_cohort,
    n_matched = sum(w$matched),
    n_den_ge1 = sum(w$elig),
    share_matched = mean(w$matched),
    share_den_ge1 = mean(w$elig)
  )
}))
print(coverage, row.names = FALSE)

## ---- candidates: outcome and analysis rows ---------------------------------------
candidate <- function(name, win, kind) {
  y <- switch(kind,
    any_discharge = as.integer(win$matched & win$elig),
    achieved      = as.integer(win$achieved)
  )
  keep <- switch(kind, any_discharge = rep(TRUE, n_cohort), achieved = win$elig %in% TRUE)
  list(name = name, y = y, keep = keep, kind = kind, share_in = mean(keep))
}
candidates <- list(
  candidate("hosp_any_fuh30", windows$fuh30, "any_discharge"),
  candidate("fuh30_followup", windows$fuh30, "achieved"),
  candidate("fuh07_followup", windows$fuh07, "achieved"),
  candidate("aap_access",     windows$aap,   "achieved")
)

status_of <- function(cd) {
  if (cd$name == "hosp_any_fuh30") {
    "unconditional IF a missing row means no discharge (see coverage); not checkable from this table alone"
  } else if (grepl("^fuh", cd$name)) {
    "SELECTED on a post-exposure event (MH discharge); prediction screen only, NOT an ATT target"
  } else if (cd$share_in >= 0.99) {
    "near-unconditional (DEN >= 1 for >= 99% of cohort)"
  } else {
    sprintf("SELECTED on measure eligibility (DEN >= 1 for %.0f%% of cohort); not an ATT target as is",
            100 * cd$share_in)
  }
}

## ---- m0 learners, scored out of fold on controls ----------------------------------
brier <- function(y, p) mean((y - p)^2)
logloss <- function(y, p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}
xgb_params <- list(max_depth = 3L, eta = 0.1, subsample = 0.8, colsample_bytree = 1,
                   min_child_weight = 1, lambda = 1)

predict_m0 <- function(X, A, Y, folds, seed) {
  n <- nrow(X)
  M1 <- ols_design(X, FALSE)
  M2 <- ols_design(X, TRUE)
  Xm <- as.matrix(X); storage.mode(Xm) <- "double"
  out <- list(constant = numeric(n), logit_main = rep(NA_real_, n),
              logit_pairwise = rep(NA_real_, n), xgboost = numeric(n))
  main_ok <- TRUE
  pairwise_ok <- TRUE
  for (k in sort(unique(folds))) {
    te <- which(folds == k)
    trc <- which(folds != k & A == 0L)
    out$constant[te] <- mean(Y[trc])
    ## A logit can (quasi-)separate on a rare outcome (a covariate cell with no
    ## events); logit_fit_predict() errors rather than return garbage. That is a
    ## finding, not a reason to drop the candidate: report the learner as
    ## unavailable (NA in every fold, so averages stay comparable) and go on.
    if (main_ok) {
      mn <- tryCatch(logit_fit_predict(M1, Y, trc, te), error = function(e) {
        cli::cli_alert_warning("additive logit unavailable: {conditionMessage(e)}")
        NULL
      })
      if (is.null(mn)) main_ok <- FALSE else out$logit_main[te] <- mn
    }
    if (pairwise_ok) {
      pw <- tryCatch(logit_fit_predict(M2, Y, trc, te), error = function(e) NULL)
      if (is.null(pw)) pairwise_ok <- FALSE else out$logit_pairwise[te] <- pw
    }
    fx <- doubletree:::.xgb_fit_one(
      Xm[trc, , drop = FALSE], Y[trc], objective = "binary:logistic",
      params = xgb_params, nrounds_max = 500L, early_stopping_rounds = 20L,
      valid_frac = 0.2, nthread = 1L, seed = seed + 1000L + k, strata = NULL
    )
    out$xgboost[te] <- doubletree:::.xgb_predict(fx, Xm[te, , drop = FALSE])
  }
  if (!main_ok) out$logit_main[] <- NA_real_
  if (!pairwise_ok) out$logit_pairwise[] <- NA_real_
  out
}

screen_one <- function(cd) {
  keep <- cd$keep
  X <- X_all[keep, , drop = FALSE]; A <- A_all[keep]; Y <- cd$y[keep]
  n <- nrow(X)
  ctrl <- A == 0L
  cnt <- c(n = n, n_treated = sum(A == 1L), n_control = sum(ctrl),
           events_treated = sum(Y[A == 1L]), events_control = sum(Y[ctrl]),
           nonevents_treated = sum(1L - Y[A == 1L]), nonevents_control = sum(1L - Y[ctrl]))
  cli::cli_h2("{cd$name}: n = {n}, control prevalence {round(cnt[['events_control']] / cnt[['n_control']], 4)}")
  if (min(cnt[c("events_control", "nonevents_control")]) < 50L) {
    cli::cli_alert_warning("Fewer than 50 control events or non-events: skipping the learners.")
    return(list(counts = cnt, metrics = NULL, gap = c(gap = NA_real_, se = NA_real_)))
  }
  per_fold <- list()
  for (r in seq_len(config_oof_repeats)) {
    folds <- doubletree::create_folds(n, 5L, strata = A, seed = config_seed + r)
    pr <- predict_m0(X, A, Y, folds, seed = config_seed + 100L * r)
    for (k in sort(unique(folds))) {
      te <- which(folds == k & ctrl)
      for (lr in names(pr)) {
        if (anyNA(pr[[lr]][te])) next
        per_fold[[length(per_fold) + 1L]] <- data.frame(
          rep = r, fold = k, learner = lr,
          brier = brier(Y[te], pr[[lr]][te]), logloss = logloss(Y[te], pr[[lr]][te])
        )
      }
    }
  }
  pf <- do.call(rbind, per_fold)
  agg <- stats::aggregate(cbind(brier, logloss) ~ learner, data = pf, FUN = mean)
  base <- agg$brier[agg$learner == "constant"]
  agg$skill <- 1 - agg$brier / base
  ## paired per-fold Brier difference, additive logit minus xgboost (positive:
  ## xgboost is better). SE from repetition 1's K fold values.
  w <- stats::reshape(pf[, c("rep", "fold", "learner", "brier")],
                      idvar = c("rep", "fold"), timevar = "learner", direction = "wide")
  print(agg, digits = 4, row.names = FALSE)
  if (!all(c("brier.logit_main", "brier.xgboost") %in% names(w))) {
    cli::cli_alert_warning("Additive logit unavailable (separation): no gap for this candidate.")
    return(list(counts = cnt, metrics = agg, gap = c(gap = NA_real_, se = NA_real_)))
  }
  d <- w$brier.logit_main - w$brier.xgboost
  d1 <- d[w$rep == 1L]
  gap_skill <- agg$skill[agg$learner == "xgboost"] - agg$skill[agg$learner == "logit_main"]
  gap_se <- stats::sd(d1) / sqrt(length(d1)) / base
  cli::cli_alert_info("gap = skill(xgboost) - skill(additive logit) = {signif(gap_skill, 3)} (rough SE {signif(gap_se, 2)}).")
  list(counts = cnt, metrics = agg, gap = c(gap = gap_skill, se = gap_se))
}

results <- lapply(candidates, screen_one)
names(results) <- vapply(candidates, `[[`, "", "name")

## ---- assemble outputs, suppressing WHOLE rows with any small count -----------------
count_cols <- c("n", "n_treated", "n_control", "events_treated", "events_control",
                "nonevents_treated", "nonevents_control")
summary_tbl <- do.call(rbind, lapply(candidates, function(cd) {
  r <- results[[cd$name]]
  m <- r$metrics
  get <- function(lr, col) if (is.null(m) || !lr %in% m$learner) NA_real_ else m[[col]][m$learner == lr]
  data.frame(
    candidate = cd$name,
    estimand_status = status_of(cd),
    share_of_cohort_analysed = cd$share_in,
    t(as.data.frame(r$counts)),
    prevalence_treated = r$counts[["events_treated"]] / r$counts[["n_treated"]],
    prevalence_control = r$counts[["events_control"]] / r$counts[["n_control"]],
    skill_logit_main = get("logit_main", "skill"),
    skill_logit_pairwise = get("logit_pairwise", "skill"),
    skill_xgboost = get("xgboost", "skill"),
    gap_xgb_minus_logit = r$gap[["gap"]],
    gap_rough_se = r$gap[["se"]],
    row.names = NULL
  )
}))
small <- Reduce(`|`, lapply(count_cols, function(cl) summary_tbl[[cl]] < config_report_min_cell))
num_cols <- setdiff(names(summary_tbl), c("candidate", "estimand_status"))
summary_tbl[small, num_cols] <- NA
summary_tbl$suppressed <- small

metrics_tbl <- do.call(rbind, lapply(names(results), function(nm) {
  m <- results[[nm]]$metrics
  if (is.null(m)) return(NULL)
  cbind(candidate = nm, m, row.names = NULL)
}))
if (!is.null(metrics_tbl)) {
  metrics_tbl <- metrics_tbl[!metrics_tbl$candidate %in% summary_tbl$candidate[small], , drop = FALSE]
}

cov_small <- coverage$n_matched < config_report_min_cell | coverage$n_den_ge1 < config_report_min_cell |
  (coverage$n_cohort - coverage$n_matched) < config_report_min_cell
cov_num <- c("rows_in_cohort", "n_matched", "n_den_ge1", "share_matched", "share_den_ge1")
coverage[cov_small, cov_num] <- NA

cli::cli_h2("Screen summary")
print(summary_tbl[, c("candidate", "n", "prevalence_control", "skill_logit_main",
                      "skill_xgboost", "gap_xgb_minus_logit", "gap_rough_se")],
      digits = 3, row.names = FALSE)

arrow::write_parquet(coverage, file.path(out_dir, "outcome_screen_coverage.parquet"))
utils::write.csv(summary_tbl, file.path(out_dir, "outcome_screen_summary.csv"), row.names = FALSE)
if (!is.null(metrics_tbl) && nrow(metrics_tbl) > 0L) {
  arrow::write_parquet(metrics_tbl, file.path(out_dir, "outcome_screen_metrics.parquet"))
}
cli::cli_alert_success("Wrote the outcome_screen_* files to {.file {out_dir}}.")
cli::cli_alert_warning(
  "Screen only. Before choosing an outcome, read {.field estimand_status}: the FUH
   follow-up rows condition on a post-exposure discharge and are not ATT targets."
)
