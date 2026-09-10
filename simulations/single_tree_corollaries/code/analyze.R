# ============================================================
# analyze.R
# Study: single_tree_corollaries  (doubletree)
# Spec:  quality_reports/specs/2026-09-09_single-outcome-tree-corollaries.md  (§5, §4)
#
# Per-regime metrics and the gate-by-gate reading. Sourced by run_pilot.R; safe to
# source standalone after common.R.
#
# ------------------------------------------------------------
# THE ONE METHODOLOGICAL CHOICE SPEC §5 IS EXPLICIT ABOUT
# ------------------------------------------------------------
#
# Super-/sub-efficiency is judged as EMPIRICAL-TO-ANALYTIC, never
# empirical-to-empirical:
#
#   ratio_emp_Vtau = n * Var-hat(theta_tie) / V_tau_full      (should be ~1)
#   ratio_emp_V    = n * Var-hat(theta_tie) / V_full          (~ V_tau/V)
#
# Spec §5 (Oracle's reason, recorded here so it is not re-litigated): the full
# two-tree estimator's OWN finite-sample variance can exceed the efficiency bound
# V because of nuisance-estimation noise, so a ratio of the tied estimator's
# empirical variance to the two-tree estimator's empirical variance could show
# "super-efficiency" for entirely the wrong reason. Both analytic denominators
# come from code/dgps.R's exact cell sums.
#
# MC SEs. A variance RATIO across R independent replications has relative MC SE
# ~sqrt(2/(R-1)) (normal approximation), which spec §4 pre-registers as ~6.3% at
# R = 1000 and which is ~20% at the pilot's R = 50. Coverage carries the
# empirical-proportion SE sqrt(p(1-p)/R). Both are printed next to every number:
# at pilot scale most of these gates are NOT resolvable, and the tables say so
# rather than inviting a reader to over-read them.
# ============================================================

#' Monte Carlo SE of a mean
st_mcse_mean <- function(x) stats::sd(x) / sqrt(length(x))

#' Monte Carlo SE of a proportion
st_mcse_prop <- function(p, R) sqrt(p * (1 - p) / R)

#' Relative MC SE of an empirical variance (hence of a variance ratio)
st_mcse_var_rel <- function(R) sqrt(2 / (R - 1))

#' Per-(regime, n, arm) summary of the spec §5 metrics
#'
#' \code{arm == "two_tree"} rows carry no naive/corrected SE by construction (the
#' split is a tied-estimator concept), so those columns come back NA -- meaningful,
#' not missing. Their variance ratios ARE computed, against the same analytic
#' denominators, because the comparator's realised variance against the efficiency
#' bound V is the reference the tied estimator's ratio is read beside.
#'
#' @param res Long results (\code{st_run_cell()} output, possibly rbound).
#' @param partition Which candidate partition's analytic variance to use as the
#'   denominator. \code{NULL} (default) uses the FIRST partition of each DGP,
#'   which is the only one in Regimes A-C. Regime E has two and must be read via
#'   \code{st_summarise_regime_E()} instead.
#' @return A data.frame, one row per (regime, n, arm).
st_summarise <- function(res, partition = NULL) {
  ok <- res[is.na(res$error_message), , drop = FALSE]
  if (!nrow(ok)) stop("Every replication failed; nothing to summarise.", call. = FALSE)
  keys <- unique(ok[, c("regime", "n", "arm")])
  keys <- keys[order(keys$regime, keys$n, keys$arm), , drop = FALSE]

  rows <- lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, ]
    sub <- ok[ok$regime == k$regime & ok$n == k$n & ok$arm == k$arm, , drop = FALSE]
    R <- nrow(sub)
    dgp <- make_st_dgp(ST_REGIMES[[k$regime]]$dgp)
    pop <- st_population_table(dgp)
    pr <- if (is.null(partition)) pop$partition[[1L]] else partition
    prow <- pop[pop$partition == pr, , drop = FALSE]

    theta0 <- sub$theta0[[1L]]
    bias <- mean(sub$theta) - theta0
    emp_var <- stats::var(sub$theta)
    emp_var_scaled <- k$n * emp_var

    covr <- function(v) {
      if (all(is.na(v))) return(c(NA_real_, NA_real_))
      p <- mean(v, na.rm = TRUE)
      c(p, st_mcse_prop(p, sum(!is.na(v))))
    }
    cn <- covr(sub$cov_naive); cc <- covr(sub$cov_corrected); ce <- covr(sub$cov_eif)
    rec <- if (all(is.na(sub$recovered))) c(NA_real_, NA_real_) else
      covr(sub$recovered)

    data.frame(
      regime = k$regime, dgp = sub$dgp[[1L]], n = k$n, arm = k$arm, R = R,
      theta0 = theta0, mean_theta = mean(sub$theta),
      bias = bias, mcse_bias = st_mcse_mean(sub$theta),
      rel_bias = bias / theta0,
      rmse = sqrt(mean((sub$theta - theta0)^2)),
      emp_var = emp_var, emp_var_scaled = emp_var_scaled,
      # Spec §5: BOTH analytic ratios, never empirical-to-empirical.
      V_tau_full = prow$V_tau_full, V_full = prow$V_full,
      ratio_emp_Vtau = emp_var_scaled / prow$V_tau_full,
      ratio_emp_V = emp_var_scaled / prow$V_full,
      mcse_ratio_rel = st_mcse_var_rel(R),
      analytic_ratio_Vtau_V = prow$ratio_Vtau_V,
      # rem:single-tree-se (Regime D, folded into A-C per spec §3).
      naive_over_correct = mean(sub$naive_over_correct, na.rm = TRUE),
      analytic_naive_over_correct = prow$naive_over_correct,
      mean_se_naive = mean(sub$se_naive, na.rm = TRUE),
      mean_se_corrected = mean(sub$se_corrected, na.rm = TRUE),
      mean_se_eif = mean(sub$se_eif, na.rm = TRUE),
      emp_sd_theta = sqrt(emp_var),
      cov_naive = cn[[1L]], mcse_cov_naive = cn[[2L]],
      cov_corrected = cc[[1L]], mcse_cov_corrected = cc[[2L]],
      cov_eif = ce[[1L]], mcse_cov_eif = ce[[2L]],
      implied_cov_naive_analytic = prow$implied_coverage_naive,
      mean_width_naive = mean(sub$ci_width_naive, na.rm = TRUE),
      mean_width_corrected = mean(sub$ci_width_corrected, na.rm = TRUE),
      recovery = rec[[1L]], mcse_recovery = rec[[2L]],
      max_identity_gap = suppressWarnings(max(sub$identity_gap, na.rm = TRUE)),
      max_mu_leaf_gap = suppressWarnings(max(sub$mu_leaf_gap, na.rm = TRUE)),
      mean_clip_frac = mean(sub$clip_frac, na.rm = TRUE),
      mean_n_leaves_m0 = mean(sub$n_leaves_m0, na.rm = TRUE),
      secs_per_rep = mean(sub$secs, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  # max() over an all-NA column returns -Inf with a warning; make that explicit.
  out$max_identity_gap[!is.finite(out$max_identity_gap)] <- NA_real_
  out$max_mu_leaf_gap[!is.finite(out$max_mu_leaf_gap)] <- NA_real_
  rownames(out) <- NULL
  out
}

#' Regime E: variance STRATIFIED by the realised partition (spec §3, §5)
#'
#' The demonstration `lem:single-tree-linear`'s caveat asks for is that the
#' stratum-specific variances DIFFER -- direct evidence that no single fixed
#' variance describes the tied estimator absent partition stabilisation. Each
#' stratum's empirical variance is therefore reported against ITS OWN analytic
#' V_tau (tree1's against V_tau(tree1), tree2's against V_tau(tree2)), which is
#' the only comparison that distinguishes "the two ties really do have different
#' variances" from "one shared variance, estimated twice".
#'
#' @param res Long results containing regime E rows.
#' @param arm Which arm to stratify (only \code{fitted_tie} is meaningful).
#' @return A data.frame, one row per realised partition.
st_summarise_regime_E <- function(res, arm = "fitted_tie") {
  sub <- res[res$regime == "E" & res$arm == arm & is.na(res$error_message), ,
             drop = FALSE]
  if (!nrow(sub)) return(NULL)
  dgp <- make_st_dgp("E")
  pop <- st_population_table(dgp)

  rows <- lapply(sort(unique(sub$realized_partition)), function(pr) {
    s <- sub[sub$realized_partition == pr, , drop = FALSE]
    prow <- pop[pop$partition == pr, , drop = FALSE]
    ev <- if (nrow(s) >= 2L) stats::var(s$theta) else NA_real_
    data.frame(
      regime = "E", arm = arm, realized_partition = pr,
      n = s$n[[1L]], R = nrow(s), share = nrow(s) / nrow(sub),
      mean_theta = mean(s$theta), theta0 = s$theta0[[1L]],
      bias = mean(s$theta) - s$theta0[[1L]],
      emp_var = ev, emp_var_scaled = s$n[[1L]] * ev,
      V_tau_full = if (nrow(prow)) prow$V_tau_full else NA_real_,
      ratio_emp_Vtau = if (nrow(prow)) s$n[[1L]] * ev / prow$V_tau_full else NA_real_,
      mcse_ratio_rel = if (nrow(s) >= 2L) st_mcse_var_rel(nrow(s)) else NA_real_,
      cov_naive = mean(s$cov_naive), cov_corrected = mean(s$cov_corrected),
      mean_se_corrected = mean(s$se_corrected),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# ---- spec §4's gates, read against the run ------------------------------

#' Gate-by-gate reading of a completed run against spec §4's table
#'
#' Every gate's target is the SPEC's pre-computed number. Nothing here adapts a
#' threshold to the observed result; a failure is reported as a failure, with the
#' observed value, so the decision about what it means is the reader's.
#'
#' Gates 2-4 are stated on V_tau/V. That ratio has TWO readings: the ANALYTIC one
#' (a property of the DGP constants, checked to machine precision by
#' \code{verify_population_gates()}) and the EMPIRICAL one
#' (\code{ratio_emp_V} = n Var-hat / V_full, which estimates it). Both are
#' reported; the analytic reading is what spec §4 computed, and the empirical one
#' is what the simulation adds.
#'
#' @param smry \code{st_summarise()} output.
#' @param arm Which arm the gates are read on. Spec §2's oracle arm is the
#'   mathematically clean one (no tree-search noise), so it is the default for
#'   Regimes A-C; \code{"fitted_tie"} is reported alongside by
#'   \code{st_print_gates()}.
#' @param band_A Spec §4's MC-SE pre-registered band for Regime A's "~1" target.
#'   \strong{Registered for n_sim = 1000}, where the ratio's MC SE is ~6.3%.
#' @param recovery_min Spec §6's fitted-arm partition-recovery floor.
#' @return A data.frame of gate rows.
st_gate_table <- function(smry, arm = "oracle_tie", band_A = c(0.88, 1.14),
                          recovery_min = 0.95) {
  pick <- function(regime, a = arm) {
    r <- smry[smry$regime == regime & smry$arm == a, , drop = FALSE]
    if (nrow(r) != 1L) NULL else r
  }
  rows <- list()
  add <- function(gate, regime, arm_used, check, observed, target, verdict) {
    rows[[length(rows) + 1L]] <<- data.frame(
      gate = gate, regime = regime, arm = arm_used, check = check,
      observed = observed, target = target, pass = verdict,
      stringsAsFactors = FALSE
    )
  }

  A <- pick("A")
  if (!is.null(A)) {
    add("2 (empirical)", "A", arm, "n*Var/V in [0.88, 1.14]",
        A$ratio_emp_V, "0.88-1.14 (registered at nsim=1000)",
        A$ratio_emp_V >= band_A[[1L]] && A$ratio_emp_V <= band_A[[2L]])
    add("2 (analytic)", "A", "-", "|V_tau/V - 1| <= 1e-12",
        abs(A$analytic_ratio_Vtau_V - 1), "<= 1e-12",
        abs(A$analytic_ratio_Vtau_V - 1) <= 1e-12)
  }
  B <- pick("B")
  if (!is.null(B)) {
    add("3 (analytic)", "B", "-", "V_tau/V <= 0.7", B$analytic_ratio_Vtau_V,
        "<= 0.7 (spec: 0.539)", B$analytic_ratio_Vtau_V <= 0.7)
    add("3 (empirical)", "B", arm, "n*Var/V <= 0.7", B$ratio_emp_V,
        "<= 0.7 (spec: 0.539)", B$ratio_emp_V <= 0.7)
  }
  C <- pick("C")
  if (!is.null(C)) {
    add("4 (analytic)", "C", "-", "V_tau/V >= 1.5", C$analytic_ratio_Vtau_V,
        ">= 1.5 (spec: 1.979)", C$analytic_ratio_Vtau_V >= 1.5)
    add("4 (empirical)", "C", arm, "n*Var/V >= 1.5", C$ratio_emp_V,
        ">= 1.5 (spec: 1.979)", C$ratio_emp_V >= 1.5)
  }
  # Gate 5 (Regime D, folded into A-C). Spec §4's verdict row discharges it on
  # Regime C ("A/B give smaller but nonzero gaps"), so C is the gated cell and
  # A/B are reported as observations at the same threshold.
  for (rg in c("A", "B", "C")) {
    x <- pick(rg)
    if (is.null(x)) next
    gate <- if (rg == "C") "5 (gated cell)" else "5 (observation)"
    add(gate, rg, arm, "naive/correct variance <= 0.6", x$naive_over_correct,
        if (rg == "C") "<= 0.6 (spec: 0.311)" else "reported, not gated",
        if (rg == "C") x$naive_over_correct <= 0.6 else NA)
  }
  # Gate 6: bias ~ 0, judged against its own MC SE rather than against zero.
  for (rg in unique(smry$regime)) {
    x <- pick(rg, if (rg == "E") "fitted_tie" else arm)
    if (is.null(x)) next
    add("6", rg, x$arm, "|bias| <= 2 * MCSE(bias)", x$bias,
        sprintf("<= %.5g", 2 * x$mcse_bias), abs(x$bias) <= 2 * x$mcse_bias)
  }
  # Spec §6(ii): fitted-arm partition recovery.
  for (rg in unique(smry$regime)) {
    x <- pick(rg, "fitted_tie")
    if (is.null(x) || is.na(x$recovery)) next
    add("§6(ii)", rg, "fitted_tie", "partition recovery >= 0.95", x$recovery,
        ">= 0.95", x$recovery >= recovery_min)
  }
  # Spec §1/§6(i): the identity, as a smoke check across every arm and regime.
  gaps <- smry$max_identity_gap[!is.na(smry$max_identity_gap)]
  if (length(gaps)) {
    add("§1/§6(i)", "all", "tie arms", "max |theta_tie - theta_forced| <= 1e-10",
        max(gaps), "<= 1e-10", max(gaps) <= IDENTITY_TOL)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# ---- printing -----------------------------------------------------------

#' Print one regime's spec §5 metrics
st_print_regime <- function(smry, regime) {
  sub <- smry[smry$regime == regime, , drop = FALSE]
  if (!nrow(sub)) {
    cli::cli_alert_warning("No rows for regime {.val {regime}}.")
    return(invisible(NULL))
  }
  cols_point <- c("arm", "R", "theta0", "mean_theta", "bias", "mcse_bias",
                  "rel_bias", "rmse")
  cols_var <- c("arm", "emp_var_scaled", "V_tau_full", "ratio_emp_Vtau",
                "V_full", "ratio_emp_V", "mcse_ratio_rel",
                "analytic_ratio_Vtau_V")
  cols_se <- c("arm", "emp_sd_theta", "mean_se_naive", "mean_se_corrected",
               "mean_se_eif", "naive_over_correct", "analytic_naive_over_correct")
  cols_cov <- c("arm", "cov_naive", "mcse_cov_naive", "implied_cov_naive_analytic",
                "cov_corrected", "mcse_cov_corrected", "cov_eif",
                "mean_width_naive", "mean_width_corrected")
  cols_diag <- c("arm", "recovery", "mcse_recovery", "mean_n_leaves_m0",
                 "mean_clip_frac", "max_identity_gap", "max_mu_leaf_gap",
                 "secs_per_rep")

  cli::cli_text("{.strong point estimate}")
  print(knitr::kable(sub[, cols_point], digits = 5, format = "simple",
                     row.names = FALSE))
  cli::cli_text("{.strong variance vs ANALYTIC V_tau and V (spec §5: never empirical-to-empirical)}")
  print(knitr::kable(sub[, cols_var], digits = 5, format = "simple",
                     row.names = FALSE))
  cli::cli_text("{.strong standard errors (naive = treated-arm term only; corrected = + control-arm term)}")
  print(knitr::kable(sub[, cols_se], digits = 5, format = "simple",
                     row.names = FALSE))
  cli::cli_text("{.strong coverage at nominal 95%}")
  print(knitr::kable(sub[, cols_cov], digits = 4, format = "simple",
                     row.names = FALSE))
  cli::cli_text("{.strong diagnostics}")
  # identity_gap and mu_leaf_gap are ~1e-15 and ~1e-7; `digits` would round both
  # to 0 and hide the very quantities spec §1's 1e-10 target is about.
  diag_tab <- sub[, cols_diag]
  for (fld in c("max_identity_gap", "max_mu_leaf_gap")) {
    diag_tab[[fld]] <- formatC(diag_tab[[fld]], format = "e", digits = 3)
  }
  print(knitr::kable(diag_tab, digits = 5, format = "simple", row.names = FALSE))
  invisible(sub)
}

#' Print the gate table for both arms
st_print_gates <- function(smry) {
  for (a in intersect(c("oracle_tie", "fitted_tie"), unique(smry$arm))) {
    cli::cli_h2("Gates read on arm {.val {a}}")
    g <- st_gate_table(smry, arm = a)
    print(knitr::kable(g, digits = 6, format = "simple", row.names = FALSE))
  }
  invisible(NULL)
}
