# ============================================================
# analyze.R
# Study: head_to_head_comparison  (doubletree, S5)
# Spec:  quality_reports/specs/2026-09-08_head-to-head-comparison.md  (see §5, §6)
#
# Metrics, with a Monte Carlo standard error on EVERY reported number. Spec §6 is
# explicit that MC SEs are the "no impress-only design" artefact: a coverage of
# 0.93 at R = 20 and a coverage of 0.93 at R = 1000 are different findings, and a
# table that prints only the point estimate makes them look identical.
#
# THE HEADLINE STATISTIC IS A RATIO, NOT A RANKING (spec §0, §5). C1 is
# "doubletree's RMSE ratio to the ORACLE-nuisance benchmark is ~1", not
# "doubletree beat the GLM". The ratio is computed PAIRED -- both arms on the same
# replication's data -- and its MC SE comes from a paired bootstrap over
# replications, because the two arms' errors are strongly positively correlated
# (same draw) and a delta-method or independent-samples SE would badly overstate
# the uncertainty in the ratio.
#
# Sourced by run_pilot.R / run_sweep.R. Sourcing has no side effects.
# ============================================================

# ---- 1 per-cell, per-arm metrics ----------------------------------------

#' Monte Carlo SE of an empirical proportion
#'
#' Uses the EMPIRICAL proportion, not spec §6's \eqn{\sqrt{0.95 \times 0.05/R}}.
#' That formula is the right thing for SIZING a study in advance (it is the SE
#' under the null that coverage is nominal); for REPORTING a realised coverage of,
#' say, 0.62 it understates the SE by a factor of ~2.2. Both are printed by
#' \code{summarise_cell()} so the sizing target stays visible.
mc_se_prop <- function(p, R) sqrt(p * (1 - p) / R)

#' Per-(regime, dgp, n, arm) metrics with Monte Carlo SEs
#'
#' @param res A results data.frame from \code{run_cell()} (or an rbind of several).
#' @return A tibble-shaped data.frame, one row per (regime, dgp, n, arm).
summarise_cell <- function(res) {
  keys <- unique(res[, c("regime", "dgp", "n", "arm")])
  keys <- keys[order(keys$regime, keys$dgp, keys$n, keys$arm), , drop = FALSE]
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, ]
    sub <- res[res$regime == k$regime & res$dgp == k$dgp &
                 res$n == k$n & res$arm == k$arm, , drop = FALSE]
    ok <- is.na(sub$error_message) & is.finite(sub$theta)
    d <- sub[ok, , drop = FALSE]
    R <- nrow(d)
    if (R == 0L) {
      return(data.frame(regime = k$regime, dgp = k$dgp, n = k$n, arm = k$arm,
                        reps = 0L, n_fail = nrow(sub),
                        stringsAsFactors = FALSE))
    }
    theta0 <- d$theta0[[1]]
    err <- d$theta - theta0
    bias <- mean(err)
    cov <- mean(d$covered)
    data.frame(
      regime = k$regime, dgp = k$dgp, n = k$n, arm = k$arm,
      reps = R, n_fail = sum(!ok),
      theta0 = theta0,
      # Spec §5: report n_treated alongside n in EVERY table -- effective
      # precision for the ATT scales with n_treated, not n.
      n_treated_mean = mean(d$n_treated),
      bias = bias,
      bias_mc_se = stats::sd(d$theta) / sqrt(R),
      rel_bias = bias / theta0,
      rmse = sqrt(mean(err^2)),
      # MC SE of the RMSE via the delta method on mean(err^2):
      # se(rmse) = se(mse) / (2 * rmse).
      rmse_mc_se = (stats::sd(err^2) / sqrt(R)) / (2 * sqrt(mean(err^2))),
      sd_theta = stats::sd(d$theta),
      mean_sigma = mean(d$sigma),
      coverage = cov,
      coverage_mc_se = mc_se_prop(cov, R),
      # The sizing-target SE of spec §6, kept visible next to the realised one.
      coverage_mc_se_nominal = sqrt(0.95 * 0.05 / R),
      ci_width = mean(d$ci_width),
      ci_width_mc_se = stats::sd(d$ci_width) / sqrt(R),
      # Interpretability proxy and certification (doubletree arms only; NA
      # elsewhere, which is meaningful, not missing).
      n_leaves_e = mean(d$n_leaves_e),
      n_leaves_m0 = mean(d$n_leaves_m0),
      certified_frac = mean(d$certified_e & d$certified_m0),
      # Spec §5 / Oracle's Rashomon-multiplicity point: coverage CONDITIONAL on a
      # certified global optimum. If several near-optimal partitions differ, the
      # "auditable" claim is weaker than a single certified fit implies, so the
      # conditional number must be reported, not just the marginal one.
      coverage_certified = if (any(!is.na(d$certified_e))) {
        cc <- d$certified_e & d$certified_m0
        if (any(cc)) mean(d$covered[cc]) else NA_real_
      } else NA_real_,
      n_certified = if (any(!is.na(d$certified_e))) {
        sum(d$certified_e & d$certified_m0)
      } else NA_integer_,
      clip_frac = mean(d$clip_frac),
      max_e_hat = max(d$max_e_hat),
      min_e_hat = min(d$min_e_hat),
      secs_mean = mean(d$secs),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, lapply(rows, function(r) {
    # Failure-only cells return a short row; pad so the rbind cannot silently
    # drop columns (the same discipline as rep_row_template()).
    template <- rows[[which.max(vapply(rows, ncol, integer(1)))]][0, , drop = FALSE]
    miss <- setdiff(names(template), names(r))
    for (m in miss) r[[m]] <- template[[m]][NA_integer_]
    r[, names(template), drop = FALSE]
  }))
  rownames(out) <- NULL
  out
}

# ---- 2 the C1/C3 statistic: a PAIRED RMSE ratio with its own MC SE -------

#' Paired RMSE ratio of every arm to a reference arm, with a bootstrap MC SE
#'
#' \eqn{RMSE(arm) / RMSE(ref)} computed on the SAME replications (same draws), and
#' its MC SE from a paired bootstrap that resamples REPLICATIONS (not arms), so
#' the strong positive correlation between two arms' errors on a shared draw is
#' preserved. Spec §5 names this -- "RMSE ratio to oracle-AIPW, with Monte Carlo
#' SE on the ratio" -- as THE C1 statistic, in place of a raw "beats GLM"
#' comparison the paper makes no claim to.
#'
#' A ratio near 1 is the C1 finding. A ratio BELOW 1 is not a doubletree victory
#' over the oracle -- it is Monte Carlo noise, and the SE is what says so.
#'
#' @param res Results for ONE (regime, dgp, n) cell.
#' @param ref Reference arm name (e.g. \code{"oracle_aipw"}).
#' @param B Bootstrap resamples.
#' @param seed Bootstrap seed. Fixed, so the reported SE is reproducible.
#' @return One row per non-reference arm, or \code{NULL} if \code{ref} is absent.
rmse_ratio_table <- function(res, ref, B = 2000L, seed = 20260908L) {
  if (is.null(ref) || !ref %in% res$arm) return(NULL)
  cells <- unique(res[, c("regime", "dgp", "n")])
  out <- lapply(seq_len(nrow(cells)), function(i) {
    c_i <- cells[i, ]
    sub <- res[res$regime == c_i$regime & res$dgp == c_i$dgp & res$n == c_i$n, ,
               drop = FALSE]
    # Only replications where EVERY arm succeeded, so the pairing is genuine.
    bad <- unique(sub$rep[!is.na(sub$error_message) | !is.finite(sub$theta)])
    sub <- sub[!sub$rep %in% bad, , drop = FALSE]
    if (!nrow(sub)) return(NULL)
    arms <- unique(sub$arm)
    reps <- sort(unique(sub$rep))
    theta0 <- sub$theta0[[1]]
    # err[rep, arm]
    err <- vapply(arms, function(a) {
      s <- sub[sub$arm == a, , drop = FALSE]
      s$theta[match(reps, s$rep)] - theta0
    }, numeric(length(reps)))
    if (length(reps) == 1L) err <- matrix(err, nrow = 1L, dimnames = list(NULL, arms))
    rmse_of <- function(idx) sqrt(colMeans(err[idx, , drop = FALSE]^2))
    point <- rmse_of(seq_along(reps))
    set.seed(seed)
    boot <- vapply(seq_len(B), function(b) {
      idx <- sample.int(length(reps), replace = TRUE)
      r <- rmse_of(idx)
      r / r[[ref]]
    }, numeric(length(arms)))
    if (length(arms) == 1L) boot <- matrix(boot, nrow = 1L, dimnames = list(arms, NULL))
    ratio <- point / point[[ref]]
    data.frame(
      regime = c_i$regime, dgp = c_i$dgp, n = c_i$n,
      arm = arms, ref = ref, reps_paired = length(reps),
      rmse = as.numeric(point),
      rmse_ratio = as.numeric(ratio),
      rmse_ratio_mc_se = apply(boot, 1, stats::sd),
      rmse_ratio_lo = apply(boot, 1, stats::quantile, probs = 0.025),
      rmse_ratio_hi = apply(boot, 1, stats::quantile, probs = 0.975),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  if (is.null(out)) return(NULL)
  rownames(out) <- NULL
  out[out$arm != ref, , drop = FALSE]
}

# ---- 3 compact printers -------------------------------------------------

#' The columns worth looking at first, per regime
#'
#' Deliberately narrow. The full \code{summarise_cell()} output has ~28 columns;
#' a 28-column console table is not read, it is scrolled past.
print_regime_summary <- function(smry, regime) {
  s <- smry[smry$regime == regime, , drop = FALSE]
  if (!nrow(s)) {
    cli::cli_alert_warning("No results for regime {.val {regime}}.")
    return(invisible(NULL))
  }
  cols <- c("dgp", "n", "arm", "reps", "n_treated_mean", "bias", "bias_mc_se",
            "rel_bias", "rmse", "coverage", "coverage_mc_se", "ci_width",
            "secs_mean")
  print(knitr::kable(s[, cols], digits = 4, format = "simple", row.names = FALSE))
  invisible(s)
}

#' Print the paired RMSE-ratio table for a regime, if it has a reference arm
print_ratio_summary <- function(res, regime) {
  ref <- REGIMES[[regime]]$ratio_ref
  if (is.null(ref)) {
    cli::cli_inform("{regime}: no RMSE-ratio reference arm (by design).")
    return(invisible(NULL))
  }
  tab <- rmse_ratio_table(res[res$regime == regime, , drop = FALSE], ref = ref)
  if (is.null(tab) || !nrow(tab)) {
    cli::cli_alert_warning("{regime}: reference arm {.val {ref}} produced no paired replications.")
    return(invisible(NULL))
  }
  cli::cli_alert_info("{regime}: RMSE ratios to {.val {ref}} (paired; bootstrap MC SE)")
  print(knitr::kable(tab, digits = 4, format = "simple", row.names = FALSE))
  invisible(tab)
}

#' Load every checkpoint on disk for a set of regimes
load_results <- function(regimes = REGIME_IDS, reps = NULL) {
  files <- list.files(DIR_RESULTS, pattern = "^cell_.*\\.rds$", full.names = TRUE)
  if (!length(files)) {
    cli::cli_abort("No checkpoints in {.path {DIR_RESULTS}}; run the pilot or sweep first.")
  }
  payloads <- lapply(files, readRDS)
  res <- do.call(rbind, lapply(payloads, `[[`, "results"))
  res <- res[res$regime %in% regimes, , drop = FALSE]
  if (!is.null(reps)) res <- res[res$rep <= reps, , drop = FALSE]
  rownames(res) <- NULL
  res
}
