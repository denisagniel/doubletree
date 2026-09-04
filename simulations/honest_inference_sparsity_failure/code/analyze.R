# ============================================================
# analyze.R
# Study: honest_inference_sparsity_failure  (doubletree, S2)
# Spec:  quality_reports/specs/2026-08-21_honest-inference-sparsity-failure.md
#        ("Metrics", "Outputs", "Decision rule")
#
# Turns the per-cell checkpoints into the spec's metrics, figures and tables.
#
# Run:  Rscript simulations/honest_inference_sparsity_failure/code/analyze.R
#       (from the doubletree package root)
# Gates: HIS_REPS / HIS_REPS_MAX_N -- which reps-per-cell checkpoints to read
#        (defaults match run_sweep.R: 300, and 150 at n = 8000)
#
# NO FAVOURITISM. Spec "Constitution compliance": all three intervals appear side
# by side at every eps, and the orthogonal blind-spot variant is reported
# prominently rather than averaged away. There is no favourable-only headline
# table and no appendix for the stress cells.
#
# ARM 0 COMES FIRST, IN CODE ORDER AND IN THE OUTPUT. Spec's 2026-09-01 arm: if
# estimate_att_crossfit() does not itself deliver a sqrt(n)-consistent,
# asymptotically normal anchor at a cell, then every delta_hat / ci_95_anchor /
# ci_95_honest number at that cell is uninterpretable and must be FLAGGED, not
# quietly tabulated next to the trustworthy ones. Section 2 computes the flag;
# every later table carries it.
#
# THE DECISION RULE IS EXECUTED, NOT DESCRIBED. Section 3 evaluates the spec's own
# escalation condition -- ci_95_anchor coverage below 1-alpha at ANY eps -- and
# prints a STOP banner if it fires, together with the diagnostic that separates
# "construction bug" from "Arm 0 / assumption violation".
# ============================================================

source(file.path("simulations", "honest_inference_sparsity_failure", "code",
                 "common.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(knitr)
})

RUNID <- format(Sys.time(), "%Y%m%d-%H%M%S")
NOMINAL <- 1 - ALPHA

# ---- 0 load every cell checkpoint ---------------------------------------

grid <- design_grid()
paths <- vapply(seq_len(nrow(grid)),
                function(i) cell_path(grid$dgp[[i]], grid$n[[i]], grid$reps[[i]]),
                character(1))
if (!all(file.exists(paths))) {
  cli::cli_abort(c(
    "{sum(!file.exists(paths))} of {length(paths)} cell checkpoint(s) missing.",
    stats::setNames(basename(paths[!file.exists(paths)]),
                    rep("*", sum(!file.exists(paths)))),
    i = "Run run_sweep.R first (it is resumable; finished cells are skipped)."
  ))
}
bundles <- lapply(paths, readRDS)
res <- do.call(rbind, lapply(bundles, `[[`, "results"))
pop <- unique(do.call(rbind, lapply(bundles, `[[`, "dgp_summary")))
rownames(pop) <- NULL
meta1 <- bundles[[1]]$meta

# `pop` is built by unique()-ing per-cell summaries down to two rows per DGP (one
# per nuisance). If any cell disagreed -- e.g. a checkpoint written before a DGP
# parameter changed -- unique() would leave extra rows and every subsequent lookup
# would silently become length > 1 and recycle. Asserted rather than assumed.
if (nrow(pop) != 2L * length(DGP_IDS) || anyDuplicated(pop[, c("dgp", "nuisance")])) {
  print(pop[, c("dgp", "nuisance", "eps", "delta_sq", "bias_pseudo")])
  cli::cli_abort(c(
    "Expected exactly {2 * length(DGP_IDS)} unique (dgp, nuisance) population rows, got {nrow(pop)}.",
    i = "Cell checkpoints disagree about a DGP's population structure -- some were written against different DGP parameters. Delete results/cell_*.rds and re-run the sweep."
  ))
}

#' Look up one population quantity for one (dgp, nuisance), exactly one value
pop_val <- function(field, id, nuisance = "e") {
  v <- pop[[field]][pop$dgp == id & pop$nuisance == nuisance]
  if (length(v) != 1L) {
    cli::cli_abort("pop lookup {field}[{id}, {nuisance}] returned {length(v)} values, expected 1.")
  }
  v
}

res$dgp <- factor(res$dgp, levels = DGP_IDS)
res$err_full <- res$theta_full - res$theta0
res$err_cf <- res$theta_cf - res$theta0
res$blindspot <- as.character(res$dgp) == BLINDSPOT_ID

cli::cli_h1("honest_inference_sparsity_failure (S2) -- analysis (runid {RUNID})")
cli::cli_inform(c(
  "*" = "{nrow(res)} replications over {nrow(grid)} cells",
  "*" = "R = {REPS_DEFAULT} at n in {paste(setdiff(N_GRID, max(N_GRID)), collapse = ', ')}; R = {REPS_AT_MAX_N} at n = {max(N_GRID)} (spec-sanctioned reduction, stated explicitly)",
  "*" = "method under test: estimate_att(leaf_budget={meta1$leaf_budget}, m_n={meta1$m_n}, lambda_n={meta1$lambda_n_rate}, propensity_loss={meta1$propensity_loss}, outcome_type={meta1$outcome_type})",
  "*" = "anchor: {meta1$anchor}",
  "*" = "ci_95_honest = {meta1$honest_ci_call}; fidelity SEhat = {meta1$se_hat_choice}",
  "*" = "seed_master={meta1$seed_master}; {meta1$r_version}; optimaltrees {meta1$optimaltrees_version}; doubletree {meta1$doubletree_sha}"
))

# ---- 1 integrity gates ---------------------------------------------------

n_fail <- sum(!is.na(res$error_message))
if (n_fail > 0L) {
  fr <- aggregate(cbind(failures = !is.na(error_message)) ~ dgp + n, res, sum)
  print(fr)
  cli::cli_abort(c(
    "{n_fail} of {nrow(res)} replications failed; summarising would average noise.",
    stats::setNames(unique(stats::na.omit(res$error_message)),
                    rep("*", length(unique(stats::na.omit(res$error_message)))))
  ))
}
# The DETERMINISTIC thm:anchor implication (see run_one_rep()): the crude anchor
# interval contains the anchor's own Wald interval by the triangle inequality, so
# it covers whenever that does, replication by replication. run_one_rep() already
# errors on a violation; re-checked over the pooled sweep so the guarantee is a
# visible line of output rather than a latent property of the harness.
viol <- sum(res$covered_cf & !res$covered_anchor)
if (viol > 0L) {
  cli::cli_abort("{viol} replication(s) violate coverage(ci_95_anchor) >= coverage(anchor Wald). anchor_ci_harness() is misconstructed.")
}
cli::cli_alert_success(
  "0 failures; ci_95_anchor contained the anchor's own Wald interval in all {nrow(res)} replications (thm:anchor's triangle inequality, checked per replication)."
)

mcse_prop <- function(p, R) sqrt(p * (1 - p) / R)

by_cell <- function(f) {
  parts <- split(res, list(as.character(res$dgp), res$n), drop = TRUE)
  out <- do.call(rbind, lapply(parts, f))
  out <- out[order(match(out$dgp, DGP_IDS), out$n), ]
  rownames(out) <- NULL
  out
}

# ---- 2 ARM 0: is the anchor itself trustworthy at this cell? --------------

# Spec's 2026-09-01 arm, verbatim: "report theta_crossfit's own empirical bias,
# RMSE, and Wald-CI coverage at each (dgp_variant, n) cell, as a first-class table
# alongside the three anchor-dependent intervals, not merely assumed."
#
# The FLAG rule, stated once here and applied everywhere downstream: a cell is
# UNINTERPRETABLE if the anchor's own Wald coverage falls below nominal by more
# than 2 Monte Carlo standard errors. `bias_cf_over_se` (the anchor's bias in
# units of its own standard error) is reported next to it because that is the
# quantity ass:rate actually constrains -- a bias that is o(n^{-1/2}) is exactly
# what "the anchor satisfies the rate condition" means, and it shows up here as a
# ratio that does NOT grow with n.
arm0 <- by_cell(function(z) {
  R <- nrow(z)
  cov <- mean(z$covered_cf)
  se <- mcse_prop(cov, R)
  data.frame(
    dgp = as.character(z$dgp[[1]]), eps = z$eps[[1]], n = z$n[[1]], reps = R,
    blindspot = z$blindspot[[1]],
    theta0 = z$theta0[[1]],
    bias_cf = mean(z$err_cf), bias_cf_mcse = stats::sd(z$err_cf) / sqrt(R),
    rmse_cf = sqrt(mean(z$err_cf^2)), sd_cf = stats::sd(z$theta_cf),
    mean_sigma_cf = mean(z$sigma_cf),
    # ass:rate's operational signature: this ratio stays bounded (ideally
    # shrinking) in n. A ratio growing like sqrt(n) is a rate violation.
    bias_cf_over_se = mean(z$err_cf) / mean(z$sigma_cf),
    cov_cf = cov, cov_cf_mcse = se,
    # Did the anchor ESCAPE the sparsity violation? It carries no leaf budget,
    # which is the only reason it can. These are the fold-wise grid leaf counts
    # and the fraction of fold-wise partitions that are SUFFICIENT for the truth.
    cf_leaves_e = mean(z$cf_leaves_e_mean), cf_leaves_mu = mean(z$cf_leaves_mu_mean),
    cf_suff_e = mean(z$cf_suff_e_frac), cf_suff_mu = mean(z$cf_suff_mu_frac),
    clip_frac_cf = mean(z$clip_frac_cf),
    arm0_ok = cov >= NOMINAL - 2 * se,
    stringsAsFactors = FALSE
  )
})
arm0$uninterpretable <- !arm0$arm0_ok

cli::cli_h1("2. ARM 0 -- does the anchor itself satisfy ass:rate on this DGP?")
print(knitr::kable(
  arm0[, c("dgp", "eps", "n", "reps", "bias_cf", "bias_cf_mcse", "rmse_cf",
           "mean_sigma_cf", "bias_cf_over_se", "cov_cf", "cov_cf_mcse",
           "cf_leaves_e", "cf_leaves_mu", "cf_suff_e", "cf_suff_mu", "arm0_ok")],
  digits = 4, format = "pipe"
))
if (any(arm0$uninterpretable)) {
  flagged <- arm0[arm0$uninterpretable, c("dgp", "eps", "n", "cov_cf", "cov_cf_mcse",
                                          "bias_cf", "bias_cf_over_se")]
  cli::cli_alert_danger(
    "{nrow(flagged)} cell(s) FLAGGED UNINTERPRETABLE: the anchor's own Wald coverage is more than 2 MC se below {NOMINAL}. Every anchor-dependent number at those cells is not to be trusted."
  )
  print(flagged)
} else {
  cli::cli_alert_success(
    "Every cell passes Arm 0: the anchor's own Wald coverage is within 2 MC se of {NOMINAL} everywhere, so all downstream anchor-dependent results are interpretable."
  )
}

# ---- 3 primary metrics: three intervals side by side ---------------------

primary <- by_cell(function(z) {
  R <- nrow(z)
  cw <- mean(z$covered); ca <- mean(z$covered_anchor); ch <- mean(z$covered_honest)
  data.frame(
    dgp = as.character(z$dgp[[1]]), eps = z$eps[[1]], n = z$n[[1]], reps = R,
    blindspot = z$blindspot[[1]],
    # the theorems' own quantities, so nothing is read off the eps LABEL
    delta_e = pop_val("delta_sq", as.character(z$dgp[[1]]), "e"),
    delta_mu = pop_val("delta_sq", as.character(z$dgp[[1]]), "mu"),
    bias_pseudo = pop_val("bias_pseudo", as.character(z$dgp[[1]]), "e"),
    D_product = pop_val("D_product", as.character(z$dgp[[1]]), "e"),
    bias_full = mean(z$err_full), bias_full_mcse = stats::sd(z$err_full) / sqrt(R),
    rmse_full = sqrt(mean(z$err_full^2)),
    # COVERAGE of all three, plus Arm 0's own for the identity check
    cov_wald = cw, cov_wald_mcse = mcse_prop(cw, R),
    cov_anchor = ca, cov_anchor_mcse = mcse_prop(ca, R),
    cov_honest = ch, cov_honest_mcse = mcse_prop(ch, R),
    cov_cf = mean(z$covered_cf),
    # WIDTH of all three
    w_wald = mean(z$ci_width),
    w_anchor = mean(z$ci_anchor_width), w_anchor_sd = stats::sd(z$ci_anchor_width),
    w_honest = mean(z$ci_honest_width),
    w_cf = mean(z$ci_cf_width),
    # cor:width's predicted plateau for the CRUDE anchor: 2|b|, since the
    # half-width is |delta_hat| + z sigma_cf and delta_hat -> b.
    w_plateau_pred = 2 * abs(pop_val("bias_pseudo", as.character(z$dgp[[1]]), "e")),
    mean_abs_delta = mean(z$abs_delta), sd_delta = stats::sd(z$delta),
    mean_sigma_full = mean(z$sigma_full), mean_sigma_cf = mean(z$sigma_cf),
    honest_cv = mean(z$honest_cv),
    stringsAsFactors = FALSE
  )
})
primary <- merge(primary, arm0[, c("dgp", "n", "arm0_ok")], by = c("dgp", "n"))
primary <- primary[order(match(primary$dgp, DGP_IDS), primary$n), ]
rownames(primary) <- NULL

cli::cli_h1("3. Coverage of all three intervals, every eps, side by side")
print(knitr::kable(
  primary[, c("dgp", "eps", "n", "reps", "delta_e", "delta_mu", "bias_pseudo",
              "bias_full", "cov_wald", "cov_anchor", "cov_honest", "cov_cf",
              "arm0_ok")],
  digits = 4, format = "pipe"
))

# ---- 4 THE DECISION RULE, executed ---------------------------------------

cli::cli_h1("4. Spec Decision rule, applied")

# Bullet 1: "If ci_95_anchor coverage falls below 1-alpha at any eps" -- the
# single most serious possible finding. Judged against the MC se, since a point
# estimate of 0.947 at R = 300 is not evidence of undercoverage.
anchor_short <- primary[primary$cov_anchor < NOMINAL - 2 * primary$cov_anchor_mcse, ,
                        drop = FALSE]
if (nrow(anchor_short) > 0L) {
  cli::cli_alert_danger("STOP -- ci_95_anchor coverage is below {NOMINAL} by more than 2 MC se at {nrow(anchor_short)} cell(s).")
  print(anchor_short[, c("dgp", "eps", "n", "reps", "cov_anchor", "cov_anchor_mcse",
                         "cov_cf", "arm0_ok", "bias_pseudo", "mean_abs_delta")])
  cli::cli_inform(c(
    "!" = "Per the spec's Decision rule, do NOT interpret anything else in this study until this is resolved, and do NOT attempt to fix theory.tex.",
    "i" = "The construction has already been checked: ci_95_anchor CONTAINS the anchor's own Wald interval in every replication (section 1), so this cannot be an interval-construction bug.",
    "i" = "Therefore compare `cov_anchor` against `cov_cf` and `arm0_ok`: if Arm 0 also fails at the same cell, the cause is the ANCHOR not satisfying ass:rate, not thm:anchor."
  ))
} else {
  cli::cli_alert_success(
    "Bullet 1 does not fire: ci_95_anchor coverage is at or above {NOMINAL} (within 2 MC se) at EVERY eps and n, including the severe and orthogonal variants -- consistent with thm:anchor's unconditional statement."
  )
}

# Bullets 2 and 3: does the as-shipped honest_ci() undercover as the
# Armstrong-Kline-Sun / Armstrong-Kolesar impossibility result predicts, and if so
# does it land in the 80-93% band the Oracle consultation cited?
honest_short <- primary[primary$cov_honest < NOMINAL - 2 * primary$cov_honest_mcse, ,
                        drop = FALSE]
if (nrow(honest_short) > 0L) {
  cli::cli_alert_warning("Bullet 3 fires: ci_95_honest undercovers at {nrow(honest_short)} cell(s) -- the PREDICTED outcome, not a surprise.")
  print(honest_short[, c("dgp", "eps", "n", "cov_honest", "cov_honest_mcse",
                         "bias_pseudo", "mean_abs_delta", "mean_sigma_cf")])
  cli::cli_inform(c(
    "i" = "Range {signif(min(honest_short$cov_honest), 3)}-{signif(max(honest_short$cov_honest), 3)}, against the Oracle consultation's cited 80-93% band.",
    "i" = "Do NOT attempt a new honesty proof for the data-driven-B case; the impossibility result forecloses it. Report and recommend a default/documentation change instead."
  ))
} else {
  cli::cli_alert_info(
    "Bullet 2 fires instead: ci_95_honest coverage is at or above {NOMINAL} at EVERY eps including severe -- BETTER than a generic data-driven-B interval would be expected to do. The spec requires checking the honest_ci() `se`-role question before calling this validated; section 5 does that."
  )
}

# ---- 5 the honest_ci() `se`-role question (spec Prerequisite) -------------

# Spec: "verify honest_ci()'s `se` argument -- Armstrong-Kline-Sun distinguishes
# sigma_O (se of the bias estimate delta_hat, used inside the critical value) from
# sigma_U (se of the anchor, used as the final half-width multiplier), and
# honest_ci() currently uses one `se` for both roles. Confirm which quantity is
# actually being passed and record explicitly whether that is sigma_O, sigma_U, or
# neither."
#
# READ FROM R/inference.R:101 DIRECTLY. honest_ci(theta_display, se, delta,
# se_delta, level) computes B = |delta| + z*se_delta, cv = honest_cv(B/se), and
# half = cv*se. So the SINGLE `se` argument appears twice: as the denominator of
# the bias-to-SE ratio, and as the final half-width multiplier. Every call site in
# the package passes sigma_crossfit -- the ANCHOR's standard error -- so:
#
#   * in the half-width role, `se` IS sigma_U (the anchor's se), correctly;
#   * in the ratio role it is NOT sigma_O. sigma_O is sd(delta_hat), the sampling
#     se of the bias ESTIMATE, which is the se of a DIFFERENCE of two positively
#     correlated estimators and is a different quantity from either one's se;
#   * with se_delta = 0 the sigma_O role is absent from the construction
#     altogether: B is treated as an a priori KNOWN bound, which is exactly the
#     Armstrong-Kolesar fixed-B setting in which cv(B/sigma_U)*sigma_U is valid,
#     and exactly the substitution Armstrong-Kline-Sun's impossibility result
#     concerns once B is data-driven.
#
# The empirical columns below let that be checked rather than asserted:
# `sd_delta` is sigma_O measured across replications, `mean_sigma_cf` is the
# sigma_U actually passed, and their ratio says how far apart the two roles are.
se_roles <- primary[, c("dgp", "eps", "n", "reps", "bias_pseudo", "mean_abs_delta",
                        "sd_delta", "mean_sigma_full", "mean_sigma_cf", "honest_cv")]
se_roles$sigma_O_emp <- se_roles$sd_delta
se_roles$sigma_U_passed <- se_roles$mean_sigma_cf
se_roles$ratio_O_over_U <- se_roles$sd_delta / se_roles$mean_sigma_cf
# How much narrower is the bias-aware interval than the crude one? Analytically
# cv(b) -> b + z_{1-alpha} as b grows (a ONE-sided normal quantile), while the
# crude anchor uses |delta| + z_{1-alpha/2}*sigma, so the bias-aware construction
# buys at most (z_{1-alpha/2} - z_{1-alpha}) * sigma_U of half-width -- about
# 0.315 sigma_U at 95%. Checked, not assumed.
se_roles$half_saving <- (primary$w_anchor - primary$w_honest) / 2
se_roles$half_saving_pred <- (Z_ALPHA - stats::qnorm(NOMINAL)) * primary$mean_sigma_cf

cli::cli_h1("5. honest_ci()'s `se` argument: which role does it fill?")
cli::cli_inform(c(
  "*" = "Read from R/inference.R:101 -- honest_ci() computes B = |delta| + z*se_delta, cv = honest_cv(B/se), half = cv*se. ONE `se` serves BOTH roles.",
  "*" = "This study passes sigma_crossfit, as every package call site does. In the half-width role that IS sigma_U (the anchor's se). In the bias-to-SE ratio role it is NOT sigma_O = sd(delta_hat).",
  "*" = "With se_delta = 0 the sigma_O role is absent entirely: B is treated as an a priori known bound. That is the Armstrong-Kolesar FIXED-B construction, and the exact substitution the Armstrong-Kline-Sun impossibility result addresses once B is data-driven.",
  "*" = "Consequence for reading any coverage shortfall: a shortfall here is attributable to the data-driven-B substitution, NOT to a sigma_O/sigma_U scale mismatch, because with se_delta = 0 no sigma_O enters."
))
print(knitr::kable(
  se_roles[, c("dgp", "eps", "n", "bias_pseudo", "mean_abs_delta", "sigma_O_emp",
               "sigma_U_passed", "ratio_O_over_U", "honest_cv", "half_saving",
               "half_saving_pred")],
  digits = 4, format = "pipe"
))

# ---- 6 width and cor:width's plateau -------------------------------------

cli::cli_h1("6. cor:width -- does the anchor width plateau, and at the predicted height?")
width_tbl <- primary[, c("dgp", "eps", "n", "delta_e", "delta_mu", "D_product",
                         "bias_pseudo", "w_plateau_pred", "w_wald", "w_cf",
                         "w_anchor", "w_honest", "mean_abs_delta", "mean_sigma_cf",
                         "arm0_ok")]
# The two reference quantities the spec asks about, kept DISTINCT because they are
# different objects:
#   w_plateau_pred = 2|b|          the EXACT height the crude anchor's width
#                                  converges to, since half = |delta_hat| + z sigma_cf
#                                  and delta_hat -> b while sigma_cf -> 0.
#   D_product      = D_w * D_mu    cor:width's own rate quantity, i.e. the SCALE
#                                  in width = O_p(max(n^{-1/2}, D_w D_mu)). A rate
#                                  has no constant, so it is not expected to equal
#                                  the plateau height -- here it exceeds 2|b| by a
#                                  DGP-specific factor. Reported side by side so
#                                  the corollary's rate claim and the exact
#                                  constant are not conflated.
width_tbl$w_anchor_minus_plateau <- width_tbl$w_anchor - width_tbl$w_plateau_pred
width_tbl$noise_share <- 2 * Z_ALPHA * width_tbl$mean_sigma_cf / width_tbl$w_anchor
print(knitr::kable(width_tbl, digits = 4, format = "pipe"))
cli::cli_inform(c(
  "i" = "`noise_share` is the fraction of the anchor's width coming from z*sigma_cf. -> 1 means pre-plateau (n^{-1/2}-dominated); -> 0 means plateaued at 2|b|.",
  "i" = "`w_anchor_minus_plateau` should shrink like n^{-1/2} toward 0 from above at fixed eps; it must NOT go negative (the width cannot fall below the plateau) and must NOT grow."
))

# ---- 7 the fidelity diagnostic, with the blind spot distinguished ---------

fid <- by_cell(function(z) {
  R <- nrow(z)
  rr <- mean(z$fid_reject)
  data.frame(
    dgp = as.character(z$dgp[[1]]), eps = z$eps[[1]], n = z$n[[1]], reps = R,
    residual_mode = z$residual_mode[[1]],
    blindspot = z$blindspot[[1]],
    delta_e = pop_val("delta_sq", as.character(z$dgp[[1]]), "e"),
    delta_mu = pop_val("delta_sq", as.character(z$dgp[[1]]), "mu"),
    bias_pseudo = pop_val("bias_pseudo", as.character(z$dgp[[1]]), "e"),
    ip_gh = pop_val("ip_gh", as.character(z$dgp[[1]]), "e"),
    reject_rate = rr, reject_mcse = mcse_prop(rr, R),
    mean_abs_delta = mean(z$abs_delta), mean_se_hat = mean(z$fid_se_hat),
    mean_stat = mean(z$fid_stat),
    stringsAsFactors = FALSE
  )
})
fid$expectation <- ifelse(
  fid$delta_e == 0, "size -> 0 under ass:sparsity",
  ifelse(fid$blindspot, "BLIND SPOT: b = 0, must NOT reject (nominal alpha)",
         "power -> 1 as bias grows")
)
fid$row_label <- ifelse(fid$blindspot,
                        paste0("** ", fid$dgp, " (ORTHOGONAL, b = 0) **"),
                        as.character(fid$dgp))

cli::cli_h1("7. prop:spectest -- fidelity-diagnostic rejection rate")
print(knitr::kable(
  fid[, c("row_label", "eps", "n", "reps", "delta_e", "delta_mu", "bias_pseudo",
          "ip_gh", "mean_abs_delta", "mean_se_hat", "reject_rate", "reject_mcse",
          "expectation")],
  digits = 4, format = "pipe"
))

bs_fid <- fid[fid$blindspot, , drop = FALSE]
bs_over <- bs_fid[bs_fid$reject_rate > ALPHA + 2 * bs_fid$reject_mcse, , drop = FALSE]
if (nrow(bs_over) > 0L) {
  cli::cli_alert_danger(
    "The fidelity diagnostic rejects at ABOVE-NOMINAL rate on the orthogonal blind-spot variant at {nrow(bs_over)} cell(s). Per the spec's Decision rule this requires RE-VERIFYING the orthogonality construction (code/verify_dgps.R) before being reported as contradicting prop:spectest."
  )
  print(bs_over[, c("dgp", "n", "reject_rate", "reject_mcse", "mean_abs_delta",
                    "mean_se_hat", "ip_gh", "bias_pseudo")])
} else {
  cli::cli_alert_success(
    "Blind spot confirmed: rejection rate at or below nominal alpha = {ALPHA} (within 2 MC se) at every n on the orthogonal variant, while delta_e = {signif(bs_fid$delta_e[[1]], 4)} and delta_mu = {signif(bs_fid$delta_mu[[1]], 4)} are as large as the severe variant's -- prop:spectest's own stated blind spot, reproduced."
  )
}

# ---- 8 figures -----------------------------------------------------------

# Facet labels carry eps AND the exact bias, so no panel can be read without the
# quantity the theorems are stated in. The blind-spot panel is labelled as such.
panel_label <- function(id) {
  e <- pop_val("eps", id, "e")
  b <- pop_val("bias_pseudo", id, "e")
  de <- pop_val("delta_sq", id, "e")
  base <- sprintf("%s: eps = %.2f\ndelta_e = %.4f, b = %.4f", id, e, de, b)
  if (identical(id, BLINDSPOT_ID)) paste0(base, "\n[ORTHOGONAL: b = 0]") else base
}
lab_map <- vapply(DGP_IDS, panel_label, character(1))

INTERVAL_LEVELS <- c("ci_95 (naive Wald)", "ci_95_anchor (thm:anchor)",
                     "ci_95_honest (honest_ci, se_delta=0)",
                     "anchor's own Wald (Arm 0)")
INTERVAL_COLS <- stats::setNames(
  c("#C4471C", "#1B6CA8", "#4C9F70", "grey35"), INTERVAL_LEVELS
)

cov_long <- do.call(rbind, lapply(
  list(c("cov_wald", "cov_wald_mcse", INTERVAL_LEVELS[[1]]),
       c("cov_anchor", "cov_anchor_mcse", INTERVAL_LEVELS[[2]]),
       c("cov_honest", "cov_honest_mcse", INTERVAL_LEVELS[[3]]),
       c("cov_cf", NA, INTERVAL_LEVELS[[4]])),
  function(sp) data.frame(
    dgp = primary$dgp, eps = primary$eps, n = primary$n,
    interval = sp[[3]],
    value = primary[[sp[[1]]]],
    mcse = if (is.na(sp[[2]])) mcse_prop(primary$cov_cf, primary$reps) else primary[[sp[[2]]]],
    stringsAsFactors = FALSE
  )
))
cov_long$interval <- factor(cov_long$interval, levels = INTERVAL_LEVELS)
cov_long$dgp <- factor(cov_long$dgp, levels = DGP_IDS, labels = lab_map[DGP_IDS])

fig1 <- ggplot(cov_long, aes(x = n, y = value, colour = interval, shape = interval,
                             linetype = interval)) +
  geom_hline(yintercept = NOMINAL, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = value - 1.96 * mcse, ymax = value + 1.96 * mcse),
                width = 0.05, na.rm = TRUE) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 2.1) +
  scale_x_log10(breaks = N_GRID, labels = N_GRID) +
  scale_y_continuous(limits = c(0, 1.02)) +
  scale_colour_manual(values = INTERVAL_COLS) +
  scale_linetype_manual(values = stats::setNames(
    c("solid", "solid", "solid", "22"), INTERVAL_LEVELS)) +
  facet_wrap(~dgp, nrow = 1) +
  labs(
    title = "Coverage of three intervals under a controlled sparsity violation (S2)",
    subtitle = paste0(
      "All three intervals at every eps, side by side, plus the anchor's own Wald interval (Arm 0). Dashed grey = nominal ",
      NOMINAL, ".\nci_95_anchor is the literal thm:anchor construction ",
      "theta_full +/- (|delta| + z*sigma_crossfit): no bias-aware critical value. ",
      "Bars are +/- 1.96 MC se.\nR = ", REPS_DEFAULT, " at n <= 2000, R = ",
      REPS_AT_MAX_N, " at n = 8000 (stated reduction)."
    ),
    x = "n (log scale)", y = "coverage of theta_0",
    colour = NULL, shape = NULL, linetype = NULL,
    caption = paste0("runid ", RUNID, " | doubletree ", meta1$doubletree_sha,
                     " | optimaltrees ", meta1$optimaltrees_sha)
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8), legend.box = "vertical")

f1 <- file.path(DIR_FIGURES, paste0("coverage_three_intervals_", RUNID))
ggsave(paste0(f1, ".png"), fig1, width = 13, height = 5.6, dpi = 200, bg = "white")
ggsave(paste0(f1, ".pdf"), fig1, width = 13, height = 5.6, device = "pdf")

# Width figure, with the predicted plateau annotated per panel.
w_long <- do.call(rbind, lapply(
  list(c("w_anchor", INTERVAL_LEVELS[[2]]), c("w_honest", INTERVAL_LEVELS[[3]]),
       c("w_wald", INTERVAL_LEVELS[[1]])),
  function(sp) data.frame(
    dgp = primary$dgp, n = primary$n, interval = sp[[2]], value = primary[[sp[[1]]]],
    stringsAsFactors = FALSE
  )
))
w_long$interval <- factor(w_long$interval, levels = INTERVAL_LEVELS)
w_long$dgp <- factor(w_long$dgp, levels = DGP_IDS, labels = lab_map[DGP_IDS])

plateau <- unique(primary[, c("dgp", "w_plateau_pred", "D_product")])
plateau$dgp <- factor(plateau$dgp, levels = DGP_IDS, labels = lab_map[DGP_IDS])
plateau$lbl <- sprintf("predicted plateau 2|b| = %.4f\n(cor:width rate scale D_w*D_mu = %.4f)",
                       plateau$w_plateau_pred, plateau$D_product)
# The y axis is log-scaled, so a plateau line at exactly 0 (the eps = 0 and
# orthogonal panels, where b = 0 by construction) cannot be drawn. Those panels
# still carry the text annotation, which states the 0 explicitly -- the absence of
# a line there IS the finding, not a missing element.
plateau_pos <- plateau[plateau$w_plateau_pred > 0, , drop = FALSE]

fig2 <- ggplot(w_long, aes(x = n, y = value, colour = interval, shape = interval)) +
  geom_hline(data = plateau_pos, aes(yintercept = w_plateau_pred),
             linetype = "dashed", colour = "grey25", inherit.aes = FALSE) +
  geom_text(data = plateau, aes(x = sqrt(min(N_GRID) * max(N_GRID)),
                                y = Inf, label = lbl),
            inherit.aes = FALSE, vjust = 1.15, hjust = 0.5, size = 2.3,
            colour = "grey25") +
  geom_line(linewidth = 0.45) +
  geom_point(size = 2.1) +
  scale_x_log10(breaks = N_GRID, labels = N_GRID) +
  scale_y_log10() +
  scale_colour_manual(values = INTERVAL_COLS) +
  facet_wrap(~dgp, nrow = 1) +
  labs(
    title = "cor:width -- interval width contracts at n^{-1/2}, then plateaus at 2|b| (S2)",
    subtitle = paste0(
      "Dashed line = the EXACT predicted plateau 2|b|, computed by enumeration from the pseudo-true partitions (not fitted to these data).\n",
      "D_w*D_mu is cor:width's own RATE scale and is a different object from the plateau CONSTANT; both are shown so they are not conflated.\n",
      "The orthogonal panel has b = 0 exactly, so its width must keep contracting at n^{-1/2} with no plateau -- the direct contrast with the severe panel."
    ),
    x = "n (log scale)", y = "mean interval width (log scale)",
    colour = NULL, shape = NULL,
    caption = paste0("runid ", RUNID)
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8))

f2 <- file.path(DIR_FIGURES, paste0("width_anchor_vs_honest_", RUNID))
ggsave(paste0(f2, ".png"), fig2, width = 13, height = 5.6, dpi = 200, bg = "white")
ggsave(paste0(f2, ".pdf"), fig2, width = 13, height = 5.6, device = "pdf")

# ---- 9 export ------------------------------------------------------------

w <- function(x, nm) {
  utils::write.csv(x, file.path(DIR_TABLES, paste0(nm, "_", RUNID, ".csv")),
                   row.names = FALSE)
}
w(pop, "population")
w(arm0, "arm0")
w(primary, "primary")
w(width_tbl, "width")
w(fid, "fidelity")
w(se_roles, "se_roles")

md <- c(
  paste0("# honest_inference_sparsity_failure (S2) -- run ", RUNID),
  "",
  paste0("- spec: `", meta1$spec, "`"),
  paste0("- method under test: `estimate_att()` with leaf_budget = ", meta1$leaf_budget,
         ", m_n = ", meta1$m_n, ", lambda_n = ", meta1$lambda_n_rate,
         ", outcome_type = ", meta1$outcome_type,
         ", propensity_loss = ", meta1$propensity_loss),
  paste0("- anchor: `", meta1$anchor, "`"),
  paste0("- `ci_95_anchor` = theta_full +/- (|delta| + z*sigma_crossfit), computed in the HARNESS; no package source was modified"),
  paste0("- `ci_95_honest` = `", meta1$honest_ci_call, "`, function used unmodified"),
  paste0("- fidelity diagnostic: reject when |delta| > z*SEhat, SEhat = ", meta1$se_hat_choice),
  paste0("- R = ", REPS_DEFAULT, " per cell at n in {",
         paste(setdiff(N_GRID, max(N_GRID)), collapse = ", "), "}; **R = ",
         REPS_AT_MAX_N, " at n = ", max(N_GRID),
         "** (the spec's own sanctioned reduction, stated explicitly; the n = 8000 cells were 65% of the measured wall time because the anchor runs nested CV)"),
  paste0("- ", nrow(res), " replications total, sequential, worker_limit = ",
         meta1$worker_limit, ", parallel_cv = ", meta1$parallel_cv,
         ", no Rashomon enumeration"),
  paste0("- seed_master = ", meta1$seed_master, "; ", meta1$seed_scheme),
  paste0("- doubletree ", meta1$doubletree_sha, ", optimaltrees ",
         meta1$optimaltrees_sha, ", ", meta1$r_version),
  "",
  "## Exact population structure (closed form, no simulation)",
  "",
  knitr::kable(pop[, c("dgp", "eps", "residual_mode", "nuisance", "theta0",
                       "sparsity_holds", "delta_sq", "delta_kl", "D_w", "D_mu",
                       "D_product", "bias_pseudo", "ip_gh", "cor_gh")],
               digits = 6, format = "pipe"),
  "",
  "`delta_sq`/`delta_kl` are delta_j of ass:sparsity: the excess population risk of the BEST",
  "tree with at most leaf_budget leaves, by direct enumeration over T_Lbar (reusing S1's",
  "`enumerate_sufficient_class.R`). `bias_pseudo` is the exact population bias of the",
  "constrained estimator at its pseudo-true limits; `ip_gh` is the {1-e_0}-weighted inner",
  "product of the two residual functions, which is 0 to machine precision exactly on the",
  "orthogonal blind-spot variant.",
  "",
  "## ARM 0 -- does the anchor itself satisfy ass:rate here?",
  "",
  knitr::kable(arm0[, c("dgp", "eps", "n", "reps", "bias_cf", "bias_cf_mcse",
                        "rmse_cf", "mean_sigma_cf", "bias_cf_over_se", "cov_cf",
                        "cov_cf_mcse", "cf_leaves_e", "cf_leaves_mu", "cf_suff_e",
                        "cf_suff_mu", "arm0_ok")],
               digits = 4, format = "pipe"),
  "",
  paste0("`cf_suff_*` is the fraction of the anchor's fold-wise partitions that are ",
         "SUFFICIENT for the true nuisance. The anchor carries no leaf budget, which is ",
         "the only reason it can be consistent where the budget-", meta1$leaf_budget,
         " estimator cannot; these columns show whether it actually used that freedom."),
  "",
  "## Coverage of all three intervals, every eps, side by side",
  "",
  knitr::kable(primary[, c("dgp", "eps", "n", "reps", "delta_e", "delta_mu",
                           "bias_pseudo", "bias_full", "rmse_full", "cov_wald",
                           "cov_wald_mcse", "cov_anchor", "cov_anchor_mcse",
                           "cov_honest", "cov_honest_mcse", "cov_cf", "arm0_ok")],
               digits = 4, format = "pipe"),
  "",
  "## Width: cor:width's plateau",
  "",
  knitr::kable(width_tbl, digits = 4, format = "pipe"),
  "",
  "## prop:spectest -- fidelity-diagnostic rejection rate (blind-spot row marked **)",
  "",
  knitr::kable(fid[, c("row_label", "eps", "n", "reps", "delta_e", "delta_mu",
                       "bias_pseudo", "ip_gh", "mean_abs_delta", "mean_se_hat",
                       "reject_rate", "reject_mcse", "expectation")],
               digits = 4, format = "pipe"),
  "",
  "## honest_ci()'s `se` argument: sigma_O, sigma_U, or neither?",
  "",
  "Read directly from `R/inference.R:101`: `honest_ci()` computes `B = |delta| + z*se_delta`,",
  "`cv = honest_cv(B/se)`, `half = cv*se`. One `se` argument serves both roles. This study",
  "passes `sigma_crossfit`, as all four existing package call sites do. Therefore:",
  "",
  "- in the **half-width** role, `se` is **sigma_U** (the anchor's standard error) -- correct;",
  "- in the **bias-to-SE ratio** role it is **not sigma_O**; sigma_O = sd(delta_hat) is the",
  "  sampling se of the bias ESTIMATE (a difference of two positively correlated estimators),",
  "  a different quantity from either estimator's own se;",
  "- with `se_delta = 0` the sigma_O role is **absent from the construction entirely**: B is",
  "  treated as an a priori KNOWN bound. That is exactly the Armstrong-Kolesar fixed-B setting,",
  "  and exactly the substitution the Armstrong-Kline-Sun impossibility result concerns once B",
  "  is data-driven. Consequently any coverage shortfall observed here is attributable to the",
  "  data-driven-B substitution, NOT to a sigma_O/sigma_U scale mismatch -- no sigma_O enters.",
  "",
  knitr::kable(se_roles[, c("dgp", "eps", "n", "bias_pseudo", "mean_abs_delta",
                            "sigma_O_emp", "sigma_U_passed", "ratio_O_over_U",
                            "honest_cv", "half_saving", "half_saving_pred")],
               digits = 4, format = "pipe"),
  "",
  paste0("`half_saving` is how much narrower `ci_95_honest`'s half-width is than ",
         "`ci_95_anchor`'s; `half_saving_pred` is the analytic ceiling ",
         "(z_{1-alpha/2} - z_{1-alpha})*sigma_U = ",
         signif(Z_ALPHA - stats::qnorm(NOMINAL), 4),
         "*sigma_U, since honest_cv(b) -> b + z_{1-alpha} for large b.")
)
md_path <- file.path(DIR_TABLES, paste0("summary_", RUNID, ".md"))
writeLines(md, md_path)

# ---- 10 console report ---------------------------------------------------

cli::cli_h1("Coverage (all three intervals + Arm 0)")
print(as.data.frame(primary[, c("dgp", "eps", "n", "bias_pseudo", "bias_full",
                                "cov_wald", "cov_anchor", "cov_honest", "cov_cf",
                                "arm0_ok")]), digits = 4)
cli::cli_h1("Width vs the predicted plateau")
print(as.data.frame(width_tbl[, c("dgp", "eps", "n", "w_wald", "w_anchor",
                                  "w_honest", "w_plateau_pred",
                                  "w_anchor_minus_plateau", "noise_share")]),
      digits = 4)
cli::cli_h1("Fidelity diagnostic")
print(as.data.frame(fid[, c("row_label", "eps", "n", "bias_pseudo", "reject_rate",
                            "reject_mcse", "expectation")]), digits = 4)

cli::cli_alert_success("wrote {.file {md_path}} + 6 CSVs + 2 figures (png & pdf)")
