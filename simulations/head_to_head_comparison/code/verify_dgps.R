# ============================================================
# verify_dgps.R
# Study: head_to_head_comparison  (doubletree, S5)
# Spec:  quality_reports/specs/2026-09-08_head-to-head-comparison.md  (§3)
#
# Run:  Rscript simulations/head_to_head_comparison/code/verify_dgps.R
#       (from the doubletree package root)
#
# STANDS ALONE AND MUST PASS BEFORE ANY REPLICATION RUNS. It fits nothing and
# runs no estimator; every number below is either a closed-form population
# quantity or a deterministic property of a draw.
#
# WHY THIS FILE IS NOT OPTIONAL
#
# Spec §3 records that TWO DGPs previously proposed for claim C1 -- S1's F1 reused
# verbatim, and single_tree_coverage's "complex" DGP -- each turned out to have a
# population bilinear remainder of EXACTLY ZERO (1.7e-18 and -5e-20) under
# main-effects-only misspecification, for a structural reason (OLS orthogonality
# plus conditional independence given the shared confounder) that no amount of
# effect-size tuning repairs. Either would have reported "GLM main-effects:
# unbiased" for C1 -- a null result caused entirely by the DGP -- discoverable
# only after a full implementation and a full 1000-replication run.
#
# Section 1 below is the guard against a third instance of that. If DGP-A's
# remainder comes back near zero, the gate-swap is implemented wrong, NOTHING
# downstream is informative, and this script stops.
#
# WHAT IT CHECKS
#   1. DGP-A's population bilinear remainder is the verified 0.01445 (NOT zero),
#      and prop:bilinear's identity pi*bias == remainder holds.
#   2. DGP-A is genuinely confounded (cor(A, Y0) from a draw, plus the four
#      conditional means spec §3 verified by hand).
#   3. Every DGP's population spec agrees with its own sample draw, per cell, at
#      1e-12 -- including the two REUSED families, since a transcription drift in
#      a file this study only reads would silently invalidate every number.
#   4. Every DGP's declared coordinate roles hold on the exact cell grid.
#   5. R4 CALIBRATION PROBE (user decision 4, 2026-09-09): whether ST3's shipped
#      calibration reliably produces a truly ILL-DEFINED propensity weight (a leaf
#      with zero controls) at the n values this study uses, as opposed to merely a
#      large one. Not pre-emptively escalated; measured, and the observation
#      recorded in README.md either way.
# ============================================================

source(file.path("simulations", "head_to_head_comparison", "code", "common.R"))

suppressPackageStartupMessages(library(knitr))

TABLES <- list()
STAMP <- format(Sys.time(), "%Y%m%d-%H%M%S")

# ---- 1 THE MANDATORY GUARD: DGP-A's remainder is nonzero ------------------

cli::cli_h1("1. DGP-A: population bilinear remainder under main-effects misspecification")

spec_A <- dgp_spec_shared_interaction()
rem <- verify_main_effects_remainder(spec_A)   # ERRORS if ~0 or off spec's value
TABLES$remainder <- rem

cli::cli_alert_success(
  "Remainder = {signif(rem$bilinear_remainder_P, 6)} (spec §3 verified {rem$expected}); NOT cancelled."
)
cli::cli_inform(c(
  "*" = "pi = P(A=1) = {signif(rem$pi_pop, 5)}; theta_0 = {signif(rem$theta0, 5)}",
  "*" = "implied asymptotic bias = remainder/pi = {signif(rem$implied_bias_P, 5)}",
  "*" = "relative bias against theta_0 = {signif(100 * rem$relative_bias_P, 4)}%  (spec §3 predicted ~45%)",
  "i" = "control-weighted (nu) mu-projection variant: remainder {signif(rem$bilinear_remainder_nu, 5)}, relative bias {signif(100 * rem$relative_bias_nu, 4)}% -- reported as a sensitivity; the assertion is on the P-weighted value spec §3 states.",
  "!" = "This remainder is a LINEAR-PROBABILITY population diagnostic establishing that the misspecification is non-degenerate. att_linear() fits a LOGIT main-effects model, whose population limit differs slightly, so the realised GLM bias in R1 will be near but not equal to {signif(rem$implied_bias_P, 4)}."
))

# The four cell values, printed so a reader can check the construction by eye
# against spec §3's formula rather than trusting a function name.
cli::cli_h2("DGP-A's (X1, X2) cells (X3-X5 are pure noise and do not enter)")
grid_A <- unique(data.frame(
  X1 = spec_A$cells$X1, X2 = spec_A$cells$X2,
  e0 = spec_A$e0, mu0 = spec_A$mu0, mu1 = spec_A$mu1
))
grid_A <- grid_A[order(grid_A$X1, grid_A$X2), ]
TABLES$cells_A <- grid_A
print(knitr::kable(grid_A, digits = 4, format = "simple", row.names = FALSE))

cli::cli_h2("Exact tree-sufficiency of DGP-A at Lbar = 4")
cli::cli_inform(c(
  "*" = "e_0 takes {length(unique(spec_A$e0))} distinct values (3-leaf tree: split X1, then X2 within X1 = 1)",
  "*" = "mu_0 takes {length(unique(spec_A$mu0))} distinct values (3-leaf tree: split X2, then X1 within X2 = 1)",
  "*" = "Lbar = {spec_A$leaf_budget} (the full 2x2 grid), so the global optimiser needs no particular split ORDER; ass:sparsity holds EXACTLY, which is what makes R1 a test of C1 and not of C2."
))
stopifnot(length(unique(spec_A$e0)) == 3L, length(unique(spec_A$mu0)) == 3L)

# ---- 2 DGP-A is genuinely confounded ------------------------------------

cli::cli_h1("2. DGP-A: confounding, from an actual draw (not from a comment)")

conf <- confounding_check_shared_interaction(spec_A, n = 200000L)
TABLES$confounding <- conf
print(knitr::kable(t(conf), format = "simple"))

if (abs(conf$emp_cor_A_Y0) < 0.02) {
  cli::cli_abort(c(
    "cor(A, Y0) = {signif(conf$emp_cor_A_Y0, 4)}: DGP-A is effectively UNCONFOUNDED.",
    x = "No amount of nuisance misfit could then move theta-hat, and C1 would be vacuous.",
    i = "This is the exact defect S1's design audit caught on its own DGP."
  ))
}
cli::cli_alert_success(
  "Confounded: cor(A, Y0) = {signif(conf$emp_cor_A_Y0, 4)}; E[e_0|X1] = {signif(conf$pop_Ee0_x1_0, 3)} -> {signif(conf$pop_Ee0_x1_1, 3)}; E[mu_0|X1] = {signif(conf$pop_Emu0_x1_0, 3)} -> {signif(conf$pop_Emu0_x1_1, 3)} (spec §3: 0.20/0.60 and 0.30/0.50)."
)

# ---- 3 & 4 every DGP: spec vs draw, and declared coordinate roles ---------

cli::cli_h1("3. Every DGP: population spec vs sample draw (per cell, tol 1e-12) and declared roles")

# make_dgp() runs check_variable_roles() (or its _dial sibling) unconditionally,
# so building each DGP IS check 4.
dgps <- lapply(DGP_IDS, make_dgp)
names(dgps) <- DGP_IDS

role_rows <- lapply(DGP_IDS, function(id) {
  d <- dgps[[id]]
  cmp <- check_spec_vs_draw(d, n = 40000L)
  cli::cli_alert_success(
    "{id}: spec matches its draw on all {nrow(cmp)} cells; roles e={paste(d$e_vars, collapse='+')}, mu={paste(d$mu_vars, collapse='+')} verified."
  )
  dgp_summary_row(d)
})
TABLES$dgps <- do.call(rbind, role_rows)
cli::cli_h2("Population summary of every DGP this study runs")
print(knitr::kable(TABLES$dgps, digits = 5, format = "simple", row.names = FALSE))

# The three reused families must be UNMODIFIED. Their key population constants are
# restated here as an assertion, so a drift in a file this study only READS
# surfaces as an error rather than as quietly different results.
cli::cli_h2("Reused DGPs: their spec-stated population constants")
local({
  a2 <- dgps$A2
  if (abs(a2$theta0 - 0.08) > 1e-12) {
    cli::cli_abort("F1 (DGP-A2) theta_0 = {signif(a2$theta0, 8)}, expected 0.08 (its post-2026-09-01 favourable calibration).")
  }
  cc <- dgps$C
  # Spec §3 verified e_0(X1=1, X2=1) = 0.9820 by direct arithmetic.
  if (abs(max(cc$e0) - stats::plogis(-2.5 + 4.0 + 2.5)) > 1e-12) {
    cli::cli_abort("ST3 (DGP-C) max e_0 = {signif(max(cc$e0), 8)}, not plogis(-2.5+4+2.5).")
  }
  cli::cli_alert_success(
    "F1 theta_0 = 0.08; ST3 max e_0 = {signif(max(cc$e0), 6)} (1 - e_0 = {signif(1 - max(cc$e0), 4)}); DGP-B bias predictions present for all four eps."
  )
  b_bias <- vapply(dgps[c("B_eps0", "B_eps05", "B_eps10", "B_eps15")],
                   function(z) z$pseudo$bias, numeric(1))
  print(knitr::kable(data.frame(
    dgp = names(b_bias),
    eps = vapply(dgps[names(b_bias)], `[[`, numeric(1), "eps"),
    delta_e_kl = vapply(dgps[names(b_bias)], function(z) z$class_e$delta_kl, numeric(1)),
    delta_mu_kl = vapply(dgps[names(b_bias)], function(z) z$class_mu$delta_kl, numeric(1)),
    predicted_bias = as.numeric(b_bias),
    row.names = NULL
  ), digits = 5, format = "simple", row.names = FALSE))
})

# ---- 5 R4 calibration probe (user decision 4) ----------------------------

cli::cli_h1("5. R4 / DGP-C calibration probe: is the weak-overlap weight ILL-DEFINED or merely large?")

cli::cli_inform(c(
  "i" = "Decision (user, 2026-09-09): do NOT pre-emptively escalate ST3's calibration. Use it as-is, MEASURE whether it produces a truly ill-defined weight at the n values used, and escalate only if it does not.",
  "*" = "A doubletree leaf that isolates the extreme cell has an UNDEFINED m_0 (and e-hat driven to the clip) exactly when that cell contains ZERO CONTROLS. That event -- not 'the weight is big' -- is what this probe counts."
))

probe_c <- local({
  spec <- dgps$C
  x_extreme <- spec$cells$X1 == 1L & spec$cells$X2 == 1L
  e_extreme <- unique(spec$e0[x_extreme])
  p_extreme <- sum(spec$cell_prob[x_extreme])
  ns <- sort(unique(c(REGIMES$R4$n, N_GRID_FULL$R4)))
  reps <- 400L
  rows <- lapply(ns, function(n) {
    set.seed(20260908L + n)
    counts <- vapply(seq_len(reps), function(r) {
      d <- spec$draw(n)
      in_cell <- d$X$X1 == 1L & d$X$X2 == 1L
      c(n_cell = sum(in_cell), n_control = sum(in_cell & d$A == 0L))
    }, numeric(2))
    data.frame(
      n = n, reps = reps,
      e0_extreme = e_extreme, one_minus_e0 = 1 - e_extreme,
      p_cell = p_extreme,
      expected_cell_n = n * p_extreme,
      expected_controls = n * p_extreme * (1 - e_extreme),
      mean_controls = mean(counts["n_control", ]),
      min_controls = min(counts["n_control", ]),
      frac_zero_controls = mean(counts["n_control", ] == 0),
      frac_le2_controls = mean(counts["n_control", ] <= 2),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
})
TABLES$probe_C <- probe_c
print(knitr::kable(probe_c, digits = 4, format = "simple", row.names = FALSE))

pilot_ns <- REGIMES$R4$n
p_zero <- probe_c$frac_zero_controls[probe_c$n %in% pilot_ns]
if (all(p_zero < 0.01)) {
  cli::cli_alert_warning(c(
    "ST3 as shipped produces a ZERO-control extreme cell in {signif(100 * max(p_zero), 3)}% of draws at the pilot's n ({paste(pilot_ns, collapse = ', ')}).",
    "i" = "That is a LARGE weight regime, not an ILL-DEFINED one. If R4's pilot shows no inference degradation, escalation (steeper propensity coefficients, or a smaller n) is the next step -- see README.md's R4 note. NOT escalated here, per the user's decision."
  ))
} else {
  cli::cli_alert_success(
    "ST3 as shipped hits a ZERO-control extreme cell in up to {signif(100 * max(p_zero), 3)}% of draws at the pilot's n: genuinely ill-defined, no escalation needed."
  )
}

# ---- 6 write the tables --------------------------------------------------

for (nm in names(TABLES)) {
  utils::write.csv(TABLES[[nm]],
                   file.path(DIR_TABLES, sprintf("verify_%s_%s.csv", nm, STAMP)),
                   row.names = FALSE)
}
saveRDS(TABLES, file.path(DIR_RESULTS, sprintf("verify_dgps_%s.rds", STAMP)))
cli::cli_alert_success("Wrote {length(TABLES)} verification table(s) to {.path {DIR_TABLES}} (stamp {STAMP}).")
cli::cli_h1("All DGP verifications passed. Safe to run replications.")
