# ============================================================
# verify_dgps.R
# Study: honest_inference_sparsity_failure  (doubletree, S2)
#
# VERIFY THE DGP FAMILY'S MATHEMATICAL PROPERTIES BY COMPUTATION, BEFORE ANY
# REPLICATION RUNS. Nothing here is simulated except where a simulation is used
# as an INDEPENDENT check on a closed-form quantity.
#
# Run:  Rscript simulations/honest_inference_sparsity_failure/code/verify_dgps.R
#       (from the doubletree package root)
#
# WHY THIS FILE EXISTS AS A SEPARATE, MANDATORY STEP
#
# S1 (partition_recovery_clt) recorded a defect worth repeating the lesson from:
# a documented design-audit fix (add confounding) never reached the code that
# actually ran, and the resulting DGP -- A independent of Y(0) by construction --
# made that study's central metric pairing VACUOUS while producing no error, no
# warning and entirely plausible-looking output. It was caught only by computing
# cor(A, Y0) and E[Y0|A=1] - E[Y0|A=0] empirically rather than trusting the
# code's own comments.
#
# This study's three central claims each rest on a specific mathematical property
# of the new DGP family, and each property is exactly the kind that fails
# silently:
#
#   thm:anchor      needs the sparsity violation to be REAL, i.e. no <=Lbar-leaf
#                   tree absorbs g or h, with delta_e, delta_mu > 0 and scaling
#                   sensibly in eps. If delta were secretly 0 the "unconditional"
#                   claim would never be tested away from the boundary.
#   cor:width       needs the plateau height to be a PREDICTION, i.e. the exact
#                   bias at the pseudo-true limits must be computable and the
#                   pseudo-true partition must be unique.
#   prop:spectest   needs the blind-spot variant's g, h to be GENUINELY orthogonal
#                   under the {1-e_0}-weighted inner product, with the true bias
#                   exactly 0 while both deltas stay substantial. The spec's own
#                   Decision rule says that if the diagnostic rejects at
#                   above-nominal rate there, the FIRST suspect is this
#                   construction, not the proposition -- so it must be verified
#                   before the sweep, not after.
#
# Sections:
#   1 spec-vs-draw and declared coordinate roles, per DGP
#   2 exact population structure (closed form): delta_j, argmin, D_w, D_mu, bias
#   3 INEXPRESSIBILITY: S_j empty, delta_j > 0, delta_j = eps^2, and the smallest
#     leaf budget at which the truth IS representable
#   4 ORTHOGONALITY: <g,h>_nu, the residual-level inner product, and bias == 0
#   5 CONFOUNDING: empirical cor(A, Y0) and E[Y0|A=1] - E[Y0|A=0] vs closed form
#   6 INDEPENDENT CHECK of population_bias(): the package's own psi_att() solved
#     at the pseudo-true nuisances on a huge draw, against the closed form
# ============================================================

source(file.path("simulations", "honest_inference_sparsity_failure", "code",
                 "common.R"))

suppressPackageStartupMessages(library(knitr))

# Gates, all defaulted so a bare invocation does the real thing (and so the file
# can be sourced in a live session and stepped through).
N_CHECK_DRAW <- as.integer(Sys.getenv("HIS_N_CHECK", "40000"))
N_CONFOUND   <- as.integer(Sys.getenv("HIS_N_CONFOUND", "200000"))
N_ORACLE     <- as.integer(Sys.getenv("HIS_N_ORACLE", "2000000"))
MAX_BUDGET_PROBE <- as.integer(Sys.getenv("HIS_MAX_BUDGET", "5"))
cli::cli_inform(c("*" = "gates: HIS_N_CHECK={N_CHECK_DRAW} HIS_N_CONFOUND={N_CONFOUND} HIS_N_ORACLE={N_ORACLE} HIS_MAX_BUDGET={MAX_BUDGET_PROBE}"))

set.seed(SEED_MASTER)

dgps <- lapply(DGP_IDS, make_dgp)   # build_dgp_dial() runs check_variable_roles_dial()
names(dgps) <- DGP_IDS

# ---- 1 population spec vs sample draw, and declared coordinate roles -----

cli::cli_h1("1. Population spec vs sample draw; declared coordinate roles")
for (id in DGP_IDS) {
  cmp <- check_spec_vs_draw_dial(dgps[[id]], n = N_CHECK_DRAW)
  cli::cli_inform(paste0(
    "  OK ", id, ": spec matches its draw on all ", nrow(cmp),
    " cells (tol 1e-12); roles e = ",
    paste(dgps[[id]]$e_vars, collapse = "+"), ", mu = ",
    paste(dgps[[id]]$mu_vars, collapse = "+"), "; X1 in both (confounded)."
  ))
}

# ---- 2 exact population structure ---------------------------------------

cli::cli_h1("2. Exact population structure (closed form, no simulation)")
pop <- do.call(rbind, lapply(dgps, dgp_population_summary_dial))
rownames(pop) <- NULL
print(knitr::kable(
  pop[, c("dgp", "eps", "residual_mode", "nuisance", "theta0", "card_T",
          "sparsity_holds", "n_sufficient", "delta_sq", "delta_kl", "margin_kl",
          "n_argmin", "argmin_key", "argmin_leaves", "n_distinct_values")],
  digits = 7, format = "simple"
))
cli::cli_h2("cor:width's own quantities, and the exact bias at the pseudo-true limits")
print(knitr::kable(
  unique(pop[, c("dgp", "eps", "residual_mode", "D_w", "D_mu", "D_product",
                 "bias_pseudo", "bilinear_form", "ip_gh", "cor_gh",
                 "ip_resid", "cor_resid", "e0_min", "e0_max")]),
  digits = 7, format = "simple"
))

# ---- 3 INEXPRESSIBILITY --------------------------------------------------

cli::cli_h1("3. Inexpressibility of g/h by any tree with <= {LEAF_BUDGET} leaves")

# (a) For eps != 0, S_j must be EMPTY and delta_j strictly positive. For eps = 0
#     the reverse: S_j non-empty and delta_j exactly 0. Both are gates, not notes.
inexpr <- do.call(rbind, lapply(DGP_IDS, function(id) {
  d <- dgps[[id]]
  data.frame(
    dgp = id, eps = d$eps, residual_mode = d$residual_mode,
    n_suff_e = d$class_e$n_sufficient, n_suff_mu = d$class_mu$n_sufficient,
    delta_e_sq = d$class_e$delta_sq, delta_mu_sq = d$class_mu$delta_sq,
    delta_e_kl = d$class_e$delta_kl, delta_mu_kl = d$class_mu$delta_kl,
    # The scaling claim: delta_sq should equal eps^2 for the aligned variants
    # (the residual is exactly orthogonal to every <=2-leaf partition's leaf
    # space, so the whole residual variance survives projection).
    ratio_e = if (d$eps > 0) d$class_e$delta_sq / d$eps^2 else NA_real_,
    ratio_mu = if (d$eps > 0) d$class_mu$delta_sq / d$eps^2 else NA_real_,
    stringsAsFactors = FALSE
  )
}))
print(knitr::kable(inexpr, digits = 7, format = "simple"))

bad <- character(0)
for (id in DGP_IDS) {
  d <- dgps[[id]]
  if (d$eps == 0) {
    if (!d$class_e$sparsity_holds || !d$class_mu$sparsity_holds) {
      bad <- c(bad, paste0(id, ": eps = 0 but ass:sparsity does NOT hold"))
    }
    if (d$class_e$delta_sq > 1e-12 || d$class_mu$delta_sq > 1e-12) {
      bad <- c(bad, paste0(id, ": eps = 0 but delta_j > 0"))
    }
  } else {
    if (d$class_e$sparsity_holds || d$class_mu$sparsity_holds) {
      bad <- c(bad, paste0(id, ": eps != 0 but SOME <=Lbar-leaf partition is ",
                           "sufficient -- the violation is not real"))
    }
    if (d$class_e$delta_sq <= 0 || d$class_mu$delta_sq <= 0) {
      bad <- c(bad, paste0(id, ": eps != 0 but delta_e or delta_mu is not > 0"))
    }
    if (d$class_e$n_argmin != 1L || d$class_mu$n_argmin != 1L) {
      bad <- c(bad, paste0(id, ": the pseudo-true partition is NOT unique (",
                           d$class_e$n_argmin, " argmins for e, ",
                           d$class_mu$n_argmin, " for mu) -- 'the exact bias at ",
                           "the pseudo-true limit' would be a set, not a number"))
    }
    # delta_j must be essentially the WHOLE residual variance, i.e. the <=Lbar-leaf
    # class absorbs (almost) none of it. A ratio well below 1 means the
    # pseudo-true partition has started trading step structure for residual
    # structure -- the exact failure design note (1) records; build_dgp_dial()
    # already gates on the argmin, and this is the quantitative restatement.
    for (j in c("e", "mu")) {
      cls <- if (j == "e") d$class_e else d$class_mu
      ratio <- cls$delta_sq / d$eps^2
      if (ratio < 0.9 || ratio > 1.0000001) {
        bad <- c(bad, paste0(id, ": delta_", j, "_sq / eps^2 = ",
                             signif(ratio, 4), ", not in [0.9, 1]"))
      }
    }
  }
}
if (length(bad)) {
  cli::cli_abort(c("Inexpressibility check FAILED.",
                   stats::setNames(bad, rep("x", length(bad)))))
}
cli::cli_alert_success(
  "S_j empty and delta_j > 0 for every eps != 0; S_j non-empty and delta_j = 0 at eps = 0; pseudo-true partition unique everywhere; delta_j = eps^2 to within 10%."
)

# (a2) MONOTONICITY of the exact bias in eps, across the aligned variants. This is
#      not decoration: prop:spectest's power claim ("consistent against any
#      alternative with non-vanishing bias") and cor:width's "width grows with
#      eps, then plateaus" both read the aligned variants as an ORDERED sequence.
#      A first version of this DGP family, with a wider eps range, produced
#      b(0.10) = 0.049 > b(0.18) = 0.022 -- non-monotone, because mu's pseudo-true
#      partition had flipped (design note (1) in dgps_delta_dial.R). Every number
#      in that run still looked plausible, which is precisely why this is a gate.
cli::cli_h2("Monotonicity of the exact bias along the aligned eps sequence")
aligned_ids <- DGP_IDS[vapply(DGP_IDS, function(i) DGP_SPECS[[i]]$mode == "aligned",
                              logical(1))]
mono <- data.frame(
  dgp = aligned_ids,
  eps = vapply(aligned_ids, function(i) dgps[[i]]$eps, numeric(1)),
  abs_bias = vapply(aligned_ids, function(i) abs(dgps[[i]]$pseudo$bias), numeric(1)),
  bias_over_eps2 = vapply(aligned_ids, function(i) {
    if (dgps[[i]]$eps == 0) NA_real_ else abs(dgps[[i]]$pseudo$bias) / dgps[[i]]$eps^2
  }, numeric(1)),
  D_product = vapply(aligned_ids, function(i) {
    dgps[[i]]$pseudo$D_w * dgps[[i]]$pseudo$D_mu
  }, numeric(1)),
  stringsAsFactors = FALSE
)
mono <- mono[order(mono$eps), ]
rownames(mono) <- NULL
print(knitr::kable(mono, digits = 6, format = "simple"))
if (any(diff(mono$abs_bias) <= 0)) {
  cli::cli_abort(c(
    "|bias| is NOT strictly increasing along the aligned eps sequence.",
    i = "cor:width's plateau reading and prop:spectest's power reading both treat these variants as ordered. Do NOT run the sweep."
  ))
}
if (any(diff(mono$D_product) <= 0)) {
  cli::cli_abort("D_w * D_mu is not increasing in eps either; the realised-error product has stopped tracking the dial.")
}
cli::cli_alert_success(
  "|bias| and D_w*D_mu both strictly increase in eps; b/eps^2 is stable at ~{signif(mean(mono$bias_over_eps2, na.rm = TRUE), 4)}."
)

# (b) The decisive, budget-level statement: what IS the smallest leaf budget at
#     which some tree represents the truth exactly? Lbar = 2 must be strictly
#     below it. Computed by re-enumerating at successively larger budgets --
#     the same enumeration machinery, so the answer is in the same currency.
cli::cli_h2("Smallest leaf budget at which the truth IS exactly representable")
budget_probe <- do.call(rbind, lapply(DGP_IDS, function(id) {
  d <- dgps[[id]]
  min_L <- function(gamma0) {
    for (L in seq_len(MAX_BUDGET_PROBE)) {
      en <- enumerate_tree_partitions(d$p, L, d$var_names)
      ok <- vapply(seq_len(nrow(en$labels)),
                   function(i) is_sufficient_partition(en$labels[i, ], gamma0),
                   logical(1))
      if (any(ok)) return(L)
    }
    NA_integer_
  }
  data.frame(
    dgp = id, eps = d$eps, leaf_budget_used = d$leaf_budget,
    min_leaves_exact_e = min_L(d$e0), min_leaves_exact_mu = min_L(d$mu0),
    n_distinct_e = length(unique(d$e0)), n_distinct_mu = length(unique(d$mu0)),
    stringsAsFactors = FALSE
  )
}))
print(knitr::kable(budget_probe, format = "simple"))
if (any(budget_probe$eps > 0 &
        (budget_probe$min_leaves_exact_e <= LEAF_BUDGET |
         budget_probe$min_leaves_exact_mu <= LEAF_BUDGET))) {
  cli::cli_abort("Some eps != 0 DGP is exactly representable within Lbar = {LEAF_BUDGET}.")
}
if (any(budget_probe$eps > 0 & (is.na(budget_probe$min_leaves_exact_e) |
                                is.na(budget_probe$min_leaves_exact_mu)))) {
  cli::cli_alert_warning(
    "For some variant the truth is not representable within {MAX_BUDGET_PROBE} leaves either; raise HIS_MAX_BUDGET to locate it."
  )
}
cli::cli_alert_success(
  "Every eps != 0 variant needs strictly more than Lbar = {LEAF_BUDGET} leaves; the anchor, which carries NO leaf budget, is the only estimator that can reach the truth."
)

# ---- 4 ORTHOGONALITY of the blind-spot construction ----------------------

cli::cli_h1("4. Orthogonality: <g,h>_nu, the residual-level inner product, and b")

orth <- do.call(rbind, lapply(DGP_IDS, function(id) {
  d <- dgps[[id]]
  data.frame(
    dgp = id, eps = d$eps, residual_mode = d$residual_mode,
    # prop:bilinear's inner product on the RESIDUAL FUNCTIONS g, h themselves --
    # the object the spec's blind-spot instruction names.
    ip_gh = d$ip_gh$ip, cor_gh = d$ip_gh$cor,
    norm_g = d$ip_gh$norm_u, norm_h = d$ip_gh$norm_v,
    # The same inner product on the REALISED pseudo-true error functions
    # (odds(e_hat) - odds(e_0), m0_hat - mu_0) -- the object the bias actually
    # depends on. Reported alongside because orthogonality of g and h does not
    # by itself imply orthogonality of the projection residuals; agreement
    # between the two columns is the substantive content of the construction.
    ip_resid = d$ip_resid$ip, cor_resid = d$ip_resid$cor,
    bilinear_form = d$pseudo$bilinear_form,
    bias_pseudo = d$pseudo$bias,
    # Both deltas must stay SUBSTANTIAL in the blind-spot variant: "both
    # nuisances badly wrong, bias exactly zero" is the claim, not "one nuisance
    # is fine".
    delta_e_sq = d$class_e$delta_sq, delta_mu_sq = d$class_mu$delta_sq,
    stringsAsFactors = FALSE
  )
}))
print(knitr::kable(orth, digits = 8, format = "simple"))

bs <- dgps[[BLINDSPOT_ID]]
sev <- dgps[["eps_sev"]]
cli::cli_h2("The blind-spot contrast, stated as numbers")
cli::cli_inform(c(
  "*" = "eps_sev    (g = h = r4): <g,h>_nu = {signif(sev$ip_gh$ip, 8)}, cor = {signif(sev$ip_gh$cor, 6)}, bias = {signif(sev$pseudo$bias, 6)}",
  "*" = "orth_blind (g = r4, h = r5): <g,h>_nu = {format(bs$ip_gh$ip, scientific = TRUE, digits = 3)}, cor = {format(bs$ip_gh$cor, scientific = TRUE, digits = 3)}, bias = {format(bs$pseudo$bias, scientific = TRUE, digits = 3)}",
  "*" = "same eps ({bs$eps}), delta_e {signif(sev$class_e$delta_sq, 6)} vs {signif(bs$class_e$delta_sq, 6)}; delta_mu {signif(sev$class_mu$delta_sq, 6)} vs {signif(bs$class_mu$delta_sq, 6)}"
))

TOL_ORTH <- 1e-12
fail <- character(0)
if (abs(bs$ip_gh$ip) > TOL_ORTH) {
  fail <- c(fail, paste0("<g,h>_nu = ", signif(bs$ip_gh$ip, 6),
                         " is not 0 to ", TOL_ORTH))
}
if (abs(bs$ip_resid$ip) > TOL_ORTH) {
  fail <- c(fail, paste0("residual-level inner product = ",
                         signif(bs$ip_resid$ip, 6), " is not 0 to ", TOL_ORTH))
}
if (abs(bs$pseudo$bias) > 1e-12) {
  fail <- c(fail, paste0("exact bias at the pseudo-true limits = ",
                         signif(bs$pseudo$bias, 6), " is not 0"))
}
if (bs$class_e$delta_sq <= 0 || bs$class_mu$delta_sq <= 0) {
  fail <- c(fail, "the blind-spot variant does not have delta_e, delta_mu > 0")
}
# The contrast is only sharp if the deltas MATCH the severe variant's: otherwise
# a difference in rejection rate could be a magnitude effect rather than an
# orthogonality effect.
if (abs(bs$class_e$delta_sq - sev$class_e$delta_sq) > 1e-9) {
  fail <- c(fail, paste0("delta_e differs between eps_sev (",
                         signif(sev$class_e$delta_sq, 6), ") and orth_blind (",
                         signif(bs$class_e$delta_sq, 6),
                         "), so the blind-spot contrast is confounded with magnitude"))
}
# delta_mu cannot be EXACTLY equal across the two variants, and the direction of
# the difference is what matters. nu_x = P_x{1 - e_0(x)} depends on X4 through
# e_0, so the aligned h = r4 is very slightly NON-centred under nu
# (E_nu[r4] = -eps/(1-e_step) within each X1 leaf) and loses that much of its
# variance to the projection, whereas the orthogonal h = r5 is exactly centred and
# keeps all of it. Hence delta_mu(orth_blind) = eps^2 exactly while
# delta_mu(eps_sev) = 0.90 * eps^2 at eps = 0.15. The blind-spot variant is
# therefore, if anything, MORE misspecified in mu than the severe one -- which
# strengthens rather than weakens "both nuisances substantially wrong, bias
# exactly zero". The gate is that directional statement, not numerical equality:
# an 11% delta_mu difference could not plausibly explain a rejection rate moving
# from ~1 to ~alpha, and delta_e is identical to machine precision.
if (bs$class_mu$delta_sq < sev$class_mu$delta_sq - 1e-12) {
  fail <- c(fail, paste0("orth_blind is LESS misspecified in mu than eps_sev ",
                         "(delta_mu ", signif(bs$class_mu$delta_sq, 6), " < ",
                         signif(sev$class_mu$delta_sq, 6),
                         "), so a lower rejection rate there could be a ",
                         "magnitude effect rather than an orthogonality effect"))
}
if (bs$class_e$delta_sq < sev$class_e$delta_sq - 1e-12) {
  fail <- c(fail, "orth_blind is less misspecified in e than eps_sev")
}
if (abs(bs$pseudo$bias - bs$pseudo$bias_from_bilinear) > 1e-12) {
  fail <- c(fail, "population_bias()'s two routes to b disagree")
}
if (length(fail)) {
  cli::cli_abort(c("Orthogonality check FAILED -- do NOT run the sweep.",
                   stats::setNames(fail, rep("x", length(fail)))))
}
cli::cli_alert_success(paste0(
  "orth_blind: <g,h>_nu and the residual-level inner product are 0 to machine ",
  "precision, the exact bias is 0, delta_e matches eps_sev's exactly and ",
  "delta_mu is >= eps_sev's -- so the only difference between the two variants ",
  "that could explain a behavioural gap is orthogonality."
))

# Consistency of population_bias()'s two algebraic routes, for EVERY variant.
route_gap <- vapply(dgps, function(d) abs(d$pseudo$bias - d$pseudo$bias_from_bilinear),
                    numeric(1))
if (any(route_gap > 1e-12)) {
  print(route_gap)
  cli::cli_abort("population_bias()'s direct and bilinear routes disagree somewhere.")
}
cli::cli_alert_success("population_bias(): direct and bilinear routes agree to 1e-12 for all {length(dgps)} variants.")

# ---- 5 CONFOUNDING, empirically (the S1 lesson) --------------------------

cli::cli_h1("5. Confounding: computed from an actual draw, not read off the comments")
conf <- do.call(rbind, lapply(dgps, confounding_check, n = N_CONFOUND))
rownames(conf) <- NULL
print(knitr::kable(conf, digits = 5, format = "simple"))

if (any(abs(conf$emp_gap) < 0.02)) {
  cli::cli_abort(c(
    "E[Y0|A=1] - E[Y0|A=0] is essentially 0 for some variant: the DGP is ",
    "effectively UNCONFOUNDED, so nuisance misfit cannot bias theta_hat and ",
    "this study's central claims are untestable on it. This is the exact defect ",
    "S1's design audit caught."
  ))
}
if (any(abs(conf$emp_gap - conf$pop_gap) > 0.01)) {
  cli::cli_abort("Empirical and closed-form confounding gaps disagree by > 0.01.")
}
cli::cli_alert_success(
  "All variants confounded: |E[Y0|A=1] - E[Y0|A=0]| in [{signif(min(abs(conf$emp_gap)), 3)}, {signif(max(abs(conf$emp_gap)), 3)}], matching the closed form."
)

# ---- 6 INDEPENDENT CHECK of the exact bias -------------------------------

cli::cli_h1("6. Independent check: the package's own psi_att() at the pseudo-true nuisances")

# Solve doubletree's OWN estimating equation on a huge draw, with the ORACLE
# pseudo-true nuisance functions plugged in cell-wise. If population_bias()'s
# closed form is right, this converges to theta_0 + bias_pseudo. This is the one
# place in this file where simulation is used -- as an independent check on the
# algebra, not as the source of the number.
oracle_bias <- do.call(rbind, lapply(DGP_IDS, function(id) {
  d <- dgps[[id]]
  set.seed(SEED_MASTER + 7L)
  dat <- d$draw(N_ORACLE)
  key_of <- function(df) {
    do.call(paste, c(lapply(d$var_names, function(v) df[[v]]), sep = "-"))
  }
  idx <- match(key_of(dat$X), key_of(d$cells))
  if (anyNA(idx)) stop("Unmatched covariate pattern in the oracle draw.", call. = FALSE)
  e_hat <- pmax(CLIP_LO, pmin(CLIP_HI, d$class_e$pseudo_true[idx]))
  m0_hat <- d$class_mu$pseudo_true[idx]
  pi_hat <- mean(dat$A)
  eta <- list(e = e_hat, m0 = m0_hat, m1 = NULL)
  s0 <- doubletree::psi_att(dat$Y, dat$A, theta = 0, eta, pi_hat)
  theta_hat <- sum(s0) / sum(dat$A / pi_hat)
  data.frame(
    dgp = id, eps = d$eps, residual_mode = d$residual_mode, n_oracle = N_ORACLE,
    theta0 = d$theta0,
    bias_exact = d$pseudo$bias,
    bias_oracle_mc = theta_hat - d$theta0,
    abs_gap = abs((theta_hat - d$theta0) - d$pseudo$bias),
    # Monte Carlo se of the oracle estimate, so "agree" has a scale.
    mc_se = sqrt(mean(doubletree::psi_att(dat$Y, dat$A, theta_hat, eta, pi_hat)^2) / N_ORACLE),
    stringsAsFactors = FALSE
  )
}))
print(knitr::kable(oracle_bias, digits = 7, format = "simple"))
if (any(oracle_bias$abs_gap > 4 * oracle_bias$mc_se)) {
  print(oracle_bias[oracle_bias$abs_gap > 4 * oracle_bias$mc_se, ])
  cli::cli_abort(c(
    "The closed-form bias and the oracle Monte Carlo disagree by more than 4 MC se.",
    i = "population_bias() is wrong, or the pseudo-true projections are not what the estimator converges to. Do NOT run the sweep."
  ))
}
cli::cli_alert_success(
  "Closed-form bias agrees with the oracle Monte Carlo within 4 MC se for all {nrow(oracle_bias)} variants (max gap {signif(max(oracle_bias$abs_gap), 3)}, max MC se {signif(max(oracle_bias$mc_se), 3)})."
)

# ---- 7 what the bias means at the sampled n ------------------------------

cli::cli_h1("7. Is the plateau TRANSITION inside the sweep?")
# cor:width: width contracts at n^{-1/2} while the approximation error is small
# relative to it, then plateaus at a height set by the approximation-error
# product. The transition is observable only if |b| crosses the Wald scale inside
# {500, 2000, 8000}. sigma is not known in closed form here, so the crude
# reference is the EIF scale sqrt(Var(psi_0)/n) at the ORACLE nuisances, computed
# once per DGP from the same huge draw used above -- an order-of-magnitude guide
# printed for sizing, not a number any claim rests on.
scale_tbl <- do.call(rbind, lapply(DGP_IDS, function(id) {
  d <- dgps[[id]]
  set.seed(SEED_MASTER + 11L)
  dat <- d$draw(200000L)
  key_of <- function(df) do.call(paste, c(lapply(d$var_names, function(v) df[[v]]), sep = "-"))
  idx <- match(key_of(dat$X), key_of(d$cells))
  eta <- list(e = pmax(CLIP_LO, pmin(CLIP_HI, d$e0[idx])), m0 = d$mu0[idx], m1 = NULL)
  pi_hat <- mean(dat$A)
  sd_psi <- sqrt(mean(doubletree::psi_att(dat$Y, dat$A, d$theta0, eta, pi_hat)^2))
  do.call(rbind, lapply(N_GRID, function(n) data.frame(
    dgp = id, eps = d$eps, n = n,
    abs_bias = abs(d$pseudo$bias),
    sd_psi = sd_psi,
    z_sigma = Z_ALPHA * sd_psi / sqrt(n),
    bias_over_z_sigma = abs(d$pseudo$bias) / (Z_ALPHA * sd_psi / sqrt(n)),
    stringsAsFactors = FALSE
  )))
}))
print(knitr::kable(scale_tbl, digits = 5, format = "simple"))
cli::cli_inform(c(
  "i" = "bias_over_z_sigma < 1 means the anchor half-width is still noise-dominated (pre-plateau); > 1 means bias-dominated (plateaued).",
  "i" = "The aligned variants should span both sides across (eps, n) -- otherwise cor:width's TRANSITION is not inside the sweep."
))

# ---- 8 export -----------------------------------------------------------

out <- file.path(DIR_RESULTS, sprintf("verify_dgps_%s.rds",
                                      format(Sys.time(), "%Y%m%d-%H%M%S")))
saveRDS(list(population = pop, inexpressibility = inexpr,
             budget_probe = budget_probe, orthogonality = orth,
             confounding = conf, oracle_bias = oracle_bias,
             scale = scale_tbl,
             gates = list(n_check = N_CHECK_DRAW, n_confound = N_CONFOUND,
                          n_oracle = N_ORACLE, max_budget = MAX_BUDGET_PROBE),
             meta = cell_metadata(dgps[[1]], NA_integer_, 0L, NA_real_,
                                  rep_row_template())),
        out)
for (nm in c("population", "inexpressibility", "budget_probe", "orthogonality",
             "confounding", "oracle_bias")) {
  tab <- switch(nm, population = pop, inexpressibility = inexpr,
                budget_probe = budget_probe, orthogonality = orth,
                confounding = conf, oracle_bias = oracle_bias)
  utils::write.csv(tab, file.path(DIR_TABLES, paste0("verify_", nm, ".csv")),
                   row.names = FALSE)
}
cli::cli_alert_success("wrote {.file {out}} + 6 verification CSVs")
cli::cli_alert_success("ALL DGP PROPERTY CHECKS PASSED. Next: Rscript .../code/run_pilot.R")
