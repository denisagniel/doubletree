# ============================================================
# analyze.R
# Study: partition_recovery_clt  (doubletree, S1)
# Spec:  quality_reports/specs/2026-09-01_partition-recovery-clt.md  (§4, §6)
#
# Turns the per-cell checkpoints into the spec's metrics, figures and tables.
#
# Run:  Rscript simulations/partition_recovery_clt/code/analyze.R
#       (from the doubletree package root)
# Gate: PRC_REPS  which reps-per-cell checkpoints to read  [1000]
#
# NO FAVOURITISM (spec §6). Every table and every figure panel carries all four
# regimes (F1, ST1, ST2, ST3) side by side. There is no favourable-only headline
# table and no appendix for the stress cells.
#
# THE ATTRIBUTION TABLE IS THE POINT. Spec §3's ST3 requirement -- and the
# reason this study runs before S2 -- is that a coverage number alone cannot say
# WHY. Section 5 below crosses coverage against the two candidate causes
# (partition recovery, propensity clipping) and section 6 adds the third that the
# design surfaced: whether lambda_n = log(n)/n even permits the necessary split.
# ============================================================

source(file.path("simulations", "partition_recovery_clt", "code", "common.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(knitr)
})

REPS <- as.integer(Sys.getenv("PRC_REPS", "1000"))
RUNID <- format(Sys.time(), "%Y%m%d-%H%M%S")

# ---- 0 load every cell checkpoint ---------------------------------------

grid <- design_grid()
paths <- vapply(seq_len(nrow(grid)),
                function(i) cell_path(grid$dgp[[i]], grid$n[[i]], REPS),
                character(1))
if (!all(file.exists(paths))) {
  cli::cli_abort(c(
    "{sum(!file.exists(paths))} of {length(paths)} cell checkpoint(s) missing at PRC_REPS={REPS}.",
    stats::setNames(basename(paths[!file.exists(paths)]),
                    rep("*", sum(!file.exists(paths)))),
    i = "Run run_sweep.R first."
  ))
}
bundles <- lapply(paths, readRDS)
res <- do.call(rbind, lapply(bundles, `[[`, "results"))
pop <- unique(do.call(rbind, lapply(bundles, `[[`, "dgp_summary")))
rownames(pop) <- NULL
meta1 <- bundles[[1]]$meta

# `pop` is built by unique()-ing 18 identical-per-DGP summaries down to one row
# per (dgp, nuisance). If any cell disagreed -- e.g. a checkpoint written before a
# DGP parameter changed -- unique() would leave extra rows and every subsequent
# `pop$delta_kl[pop$dgp == id & pop$nuisance == j]` lookup would silently become
# length > 1 and recycle. Asserted rather than assumed.
if (nrow(pop) != 2L * length(DGP_IDS) || anyDuplicated(pop[, c("dgp", "nuisance")])) {
  print(pop[, c("dgp", "nuisance", "delta_sq", "delta_kl", "card_T")])
  cli::cli_abort(c(
    "Expected exactly {2 * length(DGP_IDS)} unique (dgp, nuisance) population rows, got {nrow(pop)}.",
    i = "Cell checkpoints disagree about a DGP's population structure -- some were written against different DGP parameters. Delete results/ and re-run the sweep."
  ))
}

#' Look up one population quantity for one (dgp, nuisance), exactly one value
pop_val <- function(field, id, nuisance) {
  v <- pop[[field]][pop$dgp == id & pop$nuisance == nuisance]
  if (length(v) != 1L) {
    cli::cli_abort("pop lookup {field}[{id}, {nuisance}] returned {length(v)} values, expected 1.")
  }
  v
}

res$regime <- factor(res$regime, levels = REGIME_LEVELS)
res$err <- res$theta - res$theta0

cli::cli_h1("partition_recovery_clt -- analysis (runid {RUNID})")
cli::cli_inform(c(
  "*" = "{nrow(res)} replications over {nrow(grid)} cells at {REPS} reps/cell",
  "*" = "leaf_budget={meta1$leaf_budget}, m_n={meta1$m_n}, lambda_n={meta1$lambda_n_rate}, propensity_loss={meta1$propensity_loss}",
  "*" = "seed_master={meta1$seed_master}; {meta1$r_version}; optimaltrees {meta1$optimaltrees_version}"
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
# The margin was computed over T_Lbar; a fit outside it invalidates the bound
# comparison in section 4, so this is a gate, not a note.
if (!all(res$in_class_e) || !all(res$in_class_mu)) {
  cli::cli_abort("Fitted partitions outside the enumerated T_Lbar: e={sum(!res$in_class_e)}, mu={sum(!res$in_class_mu)}.")
}
cli::cli_alert_success("0 failures; all fitted partitions inside T_Lbar.")

# ---- 2 primary metrics, all regimes in one table -------------------------

mcse_prop <- function(p, n) sqrt(p * (1 - p) / n)

by_cell <- function(f) {
  parts <- split(res, list(res$dgp, res$n), drop = TRUE)
  out <- do.call(rbind, lapply(parts, f))
  out <- out[order(match(out$regime, REGIME_LEVELS), out$n), ]
  rownames(out) <- NULL
  out
}

primary <- by_cell(function(z) {
  R <- nrow(z)
  cov <- mean(z$covered)
  re <- mean(z$recovered_e)
  rm <- mean(z$recovered_mu)
  rb <- mean(z$recovered_both)
  data.frame(
    regime = as.character(z$regime[[1]]), dgp = z$dgp[[1]], n = z$n[[1]], reps = R,
    theta0 = z$theta0[[1]],
    bias = mean(z$err), bias_mcse = stats::sd(z$err) / sqrt(R),
    rmse = sqrt(mean(z$err^2)), sd_theta = stats::sd(z$theta),
    coverage = cov, cov_mcse = mcse_prop(cov, R),
    ci_width = mean(z$ci_width), mean_sigma = mean(z$sigma),
    # Condition (P1), the enumeration-based indicator: reported ALONGSIDE
    # coverage in the same table, per spec §4, never instead of it.
    rec_e = re, rec_e_mcse = mcse_prop(re, R),
    rec_mu = rm, rec_mu_mcse = mcse_prop(rm, R),
    rec_both = rb, rec_both_mcse = mcse_prop(rb, R),
    clip_frac = mean(z$clip_frac),
    pct_reps_any_clip = mean(z$clip_frac > 0),
    stringsAsFactors = FALSE
  )
})

# ---- 3 lem:uniform's empirical signature ---------------------------------

# "coverage holds even though tauhat varies across replications WITHIN S_j".
# Spec §7 says to check this directly rather than assume it, so: how many
# DISTINCT sufficient partitions does the fit actually visit, and does coverage
# hold conditional on recovery?
uniformity <- by_cell(function(z) {
  ok <- z[z$recovered_both, , drop = FALSE]
  data.frame(
    regime = as.character(z$regime[[1]]), dgp = z$dgp[[1]], n = z$n[[1]],
    reps_recovered_both = nrow(ok),
    n_distinct_keys_e = length(unique(z$key_e[z$recovered_e])),
    n_distinct_keys_mu = length(unique(z$key_mu[z$recovered_mu])),
    mean_leaves_grid_e = mean(z$n_leaves_grid_e),
    mean_leaves_grid_mu = mean(z$n_leaves_grid_mu),
    # Refinement rate: fraction of RECOVERED fits that used strictly more leaves
    # than the minimum sufficient partition needs. Every one of these is a case a
    # naive same-partition / same-leaf-count check would have scored as a FAILURE.
    pct_recovered_are_refinements_e =
      if (any(z$recovered_e)) {
        mean(z$n_leaves_grid_e[z$recovered_e] >
               pop_val("min_sufficient_leaves", z$dgp[[1]], "e"))
      } else NA_real_,
    pct_recovered_are_refinements_mu =
      if (any(z$recovered_mu)) {
        mean(z$n_leaves_grid_mu[z$recovered_mu] >
               pop_val("min_sufficient_leaves", z$dgp[[1]], "mu"))
      } else NA_real_,
    coverage_given_recovered = if (nrow(ok) > 0L) mean(ok$covered) else NA_real_,
    stringsAsFactors = FALSE
  )
})

# ---- 4 ST2: realized vs theoretical selection-failure rate ---------------

# prop:selection-rate: Pr(tauhat notin S_j) <= |T_Lbar| exp(-c_1 n Delta_j^2).
# c_1 is NOT identified by the proposition ("depending only on the outcome bound,
# c, and the positive leaf masses"), so the bound is reported at the STATED
# convention c_1 = 1 and, per spec §4, NOT fitted post hoc. The checkable content
# is the decay SHAPE, which `implied_c1` exposes: if the bound's form is right,
# implied_c1 is positive and does not drift toward 0 as n grows.
rate_cmp <- do.call(rbind, lapply(c("e", "mu"), function(j) {
  do.call(rbind, lapply(split(res, list(res$dgp, res$n), drop = TRUE), function(z) {
    id <- z$dgp[[1]]
    delta <- pop_val("delta_kl", id, j)
    card_T <- pop_val("card_T", id, j)
    rec <- if (j == "e") z$recovered_e else z$recovered_mu
    rate <- 1 - mean(rec)
    data.frame(
      regime = as.character(z$regime[[1]]), dgp = id, nuisance = j, n = z$n[[1]],
      delta_kl = delta, card_T = card_T,
      lambda_n = log(z$n[[1]]) / z$n[[1]],
      penalised_opt_in_S = log(z$n[[1]]) / z$n[[1]] < delta,
      emp_failure_rate = rate,
      emp_mcse = mcse_prop(rate, nrow(z)),
      bound_c1_1 = selection_bound(z$n[[1]], delta, card_T, c1 = 1),
      n_delta2 = z$n[[1]] * delta^2,
      implied_c1 = implied_c1(rate, z$n[[1]], delta, card_T),
      stringsAsFactors = FALSE
    )
  }))
}))
rate_cmp <- rate_cmp[order(rate_cmp$nuisance, match(rate_cmp$regime, REGIME_LEVELS),
                           rate_cmp$n), ]
rownames(rate_cmp) <- NULL

# ---- 5 the attribution table: WHY is coverage what it is? ----------------

# Spec §3's ST3 requirement, generalised to every regime so the comparison is
# meaningful: cross coverage against the two candidate causes. A coverage
# shortfall in the "recovered & unclipped" row cannot be blamed on either.
attribution <- by_cell(function(z) {
  cell <- function(rec, clip) {
    s <- z[z$recovered_both == rec & (z$clip_frac > 0) == clip, , drop = FALSE]
    if (nrow(s) == 0L) return(c(n = 0, cov = NA_real_, bias = NA_real_))
    c(n = nrow(s), cov = mean(s$covered), bias = mean(s$err))
  }
  a <- cell(TRUE, FALSE); b <- cell(TRUE, TRUE)
  c_ <- cell(FALSE, FALSE); d <- cell(FALSE, TRUE)
  data.frame(
    regime = as.character(z$regime[[1]]), dgp = z$dgp[[1]], n = z$n[[1]],
    coverage_all = mean(z$covered),
    n_rec_noclip = a[["n"]], cov_rec_noclip = a[["cov"]], bias_rec_noclip = a[["bias"]],
    n_rec_clip = b[["n"]], cov_rec_clip = b[["cov"]],
    n_fail_noclip = c_[["n"]], cov_fail_noclip = c_[["cov"]], bias_fail_noclip = c_[["bias"]],
    n_fail_clip = d[["n"]], cov_fail_clip = d[["cov"]],
    stringsAsFactors = FALSE
  )
})

# ---- 6 figures -----------------------------------------------------------

fig_long <- do.call(rbind, lapply(
  c("bias", "rmse", "coverage", "rec_e", "rec_mu", "ci_width"),
  function(m) data.frame(
    regime = primary$regime, dgp = primary$dgp, n = primary$n,
    metric = m, value = primary[[m]],
    mcse = switch(m,
                  bias = primary$bias_mcse,
                  coverage = primary$cov_mcse,
                  rec_e = primary$rec_e_mcse,
                  rec_mu = primary$rec_mu_mcse,
                  NA_real_),
    stringsAsFactors = FALSE
  )
))
metric_labels <- c(bias = "ATT bias", rmse = "ATT RMSE",
                   coverage = "95% CI coverage",
                   rec_e = "Pr(tau_e in S_e)", rec_mu = "Pr(tau_mu in S_mu)",
                   ci_width = "mean CI width")
fig_long$metric <- factor(metric_labels[fig_long$metric],
                          levels = unname(metric_labels))
fig_long$regime <- factor(fig_long$regime, levels = REGIME_LEVELS)

ref <- data.frame(
  metric = factor(metric_labels[c("bias", "coverage", "rec_e", "rec_mu")],
                  levels = levels(fig_long$metric)),
  yintercept = c(0, 0.95, 1, 1)
)

fig <- ggplot(fig_long, aes(x = n, y = value, colour = regime, shape = regime)) +
  geom_hline(data = ref, aes(yintercept = yintercept), linetype = "dashed",
             colour = "grey45", inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = value - 1.96 * mcse, ymax = value + 1.96 * mcse),
                width = 0.06, na.rm = TRUE) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2.2) +
  scale_x_log10(breaks = N_GRID, labels = N_GRID) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c(F1 = "#1B6CA8", ST1 = "#4C9F70",
                                 ST2 = "#C4471C", ST3 = "#7B3294")) +
  labs(
    title = "Partition recovery and the CLT under grid-exact sparsity (S1)",
    subtitle = paste0(
      "All four regimes side by side. F1/ST1 share the wide-margin DGP (ST1 = n < 2000); ",
      "ST2 is the weak-margin near-tie; ST3 is weak overlap.\n",
      "leaf_budget = ", meta1$leaf_budget, ", lambda_n = log(n)/n, ",
      REPS, " reps/cell. Bars are +/- 1.96 MC se."
    ),
    x = "n (log scale)", y = NULL, colour = NULL, shape = NULL,
    caption = paste0("runid ", RUNID, " | doubletree ", meta1$doubletree_sha,
                     " | optimaltrees ", meta1$optimaltrees_sha)
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

fig_base <- file.path(DIR_FIGURES, paste0("coverage_bias_recovery_", RUNID))
ggsave(paste0(fig_base, ".png"), fig, width = 10, height = 8, dpi = 200, bg = "white")
ggsave(paste0(fig_base, ".pdf"), fig, width = 10, height = 8, device = "pdf")

# Selection-failure decay against the bound's own scale.
rc <- rate_cmp
rc$regime <- factor(rc$regime, levels = REGIME_LEVELS)
rc$nuisance <- factor(rc$nuisance, levels = c("e", "mu"),
                      labels = c("propensity e", "control outcome mu"))
fig2 <- ggplot(rc, aes(x = n, y = pmax(emp_failure_rate, 1e-4),
                       colour = regime, shape = regime)) +
  geom_hline(yintercept = 1e-4, linetype = "dotted", colour = "grey60") +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2.2) +
  geom_line(aes(y = pmax(bound_c1_1, 1e-4)), linetype = "dashed",
            linewidth = 0.35, show.legend = FALSE) +
  scale_x_log10(breaks = N_GRID, labels = N_GRID) +
  scale_y_log10() +
  facet_wrap(~nuisance) +
  scale_colour_manual(values = c(F1 = "#1B6CA8", ST1 = "#4C9F70",
                                 ST2 = "#C4471C", ST3 = "#7B3294")) +
  labs(
    title = "Empirical selection-failure rate vs prop:selection-rate's bound",
    subtitle = paste0(
      "Solid = empirical Pr(tau_j not in S_j) (floored at 1e-4 so exact zeros are visible on a log scale).\n",
      "Dashed = |T_Lbar| exp(-c_1 n Delta_j^2) at the STATED convention c_1 = 1; c_1 is not identified by the\n",
      "proposition and is deliberately NOT fitted to the data (spec §4)."
    ),
    x = "n (log scale)", y = "failure rate (log scale)",
    colour = NULL, shape = NULL,
    caption = paste0("runid ", RUNID)
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

fig2_base <- file.path(DIR_FIGURES, paste0("selection_failure_rate_", RUNID))
ggsave(paste0(fig2_base, ".png"), fig2, width = 10, height = 5.5, dpi = 200, bg = "white")
ggsave(paste0(fig2_base, ".pdf"), fig2, width = 10, height = 5.5, device = "pdf")

# ---- 7 export ------------------------------------------------------------

w <- function(x, nm) {
  utils::write.csv(x, file.path(DIR_TABLES, paste0(nm, "_", RUNID, ".csv")),
                   row.names = FALSE)
}
w(pop, "population")
w(primary, "primary")
w(uniformity, "uniformity")
w(rate_cmp, "selection_rate")
w(attribution, "attribution")

md <- c(
  paste0("# partition_recovery_clt (S1) -- run ", RUNID),
  "",
  paste0("- spec: `", meta1$spec, "`"),
  paste0("- leaf_budget = ", meta1$leaf_budget, ", m_n = ", meta1$m_n,
         ", lambda_n = ", meta1$lambda_n_rate, ", outcome_type = ",
         meta1$outcome_type, ", propensity_loss = ", meta1$propensity_loss),
  paste0("- ", REPS, " reps per cell x ", nrow(grid), " cells = ", nrow(res),
         " estimate_att() calls; sequential, worker_limit = ", meta1$worker_limit,
         ", no Rashomon enumeration"),
  paste0("- seed_master = ", meta1$seed_master, "; ", meta1$seed_scheme),
  paste0("- doubletree ", meta1$doubletree_sha, ", optimaltrees ",
         meta1$optimaltrees_sha, ", ", meta1$r_version),
  "",
  "## Exact population structure (closed form, no simulation)",
  "",
  knitr::kable(pop[, c("dgp", "nuisance", "theta0", "p", "card_T",
                       "n_sufficient", "min_sufficient_leaves",
                       "n_distinct_values", "delta_sq", "delta_kl")],
               digits = 6, format = "pipe"),
  "",
  "`delta_sq` is manuscript.tex's own currency (R^(j) is a squared-error risk there);",
  "`delta_kl` is the currency estimate_att() actually minimises with propensity_loss =",
  "log_loss and outcome_type = binary, so it is the margin that governs these runs.",
  "",
  "## Primary metrics -- all four regimes",
  "",
  knitr::kable(primary[, c("regime", "dgp", "n", "reps", "bias", "bias_mcse",
                           "rmse", "coverage", "cov_mcse", "ci_width",
                           "rec_e", "rec_mu", "rec_both", "clip_frac",
                           "pct_reps_any_clip")],
               digits = 4, format = "pipe"),
  "",
  "## Selection failure: empirical vs prop:selection-rate at c_1 = 1",
  "",
  knitr::kable(rate_cmp[, c("regime", "dgp", "nuisance", "n", "delta_kl",
                            "lambda_n", "penalised_opt_in_S",
                            "emp_failure_rate", "emp_mcse", "bound_c1_1",
                            "n_delta2", "implied_c1")],
               digits = 5, format = "pipe"),
  "",
  "## Attribution: coverage crossed against recovery and clipping",
  "",
  knitr::kable(attribution, digits = 4, format = "pipe"),
  "",
  "## lem:uniform's signature: within-S_j variation with valid coverage",
  "",
  knitr::kable(uniformity, digits = 4, format = "pipe")
)
md_path <- file.path(DIR_TABLES, paste0("summary_", RUNID, ".md"))
writeLines(md, md_path)

# ---- 8 console report ----------------------------------------------------

cli::cli_h1("Primary metrics (all regimes)")
print(as.data.frame(primary[, c("regime", "dgp", "n", "bias", "rmse", "coverage",
                                "cov_mcse", "rec_e", "rec_mu", "rec_both",
                                "clip_frac")]), digits = 4)
cli::cli_h1("Selection failure vs bound")
print(as.data.frame(rate_cmp[, c("regime", "nuisance", "n", "delta_kl", "lambda_n",
                                 "penalised_opt_in_S", "emp_failure_rate",
                                 "bound_c1_1", "implied_c1")]), digits = 4)
cli::cli_h1("Attribution (coverage | recovery x clipping)")
print(as.data.frame(attribution), digits = 4)
cli::cli_h1("Within-S_j variation")
print(as.data.frame(uniformity), digits = 4)

cli::cli_alert_success("wrote {.file {md_path}} + 5 CSVs + 2 figures (png & pdf)")
