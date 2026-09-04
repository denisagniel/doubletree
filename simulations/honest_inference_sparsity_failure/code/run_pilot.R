# ============================================================
# run_pilot.R
# Study: honest_inference_sparsity_failure  (doubletree, S2)
#
# SIZE THE SWEEP BEFORE COMMITTING TO IT. Spec "Sample sizes and reps" is
# explicit: "Prototype timing on 10 reps first; this is larger than either
# sibling spec (2,400 calls each) and involves two estimator calls per rep rather
# than one, so expect roughly 2-4x their wall time -- budget accordingly, and
# consider dropping the n = 8000 cell to a smaller R if timing is prohibitive,
# stating the reduced R explicitly rather than silently."
#
# Run:  Rscript simulations/honest_inference_sparsity_failure/code/run_pilot.R
#       (from the doubletree package root)
#
# Gates (all defaulted; a bare invocation runs both stages):
#   HIS_PILOT_STAGE  "big" | "grid" | "all"                        [all]
#   HIS_PILOT_REPS   reps per cell in the GRID stage                 [3]
#   HIS_BIG_REPS     reps in the BIG stage                          [10]
#   HIS_BIG_DGP      which variant the BIG stage times          [eps_sev]
#   HIS_TARGET_REPS  the R the projection extrapolates to          [300]
#
# STAGE ORDER IS DELIBERATE. The LARGEST cell runs FIRST (n = 8000), not last:
# the whole point of a timing prototype is to find out what the expensive end
# costs, and a pilot that walks up from n = 500 spends its budget learning what
# was never in doubt. estimate_att_crossfit() with cv_regularization = TRUE runs
# nested CV (K folds x adaptive lambda search x 2 nuisances), so the anchor, not
# the estimator under test, is expected to dominate; the two are timed separately
# so that expectation is checked rather than assumed.
# ============================================================

source(file.path("simulations", "honest_inference_sparsity_failure", "code",
                 "common.R"))

suppressPackageStartupMessages(library(knitr))

STAGE       <- Sys.getenv("HIS_PILOT_STAGE", "all")
PILOT_REPS  <- as.integer(Sys.getenv("HIS_PILOT_REPS", "3"))
BIG_REPS    <- as.integer(Sys.getenv("HIS_BIG_REPS", "10"))
BIG_DGP     <- Sys.getenv("HIS_BIG_DGP", "eps_sev")
TARGET_REPS <- as.integer(Sys.getenv("HIS_TARGET_REPS", "300"))
if (!STAGE %in% c("big", "grid", "all")) {
  cli::cli_abort("HIS_PILOT_STAGE must be one of big / grid / all, got {.val {STAGE}}.")
}
if (!BIG_DGP %in% DGP_IDS) {
  cli::cli_abort("HIS_BIG_DGP={BIG_DGP} is not one of {.val {DGP_IDS}}.")
}
cli::cli_inform(c(
  "*" = "gates: HIS_PILOT_STAGE={STAGE} HIS_PILOT_REPS={PILOT_REPS} HIS_BIG_REPS={BIG_REPS} HIS_BIG_DGP={BIG_DGP} HIS_TARGET_REPS={TARGET_REPS}"
))

set.seed(SEED_MASTER)
N_BIG <- max(N_GRID)
RUNID <- format(Sys.time(), "%Y%m%d-%H%M%S")

#' Per-cell pilot summary, including everything Arm 0 needs to be read early
pilot_summary <- function(z) {
  data.frame(
    dgp = z$dgp[[1]], eps = z$eps[[1]], n = z$n[[1]], reps = nrow(z),
    # ARM 0 first, because nothing downstream is interpretable without it.
    bias_cf = mean(z$theta_cf - z$theta0), cov_cf = mean(z$covered_cf),
    cf_leaves_e = mean(z$cf_leaves_e_mean), cf_leaves_mu = mean(z$cf_leaves_mu_mean),
    cf_suff_e = mean(z$cf_suff_e_frac), cf_suff_mu = mean(z$cf_suff_mu_frac),
    # the method under test and the three intervals
    bias_full = mean(z$theta_full - z$theta0),
    cov_wald = mean(z$covered), cov_anchor = mean(z$covered_anchor),
    cov_honest = mean(z$covered_honest),
    w_wald = mean(z$ci_width), w_anchor = mean(z$ci_anchor_width),
    w_honest = mean(z$ci_honest_width),
    abs_delta = mean(z$abs_delta), rej = mean(z$fid_reject),
    at_pseudo_e = mean(z$at_pseudo_e), at_pseudo_mu = mean(z$at_pseudo_mu),
    leaves_e = mean(z$n_leaves_e), leaves_m0 = mean(z$n_leaves_m0),
    clip_frac = mean(z$clip_frac), clip_frac_cf = mean(z$clip_frac_cf),
    secs_full = mean(z$secs_full), secs_cf = mean(z$secs_cf),
    secs_rep = mean(z$secs_full + z$secs_cf),
    stringsAsFactors = FALSE
  )
}

big <- NULL
grid_pilot <- NULL

# ---- 1 BIG stage: the most expensive cell, first -------------------------

if (STAGE %in% c("big", "all")) {
  cli::cli_h1("1. Largest cell FIRST: {BIG_DGP} at n = {N_BIG}, {BIG_REPS} rep(s)")
  d <- make_dgp(BIG_DGP)
  cli::cli_inform(c(
    "*" = "{d$label}",
    "*" = "delta_e = {signif(d$class_e$delta_sq, 5)} (sq) / {signif(d$class_e$delta_kl, 5)} (kl); delta_mu = {signif(d$class_mu$delta_sq, 5)} / {signif(d$class_mu$delta_kl, 5)}",
    "*" = "exact bias at the pseudo-true limits = {signif(d$pseudo$bias, 5)}"
  ))
  t0 <- Sys.time()
  rows <- vector("list", BIG_REPS)
  for (r in seq_len(BIG_REPS)) {
    # on_error = "stop": in a pilot a failure must surface at its origin.
    rows[[r]] <- run_one_rep(d, N_BIG, r, on_error = "stop")
    z <- rows[[r]]
    cli::cli_inform(paste0(
      "  rep ", r, "/", BIG_REPS,
      ": full ", signif(z$secs_full, 3), "s + cf ", signif(z$secs_cf, 3), "s",
      "  theta_full=", signif(z$theta_full, 4),
      " theta_cf=", signif(z$theta_cf, 4),
      " delta=", signif(z$delta, 4),
      " cf_leaves(e,mu)=", signif(z$cf_leaves_e_mean, 3), ",",
      signif(z$cf_leaves_mu_mean, 3),
      " cf_suff(e,mu)=", signif(z$cf_suff_e_frac, 2), ",",
      signif(z$cf_suff_mu_frac, 2)
    ))
  }
  big <- do.call(rbind, rows)
  big_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  cli::cli_h2("Timing at the expensive end")
  cli::cli_inform(c(
    "*" = "wall time {signif(big_min, 3)} min for {BIG_REPS} reps = {signif(big_min * 60 / BIG_REPS, 3)} s/rep",
    "*" = "estimate_att(): {signif(mean(big$secs_full), 3)} s/rep | estimate_att_crossfit(): {signif(mean(big$secs_cf), 3)} s/rep (anchor share {signif(100 * mean(big$secs_cf) / mean(big$secs_full + big$secs_cf), 3)}%)",
    "!" = "one n = {N_BIG} cell at R = {TARGET_REPS} projects to {signif(mean(big$secs_full + big$secs_cf) * TARGET_REPS / 60, 3)} min; all {length(DGP_IDS)} of them, {signif(length(DGP_IDS) * mean(big$secs_full + big$secs_cf) * TARGET_REPS / 60, 3)} min"
  ))
  print(knitr::kable(pilot_summary(big), digits = 4, format = "simple"))
}

# ---- 2 GRID stage: every cell, few reps ----------------------------------

if (STAGE %in% c("grid", "all")) {
  grid <- design_grid()
  cli::cli_h1("2. Pilot: {PILOT_REPS} rep(s) in each of {nrow(grid)} cells")
  dgps <- lapply(DGP_IDS, make_dgp)
  names(dgps) <- DGP_IDS
  t_all <- Sys.time()
  out <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    rows <- lapply(seq_len(PILOT_REPS),
                   function(r) run_one_rep(dgps[[g$dgp]], g$n, r, on_error = "stop"))
    out[[i]] <- do.call(rbind, rows)
    gc(full = TRUE, verbose = FALSE)
    z <- out[[i]]
    cli::cli_inform(paste0(
      "  ", g$dgp, " n=", g$n,
      "  cf_cov=", signif(mean(z$covered_cf), 3),
      " cf_lv=", signif(mean(z$cf_leaves_e_mean), 3), "/",
      signif(mean(z$cf_leaves_mu_mean), 3),
      " | wald=", signif(mean(z$covered), 3),
      " anch=", signif(mean(z$covered_anchor), 3),
      " hon=", signif(mean(z$covered_honest), 3),
      " rej=", signif(mean(z$fid_reject), 3),
      "  ", signif(mean(z$secs_full + z$secs_cf), 3), " s/rep"
    ))
  }
  grid_pilot <- do.call(rbind, out)
  grid_min <- as.numeric(difftime(Sys.time(), t_all, units = "mins"))

  cli::cli_h1("3. Per-cell pilot summary (Arm 0 columns first)")
  summ <- do.call(rbind, lapply(split(grid_pilot, list(grid_pilot$dgp, grid_pilot$n),
                                      drop = TRUE), pilot_summary))
  summ <- summ[order(match(summ$dgp, DGP_IDS), summ$n), ]
  rownames(summ) <- NULL
  print(knitr::kable(
    summ[, c("dgp", "eps", "n", "bias_cf", "cov_cf", "cf_leaves_e", "cf_leaves_mu",
             "cf_suff_e", "cf_suff_mu", "bias_full", "cov_wald", "cov_anchor",
             "cov_honest", "rej")],
    digits = 4, format = "simple"
  ))
  print(knitr::kable(
    summ[, c("dgp", "eps", "n", "w_wald", "w_anchor", "w_honest", "abs_delta",
             "at_pseudo_e", "at_pseudo_mu", "leaves_e", "leaves_m0", "clip_frac",
             "clip_frac_cf", "secs_full", "secs_cf", "secs_rep")],
    digits = 4, format = "simple"
  ))

  cli::cli_h1("4. Timing projection to R = {TARGET_REPS}")
  proj <- summ[, c("dgp", "eps", "n", "secs_full", "secs_cf", "secs_rep")]
  proj$projected_min <- proj$secs_rep * TARGET_REPS / 60
  print(knitr::kable(proj, digits = 4, format = "simple"))
  total_h <- sum(proj$projected_min) / 60
  cli::cli_inform(c(
    "*" = "pilot wall time {signif(grid_min, 3)} min for {nrow(grid_pilot)} reps ({2 * nrow(grid_pilot)} estimator calls)",
    "!" = "PROJECTED FULL GRID at R = {TARGET_REPS}: {signif(total_h, 3)} hours ({nrow(grid) * TARGET_REPS} rep-pairs)"
  ))
  # The n = 8000 row is broken out because it is the one the spec explicitly
  # permits reducing, and because it is where the projection is most sensitive.
  big_h <- sum(proj$projected_min[proj$n == N_BIG]) / 60
  cli::cli_inform(c(
    "*" = "of which n = {N_BIG}: {signif(big_h, 3)} h ({signif(100 * big_h / total_h, 3)}%); n < {N_BIG}: {signif(total_h - big_h, 3)} h"
  ))
  if (total_h > 3) {
    cli::cli_alert_danger(
      "Projected {signif(total_h, 3)} h exceeds a 2-3 h budget. Per the spec, reduce R at n = {N_BIG} ONLY and state the reduced R explicitly; do not reduce R elsewhere."
    )
  } else {
    cli::cli_alert_success(
      "Projected {signif(total_h, 3)} h is within a 2-3 h budget; the full sweep can proceed at R = {TARGET_REPS} everywhere."
    )
  }
}

# ---- 3 integrity and export ----------------------------------------------

all_rows <- do.call(rbind, Filter(Negate(is.null), list(big, grid_pilot)))
if (any(!is.na(all_rows$error_message))) {
  cli::cli_abort("Pilot recorded {sum(!is.na(all_rows$error_message))} failure(s) despite on_error = 'stop'.")
}
wm <- unique(stats::na.omit(all_rows$warning_messages))
if (length(wm)) {
  cli::cli_inform(c(
    "i" = "{sum(all_rows$n_warnings)} muffled warning(s), {length(wm)} distinct text(s):",
    stats::setNames(substr(wm, 1, 200), rep("*", length(wm)))
  ))
}
# The deterministic thm:anchor implication (see run_one_rep()): the anchor covers
# whenever the anchor's own Wald interval does. run_one_rep() already errors on a
# violation; re-checked here over the pooled pilot so the guarantee is visible in
# the pilot's own output rather than only latent in the harness.
viol <- sum(all_rows$covered_cf & !all_rows$covered_anchor)
if (viol > 0L) cli::cli_abort("{viol} replication(s) violate coverage(anchor) >= coverage(crossfit Wald).")
cli::cli_alert_success(
  "0 failures; coverage(ci_95_anchor) >= coverage(anchor Wald) held in all {nrow(all_rows)} pilot replications."
)

out_path <- file.path(DIR_RESULTS, sprintf("pilot_%s.rds", RUNID))
saveRDS(list(big = big, grid = grid_pilot,
             gates = list(stage = STAGE, pilot_reps = PILOT_REPS,
                          big_reps = BIG_REPS, big_dgp = BIG_DGP,
                          target_reps = TARGET_REPS),
             meta = cell_metadata(make_dgp(BIG_DGP), N_BIG, BIG_REPS,
                                  log(N_BIG) / N_BIG, all_rows)),
        out_path)
cli::cli_alert_success("wrote {.file {out_path}}")
