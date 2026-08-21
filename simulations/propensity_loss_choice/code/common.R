# ============================================================
# Propensity-tree loss ablation -- shared harness
#
# Stem:    2026-08-20_propensity-loss-choice
# Spec:    global-scholars/quality_reports/specs/2026-08-20_propensity-loss-choice.md
# Purpose: Paired ablation of `propensity_loss` in doubletree::estimate_att().
#          Everything is held identical between the two arms except the loss
#          used to FIT/SELECT the propensity tree.
# Inputs:  the two pre-existing DGPs named in the spec (never redefined here)
# Outputs: none (sourced by run_pilot.R / run_sim.R / analyze.R)
#
# Run from the doubletree package root.
# ============================================================

# ---- 0 setup ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(cli)
  library(digest)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

PKG_ROOT <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(PKG_ROOT, "DESCRIPTION"))) {
  cli::cli_abort(c(
    "No {.file DESCRIPTION} in {.path {PKG_ROOT}}.",
    i = "Run this study from the {.pkg doubletree} package root."
  ))
}
desc_pkg <- unname(read.dcf(file.path(PKG_ROOT, "DESCRIPTION"))[1, "Package"])
if (!identical(desc_pkg, "doubletree")) {
  cli::cli_abort(c(
    "Working directory is package {desc_pkg}, not doubletree.",
    i = "Run this study from the {.pkg doubletree} package root."
  ))
}

STUDY_DIR   <- file.path(PKG_ROOT, "simulations", "propensity_loss_choice")
DIR_CODE    <- file.path(STUDY_DIR, "code")
DIR_RESULTS <- file.path(STUDY_DIR, "results")
DIR_FIGURES <- file.path(STUDY_DIR, "figures")
DIR_TABLES  <- file.path(STUDY_DIR, "tables")
for (.d in c(DIR_RESULTS, DIR_FIGURES, DIR_TABLES)) {
  if (!dir.exists(.d)) dir.create(.d, recursive = TRUE)
}

# The installed doubletree in the user library predates the 2026-08-21
# restructuring (its estimate_att() is the old cross-fit signature and has no
# `propensity_loss`), so the SOURCE tree is loaded explicitly. optimaltrees is
# normally already dev-loaded by doubletree/.Rprofile; load it here too so the
# scripts behave identically under `Rscript --vanilla` and in worker processes.
.optimaltrees_src <- file.path(dirname(PKG_ROOT), "optimaltrees")
if (!isNamespaceLoaded("optimaltrees")) {
  if (dir.exists(.optimaltrees_src)) {
    pkgload::load_all(.optimaltrees_src, quiet = TRUE, export_all = FALSE)
  } else {
    library(optimaltrees)
  }
}
pkgload::load_all(PKG_ROOT, quiet = TRUE, export_all = FALSE)

if (!"propensity_loss" %in% names(formals(doubletree::estimate_att))) {
  cli::cli_abort(c(
    "Loaded {.fun doubletree::estimate_att} has no {.arg propensity_loss} argument.",
    x = "This ablation cannot run against that version.",
    i = "A stale installed copy is shadowing the source tree."
  ))
}

# ---- 1 fixed design constants ----------------------------------------------

# Master seed (spec "Reproducibility"): set once at the top of each run script
# AND folded into rep_seed(), so the whole study is a function of this integer.
SEED_MASTER <- 20260820L

# leaf_budget (Lbar of theory.tex ass:budget). REQUIRED by estimate_att(), no
# default, and deliberately NOT swept (spec: "choose one value and state it").
# 4 is the smallest budget that represents both nuisances of both DGPs exactly:
#   generate_dgp_simple      e(x) = expit(-0.5 + .3 x1 + .3 x2)  -> 4 cells (x1,x2)
#                            mu0(x) = 0.2 + .15 x1 + .15 x3      -> 4 cells (x1,x3)
#   generate_dgp_weak_overlap e(x) = plogis(-2.5 + 4 X1 + 2.5 X2) -> 4 cells (X1,X2)
#                            p0(x) = plogis(-0.3 + .5 X1 + .4 X2) -> 4 cells (X1,X2)
# so ass:sparsity holds with delta_e = delta_mu = 0 at Lbar = 4. run_pilot.R
# verifies this empirically (certified_* / n_leaves_*) rather than assuming it.
LEAF_BUDGET <- 4L

M_N          <- 1L    # eq:feasible leaf-mass floor; estimate_att()'s own default
LAMBDA_N     <- NULL  # NULL -> log(n)/n inside estimate_att() (prop:parsimony rate)
OUTCOME_TYPE <- "binary"

N_GRID  <- c(500L, 2000L)
DGP_IDS <- c("simple", "weak_overlap")

# Arm label -> propensity_loss. The ONLY thing that varies between arms.
ARMS <- c(L = "log_loss", S = "squared_error")

DGP_LABELS <- c(
  simple       = "generate_dgp_simple (favorable)",
  weak_overlap = "generate_dgp_weak_overlap (stress)"
)
ARM_LABELS <- c(
  L = "Arm L: log_loss",
  S = "Arm S: squared_error"
)

# ---- 2 DGPs: sourced from their existing homes, never redefined -------------

DGP_SRC_SIMPLE <- file.path(
  PKG_ROOT, "simulations", "archive", "six_approach_comparison_2026-06",
  "study", "code", "dgps.R"
)
DGP_SRC_STRESS <- file.path(PKG_ROOT, "simulations", "dgps", "dgps_stress.R")

for (dgp_src in c(DGP_SRC_SIMPLE, DGP_SRC_STRESS)) {
  if (!file.exists(dgp_src)) {
    cli::cli_abort("DGP source {.file {dgp_src}} not found.")
  }
}
source(DGP_SRC_SIMPLE, local = FALSE)
source(DGP_SRC_STRESS, local = FALSE)

for (dgp_fn in c("generate_dgp_simple", "generate_dgp_weak_overlap")) {
  if (!is.function(get0(dgp_fn))) {
    cli::cli_abort("{.fun {dgp_fn}} was not defined by the sourced DGP files.")
  }
}

#' Deterministic per-replicate seed
#'
#' A pure function of (SEED_MASTER, dgp, n, rep). This is what makes the design
#' PAIRED: both arms of a replicate are fit to the byte-identical dataset, so
#' the same-partition diagnostic is meaningful and the bias/RMSE contrast is
#' free of between-arm sampling noise.
rep_seed <- function(dgp, n, rep) {
  key <- paste(SEED_MASTER, dgp, n, rep, sep = "|")
  # 7 hex digits < 2^28: always a valid positive R integer.
  strtoi(substr(digest::digest(key, algo = "xxhash64"), 1, 7), base = 16L)
}

#' Generate one replicate's dataset from an existing DGP
#'
#' The DGP bodies are used exactly as they ship; only `n` (and the RNG state)
#' is set here. `generate_dgp_simple()` takes no seed argument, so seeding is
#' done externally for BOTH DGPs for uniformity.
generate_data <- function(dgp, n, rep) {
  seed <- rep_seed(dgp, n, rep)
  set.seed(seed)
  d <- switch(
    dgp,
    simple       = generate_dgp_simple(n = n),
    weak_overlap = generate_dgp_weak_overlap(n = n, tau = 0.10, seed = NULL),
    cli::cli_abort("Unknown dgp {.val {dgp}}; expected one of {.val {DGP_IDS}}.")
  )

  required <- c("X", "A", "Y", "true_att")
  absent <- setdiff(required, names(d))
  if (length(absent)) {
    cli::cli_abort("DGP {.val {dgp}} returned no {.field {absent}}.")
  }
  # The two DGP files disagree on the name of the true propensity vector
  # (`e_true` in dgps.R, `true_e` in dgps_stress.R); accept either, demand one.
  e_true <- if (!is.null(d$e_true)) d$e_true else d$true_e
  if (is.null(e_true)) {
    cli::cli_abort("DGP {.val {dgp}} returned neither {.field e_true} nor {.field true_e}.")
  }

  c(d, list(dgp = dgp, n = n, rep = rep, seed = seed, e_true_vec = e_true,
            # Recorded on BOTH arms' rows so analyze.R can verify, rather than
            # assume, that the two arms were fit to byte-identical data.
            data_hash = digest::digest(list(d$X, d$A, d$Y), algo = "xxhash64")))
}

# ---- 3 estimation: the two arms --------------------------------------------

#' Fit one arm of the ablation
#'
#' Identical call in both arms apart from `propensity_loss`.
fit_arm <- function(d, propensity_loss) {
  doubletree::estimate_att(
    X               = d$X,
    A               = d$A,
    Y               = d$Y,
    leaf_budget     = LEAF_BUDGET,
    outcome_type    = OUTCOME_TYPE,
    lambda_n        = LAMBDA_N,
    m_n             = M_N,
    propensity_loss = propensity_loss,
    verbose         = FALSE
  )
}

#' Partition identity of a fitted tree, ignoring leaf values
#'
#' `optimaltrees::structure_hash()` hashes the sorted set of leaf regions, so
#' two fits sharing a partition but differing in leaf values hash equal --
#' exactly the "same-partition" notion the spec asks for.
partition_hash <- function(model) {
  optimaltrees::structure_hash(optimaltrees::extract_tree_structure(model))
}

#' One arm's results as a single tidy row
arm_row <- function(fit, arm, d) {
  e_hat <- fit$nuisance_fits$propensity
  ci <- fit$ci_95
  tibble::tibble(
    dgp               = d$dgp,
    n                 = d$n,
    rep               = d$rep,
    seed              = d$seed,
    data_hash         = d$data_hash,
    arm               = arm,
    propensity_loss   = unname(ARMS[[arm]]),
    true_att          = d$true_att,
    theta             = fit$theta,
    sigma             = fit$sigma,
    ci_lo             = ci[[1]],
    ci_hi             = ci[[2]],
    ci_width          = ci[[2]] - ci[[1]],
    covered           = ci[[1]] <= d$true_att && d$true_att <= ci[[2]],
    n_leaves_e        = fit$n_leaves_e,
    n_leaves_m0       = fit$n_leaves_m0,
    certified_e       = fit$certified_e,
    certified_m0      = fit$certified_m0,
    used_search_e     = fit$used_search_e,
    used_search_m0    = fit$used_search_m0,
    lambda_n_used     = fit$lambda_n,
    max_e_hat         = max(e_hat),
    min_e_hat         = min(e_hat),
    max_true_e        = max(d$e_true_vec),
    min_true_e        = min(d$e_true_vec),
    hash_e            = partition_hash(fit$nuisance_fits$e_model),
    hash_m0           = partition_hash(fit$nuisance_fits$m0_model),
    error_message     = NA_character_
  )
}

#' A failed replicate as rows shaped like arm_row()'s output
#'
#' The condition MESSAGE is retained (never silently replaced by NA estimates),
#' so a non-zero failure rate is diagnosable rather than merely counted.
failure_rows <- function(dgp, n, rep, seed, message) {
  purrr::map_dfr(names(ARMS), function(arm) {
    tibble::tibble(
      dgp = dgp, n = n, rep = rep, seed = seed, data_hash = NA_character_,
      arm = arm,
      propensity_loss = unname(ARMS[[arm]]),
      true_att = NA_real_, theta = NA_real_, sigma = NA_real_,
      ci_lo = NA_real_, ci_hi = NA_real_, ci_width = NA_real_,
      covered = NA, n_leaves_e = NA_integer_, n_leaves_m0 = NA_integer_,
      certified_e = NA, certified_m0 = NA,
      used_search_e = NA, used_search_m0 = NA, lambda_n_used = NA_real_,
      max_e_hat = NA_real_, min_e_hat = NA_real_,
      max_true_e = NA_real_, min_true_e = NA_real_,
      hash_e = NA_character_, hash_m0 = NA_character_,
      error_message = message
    )
  })
}

# ---- 4 one replicate = one dataset, both arms -------------------------------

#' Run one paired replicate
#'
#' @param on_error "stop" (default; the interactive/pilot posture -- a failure
#'   surfaces at its origin) or "record" (the long-run posture -- the condition
#'   object's message is stored so the whole cell still finishes and the caller
#'   can report a failure RATE, then abort on it).
run_one_rep <- function(dgp, n, rep, on_error = c("stop", "record")) {
  on_error <- match.arg(on_error)

  body <- function() {
    d <- generate_data(dgp, n, rep)
    rows <- lapply(names(ARMS), function(arm) arm_row(fit_arm(d, ARMS[[arm]]), arm, d))
    out <- dplyr::bind_rows(rows)
    # A non-finite estimate is a failure too, not a number to average.
    bad <- !is.finite(out$theta) | !is.finite(out$sigma)
    if (any(bad)) {
      cli::cli_abort(
        "Non-finite output in arm(s) {.val {out$arm[bad]}}: theta = {out$theta[bad]}, sigma = {out$sigma[bad]}."
      )
    }
    out
  }

  if (on_error == "stop") return(body())

  tryCatch(
    body(),
    error = function(cnd) {
      failure_rows(dgp, n, rep, rep_seed(dgp, n, rep), conditionMessage(cnd))
    }
  )
}

#' Run every replicate of one (dgp, n) cell
run_cell <- function(dgp, n, n_reps, on_error = "record", heartbeat = 25L) {
  rows <- vector("list", n_reps)
  for (r in seq_len(n_reps)) {
    rows[[r]] <- run_one_rep(dgp, n, r, on_error = on_error)
    if (heartbeat > 0L && r %% heartbeat == 0L) {
      cli::cli_inform("  {dgp} n={n}: {r}/{n_reps}")
    }
  }
  dplyr::bind_rows(rows)
}

# ---- 5 run metadata --------------------------------------------------------

#' Provenance stamped into every results file
run_metadata <- function(n_reps, label) {
  list(
    label            = label,
    runid            = format(Sys.time(), "%Y%m%d-%H%M%S"),
    spec             = "quality_reports/specs/2026-08-20_propensity-loss-choice.md",
    stem             = "2026-08-20_propensity-loss-choice",
    seed_master      = SEED_MASTER,
    leaf_budget      = LEAF_BUDGET,
    m_n              = M_N,
    lambda_n         = "NULL -> log(n)/n",
    outcome_type     = OUTCOME_TYPE,
    n_grid           = N_GRID,
    dgps             = DGP_IDS,
    arms             = ARMS,
    n_reps           = n_reps,
    doubletree_sha   = system2("git", c("-C", PKG_ROOT, "rev-parse", "--short", "HEAD"),
                               stdout = TRUE, stderr = FALSE),
    optimaltrees_sha = system2("git", c("-C", .optimaltrees_src, "rev-parse", "--short", "HEAD"),
                               stdout = TRUE, stderr = FALSE),
    r_version        = R.version.string,
    timestamp        = Sys.time()
  )
}
