# ============================================================
# common.R
# Study: honest_inference_sparsity_failure  (doubletree, S2)
# Spec:  quality_reports/specs/2026-08-21_honest-inference-sparsity-failure.md
#        (the 2026-09-01 changelog at the top of that file is authoritative)
#
# Shared harness: paths, package loading, the design grid, deterministic seeds,
# one replication, one cell. Sourced by verify_dgps.R / run_pilot.R / run_cell.R /
# run_sweep.R / analyze.R. Sourcing this file has no effect beyond defining
# objects and creating the study's output directories.
#
# Structure follows partition_recovery_clt/code/common.R (S1), which in turn
# follows propensity_loss_choice/code/common.R -- same source-tree-vs-installed
# guard, same deterministic seeding scheme, same checkpoint-per-cell posture,
# same worker_limit = 1L memory rule. Spec "Outputs" asks for exactly that reuse.
#
# WHAT IS DIFFERENT FROM S1, AND WHY
#
# (a) TWO estimator calls per replication, not one: estimate_att() (the method
#     under test) AND estimate_att_crossfit() (the anchor). Every anchor-dependent
#     quantity in this study is a function of the PAIR, on the SAME data.
#
# (b) THREE intervals per replication, all reported side by side at every eps
#     (spec "Constitution compliance", last bullet, and the Decision rule's
#     "Regardless of outcome"): ci_95 (naive Wald, negative control),
#     ci_95_anchor (the literal thm:anchor triangle inequality), ci_95_honest
#     (the as-shipped bias-aware honest_ci(), se_delta = 0).
#
# (c) NO PACKAGE CHANGE. ci_95_anchor is computed HERE, from estimate_att()'s and
#     estimate_att_crossfit()'s existing `theta`/`sigma` returns. Spec changelog
#     item 4 dropped the original plan to add it as a return field on
#     estimate_att(): that would force the paper's sample-splitting-free primary
#     estimator to couple internally to a K-fold procedure. honest_ci() is used
#     UNMODIFIED -- the point is to test it as shipped.
#
# (d) ARM 0 is first class. estimate_att_crossfit()'s OWN bias, RMSE and Wald
#     coverage are recorded per replication and tabulated per cell BEFORE any
#     anchor-dependent number is trusted (spec's new 2026-09-01 arm). See the
#     note at `covered_anchor` below for the exact reason this is not optional.
# ============================================================

suppressPackageStartupMessages({
  library(cli)
  library(digest)
})

# ---- 0 locate the study and load the packages ---------------------------

PKG_ROOT <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(PKG_ROOT, "DESCRIPTION"))) {
  cli::cli_abort(c(
    "No {.file DESCRIPTION} in {.path {PKG_ROOT}}.",
    i = "Run this study from the {.pkg doubletree} package root."
  ))
}
.desc_pkg <- unname(read.dcf(file.path(PKG_ROOT, "DESCRIPTION"))[1, "Package"])
if (!identical(.desc_pkg, "doubletree")) {
  cli::cli_abort("Working directory is package {.desc_pkg}, not doubletree.")
}

STUDY_DIR   <- file.path(PKG_ROOT, "simulations", "honest_inference_sparsity_failure")
DIR_CODE    <- file.path(STUDY_DIR, "code")
DIR_RESULTS <- file.path(STUDY_DIR, "results")
DIR_FIGURES <- file.path(STUDY_DIR, "figures")
DIR_TABLES  <- file.path(STUDY_DIR, "tables")
for (.d in c(DIR_RESULTS, DIR_FIGURES, DIR_TABLES)) {
  if (!dir.exists(.d)) dir.create(.d, recursive = TRUE)
}

# The installed doubletree in the user library can predate the current source
# tree, so load the SOURCE explicitly -- same posture as S1's and
# propensity_loss_choice's common.R, and for the same reason.
.OPTIMALTREES_SRC <- file.path(dirname(PKG_ROOT), "optimaltrees")

# TWO PACKAGE-LOADING PATHS, ONE GATE (added 2026-09-09 for the SLURM deployment;
# simulations/honest_inference_sparsity_failure/slurm/README_O2.md documents the
# cluster side). Same gate as head_to_head_comparison/code/common.R (S5), for the
# same reason and with the same DEV default.
#
#   HIS_USE_INSTALLED unset/"0"  DEV path (default, UNCHANGED behaviour):
#       pkgload::load_all() on this source tree and on the sibling
#       ../optimaltrees source tree. Correct on the dev box, where the source is
#       what is being developed and an installed copy may be stale.
#
#   HIS_USE_INSTALLED="1"        CLUSTER path: library() against packages that
#       were R CMD INSTALLed on the cluster. Module R does not carry a working
#       pkgload/devtools dev-load, and `../optimaltrees` is a RELATIVE sibling
#       path that is only correct when the checkout happens to have both repos
#       side by side -- neither assumption is safe inside a SLURM job. The rest
#       of this file is loading-agnostic: every package call below is either
#       `doubletree::`-qualified or reached through asNamespace(), so nothing
#       else has to change.
#
# .OPTIMALTREES_SRC is still defined on BOTH paths: cell_metadata() stamps its
# git SHA, and git_sha() already records a missing/failed lookup as a reason
# string rather than collapsing it to NA.
USE_INSTALLED_PKGS <- Sys.getenv("HIS_USE_INSTALLED", "0") == "1"
if (USE_INSTALLED_PKGS) {
  suppressPackageStartupMessages({
    library(optimaltrees)
    library(doubletree)
  })
} else {
  if (!isNamespaceLoaded("optimaltrees")) {
    if (dir.exists(.OPTIMALTREES_SRC)) {
      pkgload::load_all(.OPTIMALTREES_SRC, quiet = TRUE, export_all = FALSE)
    } else {
      library(optimaltrees)
    }
  }
  pkgload::load_all(PKG_ROOT, quiet = TRUE, export_all = FALSE)
}

# Fail fast rather than at replication 200: this study needs the CURRENT
# signatures, read directly from R/ rather than assumed from the spec's prose.
# NOTE: the loop variable must not be named `.arg` -- cli >= 3.4.0 reserves
# any `{.foo}` starting with a dot exclusively for its own inline styles, so
# `{.arg {.arg}}` no longer interpolates a variable named `.arg`; it is parsed
# as an (invalid) nested style directive instead.
for (arg_name in c("leaf_budget", "propensity_loss", "outcome_type", "lambda_n", "m_n")) {
  if (!arg_name %in% names(formals(doubletree::estimate_att))) {
    cli::cli_abort(c(
      "Loaded {.fun doubletree::estimate_att} has no {.arg {arg_name}} argument.",
      i = "A stale installed copy is probably shadowing the source tree."
    ))
  }
}
for (arg_name in c("K", "outcome_type", "cv_regularization", "max_depth", "seed")) {
  if (!arg_name %in% names(formals(doubletree::estimate_att_crossfit))) {
    cli::cli_abort("Loaded {.fun doubletree::estimate_att_crossfit} has no {.arg {arg_name}} argument.")
  }
}
# honest_ci() must be the generic 5-argument function the spec describes; this
# study tests it AS SHIPPED, so its contract is asserted, not assumed.
local({
  f <- names(formals(doubletree::honest_ci))
  want <- c("theta_display", "se", "delta", "se_delta", "level")
  if (!identical(f, want)) {
    cli::cli_abort(c(
      "{.fun doubletree::honest_ci} has formals {.val {f}}, expected {.val {want}}.",
      i = "This study computes ci_95_honest from it unmodified."
    ))
  }
})

source(file.path(DIR_CODE, "dgps_delta_dial.R"))
source_s1_enumeration(PKG_ROOT)   # spec changelog item 5

# ---- 1 fixed design constants -------------------------------------------

SEED_MASTER <- 20260821L   # spec "Reproducibility"

# Lbar of ass:budget. FIXED across every DGP variant, per spec "Estimand and
# method": letting the budget grow with eps would let the tree "solve around" the
# manufactured sparsity violation and defeat the sweep. 2 is not a free choice --
# see design note (1) in dgps_delta_dial.R: the step must consume the whole budget
# or delta_j stops being a clean dial and the pseudo-true partition (hence
# cor:width's predicted plateau height) stops being unique.
LEAF_BUDGET  <- 2L
M_N          <- 1L    # eq:feasible floor; estimate_att()'s own default
LAMBDA_N     <- NULL  # NULL -> log(n)/n inside estimate_att() (prop:parsimony rate)
OUTCOME_TYPE <- "binary"
PROPENSITY_LOSS <- "log_loss"   # estimate_att()'s SHIPPED default

# Anchor settings: estimate_att_crossfit()'s shipped defaults, K = 5 per spec.
# cv_regularization = TRUE is the shipped default and is kept: the anchor's
# escape from the sparsity violation is precisely that it carries NO leaf budget
# and picks its own penalty, and weakening that would sabotage Arm 0 by
# construction.
ANCHOR_K <- 5L
ANCHOR_MAX_DEPTH <- 4L   # estimate_att_crossfit()'s shipped default

ALPHA <- 0.05
Z_ALPHA <- stats::qnorm(1 - ALPHA / 2)

# Memory safety (this project's own convention; simulations/docs/
# MEMORY_SAFE_SIMULATIONS.md and MEMORY.md [LEARN:rashomon-memory]): single best
# tree only, no Rashomon enumeration anywhere, no parallelism. worker_limit = 1L
# reaches fit_tree() through both estimators' `...`; PARALLEL_CV = FALSE
# additionally switches off cv_regularization_adaptive()'s furrr/future
# fold-parallelism, which is ON by default (it degrades to sequential when no
# future plan is set, but "no plan is set" is a property of the environment, not
# of this code, so it is disabled explicitly). BLAS/OpenMP threads are pinned to
# 1 by run_sweep.R in every child process.
WORKER_LIMIT <- 1L
PARALLEL_CV <- FALSE

# The propensity clip of ass:construct(c), applied inside both estimators. Copied
# because the package constants are unexported, then CHECKED against the package
# at load: population_bias() applies the same clip, so a package change would
# otherwise silently make the predicted plateau height wrong.
CLIP_LO <- 0.01
CLIP_HI <- 0.99
local({
  pkg_lo <- get0(".PROPENSITY_LOWER_BOUND", envir = asNamespace("doubletree"))
  pkg_hi <- get0(".PROPENSITY_UPPER_BOUND", envir = asNamespace("doubletree"))
  if (is.null(pkg_lo) || is.null(pkg_hi) ||
      !isTRUE(all.equal(c(pkg_lo, pkg_hi), c(CLIP_LO, CLIP_HI)))) {
    cli::cli_abort(c(
      "doubletree's propensity clip bounds are not ({CLIP_LO}, {CLIP_HI}).",
      i = "population_bias()'s predicted plateau height would use the wrong clip."
    ))
  }
})

# ---- 2 the design grid --------------------------------------------------

# Spec "Sample sizes and reps": three n, one more than the sibling specs, because
# cor:width's plateau TRANSITION is itself a statement across n.
N_GRID <- c(500L, 2000L, 8000L)

# The five required variants (spec "Required DGP variants"). eps values are
# chosen so the implied population bias b(eps) ~ 4.9 eps^2 spans from well BELOW
# the Wald scale at every n, through it, to well ABOVE it -- so the plateau
# transition is inside the sweep rather than only its plateaued end. The exact
# realised (delta_e, delta_mu) and the exact b are computed by enumeration in
# dgps_delta_dial.R and reported next to every metric; the eps LABEL is never the
# quantity any claim is stated in.
#
# eps_max = 0.15 is a CEILING, not a taste: above roughly 0.19 the control-outcome
# tree's pseudo-true partition flips off the X1 split (dgps_delta_dial.R design
# note (1)), which silently breaks the monotonicity of the bias in eps. An earlier
# grid ending at 0.18 did exactly that. build_dgp_dial() now refuses to build such
# a variant, so this comment cannot drift away from the code.
DGP_SPECS <- list(
  eps0       = list(eps = 0.00, mode = "aligned",
                    label = "F: eps = 0, exact sparsity (favourable)"),
  eps_small  = list(eps = 0.03, mode = "aligned",
                    label = "V1: eps = 0.03, small violation"),
  eps_mod    = list(eps = 0.09, mode = "aligned",
                    label = "V2: eps = 0.09, moderate violation"),
  eps_sev    = list(eps = 0.15, mode = "aligned",
                    label = "V3: eps = 0.15, severe violation"),
  # prop:spectest's stated BLIND SPOT, required not optional. Identical to
  # eps_sev in every respect except the coordinate h lives on, so delta_e,
  # delta_mu and eps all match eps_sev while the bilinear remainder -- hence the
  # true bias -- is exactly 0.
  orth_blind = list(eps = 0.15, mode = "orthogonal",
                    label = "B: eps = 0.15, ORTHOGONAL misspecification (b = 0)")
)
DGP_IDS <- names(DGP_SPECS)
DGP_LABELS <- vapply(DGP_SPECS, `[[`, character(1), "label")
# The blind-spot variant is singled out everywhere (spec: its row must be
# "visually distinguished", and its result must never be averaged away).
BLINDSPOT_ID <- "orth_blind"

#' The full design grid, one row per cell
design_grid <- function() {
  g <- expand.grid(n = N_GRID, dgp = DGP_IDS,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  g <- g[order(match(g$dgp, DGP_IDS), g$n), c("dgp", "n"), drop = FALSE]
  g$eps <- vapply(g$dgp, function(d) DGP_SPECS[[d]]$eps, numeric(1))
  g$residual_mode <- vapply(g$dgp, function(d) DGP_SPECS[[d]]$mode, character(1))
  g$blindspot <- g$dgp == BLINDSPOT_ID
  g$reps <- reps_for_n(g$n)
  rownames(g) <- NULL
  g
}

# Replications per cell. The spec's convention is R = 300 everywhere, with ONE
# explicit exception it grants in advance: "consider dropping the n = 8000 cell to
# a smaller R if timing is prohibitive, stating the reduced R explicitly rather
# than silently."
#
# The pilot measured it. Wall time is 1.9 s/rep at n = 500, 3.8 at n = 2000 and
# 11.1 at n = 8000, of which 98% is the ANCHOR: estimate_att_crossfit() runs
# nested CV (K = 5 outer folds x an adaptive lambda search x cv_K = 5 inner folds
# x 2 nuisances), against a single pair of tree fits for estimate_att(). At
# R = 300 uniformly that is 7.05 h sequentially, 65% of it in the three n = 8000
# cells. R = 150 there brings the sweep to ~4.7 h while leaving the coverage MC se
# at 0.018 -- ample for the study's actual questions, since the pilot's anchor and
# honest coverages sit at ~1.0, not at the 0.93-vs-0.95 margin where R = 300 would
# be needed. R is NOT reduced at any other n.
REPS_DEFAULT   <- as.integer(Sys.getenv("HIS_REPS", "300"))
REPS_AT_MAX_N  <- as.integer(Sys.getenv("HIS_REPS_MAX_N", "150"))

#' Replications for a given n (see the note above REPS_DEFAULT)
reps_for_n <- function(n) {
  ifelse(as.integer(n) >= max(N_GRID), REPS_AT_MAX_N, REPS_DEFAULT)
}

#' Build one DGP by id, with its exact enumeration quantities attached
make_dgp <- function(id) {
  if (!id %in% DGP_IDS) {
    cli::cli_abort("Unknown dgp {.val {id}}; expected one of {.val {DGP_IDS}}.")
  }
  s <- DGP_SPECS[[id]]
  build_dgp_dial(dgp_spec_delta_dial(
    eps = s$eps, residual_mode = s$mode, leaf_budget = LEAF_BUDGET,
    id = id, label = s$label
  ))
}

# ---- 3 deterministic seeds ----------------------------------------------

#' Per-replicate seed: a pure function of (SEED_MASTER, dgp, n, rep)
#'
#' Same scheme as S1's and propensity_loss_choice's common.R, so the whole sweep
#' is a function of one integer. 7 hex digits < 2^28 is always a valid R integer.
#'
#' Both estimator calls of a replication run under this ONE seed stream, in a
#' fixed order (estimate_att() then estimate_att_crossfit()), so the pair is
#' reproducible as a pair. `seed = NULL` is passed to estimate_att_crossfit() on
#' purpose: its fold assignment then draws from this stream rather than resetting
#' it, which keeps the whole replication a single reproducible sequence.
rep_seed <- function(dgp, n, rep) {
  key <- paste(SEED_MASTER, dgp, n, rep, sep = "|")
  strtoi(substr(digest::digest(key, algo = "xxhash64"), 1, 7), base = 16L)
}

# ---- 4 the three intervals and the fidelity diagnostic ------------------

#' The literal thm:anchor interval, computed in the harness
#'
#' \code{theta_full +/- (|delta| + z * sigma_anchor)} with
#' \code{delta = theta_full - theta_anchor}. NO bias-aware critical value: this is
#' the triangle inequality the theorem actually proves, and nothing else.
#'
#' ON THE \code{/sqrt(n)} IN THE SPEC'S FORMULA. The spec writes the half-width as
#' \code{|delta_hat| + z * sigma_tilde / sqrt(n)}, where \code{sigma_tilde} is the
#' asymptotic STANDARD DEVIATION, so \code{sigma_tilde / sqrt(n)} is a standard
#' ERROR. \code{estimate_att_crossfit()$sigma} is ALREADY the standard error
#' (\code{att_se()} returns \code{sqrt(mean(psi^2)/n)}), so dividing by
#' \code{sqrt(n)} a second time here would shrink the interval by a factor of
#' \code{sqrt(n)} and manufacture exactly the catastrophic undercoverage the
#' spec's Decision rule tells the analyst to suspect FIRST ("most likely: check
#' the sign/construction against the theorem statement literally before
#' suspecting the theorem"). It is therefore deliberately not divided again.
#'
#' @param theta_full,theta_anchor Point estimates.
#' @param sigma_anchor The ANCHOR's standard error.
#' @param z Normal critical value.
#' @return Numeric length-2 (lower, upper).
anchor_ci_harness <- function(theta_full, theta_anchor, sigma_anchor, z = Z_ALPHA) {
  delta <- theta_full - theta_anchor
  half <- abs(delta) + z * sigma_anchor
  c(theta_full - half, theta_full + half)
}

#' prop:spectest's fidelity diagnostic, as a standalone test statistic
#'
#' Rejects ass:sparsity when \code{|delta_hat| > z * SEhat}. Distinct object from
#' both intervals: one widens crudely, one widens via a bias-aware critical
#' value, this one rejects or does not.
#'
#' \code{SEhat = sqrt(sigma_full^2 + sigma_anchor^2)} -- the proposition quantifies
#' over "any \code{SEhat ~ n^{-1/2}}", and this is the natural already-computed
#' choice (spec's own suggestion). It is CONSERVATIVE: it ignores the positive
#' covariance between the two estimators, so it OVERSTATES \code{sd(delta_hat)}
#' and the test's size under \code{ass:sparsity} will sit well below nominal
#' \code{alpha}. That is consistent with prop:spectest's "size tending to zero"
#' and is stated here rather than discovered later. It is fixed once, in the
#' harness, and NOT tuned per DGP (spec: "state whichever choice is used
#' explicitly and do not tune it per DGP").
fidelity_test <- function(delta, sigma_full, sigma_anchor, z = Z_ALPHA) {
  se_hat <- sqrt(sigma_full^2 + sigma_anchor^2)
  list(se_hat = se_hat, stat = abs(delta) / se_hat, reject = abs(delta) > z * se_hat)
}

# ---- 5 one replication --------------------------------------------------

#' The row schema every replication returns
#'
#' Defined once so success rows and failure rows are guaranteed to have the same
#' columns: a failure must never be a shorter row that silently drops out of an
#' rbind, and must never be a plausible-looking NA estimate with no message.
rep_row_template <- function() {
  data.frame(
    dgp = NA_character_, eps = NA_real_, residual_mode = NA_character_,
    n = NA_integer_, rep = NA_integer_, seed = NA_integer_,
    theta0 = NA_real_,
    # -- the method under test: estimate_att(), Part II, full sample ----------
    theta_full = NA_real_, sigma_full = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, ci_width = NA_real_, covered = NA,
    # -- ARM 0: the anchor's OWN performance. Recorded per replication so its
    #    bias/RMSE/Wald coverage is a first-class table, not an assumption.
    theta_cf = NA_real_, sigma_cf = NA_real_,
    ci_cf_lo = NA_real_, ci_cf_hi = NA_real_, ci_cf_width = NA_real_,
    covered_cf = NA,
    # -- the bias estimate both anchor-dependent intervals are built from ------
    delta = NA_real_, abs_delta = NA_real_,
    # -- thm:anchor, literal ---------------------------------------------------
    ci_anchor_lo = NA_real_, ci_anchor_hi = NA_real_,
    ci_anchor_width = NA_real_, covered_anchor = NA,
    # -- the as-shipped honest_ci(), se_delta = 0 ------------------------------
    ci_honest_lo = NA_real_, ci_honest_hi = NA_real_,
    ci_honest_width = NA_real_, covered_honest = NA,
    honest_B = NA_real_, honest_cv = NA_real_, honest_half = NA_real_,
    # -- prop:spectest's diagnostic -------------------------------------------
    fid_se_hat = NA_real_, fid_stat = NA_real_, fid_reject = NA,
    # -- full-sample fit diagnostics ------------------------------------------
    n_leaves_e = NA_integer_, n_leaves_m0 = NA_integer_,
    certified_e = NA, certified_m0 = NA,
    used_search_e = NA, used_search_m0 = NA,
    gap_e = NA_real_, gap_m0 = NA_real_, lambda_n_used = NA_real_,
    # Grid partition of the FULL-SAMPLE trees, and whether it equals the
    # pseudo-true (population-risk argmin) partition the exact bias was computed
    # at. Without this, a bias that does not match `bias_pseudo` cannot be
    # attributed to "the fit went somewhere else" versus "the prediction is wrong".
    key_e = NA_character_, key_mu = NA_character_,
    at_pseudo_e = NA, at_pseudo_mu = NA,
    suff_e = NA, suff_mu = NA,
    # -- ARM 0 diagnostics: did the ANCHOR escape the sparsity violation? -----
    # The anchor carries no leaf budget, which is the ONLY reason it can be
    # sqrt(n)-consistent on a DGP where the constrained estimator cannot. These
    # count the leaves its fold-wise trees actually used on the covariate grid,
    # and how many fold-wise partitions are SUFFICIENT for the true nuisance.
    cf_leaves_e_mean = NA_real_, cf_leaves_mu_mean = NA_real_,
    cf_leaves_e_max = NA_integer_, cf_leaves_mu_max = NA_integer_,
    cf_suff_e_frac = NA_real_, cf_suff_mu_frac = NA_real_,
    # -- clipping (ass:construct(c)) and sample composition -------------------
    clip_frac = NA_real_, max_e_hat = NA_real_, min_e_hat = NA_real_,
    clip_frac_cf = NA_real_,
    n_treated = NA_integer_, n_control = NA_integer_,
    secs_full = NA_real_, secs_cf = NA_real_,
    # Warnings are COUNTED AND KEPT, not discarded.
    n_warnings = NA_integer_, warning_messages = NA_character_,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Grid partition and sufficiency of a fitted optimaltrees model
#'
#' Thin wrapper over S1's \code{fitted_cell_labels()} /
#' \code{is_sufficient_partition()}: the partition a fitted tree induces on the
#' 2^p covariate grid, plus whether the true nuisance is constant on its leaves.
#' Evaluated on the GRID rather than the sample, for S1's design reason (2): at
#' n = 500 several cells are typically empty in a given draw and a sample-based
#' check could call a non-sufficient partition sufficient.
grid_partition_of <- function(model, cells, gamma0) {
  lab <- fitted_cell_labels(model, cells)
  list(key = paste(lab, collapse = "."),
       n_leaves = length(unique(lab)),
       sufficient = is_sufficient_partition(lab, gamma0, tol = 1e-12))
}

#' Run one replication: draw, fit BOTH estimators, build all three intervals
#'
#' @param dgp A built DGP (\code{make_dgp()} output).
#' @param n Sample size.
#' @param rep Replication index.
#' @param on_error \code{"stop"} (pilot/interactive posture: a failure surfaces at
#'   its origin) or \code{"record"} (long-run posture: the condition message is
#'   stored so the cell finishes and a failure RATE can be reported).
#' @return A one-row data.frame matching \code{rep_row_template()}.
run_one_rep <- function(dgp, n, rep, on_error = c("stop", "record")) {
  on_error <- match.arg(on_error)
  seed <- rep_seed(dgp$id, n, rep)

  body <- function() {
    set.seed(seed)
    d <- dgp$draw(n)
    warn_msgs <- character(0)

    t0 <- Sys.time()
    fit <- withCallingHandlers(
      doubletree::estimate_att(
        X = d$X, A = d$A, Y = d$Y,
        leaf_budget = dgp$leaf_budget,
        outcome_type = OUTCOME_TYPE,
        lambda_n = LAMBDA_N,
        m_n = M_N,
        propensity_loss = PROPENSITY_LOSS,
        verbose = FALSE,
        worker_limit = WORKER_LIMIT
      ),
      warning = function(w) {
        warn_msgs <<- c(warn_msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    secs_full <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    t1 <- Sys.time()
    fit_cf <- withCallingHandlers(
      doubletree::estimate_att_crossfit(
        X = d$X, A = d$A, Y = d$Y,
        K = ANCHOR_K,
        outcome_type = OUTCOME_TYPE,
        max_depth = ANCHOR_MAX_DEPTH,
        seed = NULL,              # fold assignment draws from this rep's stream
        verbose = FALSE,
        worker_limit = WORKER_LIMIT,
        parallel = PARALLEL_CV
      ),
      warning = function(w) {
        warn_msgs <<- c(warn_msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    secs_cf <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

    if (!is.finite(fit$theta) || !is.finite(fit$sigma)) {
      stop("Non-finite estimate_att() output: theta = ", fit$theta,
           ", sigma = ", fit$sigma, call. = FALSE)
    }
    if (!is.finite(fit_cf$theta) || !is.finite(fit_cf$sigma) || fit_cf$sigma <= 0) {
      stop("Non-finite or non-positive estimate_att_crossfit() output: theta = ",
           fit_cf$theta, ", sigma = ", fit_cf$sigma, call. = FALSE)
    }
    # honest_ci() needs a strictly positive `se`; a zero would be a silent
    # degenerate interval rather than an error.
    for (fld in c("theta", "sigma", "ci_95")) {
      if (!fld %in% names(fit) || !fld %in% names(fit_cf)) {
        stop("An estimator returned no `", fld, "` field; its return contract ",
             "has changed and this study's intervals are no longer valid.",
             call. = FALSE)
      }
    }

    theta0 <- dgp$theta0
    delta <- fit$theta - fit_cf$theta

    ci <- fit$ci_95
    ci_cf <- fit_cf$ci_95
    ci_anchor <- anchor_ci_harness(fit$theta, fit_cf$theta, fit_cf$sigma)
    hon <- doubletree::honest_ci(fit$theta, fit_cf$sigma, delta,
                                 se_delta = 0, level = 1 - ALPHA)
    fid <- fidelity_test(delta, fit$sigma, fit_cf$sigma)

    covers <- function(z) z[[1]] <= theta0 && theta0 <= z[[2]]
    cov_cf <- covers(ci_cf)
    cov_anchor <- covers(ci_anchor)

    # A DETERMINISTIC implication of the construction, not a statistical
    # expectation: the anchor interval is
    #   theta_full +/- (|theta_full - theta_cf| + z sigma_cf)
    # which CONTAINS [theta_cf - z sigma_cf, theta_cf + z sigma_cf] for every
    # possible draw. So the anchor covers whenever the anchor's own Wald interval
    # covers, replication by replication, and
    #   coverage(ci_95_anchor) >= coverage(anchor's Wald interval)
    # exactly. This is the whole content of thm:anchor, and it is why Arm 0 is
    # not optional: an anchor-coverage shortfall CANNOT come from the interval's
    # construction, only from the anchor itself failing ass:rate. A violation here
    # is an implementation bug, so it stops the replication rather than being
    # averaged into a coverage number.
    if (cov_cf && !cov_anchor) {
      stop("ci_95_anchor failed to cover while the anchor's own Wald interval ",
           "did. That is impossible by the triangle inequality, so ",
           "anchor_ci_harness() is misconstructed. theta_full = ", fit$theta,
           ", theta_cf = ", fit_cf$theta, ", sigma_cf = ", fit_cf$sigma,
           ", theta0 = ", theta0, call. = FALSE)
    }

    e_hat <- fit$nuisance_fits$propensity
    e_hat_cf <- fit_cf$nuisance_fits$propensity

    gp_e <- grid_partition_of(fit$nuisance_fits$e_model, dgp$cells, dgp$e0)
    gp_mu <- grid_partition_of(fit$nuisance_fits$m0_model, dgp$cells, dgp$mu0)

    # Arm 0's escape diagnostic, fold by fold.
    cf_e <- vector("list", fit_cf$K)
    cf_mu <- vector("list", fit_cf$K)
    for (k in seq_len(fit_cf$K)) {
      cf_e[[k]] <- grid_partition_of(fit_cf$nuisance_fits[[k]]$e_model,
                                     dgp$cells, dgp$e0)
      cf_mu[[k]] <- grid_partition_of(fit_cf$nuisance_fits[[k]]$m0_model,
                                      dgp$cells, dgp$mu0)
    }
    cf_leaves_e <- vapply(cf_e, `[[`, numeric(1), "n_leaves")
    cf_leaves_mu <- vapply(cf_mu, `[[`, numeric(1), "n_leaves")

    row <- rep_row_template()
    row$dgp <- dgp$id
    row$eps <- dgp$eps
    row$residual_mode <- dgp$residual_mode
    row$n <- as.integer(n)
    row$rep <- as.integer(rep)
    row$seed <- as.integer(seed)
    row$theta0 <- theta0
    row$theta_full <- fit$theta
    row$sigma_full <- fit$sigma
    row$ci_lo <- ci[[1]]; row$ci_hi <- ci[[2]]
    row$ci_width <- ci[[2]] - ci[[1]]
    row$covered <- covers(ci)
    row$theta_cf <- fit_cf$theta
    row$sigma_cf <- fit_cf$sigma
    row$ci_cf_lo <- ci_cf[[1]]; row$ci_cf_hi <- ci_cf[[2]]
    row$ci_cf_width <- ci_cf[[2]] - ci_cf[[1]]
    row$covered_cf <- cov_cf
    row$delta <- delta
    row$abs_delta <- abs(delta)
    row$ci_anchor_lo <- ci_anchor[[1]]; row$ci_anchor_hi <- ci_anchor[[2]]
    row$ci_anchor_width <- ci_anchor[[2]] - ci_anchor[[1]]
    row$covered_anchor <- cov_anchor
    row$ci_honest_lo <- hon$ci[[1]]; row$ci_honest_hi <- hon$ci[[2]]
    row$ci_honest_width <- hon$ci[[2]] - hon$ci[[1]]
    row$covered_honest <- covers(hon$ci)
    row$honest_B <- hon$B
    row$honest_cv <- hon$cv
    row$honest_half <- hon$half_width
    row$fid_se_hat <- fid$se_hat
    row$fid_stat <- fid$stat
    row$fid_reject <- fid$reject
    row$n_leaves_e <- as.integer(fit$n_leaves_e)
    row$n_leaves_m0 <- as.integer(fit$n_leaves_m0)
    row$certified_e <- fit$certified_e
    row$certified_m0 <- fit$certified_m0
    row$used_search_e <- fit$used_search_e
    row$used_search_m0 <- fit$used_search_m0
    row$gap_e <- as.numeric(fit$gap_e)
    row$gap_m0 <- as.numeric(fit$gap_m0)
    row$lambda_n_used <- fit$lambda_n
    row$key_e <- gp_e$key
    row$key_mu <- gp_mu$key
    row$at_pseudo_e <- identical(gp_e$key, dgp$class_e$argmin_key)
    row$at_pseudo_mu <- identical(gp_mu$key, dgp$class_mu$argmin_key)
    row$suff_e <- gp_e$sufficient
    row$suff_mu <- gp_mu$sufficient
    row$cf_leaves_e_mean <- mean(cf_leaves_e)
    row$cf_leaves_mu_mean <- mean(cf_leaves_mu)
    row$cf_leaves_e_max <- as.integer(max(cf_leaves_e))
    row$cf_leaves_mu_max <- as.integer(max(cf_leaves_mu))
    row$cf_suff_e_frac <- mean(vapply(cf_e, `[[`, logical(1), "sufficient"))
    row$cf_suff_mu_frac <- mean(vapply(cf_mu, `[[`, logical(1), "sufficient"))
    row$clip_frac <- mean(e_hat <= CLIP_LO + 1e-12 | e_hat >= CLIP_HI - 1e-12)
    row$clip_frac_cf <- mean(e_hat_cf <= CLIP_LO + 1e-12 | e_hat_cf >= CLIP_HI - 1e-12)
    row$max_e_hat <- max(e_hat)
    row$min_e_hat <- min(e_hat)
    row$n_treated <- sum(d$A == 1L)
    row$n_control <- sum(d$A == 0L)
    row$secs_full <- secs_full
    row$secs_cf <- secs_cf
    row$n_warnings <- length(warn_msgs)
    row$warning_messages <- if (length(warn_msgs)) {
      paste(unique(warn_msgs), collapse = " | ")
    } else {
      NA_character_
    }
    row
  }

  if (on_error == "stop") return(body())

  tryCatch(body(), error = function(cnd) {
    row <- rep_row_template()
    row$dgp <- dgp$id
    row$eps <- dgp$eps
    row$residual_mode <- dgp$residual_mode
    row$n <- as.integer(n)
    row$rep <- as.integer(rep)
    row$seed <- as.integer(seed)
    row$theta0 <- dgp$theta0
    row$error_message <- conditionMessage(cnd)
    row
  })
}

# ---- 6 one cell, sequentially, with a checkpoint -------------------------

#' Checkpoint path for a cell
cell_path <- function(dgp_id, n, reps) {
  file.path(DIR_RESULTS, sprintf("cell_%s_n%d_r%d.rds", dgp_id, n, reps))
}

#' Run one (dgp, n) cell and write its checkpoint
#'
#' Sequential by construction: one estimate_att() / estimate_att_crossfit() PAIR
#' at a time, \code{gc(full = TRUE)} on a heartbeat, no future/furrr. A kill
#' mid-sweep costs at most the cell in flight, because every completed cell is
#' already on disk.
#'
#' @param dgp_id One of \code{DGP_IDS}.
#' @param n Sample size.
#' @param reps Replications.
#' @param out Checkpoint path.
#' @param on_error Passed to \code{run_one_rep()}.
#' @param heartbeat Progress/gc interval in replications; 0 disables.
#' @return The cell's results data.frame, invisibly.
run_cell <- function(dgp_id = "eps0", n = 500L, reps = 5L,
                     out = cell_path(dgp_id, n, reps),
                     on_error = "record", heartbeat = 25L) {
  stopifnot(length(n) == 1L, n >= 50, length(reps) == 1L, reps >= 1L)
  n <- as.integer(n)
  reps <- as.integer(reps)

  dgp <- make_dgp(dgp_id)
  lambda_n <- log(n) / n

  cli::cli_h1("cell {dgp_id} n = {n}, reps = {reps}")
  cli::cli_inform(c(
    "*" = "{dgp$label}",
    "*" = "theta_0 = {signif(dgp$theta0, 6)} (exact); e_0 in [{signif(dgp$e0_range[1], 3)}, {signif(dgp$e0_range[2], 3)}]",
    "*" = "|T_Lbar| = {dgp$enumeration$n_partitions} at Lbar = {dgp$leaf_budget}; ass:sparsity holds: e={dgp$class_e$sparsity_holds}, mu={dgp$class_mu$sparsity_holds}",
    "*" = "delta_e: sq = {signif(dgp$class_e$delta_sq, 5)}, kl = {signif(dgp$class_e$delta_kl, 5)} | delta_mu: sq = {signif(dgp$class_mu$delta_sq, 5)}, kl = {signif(dgp$class_mu$delta_kl, 5)}",
    "*" = "exact bias at the pseudo-true limits: {signif(dgp$pseudo$bias, 5)}; D_w*D_mu = {signif(dgp$pseudo$D_w * dgp$pseudo$D_mu, 5)}",
    "*" = "<g,h>_nu = {signif(dgp$ip_gh$ip, 5)} (cor {signif(dgp$ip_gh$cor, 5)}); lambda_n = {signif(lambda_n, 4)}; z*sigma scale ~ {signif(1.96 / sqrt(n) * 0.6, 4)}"
  ))

  rows <- vector("list", reps)
  for (r in seq_len(reps)) {
    rows[[r]] <- run_one_rep(dgp, n, r, on_error = on_error)
    if (heartbeat > 0L && r %% heartbeat == 0L) {
      done <- do.call(rbind, rows[seq_len(r)])
      gc(full = TRUE, verbose = FALSE)
      cli::cli_inform(paste0(
        "  ", r, "/", reps,
        "  cov=", signif(mean(done$covered, na.rm = TRUE), 3),
        " cov_anch=", signif(mean(done$covered_anchor, na.rm = TRUE), 3),
        " cov_hon=", signif(mean(done$covered_honest, na.rm = TRUE), 3),
        " cov_cf=", signif(mean(done$covered_cf, na.rm = TRUE), 3),
        " rej=", signif(mean(done$fid_reject, na.rm = TRUE), 3),
        " s/rep=", signif(mean(done$secs_full + done$secs_cf, na.rm = TRUE), 3)
      ))
    }
  }
  res <- do.call(rbind, rows)

  payload <- list(
    results = res,
    dgp_summary = dgp_population_summary_dial(dgp),
    meta = cell_metadata(dgp, n, reps, lambda_n, res)
  )
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(payload, out)

  n_fail <- sum(!is.na(res$error_message))
  cli::cli_alert_success(
    "wrote {.file {basename(out)}} -- {reps} reps, {n_fail} failure(s), {signif(sum(res$secs_full + res$secs_cf, na.rm = TRUE) / 60, 3)} min"
  )
  # Distinct warning texts once per cell, not once per fit: a NEW warning must be
  # visible without hundreds of copies of an expected one drowning it out.
  wm <- unique(stats::na.omit(res$warning_messages))
  if (length(wm)) {
    cli::cli_inform(c(
      "i" = "{sum(res$n_warnings, na.rm = TRUE)} muffled warning(s), {length(wm)} distinct text(s):",
      stats::setNames(substr(wm, 1, 160), rep("*", length(wm)))
    ))
  }
  invisible(res)
}

#' Provenance stamped into every checkpoint
cell_metadata <- function(dgp, n, reps, lambda_n, res) {
  # A missing git SHA is a legitimate outcome (the path may not be a repository),
  # but the REASON is recorded rather than collapsed to NA: a bare NA is
  # indistinguishable from "git is broken".
  git_sha <- function(path) {
    out <- tryCatch(
      system2("git", c("-C", path, "rev-parse", "--short", "HEAD"),
              stdout = TRUE, stderr = FALSE),
      error = function(e) paste0("<git failed: ", conditionMessage(e), ">"),
      warning = function(w) paste0("<git warned: ", conditionMessage(w), ">")
    )
    if (length(out) == 0L) "<no git output>" else out[[1L]]
  }
  list(
    stem = "honest_inference_sparsity_failure",
    spec = "quality_reports/specs/2026-08-21_honest-inference-sparsity-failure.md",
    dgp = dgp$id, dgp_label = dgp$label,
    eps = dgp$eps, residual_mode = dgp$residual_mode,
    n = n, reps = reps,
    seed_master = SEED_MASTER,
    seed_scheme = "rep_seed(dgp, n, rep) = xxhash64(SEED_MASTER|dgp|n|rep)",
    leaf_budget = dgp$leaf_budget, m_n = M_N,
    lambda_n_rate = "NULL -> log(n)/n", lambda_n_value = lambda_n,
    propensity_loss = PROPENSITY_LOSS, outcome_type = OUTCOME_TYPE,
    anchor = "estimate_att_crossfit(K = 5, cv_regularization = TRUE, max_depth = 4L)",
    anchor_K = ANCHOR_K, anchor_max_depth = ANCHOR_MAX_DEPTH,
    alpha = ALPHA, z_alpha = Z_ALPHA,
    se_hat_choice = "sqrt(sigma_full^2 + sigma_crossfit^2), fixed, never tuned per DGP",
    honest_ci_call = "honest_ci(theta_full, sigma_crossfit, delta, se_delta = 0)",
    worker_limit = WORKER_LIMIT, parallel_cv = PARALLEL_CV, use_rashomon = FALSE,
    theta0 = dgp$theta0,
    card_T = dgp$enumeration$n_partitions,
    delta_e_sq = dgp$class_e$delta_sq, delta_e_kl = dgp$class_e$delta_kl,
    delta_mu_sq = dgp$class_mu$delta_sq, delta_mu_kl = dgp$class_mu$delta_kl,
    sparsity_holds_e = dgp$class_e$sparsity_holds,
    sparsity_holds_mu = dgp$class_mu$sparsity_holds,
    bias_pseudo = dgp$pseudo$bias,
    D_w = dgp$pseudo$D_w, D_mu = dgp$pseudo$D_mu,
    ip_gh = dgp$ip_gh$ip, cor_gh = dgp$ip_gh$cor,
    dgp_params = dgp$params,
    e_vars = dgp$e_vars, mu_vars = dgp$mu_vars,
    clip = c(lo = CLIP_LO, hi = CLIP_HI),
    secs_total = sum(res$secs_full + res$secs_cf, na.rm = TRUE),
    n_failures = sum(!is.na(res$error_message)),
    doubletree_sha = git_sha(PKG_ROOT),
    optimaltrees_sha = git_sha(.OPTIMALTREES_SRC),
    r_version = R.version.string,
    optimaltrees_version = as.character(utils::packageVersion("optimaltrees")),
    run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}
