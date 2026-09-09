# ============================================================
# common.R
# Study: partition_recovery_clt  (doubletree, S1)
# Spec:  quality_reports/specs/2026-09-01_partition-recovery-clt.md
#
# Shared harness: paths, package loading, the design grid, deterministic seeds,
# one replication, one cell. Sourced by run_pilot.R / run_cell.R / run_sweep.R /
# analyze.R. Sourcing this file has no effect beyond defining objects and
# creating the study's output directories.
#
# Run everything from the doubletree package root.
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

STUDY_DIR   <- file.path(PKG_ROOT, "simulations", "partition_recovery_clt")
DIR_CODE    <- file.path(STUDY_DIR, "code")
DIR_RESULTS <- file.path(STUDY_DIR, "results")
DIR_FIGURES <- file.path(STUDY_DIR, "figures")
DIR_TABLES  <- file.path(STUDY_DIR, "tables")
for (.d in c(DIR_RESULTS, DIR_FIGURES, DIR_TABLES)) {
  if (!dir.exists(.d)) dir.create(.d, recursive = TRUE)
}

# The installed doubletree in the user library can predate the current source
# tree, so load the SOURCE explicitly -- same posture as
# propensity_loss_choice/code/common.R, and for the same reason.
#
# TWO PACKAGE-LOADING PATHS, ONE GATE (added 2026-09-09 for SLURM deployment):
#
#   DOUBLETREE_USE_INSTALLED unset/"0"  DEV path (default, unchanged behaviour):
#       pkgload::load_all() on this source tree and on the sibling
#       ../optimaltrees source tree. Correct on the dev box.
#
#   DOUBLETREE_USE_INSTALLED="1"        CLUSTER path: library() against packages
#       that were R CMD INSTALLed on the cluster. Module R does not carry a working
#       pkgload/devtools dev-load, and `../optimaltrees` is a RELATIVE sibling path
#       that is only correct when the checkout happens to have both repos side by side.
.OPTIMALTREES_SRC <- file.path(dirname(PKG_ROOT), "optimaltrees")

USE_INSTALLED_PKGS <- Sys.getenv("DOUBLETREE_USE_INSTALLED", "0") == "1"
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

# Fail fast rather than at replication 700: this study needs the CURRENT
# estimate_att() signature, which it read directly rather than assumed.
for (.arg in c("leaf_budget", "propensity_loss", "outcome_type", "lambda_n", "m_n")) {
  if (!.arg %in% names(formals(doubletree::estimate_att))) {
    cli::cli_abort(c(
      "Loaded {.fun doubletree::estimate_att} has no {.arg {.arg}} argument.",
      i = "A stale installed copy is probably shadowing the source tree."
    ))
  }
}

source(file.path(DIR_CODE, "enumerate_sufficient_class.R"))
source(file.path(DIR_CODE, "dgps.R"))

# ST3 reuses generate_dgp_weak_overlap() unmodified from its existing home.
DGP_SRC_STRESS <- file.path(PKG_ROOT, "simulations", "dgps", "dgps_stress.R")
if (!file.exists(DGP_SRC_STRESS)) {
  cli::cli_abort("DGP source {.file {DGP_SRC_STRESS}} not found.")
}
source(DGP_SRC_STRESS, local = FALSE)
if (!is.function(get0("generate_dgp_weak_overlap"))) {
  cli::cli_abort("{.fun generate_dgp_weak_overlap} was not defined by {.file {DGP_SRC_STRESS}}.")
}

# ---- 1 fixed design constants -------------------------------------------

SEED_MASTER <- 20260901L

# Lbar of ass:budget. Not swept: it is a fixed structural parameter, and 4 is the
# smallest budget at which BOTH nuisances of BOTH grid-sparse DGPs are exactly
# representable (each needs 3 leaves) while still leaving room for a redundant
# refinement -- which is what makes S_j non-singleton and exercises the
# currency-correct recovery indicator.
LEAF_BUDGET  <- 4L
M_N          <- 1L    # eq:feasible floor; estimate_att()'s own default
LAMBDA_N     <- NULL  # NULL -> log(n)/n inside estimate_att() (prop:parsimony rate)
OUTCOME_TYPE <- "binary"

# estimate_att()'s SHIPPED default. The method under test is the estimator as it
# ships (spec §1), so the default is used rather than the squared-error variant.
# Consequence, stated because it matters for reading the margin columns:
# manuscript.tex defines R^{(j)} as a squared-error risk, so Delta_j there is
# `delta_sq`, but selection here actually minimises log-loss, so the margin that
# governs these runs is `delta_kl`. Both are reported everywhere.
PROPENSITY_LOSS <- "log_loss"

# Memory safety (this project's own convention; see
# simulations/docs/MEMORY_SAFE_SIMULATIONS.md and MEMORY.md [LEARN:rashomon-memory]):
# single best tree only -- no Rashomon enumeration anywhere on this path -- and
# no parallelism. worker_limit = 1L is forwarded to fit_tree() through
# estimate_att()'s `...`; BLAS/OpenMP threads are pinned by the runner scripts.
WORKER_LIMIT <- 1L

# The propensity clip of ass:construct(c), applied inside estimate_att(). Copied
# rather than imported because the package constants are unexported, and then
# CHECKED against the package at load so the ST3 clip diagnostic cannot silently
# measure the wrong boundary if the package changes.
CLIP_LO <- 0.01
CLIP_HI <- 0.99
local({
  pkg_lo <- get0(".PROPENSITY_LOWER_BOUND", envir = asNamespace("doubletree"))
  pkg_hi <- get0(".PROPENSITY_UPPER_BOUND", envir = asNamespace("doubletree"))
  if (is.null(pkg_lo) || is.null(pkg_hi) ||
      !isTRUE(all.equal(c(pkg_lo, pkg_hi), c(CLIP_LO, CLIP_HI)))) {
    cli::cli_abort(c(
      "doubletree's propensity clip bounds are not ({CLIP_LO}, {CLIP_HI}).",
      i = "The ST3 clip-fraction diagnostic would measure the wrong boundary."
    ))
  }
})

# ---- 2 the design grid --------------------------------------------------

# Spec §3/§5: ST1's small n and F1's large n are ONE grid, not a favourable-only
# sweep with the stress cells split off elsewhere. The regime label is a function
# of (dgp, n), so F1 and ST1 share a population DGP and differ only in n.
N_GRID_SMALL <- c(200L, 500L, 1000L)     # ST1: selection not yet converged
N_GRID_LARGE <- c(2000L, 4000L, 8000L)   # F1 proper
N_GRID       <- c(N_GRID_SMALL, N_GRID_LARGE)

DGP_IDS <- c("f1", "f1_weak", "st3_weak_overlap")

DGP_LABELS <- c(
  f1               = "F1/ST1: grid-exact sparsity, wide margin",
  f1_weak          = "ST2: grid-exact sparsity, weak margin",
  st3_weak_overlap = "ST3: weak overlap"
)

REGIME_LEVELS <- c("F1", "ST1", "ST2", "ST3")

#' Regime label for a (dgp, n) cell
regime_of <- function(dgp, n) {
  ifelse(dgp == "f1_weak", "ST2",
         ifelse(dgp == "st3_weak_overlap", "ST3",
                ifelse(n >= min(N_GRID_LARGE), "F1", "ST1")))
}

#' The full design grid, one row per cell
design_grid <- function() {
  g <- expand.grid(n = N_GRID, dgp = DGP_IDS,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  g <- g[order(match(g$dgp, DGP_IDS), g$n), c("dgp", "n"), drop = FALSE]
  g$regime <- regime_of(g$dgp, g$n)
  rownames(g) <- NULL
  g
}

# ---- 3 DGP construction (enumeration done ONCE per DGP) -----------------

#' Build one DGP by id, with its exact sufficient classes and margins attached
#'
#' The margin values quoted in the comments are what the enumeration returns for
#' these parameters; they are NOT hard-coded anywhere, and run_pilot.R prints the
#' computed values so the two can be compared.
make_dgp <- function(id) {
  spec <- switch(
    id,
    # WIDE MARGIN, CONFOUNDED. e_0 on (X1, X2), mu_0 on (X1, X3) -- X1 SHARED, so
    # A and Y(0) are dependent and nuisance misfit can actually move theta_hat
    # (dgps.R header note (c); design audit §2; spec §8 addendum item 1).
    #
    # Leaf values are NOT the widest available spread, and that is the whole
    # point. Sharing the X1 gate depletes control mass exactly where the
    # propensity is high, and mu_0's fine split lives inside {X1 = 1}, so a wide
    # e starves Delta_mu. Narrowing e to 0.20/0.45/0.75 while keeping mu wide at
    # 0.10/0.50/0.90 balances the two nuisances. All margin values here are
    # reproduced by code/calibrate_margin.R (a required deliverable of the spec's
    # §8 addendum) and are printed per cell by run_cell(); they are not asserted
    # anywhere and must never be trusted from this comment alone.
    f1 = dgp_spec_grid_sparse(
      p = 5L,
      e_lo = 0.20, e_mid = 0.45, e_gap = 0.30,
      mu_lo = 0.10, mu_mid = 0.50, mu_gap = 0.40,
      tau = 0.08, leaf_budget = LEAF_BUDGET,
      id = "f1", label = DGP_LABELS[["f1"]]
    ),
    # ST2's NEAR-TIE. Structurally identical to f1 -- same shared-X1 confounding,
    # same e_lo/e_mid, so the X1 gate is held fixed -- and ONLY the X2 (resp. X3)
    # split narrows, which is the split the margin is measured on. The gap is
    # chosen so the empirical failure rate CROSSES inside the study's own n grid:
    # a gap that is penalty-dominated at every n leaves prop:selection-rate's
    # decay-shape comparison with nothing to show (spec §8 addendum item 3), and
    # a gap clean from n = 200 shows no near-tie effect at all. See
    # code/calibrate_margin.R section 4 for the surface this was read off.
    f1_weak = dgp_spec_grid_sparse(
      p = 5L,
      e_lo = 0.20, e_mid = 0.45, e_gap = 0.20,
      mu_lo = 0.10, mu_mid = 0.50, mu_gap = 0.20,
      tau = 0.08, leaf_budget = LEAF_BUDGET,
      id = "f1_weak", label = DGP_LABELS[["f1_weak"]]
    ),
    st3_weak_overlap = dgp_spec_weak_overlap(
      tau = 0.10, leaf_budget = LEAF_BUDGET,
      id = "st3_weak_overlap", label = DGP_LABELS[["st3_weak_overlap"]]
    ),
    cli::cli_abort("Unknown dgp {.val {id}}; expected one of {.val {DGP_IDS}}.")
  )
  build_dgp(spec)
}

# ---- 4 deterministic seeds ----------------------------------------------

#' Per-replicate seed: a pure function of (SEED_MASTER, dgp, n, rep)
#'
#' Same scheme as propensity_loss_choice/code/common.R, so the whole sweep is a
#' function of one integer. 7 hex digits < 2^28 is always a valid R integer.
rep_seed <- function(dgp, n, rep) {
  key <- paste(SEED_MASTER, dgp, n, rep, sep = "|")
  strtoi(substr(digest::digest(key, algo = "xxhash64"), 1, 7), base = 16L)
}

# ---- 5 one replication --------------------------------------------------

#' The row schema every replication returns
#'
#' Defined once so success rows and failure rows are guaranteed to have the same
#' columns -- a failure must never be a shorter row that silently drops out of a
#' rbind, and must never be a plausible-looking NA estimate with no message.
rep_row_template <- function() {
  data.frame(
    dgp = NA_character_, regime = NA_character_, n = NA_integer_,
    rep = NA_integer_, seed = NA_integer_,
    theta0 = NA_real_, theta = NA_real_, sigma = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, ci_width = NA_real_, covered = NA,
    # (P1), the enumeration-based indicator. THE headline diagnostic: a coverage
    # number without it cannot distinguish "Condition (P) held and the CLT
    # delivered" from "(P) failed but the bias happened to be small at this n".
    recovered_e = NA, recovered_mu = NA, recovered_both = NA,
    excess_sq_e = NA_real_, excess_kl_e = NA_real_,
    excess_sq_mu = NA_real_, excess_kl_mu = NA_real_,
    n_leaves_grid_e = NA_integer_, n_leaves_grid_mu = NA_integer_,
    in_class_e = NA, in_class_mu = NA,
    # Partition identity, so WITHIN-S_j variation across reps is measurable.
    # That variation plus valid coverage is lem:uniform's empirical signature
    # (spec §7); it is checked, not assumed.
    key_e = NA_character_, key_mu = NA_character_,
    # estimate_att()'s own fitted-leaf and certification diagnostics.
    n_leaves_e = NA_integer_, n_leaves_m0 = NA_integer_,
    certified_e = NA, certified_m0 = NA,
    used_search_e = NA, used_search_m0 = NA,
    gap_e = NA_real_, gap_m0 = NA_real_,
    lambda_n_used = NA_real_,
    # ST3's disentangling device (spec §3): the FRACTION OF UNITS whose fitted
    # propensity sits on the clip boundary. Recorded in every regime, not just
    # ST3, so "coverage is off" can be attributed to (P1) failing, to
    # ass:construct's clip biting, to both, or to neither.
    clip_frac = NA_real_, clip_frac_hi = NA_real_, clip_frac_lo = NA_real_,
    max_e_hat = NA_real_, min_e_hat = NA_real_,
    n_treated = NA_integer_, n_control = NA_integer_,
    secs = NA_real_,
    # Warnings are COUNTED AND KEPT, not discarded. optimaltrees emits one
    # expected, benign warning per tree fit here ("High-dimensional data detected
    # (estimated 155 binary features after discretization)"): the estimate is
    # made BEFORE discretization runs and so does not know that all-binary X
    # takes discretize_features()' no-op fast path, and the model_limit it then
    # sets caps the Rashomon MODEL COUNT, which is irrelevant on this single-tree
    # path. It cannot be silenced by passing model_limit explicitly --
    # optimaltrees() then errors with "formal argument model_limit matched by
    # multiple actual arguments" -- so it is muffled for legibility and recorded
    # here instead. certified_e / certified_m0 == TRUE is the independent
    # evidence that nothing was truncated.
    n_warnings = NA_integer_, warning_messages = NA_character_,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Run one replication: draw, fit, measure
#'
#' @param dgp A built DGP (\code{make_dgp()} output).
#' @param n Sample size.
#' @param rep Replication index.
#' @param on_error \code{"stop"} (pilot/interactive posture: a failure surfaces
#'   at its origin) or \code{"record"} (long-run posture: the condition message
#'   is stored so the cell finishes and a failure RATE can be reported).
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
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (!is.finite(fit$theta) || !is.finite(fit$sigma)) {
      stop("Non-finite output: theta = ", fit$theta, ", sigma = ", fit$sigma,
           call. = FALSE)
    }

    rec_e <- partition_recovered(fit$nuisance_fits$e_model, dgp$class_e)
    rec_mu <- partition_recovered(fit$nuisance_fits$m0_model, dgp$class_mu)

    e_hat <- fit$nuisance_fits$propensity
    ci <- fit$ci_95

    row <- rep_row_template()
    row$dgp <- dgp$id
    row$regime <- regime_of(dgp$id, n)
    row$n <- as.integer(n)
    row$rep <- as.integer(rep)
    row$seed <- as.integer(seed)
    row$theta0 <- dgp$theta0
    row$theta <- fit$theta
    row$sigma <- fit$sigma
    row$ci_lo <- ci[[1]]
    row$ci_hi <- ci[[2]]
    row$ci_width <- ci[[2]] - ci[[1]]
    row$covered <- ci[[1]] <= dgp$theta0 && dgp$theta0 <= ci[[2]]
    row$recovered_e <- rec_e$recovered
    row$recovered_mu <- rec_mu$recovered
    row$recovered_both <- rec_e$recovered && rec_mu$recovered
    row$excess_sq_e <- rec_e$excess_sq
    row$excess_kl_e <- rec_e$excess_kl
    row$excess_sq_mu <- rec_mu$excess_sq
    row$excess_kl_mu <- rec_mu$excess_kl
    row$n_leaves_grid_e <- rec_e$n_leaves_grid
    row$n_leaves_grid_mu <- rec_mu$n_leaves_grid
    row$in_class_e <- rec_e$in_enumerated_class
    row$in_class_mu <- rec_mu$in_enumerated_class
    row$key_e <- rec_e$partition_key
    row$key_mu <- rec_mu$partition_key
    row$n_leaves_e <- as.integer(fit$n_leaves_e)
    row$n_leaves_m0 <- as.integer(fit$n_leaves_m0)
    row$certified_e <- fit$certified_e
    row$certified_m0 <- fit$certified_m0
    row$used_search_e <- fit$used_search_e
    row$used_search_m0 <- fit$used_search_m0
    # `gap_*` is documented to be 0 when certified and NA when infeasible, so NA
    # is a LEGITIMATE value here -- but the FIELD must exist. An absent field
    # would mean estimate_att()'s return contract changed, which is exactly the
    # kind of thing that must surface rather than be papered over with NA (the
    # earlier `if (is.null(...)) NA_real_` did the latter, and would have
    # reported a column of NAs indistinguishable from "all infeasible").
    for (fld in c("gap_e", "gap_m0", "certified_e", "certified_m0",
                  "used_search_e", "used_search_m0", "n_leaves_e", "n_leaves_m0",
                  "lambda_n")) {
      if (!fld %in% names(fit)) {
        stop("estimate_att() returned no `", fld, "` field; its return contract ",
             "has changed and this study's diagnostics are no longer valid.",
             call. = FALSE)
      }
    }
    row$gap_e <- as.numeric(fit$gap_e)
    row$gap_m0 <- as.numeric(fit$gap_m0)
    row$lambda_n_used <- fit$lambda_n
    row$clip_frac_hi <- mean(e_hat >= CLIP_HI - 1e-12)
    row$clip_frac_lo <- mean(e_hat <= CLIP_LO + 1e-12)
    row$clip_frac <- row$clip_frac_hi + row$clip_frac_lo
    row$max_e_hat <- max(e_hat)
    row$min_e_hat <- min(e_hat)
    row$n_treated <- sum(d$A == 1L)
    row$n_control <- sum(d$A == 0L)
    row$secs <- secs
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
    row$regime <- regime_of(dgp$id, n)
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
#' Sequential by construction: one \code{estimate_att()} fit at a time, with
#' \code{gc(full = TRUE)} on a heartbeat. No \code{future}/\code{furrr}. Mirrors
#' threshold_superconsistency/code/run_cell.R: a kill mid-sweep costs at most the
#' cell in flight, because every completed cell is already on disk.
#'
#' @param dgp_id One of \code{DGP_IDS}.
#' @param n Sample size.
#' @param reps Replications.
#' @param out Checkpoint path.
#' @param on_error Passed to \code{run_one_rep()}.
#' @param heartbeat Progress/gc interval in replications; 0 disables.
#' @return The cell's results data.frame, invisibly.
run_cell <- function(dgp_id = "f1", n = 200L, reps = 5L,
                     out = cell_path(dgp_id, n, reps),
                     on_error = "record", heartbeat = 100L) {
  stopifnot(length(n) == 1L, n >= 50, length(reps) == 1L, reps >= 1L)
  n <- as.integer(n)
  reps <- as.integer(reps)

  dgp <- make_dgp(dgp_id)
  lambda_n <- log(n) / n
  regime <- regime_of(dgp_id, n)

  cli::cli_h1("cell {dgp_id} (regime {regime}) n = {n}, reps = {reps}")
  cli::cli_inform(c(
    "*" = "theta_0 = {signif(dgp$theta0, 6)} (exact); e_0 in [{signif(dgp$e0_range[1], 3)}, {signif(dgp$e0_range[2], 3)}]",
    "*" = "|T_Lbar| = {dgp$enumeration$n_partitions}; |S_e| = {dgp$class_e$n_sufficient}, |S_mu| = {dgp$class_mu$n_sufficient}",
    "*" = "Delta_e: sq = {signif(dgp$class_e$delta_sq, 4)}, kl = {signif(dgp$class_e$delta_kl, 4)}",
    "*" = "Delta_mu: sq = {signif(dgp$class_mu$delta_sq, 4)}, kl = {signif(dgp$class_mu$delta_kl, 4)}",
    "*" = "lambda_n = log(n)/n = {signif(lambda_n, 4)}"
  ))
  # A penalty above the margin means the PENALISED population optimum is outside
  # S_j, so recovery must fail for a reason that has nothing to do with sampling
  # noise. Surfaced per cell rather than left to be inferred from the results.
  for (j in c("e", "mu")) {
    dkl <- if (j == "e") dgp$class_e$delta_kl else dgp$class_mu$delta_kl
    if (lambda_n >= dkl) {
      cli::cli_inform(c("!" = paste0(
        "lambda_n = ", signif(lambda_n, 4), " >= Delta_", j, "(kl) = ",
        signif(dkl, 4), ": the penalised population optimum for ", j,
        " lies OUTSIDE S_", j, ", so selection failure in this cell is ",
        "STRUCTURAL (the penalty excludes the necessary split), not noise."
      )))
    }
  }

  rows <- vector("list", reps)
  for (r in seq_len(reps)) {
    rows[[r]] <- run_one_rep(dgp, n, r, on_error = on_error)
    if (heartbeat > 0L && r %% heartbeat == 0L) {
      done <- do.call(rbind, rows[seq_len(r)])
      gc(full = TRUE, verbose = FALSE)
      cli::cli_inform(paste0(
        "  ", r, "/", reps,
        "  rec_e=", signif(mean(done$recovered_e, na.rm = TRUE), 3),
        " rec_mu=", signif(mean(done$recovered_mu, na.rm = TRUE), 3),
        " cov=", signif(mean(done$covered, na.rm = TRUE), 3),
        " s/fit=", signif(mean(done$secs, na.rm = TRUE), 3)
      ))
    }
  }
  res <- do.call(rbind, rows)

  payload <- list(
    results = res,
    dgp_summary = dgp_population_summary(dgp),
    meta = cell_metadata(dgp, n, reps, lambda_n, regime, res)
  )
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(payload, out)

  n_fail <- sum(!is.na(res$error_message))
  cli::cli_alert_success(
    "wrote {.file {basename(out)}} -- {reps} reps, {n_fail} failure(s), {signif(sum(res$secs, na.rm = TRUE) / 60, 3)} min"
  )
  # Distinct warning texts once per cell, not once per fit: a NEW warning must be
  # visible without 1000 copies of the expected one drowning it out.
  wm <- unique(stats::na.omit(res$warning_messages))
  if (length(wm)) {
    cli::cli_inform(c(
      "i" = "{sum(res$n_warnings, na.rm = TRUE)} muffled warning(s), {length(wm)} distinct text(s):",
      stats::setNames(substr(wm, 1, 150), rep("*", length(wm)))
    ))
  }
  invisible(res)
}

#' Provenance stamped into every checkpoint
cell_metadata <- function(dgp, n, reps, lambda_n, regime, res) {
  # A missing git SHA is a legitimate outcome (the path may not be a repository),
  # but the REASON must be recorded rather than collapsed to NA: a bare NA here is
  # indistinguishable from "git is broken", and provenance that silently goes
  # missing is worse than provenance that says why it is missing.
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
    stem = "partition_recovery_clt",
    spec = "quality_reports/specs/2026-09-01_partition-recovery-clt.md",
    dgp = dgp$id, dgp_label = dgp$label, regime = regime,
    n = n, reps = reps,
    seed_master = SEED_MASTER,
    seed_scheme = "rep_seed(dgp, n, rep) = xxhash64(SEED_MASTER|dgp|n|rep)",
    leaf_budget = dgp$leaf_budget, m_n = M_N,
    lambda_n_rate = "NULL -> log(n)/n", lambda_n_value = lambda_n,
    propensity_loss = PROPENSITY_LOSS, outcome_type = OUTCOME_TYPE,
    worker_limit = WORKER_LIMIT, use_rashomon = FALSE,
    theta0 = dgp$theta0,
    card_T = dgp$enumeration$n_partitions,
    delta_e_sq = dgp$class_e$delta_sq, delta_e_kl = dgp$class_e$delta_kl,
    delta_mu_sq = dgp$class_mu$delta_sq, delta_mu_kl = dgp$class_mu$delta_kl,
    n_sufficient_e = dgp$class_e$n_sufficient,
    n_sufficient_mu = dgp$class_mu$n_sufficient,
    dgp_params = dgp$params,
    e_vars = dgp$e_vars, mu_vars = dgp$mu_vars,
    clip = c(lo = CLIP_LO, hi = CLIP_HI),
    secs_total = sum(res$secs, na.rm = TRUE),
    n_failures = sum(!is.na(res$error_message)),
    doubletree_sha = git_sha(PKG_ROOT),
    optimaltrees_sha = git_sha(.OPTIMALTREES_SRC),
    r_version = R.version.string,
    optimaltrees_version = as.character(utils::packageVersion("optimaltrees")),
    run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}
