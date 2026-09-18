# ============================================================
# common.R
# Study: single_tree_corollaries  (doubletree)
# Spec:  quality_reports/specs/2026-09-09_single-outcome-tree-corollaries.md
#
# Shared harness: paths, package loading, the tied-estimator computation, the
# oracle-tau and fitted-tau arms of spec §2, the naive/corrected SEs of spec §5,
# partition recovery, one replication, one cell. Sourced by run_pilot.R and
# analyze.R. Sourcing has no effect beyond defining objects and creating the
# study's output directories.
#
# Structure follows head_to_head_comparison/code/common.R (S5): same
# source-tree-vs-installed guard, same deterministic per-replication seeding,
# same long results format, same "no silent fallback" posture.
#
# ------------------------------------------------------------
# THE OBJECT UNDER TEST (spec §1)
# ------------------------------------------------------------
#
#   theta_tie = n_1^{-1} sum_{A_i = 1} { Y_i - mu-hat(X_i) }
#
# where mu-hat is DOUBLETREE'S OWN outcome tree -- estimate_att()'s
# `nuisance_fits$m0_model` / `outcome_control`, fit on control-only data at a
# fixed leaf_budget, leaf-constant, predicted on all n. NO propensity tree is
# fitted for this estimator: it IS doubletree's two-tree estimator with the
# propensity tree forced equal to the outcome tree (`cor:single-tree`).
#
# WHY THAT ONE-LINER IS THE WHOLE ESTIMATOR. eq:score's control-arm term is
#   -n_1^{-1} sum_l w-hat(l) sum_{A=0, l} {Y - mu-hat(l)}.
# When mu-hat(l) is the within-leaf CONTROL mean, the inner sum is exactly zero in
# every leaf, so the term vanishes for ANY w-hat constant within those leaves --
# leaving the treated-arm mean above. Two preconditions, then: mu-hat's leaf values
# are the control means, and e-hat is leaf-constant on mu-hat's OWN partition. The
# VALUES of e-hat are irrelevant to the identity (spec §1: "for any
# leaf-wise-constant e-hat"); setting e-hat(l) = n_1(l)/n(l) is what makes this the
# outcome tree used AS the propensity tree, which is the object under test. This is
# an algebraic identity, exact at finite n, verified per replication here
# (`identity_gap`) and in tests/testthat/test-single-tree-identity.R (spec §1).
#
# NO EXPORTED PACKAGE FUNCTION (spec §7). theta_tie stays simulation-local: this
# construction's own corollary (`rem:single-tree-se`) states that its DEFAULT
# (naive) SE is wrong, so shipping an `estimate_att_tie()` would invite users to
# call it and report an under-covering interval.
#
# ------------------------------------------------------------
# THE TWO ARMS (spec §2), AND WHY THE SPLIT IS LOAD-BEARING
# ------------------------------------------------------------
#
#   ORACLE-TAU: tau imposed analytically. Leaf-wise sample means on the KNOWN
#     true leaves, no tree search at all. Isolates the fixed-tau mathematical
#     mechanism (the Jensen ordering, Regime C's heteroskedastic reversal, the
#     naive-SE under-coverage) with zero contamination from tree-search noise.
#   FITTED-TAU: doubletree's actual estimate_att() at leaf_budget EXACTLY equal to
#     |leaves(tau)|, with leaf assignments read via the existing internal helper
#     doubletree:::.tree_leaf_paths() -- reused, not reimplemented (spec §7).
#
# Without the split, a regime failure is ambiguous between "the corollary's math
# is wrong" and "the tree missed tau". Regime E is fitted-arm ONLY, because the
# whole point there is which tau the fitting realises.
#
# THE COMPARATOR (spec §1) is doubletree's full two-tree estimate_att() at the
# SAME leaf_budget on the SAME data -- not a rival method, the other half of the
# same identity. It runs as its own arm (`two_tree`).
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

STUDY_DIR   <- file.path(PKG_ROOT, "simulations", "single_tree_corollaries")
DIR_CODE    <- file.path(STUDY_DIR, "code")
DIR_RESULTS <- file.path(STUDY_DIR, "results")
DIR_TABLES  <- file.path(STUDY_DIR, "tables")
DIR_FIGURES <- file.path(STUDY_DIR, "figures")
for (.d in c(DIR_RESULTS, DIR_TABLES, DIR_FIGURES)) {
  if (!dir.exists(.d)) dir.create(.d, recursive = TRUE)
}

# Two loading paths, one gate -- S5's convention, unchanged. The installed
# doubletree can predate this source tree, so the DEV path loads the source.
.OPTIMALTREES_SRC <- file.path(dirname(PKG_ROOT), "optimaltrees")
USE_INSTALLED_PKGS <- Sys.getenv("STC_USE_INSTALLED", "0") == "1"
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

# Fail fast rather than at replication 40: read the contracts, do not assume them.
# NOTE: the loop variable must not be named `.arg` -- cli >= 3.4.0 reserves
# any `{.foo}` starting with a dot exclusively for its own inline styles, so
# `{.arg {.arg}}` no longer interpolates a variable named `.arg`; it is parsed
# as an (invalid) nested style directive instead.
for (arg_name in c("leaf_budget", "outcome_type", "lambda_n", "m_n")) {
  if (!arg_name %in% names(formals(doubletree::estimate_att))) {
    cli::cli_abort(c(
      "Loaded {.fun doubletree::estimate_att} has no {.arg {arg_name}} argument.",
      i = "A stale installed copy is probably shadowing the source tree."
    ))
  }
}
# The two internals this study REUSES rather than reimplements (spec §2, §7).
.tree_leaf_paths_fn <- get0(".tree_leaf_paths", envir = asNamespace("doubletree"))
.eif_att_solve_fn <- get0("eif_att_solve", envir = asNamespace("doubletree"))
if (!is.function(.tree_leaf_paths_fn) || !is.function(.eif_att_solve_fn)) {
  cli::cli_abort(c(
    "doubletree no longer provides {.fun .tree_leaf_paths} and/or {.fun eif_att_solve}.",
    i = "Spec §2 requires reusing the existing leaf-path helper, not copy-pasting it."
  ))
}

source(file.path(DIR_CODE, "dgps.R"), local = FALSE)

# ---- 1 fixed design constants -------------------------------------------

SEED_MASTER  <- 20260909L   # the spec's own date
OUTCOME_TYPE <- "continuous"  # spec §3: Y is continuous in every regime
M_N          <- 1L            # eq:feasible floor; estimate_att()'s own default
LAMBDA_N     <- NULL          # NULL -> log(n)/n inside estimate_att()
ALPHA        <- 0.05
Z_ALPHA      <- stats::qnorm(1 - ALPHA / 2)

# Memory safety (this project's convention; simulations/docs/
# MEMORY_SAFE_SIMULATIONS.md): single best tree only, no Rashomon enumeration,
# no parallelism inside a fit.
WORKER_LIMIT <- 1L

# ass:construct(c)'s propensity clip. Copied because the package constants are
# unexported, then CHECKED against the package: a package change would otherwise
# make the tied estimator's clip silently differ from estimate_att()'s, and the
# `cor:single-tree` identity holds only when the clip is INERT.
CLIP_LO <- 0.01
CLIP_HI <- 0.99
local({
  pkg_lo <- get0(".PROPENSITY_LOWER_BOUND", envir = asNamespace("doubletree"))
  pkg_hi <- get0(".PROPENSITY_UPPER_BOUND", envir = asNamespace("doubletree"))
  if (is.null(pkg_lo) || is.null(pkg_hi) ||
      !isTRUE(all.equal(c(pkg_lo, pkg_hi), c(CLIP_LO, CLIP_HI)))) {
    cli::cli_abort(c(
      "doubletree's propensity clip bounds are not ({CLIP_LO}, {CLIP_HI}).",
      i = "The tied estimator would then clip differently from estimate_att()."
    ))
  }
})

# The identity of `cor:single-tree` is exact arithmetic, so its tolerance is a
# floating-point tolerance, not a statistical one. Spec §1/§6 fix 1e-10.
IDENTITY_TOL <- 1e-10

# MEASURED PACKAGE PRECISION LIMIT (2026-09-09), and the reason the fitted arm
# refits its leaf values -- see st_arm_fitted().
#
# optimaltrees stores a leaf's fitted value to SIX significant decimal digits
# (measured: a leaf whose control mean is -0.058940614050488 comes back from
# estimate_att()'s m0_model as exactly -0.058941, a 3.9e-07 gap). The tied
# estimator's identity is exact only when the leaf value IS the within-leaf
# control mean, so using the serialised prediction verbatim caps the achievable
# identity gap at ~4e-07 -- 3.6 orders of magnitude above spec §1's 1e-10, for a
# reason that has nothing to do with cor:single-tree.
#
# MU_LEAF_TOL is therefore the tolerance at which "these leaf values ARE the
# control means" is asserted: loose enough to admit the 6-digit rounding, tight
# enough that a genuinely different leaf value (a wrong subset, a treated-inclusive
# mean, a shrunk value) cannot pass. The actual gap is REPORTED per replication as
# `mu_leaf_gap` rather than merely tolerated.
MU_LEAF_TOL <- 1e-5

# ---- 2 the regime grid --------------------------------------------------

# Spec §3's regimes. `arms` is explicit per regime because Regime E has NO oracle
# arm (spec §3: the point is which tau the fitting realises), and that is a fact
# about the design that belongs in the design table rather than in an `if`.
ST_REGIMES <- list(
  A = list(dgp = "A", claims = "cor:single-tree-saturated",
           arms = c("oracle_tie", "fitted_tie", "two_tree")),
  B = list(dgp = "B", claims = "cor:single-tree-coarsening (<= direction)",
           arms = c("oracle_tie", "fitted_tie", "two_tree")),
  C = list(dgp = "C", claims = "cor:single-tree-coarsening (reversal sub-claim)",
           arms = c("oracle_tie", "fitted_tie", "two_tree")),
  E = list(dgp = "E", claims = "lem:single-tree-linear (random-index caveat)",
           arms = c("fitted_tie", "two_tree"))
)
ST_REGIME_IDS <- names(ST_REGIMES)

# Spec §6's local pilot: nsim = 50 at a single n = 1000. The cluster stage
# (nsim = 1000, n in {1000, 4000}) is a SEPARATE decision and is deliberately not
# encoded as a runnable default here.
PILOT_NSIM <- 50L
PILOT_N    <- 1000L
CLUSTER_NSIM <- 1000L                 # recorded for sizing only; not run here
CLUSTER_N    <- c(1000L, 4000L)       # ditto

#' Per-replicate seed: a pure function of (SEED_MASTER, regime, n, rep)
#'
#' Same scheme as S1/S2/S5's common.R, so the whole study is a function of one
#' integer. 7 hex digits < 2^28 is always a valid R integer. The REGIME is in the
#' key so two regimes at the same n never run on byte-identical data.
st_rep_seed <- function(regime, n, rep) {
  key <- paste(SEED_MASTER, regime, n, rep, sep = "|")
  strtoi(substr(digest::digest(key, algo = "xxhash64"), 1, 7), base = 16L)
}

# ---- 3 leaf assignments from a fitted outcome tree ----------------------

#' Leaf index of every row under estimate_att()'s fitted outcome tree
#'
#' REUSES \code{doubletree:::.tree_leaf_paths()} (spec §2, §7). That helper walks
#' the nested-list tree against the design whose columns its 0-based feature
#' indices address, so the discretised design must be the one the model was fit
#' on -- hence \code{apply_discretization()} whenever the model carries
#' discretisation metadata, mirroring \code{estimate_att_single_tree()}'s own call
#' site exactly. With binary X (the only case estimate_att() accepts) this is a
#' no-op, but relying on that silently would break the moment the package gains
#' internal discretisation.
#'
#' @param fit An \code{estimate_att()} result.
#' @param X The design passed to \code{estimate_att()}.
#' @return Integer leaf index per row of X (1-based, arbitrary but stable order).
st_fitted_leaf <- function(fit, X) {
  model <- fit$nuisance_fits$m0_model
  tree <- model@trees[[1L]]
  meta <- model@discretization_metadata
  Xb <- if (!is.null(meta)) optimaltrees::apply_discretization(X, meta) else X
  paths <- .tree_leaf_paths_fn(tree, Xb)
  leaf <- match(paths, sort(unique(paths)))

  # The strongest available check that the leaf extraction is CORRECT rather than
  # merely plausible: estimate_att()'s own leaf-constant predictions must be
  # constant within these leaves. If the design or the path convention were
  # misaligned, they would not be, and every downstream leaf-wise quantity
  # (e-hat, sigma-hat^2, w_tie) would be quietly computed on the wrong groups.
  spread <- tapply(fit$nuisance_fits$outcome_control, leaf,
                   function(z) diff(range(z)))
  if (max(spread) > 1e-9) {
    stop("Leaf extraction disagrees with estimate_att()'s own predictions: ",
         "outcome_control varies by up to ", signif(max(spread), 4),
         " within an extracted leaf. .tree_leaf_paths() was called against the ",
         "wrong design or the path convention has changed.", call. = FALSE)
  }
  if (length(unique(leaf)) != fit$n_leaves_m0) {
    stop("Extracted ", length(unique(leaf)), " leaves but estimate_att() reports ",
         "n_leaves_m0 = ", fit$n_leaves_m0, ".", call. = FALSE)
  }
  as.integer(leaf)
}

# ---- 4 the tied estimator, and its two standard errors ------------------

#' theta_tie plus the naive and corrected SEs of spec §5
#'
#' \code{m0_hat} must be leaf-constant on \code{leaf} AND equal to the within-leaf
#' CONTROL mean -- both are what make \code{cor:single-tree}'s identity exact, and
#' both are asserted rather than assumed.
#'
#' SEs (spec §5). The naive SE keeps ONLY the treated-arm term of eq:score, which
#' is what a user who reads theta_tie as "a mean over the treated" would compute:
#'   SE_naive^2 = n^{-1} pi-hat^{-2} n^{-1} sum_i A_i {Y_i - mu-hat_i - theta}^2.
#' \code{rem:single-tree-se} says that omits the nonnegative control-arm term
#'   pi^{-2} E[(1 - e_0) w_tie^2 sigma_0^2],
#' so the corrected SE adds its leaf-wise plug-in:
#'   SE_corr^2 = SE_naive^2 + n^{-1} pi-hat^{-2}
#'               sum_l P-hat(l){1 - e-hat(l)} w-hat_tie(l)^2 sigma-hat_0^2(l).
#' Spec §5 writes that second term without the n^{-1}; it is supplied here because
#' SE_naive^2 already carries one and the two must be commensurate. Flagged rather
#' than silently fixed.
#'
#' \code{se_eif} is the ordinary empirical EIF plug-in from the package's own
#' \code{eif_att_solve()} at the tied nuisances -- an independent cross-check on
#' the corrected SE, which should agree closely (the two differ only in whether
#' the control-arm residual variance is pooled per leaf or taken pointwise).
#'
#' @param Y,A Outcome and treatment.
#' @param m0_hat Leaf-constant control-outcome predictions for all n rows.
#' @param leaf Integer leaf index per row.
#' @return A list of scalars plus the per-leaf table.
st_tied_estimator <- function(Y, A, m0_hat, leaf) {
  n <- length(Y)
  stopifnot(length(A) == n, length(m0_hat) == n, length(leaf) == n)
  n_treated <- sum(A == 1L)
  if (n_treated < 1L) stop("No treated units; the ATT is undefined.", call. = FALSE)

  lv <- sort(unique(leaf))
  n_l <- as.integer(tapply(rep(1L, n), leaf, sum)[as.character(lv)])
  n0_l <- as.integer(tapply(A == 0L, leaf, sum)[as.character(lv)])
  if (any(n0_l < 2L)) {
    stop("Leaf(es) with fewer than 2 control units (", paste(n0_l, collapse = ", "),
         "): the leaf-wise sigma-hat_0^2 the corrected SE needs is undefined. ",
         "This is a real data-scarcity event, not something to impute.",
         call. = FALSE)
  }

  # e-hat(l): the leaf's treated fraction, i.e. the propensity tree FORCED equal
  # to the outcome tree. Clipped with estimate_att()'s own bounds; the identity
  # requires the clip to be inert, so the clipped fraction is REPORTED.
  ehat_l_raw <- as.numeric(tapply(A, leaf, mean)[as.character(lv)])
  ehat_l <- pmin(pmax(ehat_l_raw, CLIP_LO), CLIP_HI)
  clip_frac <- sum(n_l[ehat_l_raw != ehat_l]) / n
  w_l <- ehat_l / (1 - ehat_l)

  idx <- match(leaf, lv)
  mu_l <- as.numeric(tapply(m0_hat, leaf, mean)[as.character(lv)])
  if (max(abs(m0_hat - mu_l[idx])) > 1e-9) {
    stop("m0_hat is not leaf-constant on `leaf`; cor:single-tree's identity ",
         "does not apply.", call. = FALSE)
  }
  ctrl_mean_l <- as.numeric(tapply(Y[A == 0L], leaf[A == 0L], mean)[as.character(lv)])
  mu_leaf_gap <- max(abs(mu_l - ctrl_mean_l))
  if (mu_leaf_gap > MU_LEAF_TOL) {
    stop("The leaf values of m0_hat are not the within-leaf CONTROL means ",
         "(max gap ", signif(mu_leaf_gap, 4), " > ", MU_LEAF_TOL, "). The ",
         "identity's second precondition fails, so theta_tie is not the tied ",
         "plug-in of cor:single-tree.", call. = FALSE)
  }

  resid <- Y - m0_hat
  theta_tie <- mean(resid[A == 1L])

  # The SAME quantity via the package's full two-tree solver with the propensity
  # tree forced equal to the outcome tree. Equality to IDENTITY_TOL IS
  # cor:single-tree, checked every replication (spec §6's smoke check).
  forced <- .eif_att_solve_fn(Y, A, ehat_l[idx], m0_hat, n)
  identity_gap <- abs(theta_tie - forced$theta)

  pi_hat <- mean(A)
  Vnum_naive <- pi_hat^(-2) * mean(A * (resid - theta_tie)^2)
  # Plug-in sigma-hat_0^2(l): within-leaf control residual variance. Denominator
  # n0(l) (not n0(l) - 1) matches the population quantity E[(Y - mu_0)^2 | l] the
  # formula is a plug-in for; with n0(l) >= 2 asserted above it is well defined.
  sig0_l <- as.numeric(tapply(resid[A == 0L]^2, leaf[A == 0L], mean)[as.character(lv)])
  P_l <- n_l / n
  ctrl_term <- pi_hat^(-2) * sum(P_l * (1 - ehat_l) * w_l^2 * sig0_l)
  Vnum_corrected <- Vnum_naive + ctrl_term

  list(
    theta = theta_tie,
    theta_forced = forced$theta,
    identity_gap = identity_gap,
    se_naive = sqrt(Vnum_naive / n),
    se_corrected = sqrt(Vnum_corrected / n),
    se_eif = forced$sigma,
    Vnum_naive = Vnum_naive,
    Vnum_corrected = Vnum_corrected,
    naive_over_correct = Vnum_naive / Vnum_corrected,
    pi_hat = pi_hat,
    clip_frac = clip_frac,
    mu_leaf_gap = mu_leaf_gap,
    n_leaves = length(lv),
    leaf_table = data.frame(leaf = lv, n = n_l, n0 = n0_l, P = P_l,
                            ehat = ehat_l, w_tie = w_l, sigma0_sq = sig0_l,
                            mu = mu_l, stringsAsFactors = FALSE)
  )
}

#' Leaf-wise control means on a GIVEN partition (the oracle arm's mu-hat)
#'
#' No tree search: this is spec §2's oracle arm verbatim. Errors rather than
#' imputing when a leaf has no controls -- an empty leaf means mu-hat is genuinely
#' undefined there, and a fallback would report a plausible number for a quantity
#' that does not exist.
st_leafwise_control_mean <- function(Y, A, leaf) {
  lv <- sort(unique(leaf))
  n0 <- as.integer(tapply(A == 0L, leaf, sum)[as.character(lv)])
  if (any(n0 < 1L)) {
    stop("Leaf(es) with zero control units: the oracle arm's leaf-wise control ",
         "mean is undefined there.", call. = FALSE)
  }
  mu_l <- as.numeric(tapply(Y[A == 0L], leaf[A == 0L], mean)[as.character(lv)])
  mu_l[match(leaf, lv)]
}

# ---- 5 partition recovery ------------------------------------------------

#' Do two labellings induce the SAME partition of the rows?
#'
#' Label VALUES are arbitrary (the fitted tree's leaf ids come from path order),
#' so equality is tested as a bijection between the two labellings' blocks: the
#' cross-tabulation must have exactly one nonzero entry per row and per column.
#' A fitted tree with the right number of leaves but a different split is caught;
#' so is a relabelled but identical partition, correctly, as a MATCH.
st_same_partition <- function(a, b) {
  tab <- table(a, b)
  nrow(tab) == ncol(tab) &&
    all(rowSums(tab > 0) == 1L) && all(colSums(tab > 0) == 1L)
}

#' Which candidate partition (if any) did the fitted tree realise?
#'
#' Regimes A-C have one candidate (\code{"true"}), so this doubles as the
#' partition-recovery indicator of spec §5. Regime E has two genuine ties, and
#' naming WHICH one was realised is the regime's primary metric (spec §3).
#'
#' @param fitted_leaf Integer leaf index per row from \code{st_fitted_leaf()}.
#' @param d A draw (\code{spec$draw()} output); its X is mapped to cells.
#' @param dgp The DGP spec.
#' @return Character: a name from \code{dgp$leaf_maps}, or \code{"other"}.
st_realized_partition <- function(fitted_leaf, d, dgp) {
  key_cells <- do.call(paste, c(dgp$cells, sep = "-"))
  key_rows <- do.call(paste, c(d$X, sep = "-"))
  ci <- match(key_rows, key_cells)
  for (nm in names(dgp$leaf_maps)) {
    cand <- dgp$leaf_maps[[nm]][ci]
    if (st_same_partition(fitted_leaf, cand)) return(nm)
  }
  "other"
}

# ---- 6 the arms ----------------------------------------------------------

#' Oracle-tau arm: tied estimator on the KNOWN true leaves, no tree fitting
st_arm_oracle_tie <- function(d, dgp) {
  if (!isTRUE(dgp$oracle_arm) || is.null(dgp$leaf_fun)) {
    stop("DGP '", dgp$id, "' has no single true partition, so spec §2's oracle ",
         "arm is not defined for it (Regime E is fitted-arm only).", call. = FALSE)
  }
  leaf <- dgp$leaf_fun(d$X)
  if (!identical(as.integer(leaf), as.integer(d$leaf_true))) {
    stop("dgp$leaf_fun() disagrees with the draw's own leaf_true; the oracle arm ",
         "would be built on a different partition from the one the population ",
         "quantities were computed for.", call. = FALSE)
  }
  m0_hat <- st_leafwise_control_mean(d$Y, d$A, leaf)
  res <- st_tied_estimator(d$Y, d$A, m0_hat, leaf)
  res$realized_partition <- names(dgp$leaf_maps)[[1L]]
  res$recovered <- TRUE   # imposed, not estimated
  res
}

#' Fitted-tau arm: doubletree's estimate_att() at leaf_budget = |leaves(tau)|
#'
#' The SAME fit supplies both the tied estimator (via its outcome tree alone) and
#' the \code{two_tree} comparator (its own theta/sigma/ci), so the two arms are on
#' byte-identical nuisance fits and any difference between them is the tie, not
#' the data or the tuning.
#'
#' DELIBERATE DEVIATION FROM A LITERAL READING OF SPEC §1, FLAGGED. Spec §1 says
#' mu-hat is "estimate_att()'s m0_model/m0_hat". The tied estimator here is built
#' on the fitted tree's PARTITION with its leaf values RE-COMPUTED as exact
#' within-leaf control means, rather than on the serialised
#' \code{nuisance_fits$outcome_control} vector, because optimaltrees rounds a leaf
#' value to six significant digits (see \code{MU_LEAF_TOL}). The two objects are
#' the same estimator -- leaf-constant control means on the fitted partition, spec
#' §1's own description -- and differ by <= 4e-07 in mu-hat and hence in
#' theta_tie, far below any Monte Carlo scale in this study. Without the refit,
#' spec §1's 1e-10 identity target is unreachable for a reason that has nothing to
#' do with \code{cor:single-tree}. \code{mu_leaf_gap} reports the discrepancy so
#' the substitution is visible in every results table.
st_arm_fitted <- function(d, dgp) {
  fit <- doubletree::estimate_att(
    X = d$X, A = d$A, Y = d$Y,
    leaf_budget = dgp$leaf_budget,
    outcome_type = OUTCOME_TYPE,
    lambda_n = LAMBDA_N, m_n = M_N,
    verbose = FALSE, worker_limit = WORKER_LIMIT
  )
  leaf <- st_fitted_leaf(fit, d$X)
  m0_hat <- st_leafwise_control_mean(d$Y, d$A, leaf)
  # The gap against the package's own serialised predictions, on the SAME leaves.
  serial_gap <- max(abs(m0_hat - fit$nuisance_fits$outcome_control))
  res <- st_tied_estimator(d$Y, d$A, m0_hat, leaf)
  res$mu_leaf_gap <- serial_gap
  res$realized_partition <- st_realized_partition(leaf, d, dgp)
  # In A-C there is one candidate, so "realised the candidate" IS recovery. In E
  # there is no single intended tau, so recovery is not a meaningful yes/no and is
  # reported as NA rather than as a coin flip between the two ties.
  res$recovered <- if (length(dgp$leaf_maps) == 1L) {
    res$realized_partition == names(dgp$leaf_maps)[[1L]]
  } else {
    NA
  }
  res$fit <- fit
  res
}

# ---- 7 one replication --------------------------------------------------

#' The row schema every (replication, arm) pair returns
#'
#' Defined once so success and failure rows have identical columns: a failure must
#' never be a shorter row that silently drops out of an rbind, and never a
#' plausible-looking NA estimate with no message.
st_row_template <- function() {
  data.frame(
    regime = NA_character_, dgp = NA_character_, n = NA_integer_,
    rep = NA_integer_, seed = NA_integer_, arm = NA_character_,
    theta0 = NA_real_, theta = NA_real_,
    se_naive = NA_real_, se_corrected = NA_real_, se_eif = NA_real_,
    cov_naive = NA, cov_corrected = NA, cov_eif = NA,
    ci_width_naive = NA_real_, ci_width_corrected = NA_real_,
    naive_over_correct = NA_real_,
    identity_gap = NA_real_,
    realized_partition = NA_character_, recovered = NA,
    n_leaves_m0 = NA_integer_, n_leaves_e = NA_integer_,
    n_treated = NA_integer_, n_control = NA_integer_,
    pi_hat = NA_real_, clip_frac = NA_real_, mu_leaf_gap = NA_real_,
    secs = NA_real_,
    n_warnings = NA_integer_, warning_messages = NA_character_,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Fill a tied-estimator row (oracle_tie / fitted_tie)
st_fill_tie_row <- function(row, res, theta0) {
  covers <- function(se) {
    abs(res$theta - theta0) <= Z_ALPHA * se
  }
  row$theta <- res$theta
  row$se_naive <- res$se_naive
  row$se_corrected <- res$se_corrected
  row$se_eif <- res$se_eif
  row$cov_naive <- covers(res$se_naive)
  row$cov_corrected <- covers(res$se_corrected)
  row$cov_eif <- covers(res$se_eif)
  row$ci_width_naive <- 2 * Z_ALPHA * res$se_naive
  row$ci_width_corrected <- 2 * Z_ALPHA * res$se_corrected
  row$naive_over_correct <- res$naive_over_correct
  row$identity_gap <- res$identity_gap
  row$realized_partition <- res$realized_partition
  row$recovered <- res$recovered
  row$n_leaves_m0 <- as.integer(res$n_leaves)
  row$pi_hat <- res$pi_hat
  row$clip_frac <- res$clip_frac
  row$mu_leaf_gap <- res$mu_leaf_gap
  row
}

#' Run one replication: draw once, run every arm of the regime on that draw
#'
#' One seed stream per replication (\code{set.seed(st_rep_seed(...))} then ONE
#' draw), so the oracle and fitted arms and the comparator are compared on
#' identical data. \code{estimate_att()} is called once and its result feeds both
#' \code{fitted_tie} and \code{two_tree}.
#'
#' The per-replication \code{identity_gap} assertion is deliberately FATAL under
#' \code{on_error = "stop"}: \code{cor:single-tree} is exact arithmetic, so a gap
#' above \code{IDENTITY_TOL} is an implementation bug in this harness, not a
#' statistical event to be averaged into a summary.
#'
#' @param regime One of \code{ST_REGIME_IDS}.
#' @param dgp A built DGP (\code{make_st_dgp()} output).
#' @param n Sample size.
#' @param rep Replication index.
#' @param on_error \code{"stop"} (interactive posture) or \code{"record"} (the
#'   condition message is stored so a failure RATE is computable).
#' @return A data.frame, one row per arm.
st_run_one_rep <- function(regime, dgp, n, rep, on_error = c("stop", "record")) {
  on_error <- match.arg(on_error)
  arms <- ST_REGIMES[[regime]]$arms
  seed <- st_rep_seed(regime, n, rep)

  blank <- function(arm) {
    row <- st_row_template()
    row$regime <- regime; row$dgp <- dgp$id; row$n <- as.integer(n)
    row$rep <- as.integer(rep); row$seed <- as.integer(seed)
    row$arm <- arm; row$theta0 <- dgp$theta0
    row
  }

  body <- function() {
    set.seed(seed)
    d <- dgp$draw(n)
    n_treated <- sum(d$A == 1L)
    n_control <- sum(d$A == 0L)

    warn_msgs <- character(0)
    t0 <- Sys.time()
    fitted <- withCallingHandlers(
      if (any(c("fitted_tie", "two_tree") %in% arms)) st_arm_fitted(d, dgp) else NULL,
      warning = function(w) {
        warn_msgs <<- c(warn_msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    secs_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    rows <- list()

    if ("oracle_tie" %in% arms) {
      t1 <- Sys.time()
      orc <- st_arm_oracle_tie(d, dgp)
      rows$oracle_tie <- st_fill_tie_row(blank("oracle_tie"), orc, dgp$theta0)
      rows$oracle_tie$secs <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    }
    if ("fitted_tie" %in% arms) {
      rows$fitted_tie <- st_fill_tie_row(blank("fitted_tie"), fitted, dgp$theta0)
      rows$fitted_tie$n_leaves_e <- as.integer(fitted$fit$n_leaves_e)
      rows$fitted_tie$secs <- secs_fit
    }
    if ("two_tree" %in% arms) {
      f <- fitted$fit
      row <- blank("two_tree")
      row$theta <- f$theta
      # The comparator's SE is the package's own EIF plug-in; the naive/corrected
      # split is a TIED-estimator concept and stays NA here, which is meaningful
      # (not applicable), not missing.
      row$se_eif <- f$sigma
      row$cov_eif <- f$ci_95[[1L]] <= dgp$theta0 && dgp$theta0 <= f$ci_95[[2L]]
      row$n_leaves_m0 <- as.integer(f$n_leaves_m0)
      row$n_leaves_e <- as.integer(f$n_leaves_e)
      row$realized_partition <- fitted$realized_partition
      row$recovered <- fitted$recovered
      row$pi_hat <- mean(d$A)
      row$secs <- secs_fit
      rows$two_tree <- row
    }

    for (arm in names(rows)) {
      g <- rows[[arm]]$identity_gap
      if (!is.na(g) && g > IDENTITY_TOL) {
        stop("cor:single-tree's identity FAILED in arm '", arm, "': |theta_tie - ",
             "theta_forced| = ", signif(g, 4), " > ", IDENTITY_TOL, ". That is ",
             "exact arithmetic, so this is a harness bug (leaf extraction, the ",
             "clip biting, or leaf values that are not within-leaf means), not a ",
             "statistical event.", call. = FALSE)
      }
    }

    out <- do.call(rbind, rows[intersect(arms, names(rows))])
    out$n_treated <- as.integer(n_treated)
    out$n_control <- as.integer(n_control)
    out$n_warnings <- length(warn_msgs)
    out$warning_messages <- if (length(warn_msgs)) {
      paste(unique(warn_msgs), collapse = " | ")
    } else {
      NA_character_
    }
    out
  }

  if (on_error == "stop") return(body())

  tryCatch(body(), error = function(cnd) {
    # One failure row PER ARM, so a failed cell keeps the shape of a successful
    # one. The whole replication fails together because every arm shares the draw.
    out <- do.call(rbind, lapply(arms, blank))
    out$error_message <- conditionMessage(cnd)
    out
  })
}

# ---- 8 one cell ----------------------------------------------------------

#' Checkpoint path for a cell
st_cell_path <- function(regime, n, reps) {
  file.path(DIR_RESULTS, sprintf("cell_%s_n%d_r%d.rds", regime, n, reps))
}

#' Run one (regime, n) cell and write its checkpoint
#'
#' Sequential by construction, \code{gc(full = TRUE)} on a heartbeat, no
#' future/furrr. A kill mid-run costs at most the cell in flight.
st_run_cell <- function(regime, n = PILOT_N, reps = PILOT_NSIM,
                        out = st_cell_path(regime, n, reps),
                        on_error = "record", heartbeat = 10L) {
  if (!regime %in% ST_REGIME_IDS) {
    cli::cli_abort("Unknown regime {.val {regime}}; expected {.val {ST_REGIME_IDS}}.")
  }
  stopifnot(length(n) == 1L, n >= 100, length(reps) == 1L, reps >= 1L)
  n <- as.integer(n); reps <- as.integer(reps)

  dgp <- make_st_dgp(ST_REGIMES[[regime]]$dgp)
  pop <- st_population_table(dgp)
  arms <- ST_REGIMES[[regime]]$arms

  cli::cli_h1("Regime {regime} / n = {n} / reps = {reps}")
  cli::cli_inform(c(
    "*" = "{dgp$label}",
    "*" = "claims: {ST_REGIMES[[regime]]$claims}",
    "*" = "theta_0 = {signif(dgp$theta0, 8)} (exact); pi = {signif(dgp$pi_pop, 6)}; Lbar = {dgp$leaf_budget} = |leaves(tau)|",
    "*" = "analytic: Tr = {signif(dgp$Tr, 6)}, C_untied = {signif(dgp$C_untied, 6)}, C_tied = {paste(sprintf('%s:%.6f', pop$partition, pop$C_tied), collapse = ' ')}",
    "*" = "analytic V_tau/V = {paste(sprintf('%s:%.6f', pop$partition, pop$ratio_Vtau_V), collapse = ' ')}; naive/correct = {paste(sprintf('%s:%.4f', pop$partition, pop$naive_over_correct), collapse = ' ')}",
    "*" = "arms: {paste(arms, collapse = ', ')}"
  ))

  rows <- vector("list", reps)
  for (r in seq_len(reps)) {
    rows[[r]] <- st_run_one_rep(regime, dgp, n, r, on_error = on_error)
    if (heartbeat > 0L && r %% heartbeat == 0L) {
      done <- do.call(rbind, rows[seq_len(r)])
      gc(full = TRUE, verbose = FALSE)
      cli::cli_inform(paste0(
        "  ", r, "/", reps, "  max identity gap = ",
        signif(max(done$identity_gap, na.rm = TRUE), 3),
        "  recovery = ",
        signif(mean(done$recovered[done$arm == "fitted_tie"], na.rm = TRUE), 3),
        "  s/rep = ", signif(sum(done$secs, na.rm = TRUE) / r, 3)
      ))
    }
  }
  res <- do.call(rbind, rows)
  rownames(res) <- NULL

  payload <- list(
    results = res, population = pop,
    meta = list(regime = regime, n = n, reps = reps, arms = arms,
                seed_master = SEED_MASTER, dgp_label = dgp$label,
                leaf_budget = dgp$leaf_budget,
                r_version = R.version.string, when = Sys.time())
  )
  saveRDS(payload, out)

  n_fail <- length(unique(res$rep[!is.na(res$error_message)]))
  cli::cli_alert_success(
    "wrote {.file {basename(out)}} -- {reps} reps x {length(arms)} arms, {n_fail} failed rep(s), {signif(sum(res$secs, na.rm = TRUE) / 60, 3)} min"
  )
  wm <- unique(stats::na.omit(res$warning_messages))
  if (length(wm)) {
    cli::cli_inform(c(
      "i" = "{sum(res$n_warnings, na.rm = TRUE)} muffled warning(s), {length(wm)} distinct:",
      stats::setNames(substr(wm, 1, 150), rep("*", length(wm)))
    ))
  }
  em <- unique(stats::na.omit(res$error_message))
  if (length(em)) {
    cli::cli_warn(c(
      "{n_fail} replication(s) failed, {length(em)} distinct message(s):",
      stats::setNames(substr(em, 1, 300), rep("x", length(em)))
    ))
  }
  invisible(res)
}
