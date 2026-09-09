# ============================================================
# common.R
# Study: head_to_head_comparison  (doubletree, S5)
# Spec:  quality_reports/specs/2026-09-08_head-to-head-comparison.md
#
# Shared harness: paths, package loading, the arm registry, the regime grid,
# deterministic seeds, one replication, one cell. Sourced by verify_dgps.R /
# run_pilot.R / run_cell.R / run_sweep.R / analyze.R. Sourcing this file has no
# effect beyond defining objects and creating the study's output directories.
#
# Structure follows partition_recovery_clt/code/common.R (S1) and
# honest_inference_sparsity_failure/code/common.R (S2): same
# source-tree-vs-installed guard, same deterministic seeding scheme, same
# checkpoint-per-cell posture, same worker_limit = 1L memory rule.
#
# WHAT IS DIFFERENT FROM S1/S2, AND WHY
#
# (a) LONG RESULTS FORMAT, one row per (regime, dgp, n, rep, ARM). S1 and S2 have
#     a fixed set of estimators per replication and can afford one wide row. Here
#     the ARM SET VARIES BY REGIME (5 arms in R1, 2 in R2, 1 in R6), so a wide
#     schema would be mostly NA and every analysis would have to know which
#     columns are meaningful in which regime. Long format makes "which arms ran
#     in this cell" a fact in the data rather than a fact in a comment.
#
# (b) AN ARM REGISTRY. Every estimator is a closure in `ARMS` with one signature,
#     so adding an arm cannot change the harness and the harness cannot special-
#     case an arm. The two DERIVED arms (thm:anchor intervals, built from a
#     flagship/anchor PAIR rather than from data) are declared in `DERIVED_ARMS`
#     and constructed after the base arms, not smuggled into the loop.
#
# (c) NO ENUMERATION for DGP-A / A2 / C. S1's build_dgp() attaches the exact
#     sufficient class and selection margin, which this study never uses -- it
#     reports coverage and RMSE ratios, not partition recovery. The enumeration at
#     Lbar = 16 (R6's largest budget) is also combinatorially prohibitive. S1's
#     check_variable_roles() IS still run on every spec, unconditionally: it is
#     cheap and it is the check that catches a mis-keyed nuisance, which has no
#     other symptom than quantities of quietly the wrong size. DGP-B keeps its
#     full build_dgp_dial() because at Lbar = 2 the enumeration is 6 partitions
#     and it yields `pseudo$bias`, the PREDICTED bias R3 is read against.
#
# ------------------------------------------------------------
# THE TUNING CHARTER (spec §2), AND THE ASYMMETRY IT DOES NOT HIDE
# ------------------------------------------------------------
#
# 1. SHARED DATA. One seed stream per replication: set.seed(rep_seed(...)) then
#    ONE draw, used by every arm in that replication. Paired comparisons (the
#    RMSE ratio to oracle-AIPW, C3's flagship-vs-crossfit contrast) are therefore
#    on identical data, not independently resampled data per method.
#
# 2. SHARED FOLDS. One `fold_seed` per replication, passed to EVERY cross-fitted
#    arm. att_linear(), att_forest() and att_lasso_dml() use byte-identical
#    stratified-fold code, so an identical seed gives them an IDENTICAL partition
#    and cross-fitting noise is not a nuisance confound between them.
#    estimate_att_crossfit() builds its folds internally by its own routine; it
#    receives the same seed, so it is reproducible, but its partition is NOT
#    guaranteed to coincide with the other three. That is a real, small residual
#    difference and is recorded here rather than claimed away.
#
# 3. SHARED CLIP. Every arm clips e-hat to [0.01, 0.99]. That is att_linear()'s
#    existing convention AND doubletree's own ass:construct(c) constant (asserted
#    against the package below), so for once there is no discrepancy to state.
#    att_oracle() carries the same clip but ASSERTS it is inert -- clipping a true
#    propensity would bias the efficiency denominator every C1 ratio is built on.
#
# 4. CV-TUNING IS ASYMMETRIC AND THAT IS STATED, NOT HIDDEN. Spec §2 asks that
#    every method be either CV-tuned or none be, per DGP, and then concedes the
#    asymmetry is unavoidable. The arms partition into:
#      CV-tuned:  estimate_att_crossfit() (cv_regularization = TRUE, shipped
#                 default), att_lasso_dml() (cv.glmnet over lambda)
#      untuned:   att_linear() x2 (a GLM has no penalty), att_forest() (ranger
#                 defaults), att_oracle() (nothing is fit)
#      neither:   estimate_att(), whose leaf_budget is FIXED PER DGP in advance at
#                 the value that makes that DGP exactly sufficient -- not tuned by
#                 any data-driven procedure -- matching the paper's own framing
#                 that the analyst sets the budget.
#    No arm's tuning was chosen after seeing results. See README.md, which repeats
#    this so a reader of the tables does not have to find it in code.
#
# 5. leaf_budget IS A DGP PROPERTY, NOT A TUNED PARAMETER. It is set once per DGP
#    in `DGP_REGISTRY` (A: 4, A2: 4, B: 2, C: 4) at the value that makes that DGP
#    exactly representable, and is only VARIED in R6, where varying it IS the
#    experiment (the (2*Lbar+1)/n coverage boundary).
#
# ------------------------------------------------------------
# WHAT THIS STUDY IS NOT (spec §0) -- read before any number is quoted
# ------------------------------------------------------------
#
# This is NOT a leaderboard. Under standard regularity a correctly-specified
# linear working model, a correctly-specified black-box AIPW/DML estimator, and
# doubletree under exact structural sparsity ALL attain the same semiparametric
# efficiency bound (thm:regular -- doubletree attains the FULL nonparametric
# bound, not something better than it). No arm can win on efficiency
# asymptotically, and the paper claims none does. The three licensed claims are:
#   C1  under exact sparsity, doubletree MATCHES an oracle-nuisance benchmark,
#       while a linear working model WITHOUT the right basis is misspecified;
#   C2  as sparsity degrades, doubletree's plain Wald CI degrades but the anchor
#       interval retains nominal coverage;
#   C3  at matched n, the full-sample flagship is at least as precise as
#       doubletree's own cross-fit fallback.
# Regimes where a competitor's working model happens to be correctly specified
# (notably DGP-B, which is additive on the probability scale at EVERY eps -- see
# spec §3) are NOT evidence against doubletree. Nothing beyond "matched X" or
# "retained coverage where Y did not" is licensed by this design.
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

STUDY_DIR   <- file.path(PKG_ROOT, "simulations", "head_to_head_comparison")
DIR_CODE    <- file.path(STUDY_DIR, "code")
DIR_RESULTS <- file.path(STUDY_DIR, "results")
DIR_FIGURES <- file.path(STUDY_DIR, "figures")
DIR_TABLES  <- file.path(STUDY_DIR, "tables")
for (.d in c(DIR_RESULTS, DIR_FIGURES, DIR_TABLES)) {
  if (!dir.exists(.d)) dir.create(.d, recursive = TRUE)
}

# The installed doubletree in the user library can predate the current source
# tree, so load the SOURCE explicitly -- same posture as S1/S2's common.R.
.OPTIMALTREES_SRC <- file.path(dirname(PKG_ROOT), "optimaltrees")

# TWO PACKAGE-LOADING PATHS, ONE GATE (added 2026-09-09 for the SLURM deployment;
# simulations/head_to_head_comparison/slurm/README_O2.md documents the cluster side).
#
#   H2H_USE_INSTALLED unset/"0"  DEV path (default, unchanged behaviour):
#       pkgload::load_all() on this source tree and on the sibling
#       ../optimaltrees source tree. Correct on the dev box, where the source is
#       what is being developed and an installed copy may be stale.
#
#   H2H_USE_INSTALLED="1"        CLUSTER path: library() against packages that
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
USE_INSTALLED_PKGS <- Sys.getenv("H2H_USE_INSTALLED", "0") == "1"
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

# Fail fast rather than at replication 40: read the signatures, do not assume.
for (.arg in c("leaf_budget", "propensity_loss", "outcome_type", "lambda_n", "m_n")) {
  if (!.arg %in% names(formals(doubletree::estimate_att))) {
    cli::cli_abort(c(
      "Loaded {.fun doubletree::estimate_att} has no {.arg {.arg}} argument.",
      i = "A stale installed copy is probably shadowing the source tree."
    ))
  }
}
for (.arg in c("K", "outcome_type", "cv_regularization", "max_depth", "seed")) {
  if (!.arg %in% names(formals(doubletree::estimate_att_crossfit))) {
    cli::cli_abort("Loaded {.fun doubletree::estimate_att_crossfit} has no {.arg {.arg}} argument.")
  }
}

# ---- 0b reused code from the other studies (READ ONLY) ------------------

#' Source another study's file read-only, naming the dependency on failure
#'
#' An opaque "object 'dgp_spec_grid_sparse' not found" 200 lines later is the
#' failure mode this exists to prevent.
#'
#' @param rel Path relative to the doubletree package root.
#' @param provides Function names the file must define.
#' @return Invisibly, the absolute path sourced.
source_reused <- function(rel, provides) {
  path <- file.path(PKG_ROOT, rel)
  if (!file.exists(path)) {
    cli::cli_abort(c(
      "Reused source {.file {rel}} not found.",
      i = "This study is required by its spec (§3) to reuse it, not reimplement it."
    ))
  }
  source(path, local = FALSE)
  absent <- provides[!vapply(provides, function(f) is.function(get0(f)), logical(1))]
  if (length(absent)) {
    cli::cli_abort("{.file {rel}} defined no {.fun {absent}}; its contract has changed.")
  }
  invisible(path)
}

# S1's enumeration helpers (make_cell_grid is needed by every DGP spec).
source_reused(
  file.path("simulations", "partition_recovery_clt", "code",
            "enumerate_sufficient_class.R"),
  c("make_cell_grid", "enumerate_tree_partitions", "project_on_partition",
    "excess_risk", "is_sufficient_partition", "fitted_cell_labels",
    "canonical_labels")
)
# ST3 delegates its sample draw to the shipped generator; must precede S1's
# dgps.R being USED (not sourced), since dgp_spec_weak_overlap()'s draw closure
# calls it.
source_reused(
  file.path("simulations", "dgps", "dgps_stress.R"),
  "generate_dgp_weak_overlap"
)
# S1's DGPs: DGP-A2 (= F1) and DGP-C (= ST3), plus new_dgp_spec()/build_dgp()/
# check_variable_roles()/check_spec_vs_draw(), which this study's own DGP-A reuses.
source_reused(
  file.path("simulations", "partition_recovery_clt", "code", "dgps.R"),
  c("dgp_spec_grid_sparse", "dgp_spec_weak_overlap", "new_dgp_spec",
    "check_variable_roles", "check_spec_vs_draw")
)
# S2's DGP-B (the eps dial) plus population_bias()/bilinear_ip(), which this
# study's DGP-A verification reuses rather than re-deriving.
source_reused(
  file.path("simulations", "honest_inference_sparsity_failure", "code",
            "dgps_delta_dial.R"),
  c("dgp_spec_delta_dial", "build_dgp_dial", "population_bias", "bilinear_ip",
    "approximation_error", "dgp_population_summary_dial")
)
# This study's ONE new DGP.
source_reused(
  file.path("simulations", "head_to_head_comparison", "code", "dgps.R"),
  c("dgp_spec_shared_interaction", "main_effects_projection",
    "main_effects_bilinear_remainder", "verify_main_effects_remainder",
    "confounding_check_shared_interaction")
)

# The competitor arms.
source_reused(file.path("simulations", "methods", "method_linear.R"), "att_linear")
source_reused(file.path("simulations", "methods", "method_forest.R"), "att_forest")
source_reused(file.path("simulations", "methods", "method_oracle_aipw.R"), "att_oracle")
source_reused(file.path("simulations", "methods", "method_lasso_dml.R"), "att_lasso_dml")

# att_linear(interactions = TRUE) guards `ncol(X) <= 5` before using `Y ~ .^2`
# (spec §1's MUST DO: the guard was `<= 4`, which would have silently produced
# MAIN-EFFECTS results under the "interactions" label on every p = 5 DGP here).
# The fix is in method_linear.R; this asserts it, because a silent revert would
# reintroduce exactly the mislabelling the spec warns about and nothing else
# would notice.
local({
  src <- paste(deparse(att_linear), collapse = " ")
  n_ok <- lengths(regmatches(src, gregexpr("interactions && ncol\\(X\\) <= 5", src)))
  if (n_ok != 2L) {
    cli::cli_abort(c(
      "{.fun att_linear} does not guard {.code interactions && ncol(X) <= 5} in both places ({n_ok} found).",
      i = "Every DGP in this study has p = 5; a {.code <= 4} guard reports MAIN-EFFECTS results under the {.val glm_int} label.",
      i = "Spec §1, MUST DO."
    ))
  }
})

# ---- 1 fixed design constants -------------------------------------------

SEED_MASTER <- 20260908L   # the spec's own date

M_N          <- 1L    # eq:feasible floor; estimate_att()'s own default
LAMBDA_N     <- NULL  # NULL -> log(n)/n inside estimate_att() (prop:parsimony rate)
OUTCOME_TYPE <- "binary"
PROPENSITY_LOSS <- "log_loss"   # estimate_att()'s SHIPPED default

# Cross-fitting folds, shared by EVERY cross-fitted arm (tuning charter item 2).
K_FOLDS <- 5L

# estimate_att_crossfit()'s shipped default depth, used wherever it acts as the
# anchor. R5 deliberately breaks this (see REGIMES).
ANCHOR_MAX_DEPTH <- 4L
# R5's deliberately UNDERPOWERED anchor. F1's nuisances are hierarchical --
# e_0 splits X1 then X2 WITHIN X1 = 1 -- so they need depth 2. A depth-1 anchor
# cannot represent either nuisance at any n, so it cannot be sqrt(n)-consistent
# and ass:rate (thm:anchor's own load-bearing assumption) FAILS by construction.
# That is the point: thm:anchor's guarantee is conditional on the anchor, and R5
# checks what happens when the condition is not met, rather than gesturing at it.
ANCHOR_MAX_DEPTH_SHALLOW <- 1L

ALPHA <- 0.05
Z_ALPHA <- stats::qnorm(1 - ALPHA / 2)

# Memory safety (this project's convention; simulations/docs/
# MEMORY_SAFE_SIMULATIONS.md and MEMORY.md [LEARN:rashomon-memory]): single best
# tree only, no Rashomon enumeration, no parallelism. worker_limit = 1L reaches
# fit_tree() through both doubletree estimators' `...`; PARALLEL_CV = FALSE
# switches off cv_regularization_adaptive()'s furrr/future fold-parallelism,
# which is ON by default.
WORKER_LIMIT <- 1L
PARALLEL_CV <- FALSE

# ranger's thread count is pinned to 1 inside att_forest() (num.threads = 1L).
# Spec §5 requires this explicitly or the timing comparison is meaningless:
# ranger is multithreaded by default and would otherwise be timed against
# single-threaded competitors. Asserted rather than trusted.
local({
  src <- paste(deparse(att_forest), collapse = " ")
  if (!grepl("num.threads = 1", src, fixed = TRUE)) {
    cli::cli_abort(c(
      "{.fun att_forest} does not pin {.code num.threads = 1}.",
      i = "Spec §5: an unpinned thread count makes cross-method timing meaningless."
    ))
  }
})

# The propensity clip of ass:construct(c), applied inside both doubletree
# estimators AND (by the same constants) inside every competitor arm. Copied
# because the package constants are unexported, then CHECKED against the package:
# a package change would otherwise make the charter's "shared clip" claim false
# without any error.
CLIP_LO <- 0.01
CLIP_HI <- 0.99
local({
  pkg_lo <- get0(".PROPENSITY_LOWER_BOUND", envir = asNamespace("doubletree"))
  pkg_hi <- get0(".PROPENSITY_UPPER_BOUND", envir = asNamespace("doubletree"))
  if (is.null(pkg_lo) || is.null(pkg_hi) ||
      !isTRUE(all.equal(c(pkg_lo, pkg_hi), c(CLIP_LO, CLIP_HI)))) {
    cli::cli_abort(c(
      "doubletree's propensity clip bounds are not ({CLIP_LO}, {CLIP_HI}).",
      i = "The tuning charter's 'same clip for every arm' claim would be false."
    ))
  }
})

# ---- 2 the DGP registry -------------------------------------------------

# One entry per DGP id. `leaf_budget` is a DGP PROPERTY here (charter item 5),
# fixed at the value that makes that DGP exactly representable -- except in R6,
# where the A2_L* variants vary it because varying it IS the experiment.
#
# `build` returns a spec; `enumerate` says whether to run the (expensive at large
# Lbar, and unused by this study) sufficient-class enumeration.
DGP_REGISTRY <- list(
  # C1's DGP. New here; see code/dgps.R for why F1 could not be reused.
  A = list(
    label = "DGP-A: shared-interaction, gate-swapped (exact sparsity, GLM-main-effects misspecified)",
    leaf_budget = 4L, enumerate = FALSE,
    build = function() dgp_spec_shared_interaction(leaf_budget = 4L, id = "A")
  ),
  # C3's and R6's DGP. S1's F1, REUSED UNMODIFIED at its post-2026-09-01
  # favourable calibration. Its population bilinear remainder under
  # main-effects misspecification is exactly zero (spec §3's correction), which
  # is why no GLM arm runs on it -- neither C3 nor R6 needs one.
  A2 = list(
    label = "DGP-A2 = F1 (reused unmodified): grid-exact sparsity, shared-X1 confounding",
    leaf_budget = 4L, enumerate = FALSE,
    build = function() dgp_spec_grid_sparse(leaf_budget = 4L, id = "A2",
                                           label = "F1 (reused)")
  ),
  # R6's Lbar sweep on the SAME DGP. Lbar = 4 is the exactly-sufficient budget;
  # 8 and 16 are deliberately over-budget, which is what pushes (2 Lbar + 1)/n
  # toward the coverage boundary at small n.
  A2_L4 = list(
    label = "F1, Lbar = 4 (exactly sufficient)",
    leaf_budget = 4L, enumerate = FALSE,
    build = function() dgp_spec_grid_sparse(leaf_budget = 4L, id = "A2_L4",
                                           label = "F1, Lbar = 4")
  ),
  A2_L8 = list(
    label = "F1, Lbar = 8 (over-budget)",
    leaf_budget = 8L, enumerate = FALSE,
    build = function() dgp_spec_grid_sparse(leaf_budget = 8L, id = "A2_L8",
                                           label = "F1, Lbar = 8")
  ),
  A2_L16 = list(
    label = "F1, Lbar = 16 (far over-budget)",
    leaf_budget = 16L, enumerate = FALSE,
    build = function() dgp_spec_grid_sparse(leaf_budget = 16L, id = "A2_L16",
                                           label = "F1, Lbar = 16")
  ),
  # C2's DGP. S2's eps dial, REUSED UNMODIFIED, aligned residual mode, Lbar = 2.
  # eps grid is spec §3/§4's {0, 0.05, 0.10, 0.15}, not S2's own
  # {0, 0.03, 0.09, 0.15}; both satisfy dgp_spec_delta_dial()'s own ceiling
  # (eps^2 < Var_nu(mu_step) = 0.0356), which the constructor asserts.
  #
  # HONEST FRAMING, FIXED IN ADVANCE (spec §3, DGP-B). e_0 = e_step(X1) + eps*g(X4)
  # is additive on the PROBABILITY scale at every eps, so a linear-probability or
  # near-linear-logit main-effects model is approximately correctly specified
  # here at every eps. The GLM arm should NOT degrade as eps grows, and a
  # favourable-to-GLM finding here is NOT a doubletree failure. What R3 is
  # evidence for is narrower and is stated before the run: doubletree's PLAIN CI
  # degrades as eps grows (Lbar = 2 cannot represent the truth, which needs 4
  # leaves) while the ANCHOR interval recovers nominal coverage.
  B_eps0 = list(
    label = "DGP-B eps = 0.00 (exact sparsity)",
    leaf_budget = 2L, enumerate = TRUE,
    build = function() dgp_spec_delta_dial(eps = 0.00, residual_mode = "aligned",
                                          leaf_budget = 2L, id = "B_eps0")
  ),
  B_eps05 = list(
    label = "DGP-B eps = 0.05 (small violation)",
    leaf_budget = 2L, enumerate = TRUE,
    build = function() dgp_spec_delta_dial(eps = 0.05, residual_mode = "aligned",
                                          leaf_budget = 2L, id = "B_eps05")
  ),
  B_eps10 = list(
    label = "DGP-B eps = 0.10 (moderate violation)",
    leaf_budget = 2L, enumerate = TRUE,
    build = function() dgp_spec_delta_dial(eps = 0.10, residual_mode = "aligned",
                                          leaf_budget = 2L, id = "B_eps10")
  ),
  B_eps15 = list(
    label = "DGP-B eps = 0.15 (severe violation)",
    leaf_budget = 2L, enumerate = TRUE,
    build = function() dgp_spec_delta_dial(eps = 0.15, residual_mode = "aligned",
                                          leaf_budget = 2L, id = "B_eps15")
  ),
  # R4's stress DGP. S1's ST3, REUSED UNMODIFIED: p = 4,
  # e_0 = plogis(-2.5 + 4 X1 + 2.5 X2), so e_0(1,1) = 0.9820 and 1 - e_0 = 0.0180
  # in that cell. At n = 500 that leaf holds ~125 rows with an expected CONTROL
  # count of ~2, so a tree leaf's empirical e-hat (or w-hat = e/(1-e)) can be
  # driven to an extreme or ill-defined value -- a failure mode a smooth e-hat
  # from GLM/RF essentially never hits at the same rate.
  #
  # CALIBRATION IS NOT PRE-EMPTIVELY ESCALATED (user decision 4, 2026-09-09):
  # ST3's shipped coefficients are used AS-IS. Whether they actually produce a
  # truly ill-defined weight (as opposed to merely a large one) at the n values
  # used here is MEASURED by the pilot -- run_pilot.R reports the clip fraction,
  # the empirical min control count in the extreme cell, and max(e-hat) -- and the
  # observation is recorded in README.md either way. Escalation (steeper
  # coefficients) happens only if the pilot shows it is needed.
  C = list(
    label = "DGP-C = ST3 (reused unmodified): weak overlap, e_0 up to 0.982",
    leaf_budget = 4L, enumerate = FALSE,
    build = function() dgp_spec_weak_overlap(leaf_budget = 4L, id = "C",
                                            label = "ST3 weak overlap (reused)")
  )
)

DGP_IDS <- names(DGP_REGISTRY)

#' Build one DGP by id
#'
#' \code{check_variable_roles()} (S1's) runs on every spec, unconditionally: it
#' verifies on the exact cell grid that each nuisance depends on exactly the
#' coordinates it declares. A mis-keyed nuisance is a perfectly self-consistent
#' DGP that errors nowhere and shows up only as a bilinear remainder of quietly
#' the wrong size -- the exact defect S1's own design audit caught.
#'
#' The full \code{build_dgp()} / \code{build_dgp_dial()} enumeration runs only
#' where \code{enumerate = TRUE} (DGP-B): this study reports coverage and RMSE
#' ratios, not partition recovery, so S1's sufficient-class machinery is not
#' needed for A/A2/C -- and at R6's \code{Lbar = 16} it would be combinatorially
#' prohibitive. DGP-B keeps it because at \code{Lbar = 2} it is 6 partitions and
#' it yields \code{pseudo$bias}, the PREDICTED bias R3's results are read against.
#'
#' @param id One of \code{DGP_IDS}.
#' @return A built DGP spec.
make_dgp <- function(id) {
  if (!id %in% DGP_IDS) {
    cli::cli_abort("Unknown dgp {.val {id}}; expected one of {.val {DGP_IDS}}.")
  }
  entry <- DGP_REGISTRY[[id]]
  spec <- entry$build()
  if (isTRUE(entry$enumerate)) {
    spec <- build_dgp_dial(spec)   # also runs check_variable_roles_dial()
  } else {
    check_variable_roles(spec)
  }
  spec$study_label <- entry$label
  spec
}

# ---- 3 the arm registry -------------------------------------------------

# Every BASE arm has the identical signature
#   function(d, dgp, fold_seed) -> list(theta, sigma, ci, extras)
# so the harness cannot special-case an arm and adding an arm cannot change the
# harness. `extras` is a named list of scalars merged into the result row.
#
# NO SILENT FALLBACKS anywhere below: an estimator failure propagates to
# run_one_rep(), which either stops (pilot posture) or records the condition
# MESSAGE (sweep posture). Never a plausible-looking NA with no explanation.

#' Diagnostics common to every arm that produces a fitted propensity
#'
#' Clip fraction is reported for EVERY arm in EVERY regime, not just R4, so
#' "coverage is off here" can be attributed to the clip biting, to
#' misspecification, or to neither.
e_hat_diagnostics <- function(e_hat) {
  list(
    clip_frac = mean(e_hat <= CLIP_LO + 1e-12 | e_hat >= CLIP_HI - 1e-12),
    max_e_hat = max(e_hat),
    min_e_hat = min(e_hat)
  )
}

ARMS <- list(

  # The paper's own estimator: full sample, NO splitting. leaf_budget comes from
  # the DGP (charter item 5), never from a data-driven search.
  doubletree = function(d, dgp, fold_seed) {
    fit <- doubletree::estimate_att(
      X = d$X, A = d$A, Y = d$Y,
      leaf_budget = dgp$leaf_budget,
      outcome_type = OUTCOME_TYPE,
      lambda_n = LAMBDA_N, m_n = M_N,
      propensity_loss = PROPENSITY_LOSS,
      verbose = FALSE, worker_limit = WORKER_LIMIT
    )
    list(theta = fit$theta, sigma = fit$sigma, ci = fit$ci_95,
         extras = c(list(
           n_leaves_e = as.integer(fit$n_leaves_e),
           n_leaves_m0 = as.integer(fit$n_leaves_m0),
           certified_e = fit$certified_e, certified_m0 = fit$certified_m0,
           leaf_budget = dgp$leaf_budget,
           lambda_n_used = fit$lambda_n
         ), e_hat_diagnostics(fit$nuisance_fits$propensity)),
         raw = fit)
  },

  # C3's comparison partner AND C2's anchor: doubletree's own cross-fit fallback.
  doubletree_crossfit = function(d, dgp, fold_seed) {
    crossfit_arm(d, dgp, fold_seed, max_depth = ANCHOR_MAX_DEPTH)
  },

  # R5's deliberately underpowered anchor (max_depth = 1: cannot represent F1's
  # depth-2 hierarchical nuisances at any n, so ass:rate fails by construction).
  doubletree_crossfit_shallow = function(d, dgp, fold_seed) {
    crossfit_arm(d, dgp, fold_seed, max_depth = ANCHOR_MAX_DEPTH_SHALLOW)
  },

  # C1's EFFICIENCY DENOMINATOR. True e_0, mu_0 in the same eq:score. No
  # estimation, hence no folds and no fold_seed.
  oracle_aipw = function(d, dgp, fold_seed) {
    if (is.null(d$e_true) || is.null(d$mu0_true)) {
      stop("DGP '", dgp$id, "' draw() returns no e_true/mu0_true, so the oracle ",
           "arm cannot be run on it.", call. = FALSE)
    }
    fit <- att_oracle(d$X, d$A, d$Y, e_true = d$e_true, m0_true = d$mu0_true)
    list(theta = fit$theta, sigma = fit$sigma, ci = fit$ci,
         extras = c(list(n_clipped = fit$hyperparams$n_clipped),
                    e_hat_diagnostics(d$e_true)),
         raw = fit)
  },

  # "What an analyst actually does."
  glm_main = function(d, dgp, fold_seed) {
    linear_arm(d, fold_seed, interactions = FALSE)
  },

  # "The diligent analyst." Requires method_linear.R's ncol(X) <= 5 guard,
  # asserted at load above.
  glm_int = function(d, dgp, fold_seed) {
    linear_arm(d, fold_seed, interactions = TRUE)
  },

  # sec:lasso-comparison's empirical instance.
  lasso_dml = function(d, dgp, fold_seed) {
    fit <- att_lasso_dml(d$X, d$A, d$Y, K = K_FOLDS, seed = fold_seed,
                         interactions = TRUE)
    list(theta = fit$theta, sigma = fit$sigma, ci = fit$ci,
         extras = list(
           n_basis = fit$hyperparams$n_basis,
           nz_e_mean = mean(fit$hyperparams$nz_e),
           nz_m0_mean = mean(fit$hyperparams$nz_m0),
           lambda_e_mean = mean(fit$hyperparams$lambda_e),
           lambda_m0_mean = mean(fit$hyperparams$lambda_m0)
         ),
         raw = fit)
  },

  # The black-box ML arm the user explicitly requested. GBM is deferred (spec §7):
  # RF already discharges "black-box ML" and GBM's tuning surface would invite an
  # unresolvable "you tuned my competitor badly" objection.
  forest = function(d, dgp, fold_seed) {
    fit <- att_forest(d$X, d$A, d$Y, K = K_FOLDS, seed = fold_seed)
    list(theta = fit$theta, sigma = fit$sigma, ci = fit$ci,
         extras = list(num_trees = fit$hyperparams$num.trees,
                       mtry = fit$hyperparams$mtry),
         raw = fit)
  }
)

#' Shared body of the two crossfit arms (they differ only in max_depth)
crossfit_arm <- function(d, dgp, fold_seed, max_depth) {
  fit <- doubletree::estimate_att_crossfit(
    X = d$X, A = d$A, Y = d$Y,
    K = K_FOLDS, outcome_type = OUTCOME_TYPE,
    max_depth = max_depth,
    seed = fold_seed,
    verbose = FALSE, worker_limit = WORKER_LIMIT, parallel = PARALLEL_CV
  )
  list(theta = fit$theta, sigma = fit$sigma, ci = fit$ci_95,
       extras = c(list(max_depth = max_depth, K = fit$K),
                  e_hat_diagnostics(fit$nuisance_fits$propensity)),
       raw = fit)
}

#' Shared body of the two GLM arms (they differ only in the basis)
linear_arm <- function(d, fold_seed, interactions) {
  fit <- att_linear(d$X, d$A, d$Y, K = K_FOLDS, seed = fold_seed,
                    interactions = interactions)
  list(theta = fit$theta, sigma = fit$sigma, ci = fit$ci,
       extras = list(interactions = interactions,
                     glm_convergence = fit$convergence),
       raw = fit)
}

# DERIVED arms: built from a PAIR of already-computed base arms, not from data.
# Declared here rather than handled inline so the harness's loop stays one loop
# over base arms and the dependency is explicit.
DERIVED_ARMS <- list(
  # thm:anchor, LITERAL: theta_full +/- (|theta_full - theta_anchor| + z*sigma_anchor).
  # No bias-aware critical value -- this is the triangle inequality the theorem
  # actually proves and nothing else.
  #
  # sigma_anchor is ALREADY a standard error (att_se() returns
  # sqrt(mean(psi^2)/n)), so it is NOT divided by sqrt(n) a second time. S2's
  # common.R records why that matters: dividing again shrinks the interval by
  # sqrt(n) and manufactures catastrophic undercoverage that looks like a
  # theorem failure.
  doubletree_anchor = list(
    from = c("doubletree", "doubletree_crossfit")
  ),
  doubletree_anchor_shallow = list(
    from = c("doubletree", "doubletree_crossfit_shallow")
  )
)

#' The literal thm:anchor interval
#'
#' @param theta_full,theta_anchor Point estimates.
#' @param sigma_anchor The ANCHOR's standard error (not an sd).
#' @param z Normal critical value.
#' @return Numeric length-2 (lower, upper).
anchor_ci_harness <- function(theta_full, theta_anchor, sigma_anchor, z = Z_ALPHA) {
  half <- abs(theta_full - theta_anchor) + z * sigma_anchor
  c(theta_full - half, theta_full + half)
}

# ---- 4 the regime grid --------------------------------------------------

# Spec §4's table, verbatim in structure. `n` here is the PILOT grid; the full
# sweep's n values and nsim are in N_GRID_FULL / NSIM_FULL below and are NOT what
# run_pilot.R uses (spec §6's nsim = 1000/2000 is a cluster/batch job).
REGIMES <- list(
  R1 = list(
    label = "R1: exact sparsity, favourable -- C1 (no interpretability-for-efficiency cost)",
    evidence = "C1",
    dgps = "A",
    n = c(500L, 1000L, 2000L),
    arms = c("doubletree", "oracle_aipw", "glm_main", "glm_int", "lasso_dml"),
    # The C1 statistic is the RMSE RATIO to this arm, with its own MC SE --
    # never a raw "beats GLM" comparison (spec §0, §5).
    ratio_ref = "oracle_aipw"
  ),
  R2 = list(
    label = "R2: full-sample flagship vs its own cross-fit fallback -- C3",
    evidence = "C3",
    dgps = "A2",
    n = c(500L, 2000L),
    arms = c("doubletree", "doubletree_crossfit"),
    # C3 is the one place a "win" IS a paper claim (the whole point of not
    # splitting), so the ratio reference is the fallback, not the oracle.
    ratio_ref = "doubletree_crossfit"
  ),
  R3 = list(
    label = "R3: degrading structural sparsity (eps dial) -- C2 (anchor retains coverage)",
    evidence = "C2",
    dgps = c("B_eps0", "B_eps05", "B_eps10", "B_eps15"),
    n = 500L,
    arms = c("doubletree", "doubletree_crossfit", "doubletree_anchor",
             "glm_main", "oracle_aipw"),
    ratio_ref = "oracle_aipw"
  ),
  R4 = list(
    label = "R4: ATT overlap/boundary stress -- adversarial (a), doubletree-specific inference failure",
    evidence = "stress (a)",
    dgps = "C",
    n = c(500L, 1000L),
    # No oracle arm: ST3's e_0 reaches 0.982, and att_oracle() refuses to run
    # with a binding clip (a clipped oracle is not an oracle). e_0 = 0.982 is
    # still inside (0.01, 0.99) so it would in fact run -- but the RMSE ratio
    # would be to a near-degenerate-weight benchmark, which is not the C1
    # statistic and would invite over-reading. Spec §4 lists flagship/GLM/RF.
    arms = c("doubletree", "glm_main", "forest"),
    ratio_ref = NULL
  ),
  R5 = list(
    label = "R5: anchor rate-condition violation -- adversarial (b), thm:anchor's own assumption",
    evidence = "stress (b)",
    # BOTH DGPs, and the pair is the point (spec §4 offers "A2 or B (reuse)"; the
    # 2026-09-09 pilot showed only the pair is informative).
    #
    # On A2 the flagship is UNBIASED (exact sparsity), so a bad anchor only makes
    # |delta_hat| = |theta_full - theta_anchor| LARGER, which WIDENS the anchor
    # interval and preserves coverage. Measured at n = 500, 30 reps: coverage 1.00
    # for both the shipped-depth and the depth-1 anchor, with the shallow anchor's
    # interval the WIDER of the two (0.231 vs 0.223). thm:anchor's own
    # conservatism absorbs the ass:rate violation, so A2 alone cannot exhibit the
    # failure mode this regime exists to look for.
    #
    # thm:anchor can only be stressed where the FLAGSHIP is biased AND the anchor
    # is bad simultaneously -- the anchor interval is then centred on a biased
    # theta_full with a |delta_hat| that no longer measures that bias. B_eps15 is
    # the flagship-biased case (predicted bias 0.1105 at Lbar = 2, and the pilot
    # measured plain-Wald coverage 0.667 there), so pairing it with the depth-1
    # anchor is the actual experiment.
    dgps = c("A2", "B_eps15"),
    n = 500L,
    # BOTH anchors run on the same data: the underpowered one (the experiment)
    # and the shipped-default one (the reference that shows the shortfall is the
    # anchor's depth and not the DGP). Spec §4 asks only for the underpowered
    # anchor; the default is added because without a same-data reference an
    # anchor-coverage shortfall cannot be attributed. Flagged as a deliberate
    # addition in README.md.
    arms = c("doubletree",
             "doubletree_crossfit_shallow", "doubletree_anchor_shallow",
             "doubletree_crossfit", "doubletree_anchor"),
    ratio_ref = NULL
  ),
  R6 = list(
    label = "R6: (2*Lbar+1)/n coverage boundary -- adversarial (c), theory-predicted under-coverage",
    evidence = "stress (c)",
    dgps = c("A2_L4", "A2_L8", "A2_L16"),
    n = c(200L, 500L),
    arms = "doubletree",
    ratio_ref = NULL
  )
)

REGIME_IDS <- names(REGIMES)

# Spec §6's targets, recorded here so the pilot's extrapolation has something to
# extrapolate TO and so nobody has to re-read the spec to size the real sweep.
# NOT used by run_pilot.R.
NSIM_FULL <- c(R1 = 1000L, R2 = 2000L, R3 = 1000L, R4 = 1000L, R5 = 1000L, R6 = 2000L)
# R2 and R6 get the full 200-8000 sweep in the real run (spec §4); the pilot uses
# the subsets in REGIMES above.
N_GRID_FULL <- list(
  R1 = c(500L, 1000L, 2000L, 4000L),
  R2 = c(200L, 500L, 1000L, 2000L, 4000L, 8000L),
  R3 = c(500L, 2000L),
  R4 = c(500L, 1000L, 2000L),
  R5 = c(500L, 2000L),
  R6 = c(200L, 500L, 1000L, 2000L)
)

#' The full pilot design grid, one row per (regime, dgp, n) cell
design_grid <- function(regimes = REGIME_IDS) {
  rows <- lapply(regimes, function(r) {
    spec <- REGIMES[[r]]
    g <- expand.grid(n = spec$n, dgp = spec$dgps,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    g <- g[order(match(g$dgp, spec$dgps), g$n), , drop = FALSE]
    data.frame(regime = r, dgp = g$dgp, n = g$n,
               n_arms = length(spec$arms), evidence = spec$evidence,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# ---- 5 deterministic seeds ----------------------------------------------

#' Per-replicate seed: a pure function of (SEED_MASTER, regime, dgp, n, rep)
#'
#' Same scheme as S1's and S2's common.R, so the whole sweep is a function of one
#' integer. 7 hex digits < 2^28 is always a valid R integer.
#'
#' The REGIME is in the key on purpose: R2 and R5 share DGP-A2 at n = 500, and if
#' they shared the seed they would run on byte-identical data, making their
#' results correlated in a way no reader would expect from two separately
#' reported regimes.
rep_seed <- function(regime, dgp, n, rep) {
  key <- paste(SEED_MASTER, regime, dgp, n, rep, sep = "|")
  strtoi(substr(digest::digest(key, algo = "xxhash64"), 1, 7), base = 16L)
}

#' The ONE fold seed shared by every cross-fitted arm in a replication
#'
#' Charter item 2. Derived from the replication seed so it is a function of the
#' same single integer, but distinct from it so the fold draw is not the data
#' draw's own stream continued.
fold_seed_for <- function(seed) {
  strtoi(substr(digest::digest(paste0("folds|", seed), algo = "xxhash64"), 1, 7),
         base = 16L)
}

# ---- 6 one replication --------------------------------------------------

#' The row schema every (replication, arm) pair returns
#'
#' Defined once so success rows and failure rows are guaranteed to have the same
#' columns: a failure must never be a shorter row that silently drops out of an
#' rbind, and must never be a plausible-looking NA estimate with no message.
#'
#' Arm-specific diagnostics (\code{n_leaves_e}, \code{nz_e_mean}, ...) live in a
#' single JSON-ish \code{extras} string rather than as sparse columns: with 8 arms
#' whose diagnostics barely overlap, sparse columns would be ~90% NA and every
#' reader would have to know which arm makes which column meaningful. The handful
#' of diagnostics that EVERY analysis needs (leaves, certification, clip
#' fraction) are promoted to real columns.
rep_row_template <- function() {
  data.frame(
    regime = NA_character_, dgp = NA_character_, n = NA_integer_,
    rep = NA_integer_, seed = NA_integer_, fold_seed = NA_integer_,
    arm = NA_character_,
    theta0 = NA_real_, theta = NA_real_, sigma = NA_real_,
    ci_lo = NA_real_, ci_hi = NA_real_, ci_width = NA_real_, covered = NA,
    # Effective precision for the ATT scales with n_treated, not n (spec §5), so
    # it is a first-class column in every row rather than something to be
    # reconstructed from n and pi.
    n_treated = NA_integer_, n_control = NA_integer_,
    # doubletree only; NA for every competitor arm, which is meaningful (they
    # have no leaves), not missing.
    n_leaves_e = NA_integer_, n_leaves_m0 = NA_integer_,
    certified_e = NA, certified_m0 = NA, leaf_budget = NA_integer_,
    # Recorded for EVERY arm in EVERY regime, so a coverage anomaly can be
    # attributed to the clip biting rather than guessed at.
    clip_frac = NA_real_, max_e_hat = NA_real_, min_e_hat = NA_real_,
    secs = NA_real_,
    extras = NA_character_,
    # Warnings are COUNTED AND KEPT, not discarded.
    n_warnings = NA_integer_, warning_messages = NA_character_,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Flatten an arm's `extras` list into one short string
#'
#' Deliberately not JSON (no new dependency) and deliberately not columns (see
#' \code{rep_row_template()}). Round-trippable by eye, which is all any reader
#' needs from a diagnostic field.
encode_extras <- function(x) {
  if (!length(x)) return(NA_character_)
  keep <- setdiff(names(x), c("n_leaves_e", "n_leaves_m0", "certified_e",
                              "certified_m0", "leaf_budget",
                              "clip_frac", "max_e_hat", "min_e_hat"))
  if (!length(keep)) return(NA_character_)
  paste(vapply(keep, function(k) {
    v <- x[[k]]
    paste0(k, "=", paste(signif_if_num(v), collapse = ","))
  }, character(1)), collapse = "; ")
}

signif_if_num <- function(v) {
  if (is.numeric(v)) format(signif(v, 5), trim = TRUE) else as.character(v)
}

#' Run one replication: draw once, run every arm of the regime on that draw
#'
#' One seed stream per replication (charter item 1): \code{set.seed(seed)} then
#' ONE draw. Every base arm then runs on that same \code{d}, and every
#' cross-fitted arm receives the same \code{fold_seed} (charter item 2).
#'
#' The seed is re-set to \code{seed} before EACH arm. Without that, an arm's own
#' internal RNG consumption would shift the stream every later arm sees, so
#' adding or reordering an arm would silently change other arms' results and the
#' regimes would stop being independently reproducible. The DATA is unaffected
#' either way (it is drawn before any arm runs); this is about the arms.
#'
#' @param regime One of \code{REGIME_IDS}.
#' @param dgp A built DGP (\code{make_dgp()} output).
#' @param n Sample size.
#' @param rep Replication index.
#' @param on_error \code{"stop"} (pilot/interactive posture: a failure surfaces at
#'   its origin) or \code{"record"} (long-run posture: the condition message is
#'   stored so the cell finishes and a failure RATE can be reported).
#' @return A data.frame with one row per arm, matching \code{rep_row_template()}.
run_one_rep <- function(regime, dgp, n, rep, on_error = c("stop", "record")) {
  on_error <- match.arg(on_error)
  arms <- REGIMES[[regime]]$arms
  base_arms <- intersect(arms, names(ARMS))
  derived <- intersect(arms, names(DERIVED_ARMS))
  unknown <- setdiff(arms, c(names(ARMS), names(DERIVED_ARMS)))
  if (length(unknown)) {
    cli::cli_abort("Regime {.val {regime}} names unknown arm(s) {.val {unknown}}.")
  }
  for (da in derived) {
    need <- DERIVED_ARMS[[da]]$from
    if (!all(need %in% base_arms)) {
      cli::cli_abort(c(
        "Derived arm {.val {da}} needs base arms {.val {need}}.",
        i = "Regime {.val {regime}} runs only {.val {base_arms}}."
      ))
    }
  }

  seed <- rep_seed(regime, dgp$id, n, rep)
  fold_seed <- fold_seed_for(seed)

  blank <- function(arm) {
    row <- rep_row_template()
    row$regime <- regime; row$dgp <- dgp$id; row$n <- as.integer(n)
    row$rep <- as.integer(rep); row$seed <- as.integer(seed)
    row$fold_seed <- as.integer(fold_seed)
    row$arm <- arm; row$theta0 <- dgp$theta0
    row
  }

  body <- function() {
    set.seed(seed)
    d <- dgp$draw(n)
    n_treated <- sum(d$A == 1L)
    n_control <- sum(d$A == 0L)

    fitted <- list()
    rows <- list()

    for (arm in base_arms) {
      warn_msgs <- character(0)
      set.seed(seed)   # see the roxygen note: arms must not shift each other's stream
      t0 <- Sys.time()
      res <- withCallingHandlers(
        ARMS[[arm]](d, dgp, fold_seed),
        warning = function(w) {
          warn_msgs <<- c(warn_msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      if (!is.finite(res$theta) || !is.finite(res$sigma) || res$sigma <= 0) {
        stop("Arm '", arm, "' returned theta = ", res$theta, ", sigma = ",
             res$sigma, "; a non-finite or non-positive SE cannot produce an ",
             "interval.", call. = FALSE)
      }
      fitted[[arm]] <- res
      rows[[arm]] <- fill_row(blank(arm), res$theta, res$sigma, res$ci,
                              dgp$theta0, n_treated, n_control, res$extras,
                              secs, warn_msgs)
    }

    for (da in derived) {
      need <- DERIVED_ARMS[[da]]$from
      full <- fitted[[need[[1]]]]
      anch <- fitted[[need[[2]]]]
      ci <- anchor_ci_harness(full$theta, anch$theta, anch$sigma)
      # A DETERMINISTIC implication of the construction, not a statistical
      # expectation: the anchor interval CONTAINS the anchor's own Wald interval
      # for every possible draw, so it must cover whenever that one does. A
      # violation is an implementation bug, so it stops the replication rather
      # than being averaged into a coverage number. (S2's common.R records the
      # sqrt(n) mistake this catches.)
      covers <- function(z) z[[1]] <= dgp$theta0 && dgp$theta0 <= z[[2]]
      if (covers(anch$ci) && !covers(ci)) {
        stop("Derived arm '", da, "' failed to cover while its anchor's own Wald ",
             "interval covered. Impossible by the triangle inequality, so ",
             "anchor_ci_harness() is misconstructed. theta_full = ", full$theta,
             ", theta_anchor = ", anch$theta, ", sigma_anchor = ", anch$sigma,
             ", theta0 = ", dgp$theta0, call. = FALSE)
      }
      rows[[da]] <- fill_row(
        blank(da), full$theta, anch$sigma, ci, dgp$theta0,
        n_treated, n_control,
        list(built_from = paste(need, collapse = "+"),
             delta = full$theta - anch$theta),
        NA_real_, character(0)
      )
    }

    do.call(rbind, rows[arms])
  }

  if (on_error == "stop") return(body())

  tryCatch(body(), error = function(cnd) {
    # One failure row PER ARM, so a failed cell keeps the same shape as a
    # successful one and a per-arm failure rate is computable. The whole
    # replication fails together because every arm shares the draw.
    out <- do.call(rbind, lapply(arms, blank))
    out$error_message <- conditionMessage(cnd)
    out
  })
}

#' Populate one result row
fill_row <- function(row, theta, sigma, ci, theta0, n_treated, n_control,
                     extras, secs, warn_msgs) {
  row$theta <- theta
  row$sigma <- sigma
  row$ci_lo <- ci[[1]]
  row$ci_hi <- ci[[2]]
  row$ci_width <- ci[[2]] - ci[[1]]
  row$covered <- ci[[1]] <= theta0 && theta0 <= ci[[2]]
  row$n_treated <- as.integer(n_treated)
  row$n_control <- as.integer(n_control)
  for (fld in c("n_leaves_e", "n_leaves_m0", "certified_e", "certified_m0",
                "leaf_budget", "clip_frac", "max_e_hat", "min_e_hat")) {
    if (!is.null(extras[[fld]])) row[[fld]] <- extras[[fld]]
  }
  row$extras <- encode_extras(extras)
  row$secs <- secs
  row$n_warnings <- length(warn_msgs)
  row$warning_messages <- if (length(warn_msgs)) {
    paste(unique(warn_msgs), collapse = " | ")
  } else {
    NA_character_
  }
  row
}

# ---- 7 one cell, sequentially, with a checkpoint -------------------------

#' Checkpoint path for a cell
cell_path <- function(regime, dgp_id, n, reps) {
  file.path(DIR_RESULTS, sprintf("cell_%s_%s_n%d_r%d.rds", regime, dgp_id, n, reps))
}

#' Run one (regime, dgp, n) cell and write its checkpoint
#'
#' Sequential by construction, \code{gc(full = TRUE)} on a heartbeat, no
#' future/furrr. A kill mid-sweep costs at most the cell in flight, because every
#' completed cell is already on disk.
#'
#' @param regime One of \code{REGIME_IDS}.
#' @param dgp_id One of the regime's DGPs.
#' @param n Sample size.
#' @param reps Replications.
#' @param out Checkpoint path.
#' @param on_error Passed to \code{run_one_rep()}.
#' @param heartbeat Progress/gc interval in replications; 0 disables.
#' @return The cell's results data.frame, invisibly.
run_cell <- function(regime, dgp_id, n, reps = 5L,
                     out = cell_path(regime, dgp_id, n, reps),
                     on_error = "record", heartbeat = 10L) {
  stopifnot(length(n) == 1L, n >= 50, length(reps) == 1L, reps >= 1L)
  n <- as.integer(n)
  reps <- as.integer(reps)
  if (!regime %in% REGIME_IDS) {
    cli::cli_abort("Unknown regime {.val {regime}}.")
  }
  if (!dgp_id %in% REGIMES[[regime]]$dgps) {
    cli::cli_abort("Regime {.val {regime}} does not use dgp {.val {dgp_id}}.")
  }

  dgp <- make_dgp(dgp_id)
  arms <- REGIMES[[regime]]$arms

  cli::cli_h1("{regime} / {dgp_id} / n = {n}, reps = {reps}, arms = {length(arms)}")
  cli::cli_inform(c(
    "*" = "{REGIMES[[regime]]$label}",
    "*" = "{dgp$study_label}",
    "*" = "theta_0 = {signif(dgp$theta0, 6)} (exact, closed form); e_0 in [{signif(dgp$e0_range[1], 4)}, {signif(dgp$e0_range[2], 4)}]; pi = {signif(dgp$p_treated, 4)}",
    "*" = "Lbar = {dgp$leaf_budget}; (2*Lbar+1)/n = {signif((2 * dgp$leaf_budget + 1) / n, 4)}; lambda_n = log(n)/n = {signif(log(n) / n, 4)}",
    "*" = "arms: {paste(arms, collapse = ', ')}"
  ))
  if (!is.null(dgp$pseudo)) {
    cli::cli_inform(c(
      "*" = "PREDICTED bias at the pseudo-true limits: {signif(dgp$pseudo$bias, 5)} (delta_e = {signif(dgp$class_e$delta_kl, 4)}, delta_mu = {signif(dgp$class_mu$delta_kl, 4)}, kl)"
    ))
  }

  rows <- vector("list", reps)
  for (r in seq_len(reps)) {
    rows[[r]] <- run_one_rep(regime, dgp, n, r, on_error = on_error)
    if (heartbeat > 0L && r %% heartbeat == 0L) {
      done <- do.call(rbind, rows[seq_len(r)])
      gc(full = TRUE, verbose = FALSE)
      cov_by_arm <- tapply(done$covered, done$arm, function(z) mean(z, na.rm = TRUE))
      cli::cli_inform(paste0(
        "  ", r, "/", reps, "  cov: ",
        paste(sprintf("%s=%.2f", names(cov_by_arm), as.numeric(cov_by_arm)),
              collapse = " "),
        "  s/rep=", signif(sum(done$secs, na.rm = TRUE) / r, 3)
      ))
    }
  }
  res <- do.call(rbind, rows)

  payload <- list(
    results = res,
    dgp_summary = dgp_summary_row(dgp),
    meta = cell_metadata(regime, dgp, n, reps, res)
  )
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(payload, out)

  n_fail <- length(unique(res$rep[!is.na(res$error_message)]))
  cli::cli_alert_success(
    "wrote {.file {basename(out)}} -- {reps} reps x {length(arms)} arms, {n_fail} failed rep(s), {signif(sum(res$secs, na.rm = TRUE) / 60, 3)} min"
  )
  # Distinct warning texts once per cell, not once per fit.
  wm <- unique(stats::na.omit(res$warning_messages))
  if (length(wm)) {
    cli::cli_inform(c(
      "i" = "{sum(res$n_warnings, na.rm = TRUE)} muffled warning(s), {length(wm)} distinct text(s):",
      stats::setNames(substr(wm, 1, 150), rep("*", length(wm)))
    ))
  }
  em <- unique(stats::na.omit(res$error_message))
  if (length(em)) {
    cli::cli_warn(c(
      "{n_fail} replication(s) failed, {length(em)} distinct message(s):",
      stats::setNames(substr(em, 1, 200), rep("x", length(em)))
    ))
  }
  invisible(res)
}

#' One-row population summary of a built DGP, common across DGP families
#'
#' The families' own summary functions (\code{dgp_population_summary()},
#' \code{dgp_population_summary_dial()}) return different columns, so this study
#' keeps its own minimal intersection plus the two things every table here needs
#' (theta_0 and the propensity range).
dgp_summary_row <- function(dgp) {
  data.frame(
    dgp = dgp$id, label = dgp$study_label,
    p = dgp$p, leaf_budget = dgp$leaf_budget,
    theta0 = dgp$theta0, pi_pop = dgp$p_treated,
    e0_min = dgp$e0_range[[1]], e0_max = dgp$e0_range[[2]],
    e_vars = paste(dgp$e_vars, collapse = "+"),
    mu_vars = paste(dgp$mu_vars, collapse = "+"),
    # Present only for DGP-B (the only family this study enumerates).
    bias_pseudo = if (is.null(dgp$pseudo)) NA_real_ else dgp$pseudo$bias,
    delta_e_kl = if (is.null(dgp$class_e)) NA_real_ else dgp$class_e$delta_kl,
    delta_mu_kl = if (is.null(dgp$class_mu)) NA_real_ else dgp$class_mu$delta_kl,
    stringsAsFactors = FALSE
  )
}

#' Provenance stamped into every checkpoint
cell_metadata <- function(regime, dgp, n, reps, res) {
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
    stem = "head_to_head_comparison",
    spec = "quality_reports/specs/2026-09-08_head-to-head-comparison.md",
    regime = regime, regime_label = REGIMES[[regime]]$label,
    evidence = REGIMES[[regime]]$evidence,
    arms = REGIMES[[regime]]$arms,
    ratio_ref = REGIMES[[regime]]$ratio_ref,
    dgp = dgp$id, dgp_label = dgp$study_label,
    n = n, reps = reps,
    seed_master = SEED_MASTER,
    seed_scheme = "rep_seed(regime, dgp, n, rep) = xxhash64(SEED_MASTER|regime|dgp|n|rep); fold_seed = xxhash64('folds|'||seed)",
    leaf_budget = dgp$leaf_budget, m_n = M_N,
    lambda_n_rate = "NULL -> log(n)/n", lambda_n_value = log(n) / n,
    propensity_loss = PROPENSITY_LOSS, outcome_type = OUTCOME_TYPE,
    K_folds = K_FOLDS,
    anchor_max_depth = ANCHOR_MAX_DEPTH,
    anchor_max_depth_shallow = ANCHOR_MAX_DEPTH_SHALLOW,
    alpha = ALPHA, z_alpha = Z_ALPHA,
    cv_tuned_arms = c("doubletree_crossfit", "doubletree_crossfit_shallow",
                      "lasso_dml"),
    untuned_arms = c("glm_main", "glm_int", "forest", "oracle_aipw"),
    budget_fixed_arms = "doubletree",
    tuning_asymmetry_note = paste(
      "CV-tuning is NOT symmetric across arms and is not claimed to be:",
      "estimate_att_crossfit() CV-tunes its regularisation (shipped default) and",
      "att_lasso_dml() CV-tunes lambda; att_linear() has no penalty to tune and",
      "att_forest() uses ranger defaults; estimate_att()'s leaf_budget is fixed",
      "per DGP in advance, not tuned. See simulations/head_to_head_comparison/README.md."
    ),
    worker_limit = WORKER_LIMIT, parallel_cv = PARALLEL_CV,
    ranger_threads = 1L, use_rashomon = FALSE,
    theta0 = dgp$theta0,
    dgp_params = dgp$params,
    e_vars = dgp$e_vars, mu_vars = dgp$mu_vars,
    clip = c(lo = CLIP_LO, hi = CLIP_HI),
    secs_total = sum(res$secs, na.rm = TRUE),
    n_failed_reps = length(unique(res$rep[!is.na(res$error_message)])),
    doubletree_sha = git_sha(PKG_ROOT),
    optimaltrees_sha = git_sha(.OPTIMALTREES_SRC),
    # Which loading path produced the code that ran (see section 0): "installed"
    # on the cluster, "source" on the dev box. A results file whose numbers are
    # argued over should say which of the two it came from.
    pkg_load_path = if (USE_INSTALLED_PKGS) "installed" else "source",
    doubletree_version = as.character(utils::packageVersion("doubletree")),
    r_version = R.version.string,
    optimaltrees_version = as.character(utils::packageVersion("optimaltrees")),
    glmnet_version = as.character(utils::packageVersion("glmnet")),
    ranger_version = as.character(utils::packageVersion("ranger")),
    run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}
