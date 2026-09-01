# ============================================================
# check_stage_b.R
# Study: threshold_superconsistency  (doubletree)
#
# Correctness self-check for the incremental Stage B scan. The scan in
# stage_b.R updates per-leaf sufficient statistics as the candidate cut moves,
# which is fast but easy to get subtly wrong (off-by-one on the migration
# boundary, wrong branch orientation, empty-leaf handling). This re-scores every
# candidate the slow, obviously-correct way (rebuild the whole partition, re-fit
# leaf means, sum squared errors) on a small problem and demands the two agree.
#
# Also checks the branch orientation end-to-end: the coordinate-space tree read
# back out of a fitted model must reproduce the fitted model's own predictions.
#
# Run directly:  Rscript code/check_stage_b.R
# ============================================================

suppressPackageStartupMessages({
  library(optimaltrees)
  library(cli)
})

STUDY <- "/Users/dagniel/RAND/rprojects/global-scholars/doubletree/simulations/threshold_superconsistency"
source(file.path(STUDY, "code", "dgp.R"))
source(file.path(STUDY, "code", "stage_b.R"))

set.seed(20260831)

n <- 1500
dat <- generate_jump_tree(n)
X <- dat$X; y <- dat$Y

fit <- fit_tree(
  X, y,
  loss_function       = "squared_error",
  regularization      = log(n) / n,
  discretize_method   = "quantiles",
  discretize_bins     = 8L,
  max_depth           = 2L,
  store_training_data = TRUE,
  verbose             = FALSE
)

md <- fit@discretization_metadata
lk <- bin_lookup(md)
ct <- as_coord_tree(fit@trees[[1L]], lk, md)

cli::cli_h2("1. Coordinate-space tree reproduces the fitted model's predictions")
lid <- leaf_index(ct, X)
K <- count_leaves(ct)
fitted_means <- tapply(y, lid, mean)
recon <- as.numeric(fitted_means[as.character(lid)])
pkg_pred <- as.numeric(predict(fit, X))

# (a) EXACT partition check, tolerance-free. If the reconstructed tree induced a
# different partition -- in particular if a branch were flipped -- some
# reconstructed leaf would straddle two of the package's predicted values. This
# is the check with teeth; the numeric comparison below cannot distinguish a
# wrong partition from float rounding on its own.
per_leaf_distinct <- tapply(pkg_pred, lid, function(v) length(unique(v)))
n_distinct_pkg <- length(unique(pkg_pred))
cli::cli_inform("leaves = {K}; distinct package predictions = {n_distinct_pkg}; max distinct per reconstructed leaf = {max(per_leaf_distinct)}")
if (max(per_leaf_distinct) != 1L || n_distinct_pkg != K) {
  cli::cli_abort("Reconstructed partition differs from the fitted model's. Branch orientation or feature mapping is wrong.")
}

# (b) Numeric agreement. Tolerance is 1e-5, not machine epsilon: the solver
# stores leaf predictions in SINGLE precision, so predict() returns float-rounded
# means while the reconstruction re-computes them in double. Observed deviation is
# ~3e-7, i.e. float rounding, not a modelling difference.
max_dev <- max(abs(recon - pkg_pred))
cli::cli_inform("max |reconstructed - package predict| = {signif(max_dev, 3)} (float storage; tol 1e-5)")
if (max_dev > 1e-5) {
  cli::cli_abort("Reconstructed leaf values disagree with the fitted model beyond float rounding.")
}
cli::cli_alert_success("Reconstruction matches (partition exactly, values to float precision).")

cli::cli_h2("2. Incremental scan matches a naive full re-scoring, at every node")
paths <- split_paths(ct)
paths <- paths[order(vapply(paths, length, integer(1)))]

for (p in paths) {
  idx  <- rows_at(ct, X, p)
  node <- node_at(ct, p)
  br   <- node_bracket(ct, p)
  Xn <- X[idx, , drop = FALSE]; yn <- y[idx]
  xj <- Xn[[node$coord]]

  n_left    <- count_leaves(node$left)
  lid_left  <- leaf_index(node$left, Xn)
  lid_right <- n_left + leaf_index(node$right, Xn)
  Kn <- n_left + count_leaves(node$right)

  cands <- bracket_candidates(xj, node, br)
  fast  <- scan_cutoff(xj, yn, lid_left, lid_right, Kn, cands)

  # Naive: for each candidate, rebuild the local partition from scratch.
  naive_sse <- vapply(cands, function(c0) {
    nd <- node; nd$cut <- c0
    sub <- list(kind = "split", coord = nd$coord, cut = c0,
                left = nd$left, right = nd$right)
    tree_sse(sub, Xn, yn)
  }, numeric(1))
  naive_best <- cands[[which.min(naive_sse)]]

  d_sse <- abs(fast$sse - min(naive_sse))
  d_cut <- abs(fast$cut - naive_best)
  cli::cli_inform(paste0(
    "node '", if (length(p)) paste(p, collapse = "/") else "root",
    "' coord=", node$coord, " n=", length(idx),
    " bracket=(", signif(br$lo, 4), ",", signif(br$hi, 4), ")",
    " cands=", fast$n_candidates,
    " |dSSE|=", signif(d_sse, 3), " |dcut|=", signif(d_cut, 3)
  ))
  if (d_sse > 1e-8 * max(1, abs(min(naive_sse))) || d_cut > 0) {
    cli::cli_abort("Incremental scan disagrees with naive re-scoring at node '{paste(p, collapse='/')}'.")
  }
}
cli::cli_alert_success("Incremental scan is exact.")

cli::cli_h2("3. Refinement never increases training risk, and moves toward the truth")
ref <- refine_tree(ct, X, y)
sse_a <- tree_sse(ct, X, y)
sse_b <- tree_sse(ref$tree, X, y)
cli::cli_inform("SSE Stage A = {signif(sse_a, 8)}; Stage B = {signif(sse_b, 8)}")
if (sse_b > sse_a + 1e-8) cli::cli_abort("Stage B increased training risk; candidate set is wrong.")
print(ref$refined)
cli::cli_inform("true t1* = {T1_STAR}, t2* = {T2_STAR}")
cli::cli_alert_success("Self-check complete.")

cli::cli_h2("4. collapse_transitions(): synthetic transition-pair fixture")

# A hand-built coordinate-space tree on ONE coordinate ("Xc") with a 3-leaf
# grid-adjacent pair: root splits Xc at grid index k=5, its right child ALSO
# splits Xc at k=6 (adjacent), sandwiching a transition leaf at (tau5, tau6].
# n_bins irrelevant here -- only relative k matters -- so cuts are synthetic.
mk_leaf <- function() list(kind = "leaf")
mk_split <- function(coord, k, cut, left, right) {
  list(kind = "split", coord = coord, k = k, cut = cut,
       grid_cut = cut, grid_cuts = cut, collapsed = FALSE,
       left = left, right = right)
}

fixture <- mk_split(
  "Xc", k = 5L, cut = 0.5,
  left  = mk_leaf(),
  right = mk_split("Xc", k = 6L, cut = 0.6,
                    left  = mk_leaf(),   # the transition leaf, (0.5, 0.6]
                    right = mk_leaf())
)

out <- collapse_transitions(fixture)
cli::cli_inform("n_collapses = {out$n_collapses} (expect 1)")
if (out$n_collapses != 1L) cli::cli_abort("Fixture: expected exactly 1 collapse.")
ct_fix <- out$tree
if (!isTRUE(ct_fix$collapsed)) cli::cli_abort("Fixture: root should be the collapsed node.")
if (count_leaves(ct_fix) != 2L) cli::cli_abort("Fixture: collapsed tree should have 2 leaves, got {count_leaves(ct_fix)}.")
if (!isTRUE(all.equal(sort(ct_fix$grid_cuts), c(0.5, 0.6)))) {
  cli::cli_abort("Fixture: grid_cuts should span both original cutpoints (0.5, 0.6).")
}
if (!isTRUE(all.equal(ct_fix$cut, 0.55))) cli::cli_abort("Fixture: starting cut should be the midpoint 0.55.")
cli::cli_alert_success("collapse_transitions() collapses a grid-adjacent pair correctly.")

# Same-side-orientation mirror: LEFT child adjacent (k=6 parent, k=5 child).
fixture_left <- mk_split(
  "Xc", k = 6L, cut = 0.6,
  left  = mk_split("Xc", k = 5L, cut = 0.5,
                    left  = mk_leaf(),
                    right = mk_leaf()),  # transition leaf, (0.5, 0.6]
  right = mk_leaf()
)
out_left <- collapse_transitions(fixture_left)
if (out_left$n_collapses != 1L || count_leaves(out_left$tree) != 2L) {
  cli::cli_abort("Fixture (left orientation): collapse did not behave symmetrically.")
}
cli::cli_alert_success("collapse_transitions() handles the left-child orientation symmetrically.")

# Non-adjacent same-coordinate splits (k differs by >=2): must NOT collapse.
fixture_far <- mk_split(
  "Xc", k = 5L, cut = 0.5,
  left  = mk_leaf(),
  right = mk_split("Xc", k = 9L, cut = 0.9, left = mk_leaf(), right = mk_leaf())
)
out_far <- collapse_transitions(fixture_far)
if (out_far$n_collapses != 0L) cli::cli_abort("Fixture: non-adjacent same-coordinate splits must not collapse.")
cli::cli_alert_success("collapse_transitions() leaves genuinely distinct boundaries untouched.")

# Degenerate orientation (right child cuts BELOW parent) must abort loudly,
# not silently mishandle.
fixture_degenerate <- mk_split(
  "Xc", k = 6L, cut = 0.6,
  left  = mk_leaf(),
  right = mk_split("Xc", k = 5L, cut = 0.5, left = mk_leaf(), right = mk_leaf())
)
degenerate_ok <- tryCatch({
  collapse_transitions(fixture_degenerate)
  FALSE
}, error = function(e) TRUE)
if (!degenerate_ok) cli::cli_abort("Fixture: degenerate right-child-cuts-below-parent must abort, did not.")
cli::cli_alert_success("collapse_transitions() aborts loudly on a degenerate (wrong-side) split.")

# Non-leaf transition child (chained collapse) must abort loudly.
fixture_chained <- mk_split(
  "Xc", k = 5L, cut = 0.5,
  left  = mk_leaf(),
  right = mk_split("Xc", k = 6L, cut = 0.6,
                    left  = mk_split("Xc", k = 5L, cut = 0.55,  # non-leaf transition child
                                      left = mk_leaf(), right = mk_leaf()),
                    right = mk_leaf())
)
chained_ok <- tryCatch({
  collapse_transitions(fixture_chained)
  FALSE
}, error = function(e) TRUE)
if (!chained_ok) cli::cli_abort("Fixture: chained (non-leaf transition child) collapse must abort, did not.")
cli::cli_alert_success("collapse_transitions() aborts loudly on a chained (three-way) collapse.")

cli::cli_h2("5. Full pipeline: collapse_transitions() wired ahead of refine_tree() on the fitted tree")
ct2 <- collapse_transitions(as_coord_tree(fit@trees[[1L]], lk, md))$tree
ref2 <- refine_tree(ct2, X, y)
sse_a2 <- tree_sse(ct2, X, y)
sse_b2 <- tree_sse(ref2$tree, X, y)
cli::cli_inform("post-collapse SSE Stage A = {signif(sse_a2, 8)}; Stage B = {signif(sse_b2, 8)}; n_leaves = {count_leaves(ct2)}")
print(ref2$refined)
cli::cli_alert_success("Full pipeline (collapse -> refine) runs end to end on the fitted tree.")
