# threshold_superconsistency

**Go/no-go check, run 2026-08-31.** Not a production simulation study.

## Question

The paper claims the estimator locates each split threshold to within
`O_p(1/n)` of the truth, with no pre-specified grid and no cross-fitting. That
claim rests on the true nuisance function having a genuine **jump** at each
split boundary (not a slope). Does it actually hold in practice?

The estimator under test is two-stage:

* **Stage A** — fit a single best tree with the existing `optimaltrees` solver at
  its coarse default discretization (`ceiling(n^(1/3))` quantile bins per
  coordinate). This recovers the topology and brackets each threshold, but its
  cut can only ever be a grid point, so its error is floored at the grid mesh.
* **Stage B** (new, R-level, `code/stage_b.R`) — for each split, take the two
  neighbouring grid points as a bracket and scan **every observed value** of that
  coordinate inside it, keeping whichever exactly minimizes empirical risk. No
  grid refinement. Attainable resolution becomes the spacing between order
  statistics, `~1/n`.

If the jump mechanism bites, the log-log slope of `mean |t_hat - t*|` against `n`
is near **-1**. If it does not, the error stays grid-limited near **-1/3**.

## Design

| | |
|---|---|
| DGP | `p = 3`, `X_j ~ U(0,1)` iid, `X3` pure noise |
| True tree | depth 2, 3 leaves: root `X1 <= 0.45`; inside the left child `X2 <= 0.62` |
| Leaf means | 0.3 / 0.0 / 0.6 — spaced exactly `kappa = 0.3`, so the **minimum** jump across each split boundary is 0.3 |
| Noise | `Y = mu(X) + N(0, 0.25^2)` |
| Loss | `squared_error` (single regression tree; no ATT/EIF machinery needed for this question) |
| `n` | 2000, 8000, 32000 |
| Reps | 8 per cell, sequential |
| `lambda` | `log(n)/n` (theory rate) |
| `max_depth` | 2, set explicitly |

`t1* = 0.45` and `t2* = 0.62` are deliberately **not** round quantiles of
`U(0,1)`: if a true threshold coincided with a grid point, Stage A alone would
localize it at the `n^{-1/2}` rate of a sample quantile and the diagnostic could
not distinguish "Stage B beat the grid" from "the grid happened to be right".

## Result

| n | grid mesh | Stage A mean err | Stage B mean err | Stage B median err | Stage B mean, bracket-valid |
|---:|---:|---:|---:|---:|---:|
| 2000 | 0.0769 | 0.0209 | 5.94e-04 | 3.53e-04 | 5.94e-04 |
| 8000 | 0.0500 | 0.0146 | 3.01e-04 | 1.85e-04 | 3.01e-04 |
| 32000 | 0.0312 | 0.0198 | 2.42e-03 | 2.28e-05 | 3.38e-05 |

log-log slopes for split 1 (root, X1):

* **-1.03** (bracket-valid reps, 22/24) — **primary**
* **-0.99** (median)
* **+0.51** (unrestricted mean — contaminated, see below)

split 2 (X2): **-1.18** (mean), **-1.18** (median).

**VERDICT: superconsistency confirmed.** Both splits give slopes near -1.
Topology was recovered in 24/24 replicates.

### The one caveat, stated plainly

In 2 of 24 replicates (n=32000, reps 3 and 4) the solver's root cut was the
*second*-nearest grid point rather than the nearest, so the ±1-grid-point bracket
did **not** contain `t1*`. Stage B was then clipped: it returned the upper bracket
edge (within 1e-4 of it) and could not reach the truth. Those two reps alone flip
the unrestricted mean slope from -1 to +0.5.

This is a limitation of the **Stage A -> Stage B handoff** (bracket too narrow),
not evidence against the jump mechanism — the bracket and the clipping are
recorded per replicate (`t1_lo`, `t1_hi`, `bracket1_contains_truth`) so the
distinction is checkable, not asserted. The clipped reps are reported in every
table rather than dropped.

Both clipped reps have 4 leaves, i.e. the solver used a second `X1` cut and
offset the root one mesh low. **Recommended fix before relying on Stage B:**
widen the bracket (±2 grid points), or bracket over the union of all cuts the
tree uses on that coordinate. Also note the Stage A control slope came out at
-0.02 rather than the expected -1/3, for the same reason: at larger n the root
cut is increasingly *not* the tree's estimate of `t1*`.

## Memory

This codebase has previously exhausted host memory during tree simulations
(`../docs/MEMORY_SAFE_SIMULATIONS.md`, `MEMORY.md` `[LEARN:rashomon-memory]`).
Rules followed here: **no Rashomon enumeration** (single best tree only), **no
parallelism** (`worker_limit = 1`, threads pinned to 1), explicit `max_depth = 2`,
`gc(full = TRUE)` after every replicate, one **process per cell**, checkpoint per
cell, and a staged rollout (micro-test -> 2000 -> 8000 -> 32000 with an
inspection between each).

Peak `phys_footprint` per cell: **170 MB** (n=2000), **223 MB** (n=8000),
**411 MB** (n=32000). Steady state ~200 MB and flat across replicates — no
growth trend. Wall clock: 10s / 29s / 294s per cell.

Memory is measured with `phys_footprint` (`code/proc_footprint.cpp`, copied from
`single_tree_coverage/`), not RSS: macOS RSS excludes compressed pages. `ps` is
unavailable to this session entirely, so subprocess-based measurement is not an
option here anyway.

## Files

```
code/dgp.R              jump-tree DGP + population split gains
code/stage_b.R          fitted-tree -> coordinate space; the exact bracket scan
code/check_stage_b.R    self-check: partition reconstruction + scan vs naive rescoring
code/run_cell.R         one cell per process; Rscript code/run_cell.R <n> <reps>
code/analyze.R          per-cell means, log-log slopes, bracket diagnostics, verdict
code/proc_footprint.cpp macOS phys_footprint probe
results/cell_n<N>_r8.rds        per-cell checkpoints
results/microtest_n2000_r2.rds  staged-rollout micro-test
results/summary_slope.rds       slopes, verdicts, bracket misses, memory
```

## Reproduce

```sh
Rscript code/check_stage_b.R                  # verify Stage B is exact
Rscript code/run_cell.R 2000  2               # micro-test first, always
Rscript code/run_cell.R 2000  8
Rscript code/run_cell.R 8000  8
Rscript code/run_cell.R 32000 8
Rscript code/analyze.R
```

Seeded once per cell (`20260831 + n`); re-running reproduces the numbers above
exactly (verified).

**Independence from the concurrent `rho` change.** `code/run_cell.R` passes
`discretize_bins = ceiling(n^(1/3))` **explicitly** rather than relying on
`discretize_bins = "adaptive"`. A separate uncommitted change raises the package's
adaptive exponent `rho` from 1/3 to 0.45; reading the default would have silently
changed what "Stage A at the coarse resolution" means, and at n=32000 would have
produced 318 binary features instead of 93.
