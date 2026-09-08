# Requirements Specification: empirical size of the value-recovery remainder terms

**Date:** 2026-09-08
**Status:** implemented
**Question:** At practical `n` (thousands) and leaf budget `Lbar` (10–30), are the
finite-sample remainders `T_a` and `T_eps` of the value-recovery inference route
actually `o_p(n^{-1/2})`, or do they behave as badly as the crude theoretical
bound suggests?

---

## Background

The value-recovery (L2-closeness, *not* exact-topology-recovery) route to valid
EIF-based ATT inference with single-tree (no cross-fitting) nuisances rests on
Step 2 of the full-sample CLT proof in
`doubletree/inst/paper/refining-grid-standalone.tex` (drafting history in
`_refining-grid-inferential-draft.tex`, lines 660–731). That step decomposes the
cross remainder exactly as `T = T_eps - T_a` and bounds each piece via a
localised empirical-process radius.

A proof audit found that the **crude** bound for these terms can be vacuous at
`Lbar >= 20`, `n = 2000`, once correctly-sourced VC-entropy constants are used.
Crude worst-case bounds are notoriously loose, so this study measures the terms
directly instead of arguing about the bound.

**Scope:** measurement only. This study does not attempt to repair the theory.

---

## Quantities under study

With `w_0 = e_0 / (1 - e_0)`, `w_hat = e_hat / (1 - e_hat)`,
`Delta_w = w_hat_{tau_e} - w_0`, `eps_i = Y_i - mu_0(Z_i)`,
`dnu(z) = (1 - e_0(z)) dP_Z(z)`, and `a_mu = Pi^nu_{tau_mu} mu_0 - mu_0`:

```
b_i    = Delta_w(Z_i) - mean{ Delta_w(Z_k) : A_k = 0, leaf_mu(k) = leaf_mu(i) }
T      = n^{-1} sum_{i : A_i = 0} b_i { Y_i - mu_hat_{tau_mu}(Z_i) }
T_eps  = n^{-1} sum_{i : A_i = 0} b_i eps_i
T_a    = n^{-1} sum_{i : A_i = 0} b_i a_mu(Z_i)
T      = T_eps - T_a                                     (exact algebraic identity)
||f||_mu^2 = E[(1 - e_0(Z)) f(Z)^2]
D_w = ||w_hat - w_0||_mu ,   D_mu = ||mu_hat - mu_0||_mu
```

Asymptotic linearity needs `sqrt(n) (T_eps - T_a) = o_p(1)`.

---

## Design (as implemented)

### DGP — well-specified, moderate SNR

| Element | Choice |
|---|---|
| Covariates | `p = 5`, `Z_j ~ iid Uniform(0, 1)` |
| Analyst grid `X_n = Q_n(Z)` | 2 thresholds per coordinate at `1/3, 2/3` → **10 binary features**, `3^5 = 243` equal-mass grid cells |
| True `e_0`, `mu_0` | axis-aligned trees with **exactly `Lbar` leaves**, all splits **on the grid** → `mu_{0,n} = mu_0`, mesh bias `h_n = 0` |
| Partition construction | deterministic: repeatedly split the largest-volume leaf on its widest coordinate at the median interior grid threshold |
| Distinct partitions | `mu` uses coordinate priority `(1,2,3,4,5)`, `e` uses `(3,4,5,1,2)`, so `tau_e != tau_mu` (a shared partition would make the within-leaf centring in `b_i` degenerate) |
| Outcome leaf values | cycle over `0.75 * {-1, -1/3, 1/3, 1}`; `sd(mu_0) ≈ 0.55–0.61` |
| Propensity leaf values | cycle over `{0.30, 0.45, 0.55, 0.70}`; `pi = 0.5`, overlap bounded away from 0/1 |
| Outcome model | `Y = mu_0(Z) + 1 * A + eps`, `eps ~ Uniform(-sqrt 3, sqrt 3)` (variance 1, **bounded**, so the theory's `B_Y` moment bound holds literally) |
| True ATT | `theta_0 = 1` exactly (constant effect) |
| Jump size | leaf-value jumps 0.5–1.5 noise SD ⇒ moderate SNR as specified |
| Min leaf mass | 11.1% (`Lbar = 5`) → 1.65% (`Lbar = 30`); at `n = 2000, Lbar = 30` that is ~33 units, ~16 controls per leaf |

### Cells

`n in {2000, 5000, 10000, 20000}` × `Lbar in {5, 10, 20, 30}` = 16 cells,
**500 replications each**.

### Fitting — no cross-fitting, well-specified regime

* `optimaltrees::fit_tree()` (exact GOSDT/TreeFARMS search), full sample.
* Propensity: all units, `loss_function = "log_loss"`.
* Control outcome: **control units only**, `loss_function = "squared_error"`
  (matching `doubletree/R/nuisance_trees.R:134-147`).
* Leaf values recomputed as within-leaf empirical means of the fitting sample.
  For both losses the leaf-optimal constant *is* that mean, so this reproduces
  the fitted values while making the exact refit identity
  `sum_{i in leaf, A_i = 0} (Y_i - mu_hat_i) = 0` hold to machine precision —
  which `T = T_eps - T_a` depends on. **Asserted every replication**
  (`IDENTITY_TOL = 1e-8`).
* Propensity clipped to `[0.01, 0.99]` (`doubletree/R/score_att.R` convention).
* Leaf budget `Lbar` is **enforced on every draw**: `lambda` is calibrated once
  per cell (geometric bisection to the largest `lambda` whose realised leaf count
  is `<= Lbar`), then raised geometrically on any individual draw that overshoots.
* Depth budget `= depth of the true partition` (3, 5, 6, 7 for
  `Lbar = 5, 10, 20, 30`), so the truth is inside the searched class.

**Documented restriction.** The searched class is
`{trees with <= Lbar leaves AND depth <= d(Lbar)}`, not the theory's full
`T_{Lbar,n}` (which allows depth up to `Lbar - 1`). This is a computational
necessity: exact optimal-tree search is superexponential in depth. Consequence
for the bound comparison: the empirical side comes from a *smaller* class than
the crude bound describes, so `r_n` (exact log-cardinality of the class actually
searched) is reported alongside `rho_n` (the audit's VC-entropy radius).

### Computed exactly (no Monte-Carlo error)

Because every fitted function is piecewise constant on the 243 grid cells and
the truth is too, all population integrals are exact finite sums over cells:
`Pi^nu_{tau_mu} mu_0`, `a_mu`, `||a_mu||_mu`, `D_w`, `D_mu`.

### Reported per cell

* mean/SD of `D_w`, `D_mu`, `||a_mu||_mu`; realised leaf counts
* mean/SD of `T_a`, `T_eps`
* mean/SD/median/95th pct of `sqrt(n)|T_a|`, `sqrt(n)|T_eps|`
* mean/95th/max of `sqrt(n)|T_eps - T_a|` — **the quantity that must vanish**
* `kappa_n`, `rho_n`, `log|T|`, `r_n`; lemma bounds at realised `D` and at
  worst case `D ~ rho_n`; ratio empirical / bound
* ATT bias, SD, mean SE, **empirical coverage** of the nominal-95% Wald interval

### Theoretical bound

```
kappa_n = v_{Lbar,p} * log(2 * e * n / v_{Lbar,p}) + log(4 / delta_n),  delta_n = 1/n
rho_n   = sqrt(kappa_n / n)
```
`v_{Lbar,p}` from the corrected `V_sub = O(Lbar log(Lbar p))`, calibrated by
least squares to the audit's anchors `v_{10,5} ~ 185`, `v_{20,5} ~ 427`,
`v_{30,5} ~ 686`, giving `v(L) = 3.122 * L * log2(5L) + 9.93` (reproduces the
anchors to within 1% and supplies `v_{5,5} ~ 82`). The Bartlett–Jordan–McAuliffe
Lemma 9 form `4 * d_tilde * L * log2(3L)` (`vc-dimension-hybrid-version.tex:138`)
is reported as a larger reference point. Constants `A = 2`, `delta_n = 1/n` are
choices, documented as such.

Lemma bounds tested (`refining-grid-standalone.tex`, Step 2):
```
|T_eps| = O_p{ (D_w + rho_n) rho_n }
|T_a|   = O_p{ D_w D_mu + (D_w + D_mu) rho_n + rho_n^2 }
```

---

## Success criteria (decision rule, fixed in advance)

The value-recovery route is **practically viable** at a cell if, at that cell,
`mean sqrt(n)|T_eps - T_a| << 1` and the quantity does **not** grow with `n` at
fixed `Lbar`, and Wald coverage is within Monte-Carlo error of 95%
(`500` reps ⇒ ±1.9pp at nominal 95%).

The audit's pessimism is **confirmed as a practical problem** if
`sqrt(n)|T_eps - T_a|` is `O(1)` or larger, or grows in `n`, or coverage
degrades as `Lbar` rises.

---

## Files

```
code/dgp.R           DGP, grid, partitions, exact grid-cell bookkeeping
code/terms.R         fitting, exact projections, realised T_a / T_eps, ATT + Wald
code/bounds.R        kappa_n / rho_n, VC calibration, exact class cardinality
code/run_study.R     driver (resumable per cell, furrr parallel)
code/analyze.R       summary tables -> results/tables.md
code/probe_timing.R  feasibility probe (leaf-budget cost surface)
results/             replications.rds, cells/, bounds.rds, tables.md, run_log.txt
REPORT.md            the answer
```

Reproduce: `Rscript code/run_study.R && Rscript code/analyze.R`
(`RTVR_REPS`, `RTVR_WORKERS`, `RTVR_N`, `RTVR_L`, `RTVR_SEED` override defaults;
seed default `20260908`).
