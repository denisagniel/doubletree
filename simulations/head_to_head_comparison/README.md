# head_to_head_comparison (doubletree, S5)

**Spec:** [`quality_reports/specs/2026-09-08_head-to-head-comparison.md`](../../../quality_reports/specs/2026-09-08_head-to-head-comparison.md)
**Stem:** `head_to_head_comparison`
**Status:** implemented 2026-09-09; **pilot only** (nsim ≈ 30). The full sweep at
spec §6's nsim = 1000/2000 is a separate cluster/batch job — see
[`code/run_sweep.R`](code/run_sweep.R).

---

## 0. This is not a leaderboard — read this before quoting any number

The paper makes **no claim** that doubletree has lower MSE or is more efficient
than a correctly-specified GLM or a correctly-specified black-box AIPW/DML
estimator. Under standard regularity, a correctly-specified linear working
model, a correctly-specified black-box AIPW/DML estimator, and doubletree under
exact structural sparsity **all attain the same semiparametric efficiency bound**
`V` (`thm:regular`: doubletree attains the *full* nonparametric bound, not a
smaller submodel bound — it does not attain something *better* than that bound).
**No arm in this study can "win" on efficiency asymptotically.**

What the paper does claim, and what this study is scoped to check:

| Claim | Statement | Regime |
|---|---|---|
| **C1** | Under structural sparsity, doubletree pays no interpretability-for-efficiency trade-off relative to an **oracle-nuisance** benchmark, while a linear working model *without the right basis* is genuinely misspecified. | R1 |
| **C2** | As structural sparsity degrades, doubletree's plain Wald CI degrades, but the **anchor** interval (`thm:anchor`) retains nominal coverage. | R3 |
| **C3** | At matched `n`, the full-sample flagship (no splitting) is at least as precise as doubletree's own cross-fit fallback. **The one place a "win" is actually a paper claim** (the whole point of not splitting). | R2 |

Consequences that are fixed **in advance**, regardless of the numbers:

- The headline C1 statistic is the **RMSE ratio to `oracle_aipw`, with a Monte
  Carlo SE on the ratio** — never a raw "beats GLM" comparison. A ratio near 1
  is the finding; a ratio *below* 1 is Monte Carlo noise, and the SE says so.
- **DGP-B (R3) is additive on the probability scale at every ε** (spec §3,
  verified numerically). So a linear-probability or near-linear-logit
  main-effects model is approximately *correctly specified* there at every ε and
  **should not degrade as ε grows**. A favourable-to-GLM finding in R3 is *not* a
  doubletree failure. The honest claim R3 licenses is narrower: *doubletree's
  plain CI degrades as ε grows (at `Lbar = 2` it cannot represent a truth needing
  4 leaves); the anchor interval recovers nominal coverage; a correctly-specified
  linear model is simply not stressed by this particular violation.*
- No comparative sentence beyond "matched X" or "retained coverage where Y did
  not" is licensed by this design.
- Spec §8's decision 1 (placement in the paper) is **deferred**, and decision 2
  (rescoping `sec:discussion`'s "left to future work" sentence to *real-data*
  applications) is deferred to a separate session after results land.
  `manuscript.tex` and `theory.tex` are **untouched** by this study.

---

## 1. Layout

```
head_to_head_comparison/
├── README.md                  # this file
├── code/
│   ├── common.R               # harness: paths, arm registry, regime grid, seeds, one cell
│   ├── dgps.R                 # the ONE new DGP (DGP-A) + the bilinear-remainder guard
│   ├── verify_dgps.R          # MUST PASS before any replication (see §3)
│   ├── analyze.R              # metrics with MC SEs; the paired RMSE-ratio bootstrap
│   ├── run_pilot.R            # small-scale run + wall-time projection
│   └── run_sweep.R            # full-scale run at spec §6's nsim, resumable
├── results/                   # per-cell .rds checkpoints, pilot_*.rds, verify_dgps_*.rds
├── tables/                    # verify_*.csv, pilot_*.csv, sweep_*.csv
└── figures/
```

Run everything **from the doubletree package root**:

```sh
Rscript simulations/head_to_head_comparison/code/verify_dgps.R
Rscript simulations/head_to_head_comparison/code/run_pilot.R 30
Rscript simulations/head_to_head_comparison/code/run_sweep.R          # cluster/batch
```

---

## 2. Arms

| Arm id | Function | Role |
|---|---|---|
| `doubletree` | `doubletree::estimate_att()` | the paper's own claim: full sample, **no splitting** |
| `doubletree_crossfit` | `doubletree::estimate_att_crossfit()` (`max_depth = 4`) | C3's comparison partner; C2's anchor |
| `doubletree_crossfit_shallow` | same, `max_depth = 1` | R5's deliberately **underpowered** anchor |
| `doubletree_anchor` | *derived* from `doubletree` + `doubletree_crossfit` | the literal `thm:anchor` interval |
| `doubletree_anchor_shallow` | *derived* from `doubletree` + `..._shallow` | R5's anchor under a violated `ass:rate` |
| `oracle_aipw` | **new**, `simulations/methods/method_oracle_aipw.R` | C1's **efficiency denominator** |
| `glm_main` | `att_linear(interactions = FALSE)` | "what an analyst actually does" |
| `glm_int` | `att_linear(interactions = TRUE)` | "the diligent analyst" |
| `lasso_dml` | **new**, `simulations/methods/method_lasso_dml.R` | `sec:lasso-comparison`'s empirical instance |
| `forest` | `att_forest()` (`ranger`) | the black-box ML arm |

Every arm computes the **identical** score, `manuscript.tex` `eq:score`:

```
psi(O; theta, eta) = pi^{-1} [ A{Y - mu(X) - theta} - (1-A) w(X){Y - mu(X)} ],  w = e/(1-e)
```

verified term-for-term against `att_linear()`'s existing implementation for both
new arms. A score difference between arms would be a second, uncontrolled source
of difference in a study whose whole subject is the *nuisance estimator*.

**Derived arms** (`doubletree_anchor*`) are built from a flagship/anchor **pair**
on the same replication, not from data:
`theta_full ± (|theta_full − theta_anchor| + z·sigma_anchor)`. No bias-aware
critical value — this is the triangle inequality `thm:anchor` actually proves.
`sigma_anchor` is already a standard error (`att_se()` returns
`sqrt(mean(psi^2)/n)`), so it is **not** divided by `sqrt(n)` again; the harness
asserts, per replication, that the anchor interval covers whenever the anchor's
own Wald interval does — a deterministic implication of the construction, so a
violation is an implementation bug and stops the run.

**GBM, the saturated-logistic reference, and the demoted "continuous" DGP are
deliberately absent** (spec §7). RF already discharges "black-box ML"; GBM's
tuning surface would invite an unresolvable "you tuned my competitor badly"
objection with no corresponding gain in what the study can conclude.

### `glmnet` is an ad-hoc dependency, by precedent

`glmnet` is used through `requireNamespace()` inside
`simulations/methods/method_lasso_dml.R` and is **not** added to doubletree's
`DESCRIPTION`. That mirrors `ranger`, which `method_forest.R` uses the same way
and which is likewise absent from `DESCRIPTION` (checked 2026-09-09). These files
live under `simulations/`, are never installed with the package, and are on no
package code path, so a `DESCRIPTION` entry would impose a dependency on package
users for code they cannot reach. (Spec §8 decision 3 left the convention open;
this is the resolution, by following the existing one.)

---

## 3. The tuning charter (spec §2) — including the asymmetry it does not hide

1. **Shared data.** One seed stream per replication:
   `set.seed(rep_seed(regime, dgp, n, rep))`, then **one** draw used by every arm
   in that replication. Paired comparisons (the RMSE ratio to `oracle_aipw`,
   C3's contrast) are therefore on identical data.
2. **Shared folds.** One `fold_seed` per replication, passed to every
   cross-fitted arm. `att_linear()`, `att_forest()` and `att_lasso_dml()` use
   **byte-identical** stratified-fold code, so an identical seed gives them an
   identical partition. `estimate_att_crossfit()` builds folds by its own
   internal routine; it receives the same seed, so it is reproducible, but its
   partition is **not guaranteed to coincide** with the other three. Small, real,
   and stated rather than claimed away.
3. **Shared clip.** Every arm clips `e_hat` to `[0.01, 0.99]`. That is
   `att_linear()`'s existing convention **and** doubletree's own
   `ass:construct(c)` constant, asserted against the package at load — so for
   once there is no discrepancy to state. `att_oracle()` carries the same clip
   but **asserts it is inert**: clipping a *true* propensity would bias the
   efficiency denominator every C1 ratio is built on.
4. **CV-tuning is asymmetric, and that is stated, not hidden.** Spec §2 asks for
   "every method CV-tuned or none, per DGP" and then concedes the asymmetry is
   unavoidable. The arms partition into:
   - **CV-tuned:** `doubletree_crossfit*` (`cv_regularization = TRUE`, shipped
     default), `lasso_dml` (`cv.glmnet` over λ)
   - **untuned:** `glm_main`, `glm_int` (a GLM has no penalty), `forest`
     (`ranger` defaults), `oracle_aipw` (nothing is fit)
   - **neither:** `doubletree`, whose `leaf_budget` is **fixed per DGP in
     advance** at the value making that DGP exactly sufficient — not chosen by
     any data-driven procedure — matching the paper's own framing that the
     analyst sets the budget.

   No arm's tuning was chosen after seeing results. `lambda_rule = "min"`
   (`lambda.min`) is the lasso default because a nuisance estimate feeding an
   orthogonal score wants prediction-optimality; `lambda.1se` deliberately
   under-fits, which for a *nuisance* is a bias source, not parsimony.
5. **`leaf_budget` is a DGP property, not a tuned parameter** — A: 4, A2: 4,
   B: 2, C: 4 — and is varied **only** in R6, where varying it *is* the
   experiment.
6. **Thread pinning.** `worker_limit = 1`, `parallel_cv = FALSE`, `ranger`
   `num.threads = 1` (asserted in `common.R`), and BLAS/OpenMP pinned to 1 in
   `run_sweep.R`. Spec §5 requires this explicitly or the per-method wall-clock
   comparison is meaningless.

---

## 4. DGPs

| id | Source | Status | Serves |
|---|---|---|---|
| `A` | **new**, `code/dgps.R` | new | R1 (C1) |
| `A2` | `partition_recovery_clt/code/dgps.R::dgp_spec_grid_sparse()` | **reused unmodified** | R2 (C3), R5, R6 |
| `A2_L4/L8/L16` | same, `leaf_budget` varied | reused unmodified | R6 |
| `B_eps0/05/10/15` | `honest_inference_sparsity_failure/code/dgps_delta_dial.R::dgp_spec_delta_dial()` | **reused unmodified** | R3 (C2) |
| `C` | `partition_recovery_clt/code/dgps.R::dgp_spec_weak_overlap()` | **reused unmodified** | R4 |

The three reused files are **read-only** for this study. `common.R` sources them
through `source_reused()`, which asserts the functions it needs exist, and
`verify_dgps.R` re-asserts their spec-stated population constants (F1's
`theta_0 = 0.08`; ST3's `max e_0 = plogis(-2.5+4+2.5) = 0.982014`) so a drift in
a file this study only reads surfaces as an error rather than as quietly
different results.

### Why DGP-A had to be written at all — and the guard that proves it

Spec §3's **Correction** records that the DGP originally planned for C1 (F1
reused verbatim) has a population **bilinear remainder of exactly zero**
(`1.7e-18`) under main-effects-only misspecification, as does the rejected
"complex" DGP from `single_tree_coverage` (`-5e-20`). Both would have reported
"GLM main-effects: unbiased" for C1 — a **null result caused entirely by the
DGP** — discoverable only after a full implementation and a full run.

The mechanism is general: a population OLS main-effects fit including the shared
confounder `X1` has a residual with `E[residual | X1] = 0` exactly. If one
nuisance's "other" interacting covariate is conditionally independent of the
other nuisance's, given `X1`, the cross term integrates to exactly zero
**regardless of effect size**. Turning up `e_gap`/`mu_gap` cannot fix it.

DGP-A defeats that by construction: both nuisances depend on the **same** pair
`(X1, X2)` with the **gate swapped** —

```
e_0(X1,X2)  = e_lo                    X1 = 0            (X2 inert here)
              e_mid                   X1 = 1, X2 = 0
              e_mid + e_gap           X1 = 1, X2 = 1

mu_0(X1,X2) = mu_lo                   X2 = 0            (X1 inert here)
              mu_mid                  X1 = 0, X2 = 1
              mu_mid + mu_gap         X1 = 1, X2 = 1

mu_1 = mu_0 + tau ;  X3,X4,X5 pure noise ;  p = 5
e_lo=0.20 e_mid=0.45 e_gap=0.30 ; mu_lo=0.10 mu_mid=0.50 mu_gap=0.40 ; tau=0.08
```

Constants are F1's own, so the two DGPs differ in **dependence structure and
nothing else** — "GLM is biased here and not there" must not be readable as an
effect-size difference.

**`verify_main_effects_remainder()` in `code/dgps.R` is a hard gate**, re-run at
the top of both `run_pilot.R` and `run_sweep.R`. It recomputes the remainder from
this study's own formulas and **errors** if it is near zero (with the
cancellation diagnosis in the message) or if it has drifted from spec §3's
verified value. Verified 2026-09-09:

| quantity | computed | spec §3 |
|---|---|---|
| bilinear remainder (P-weighted) | **0.0144457** | 0.01445 |
| `pi = P(A=1)` | 0.40000 | 0.4 |
| implied asymptotic bias | 0.036114 | 0.0361 |
| relative bias vs `theta_0 = tau = 0.08` | **45.14 %** | ~45 % |
| `E[e_0 | X1=0] / E[e_0 | X1=1]` | 0.20 / 0.60 | 0.20 / 0.60 |
| `E[mu_0 | X1=0] / E[mu_0 | X1=1]` | 0.30 / 0.50 | 0.30 / 0.50 |
| empirical `cor(A, Y0)` (n = 2e5) | 0.2077 | — (confounding confirmed) |
| distinct values of `e_0` / `mu_0` | 3 / 3 | 3 / 3 (exact sparsity at `Lbar = 4`) |

The control-weighted (`nu`) variant of the `mu` projection gives 0.01739
(relative bias 54.3 %) and is reported as a **sensitivity**; the assertion is on
the P-weighted value spec §3 states.

**Two caveats on the 0.01445, stated so a small discrepancy is not misread as a
bug:** it is a *linear-probability* population projection, whereas `att_linear()`
fits a *logit* main-effects model whose population limit differs slightly; and it
is a **DGP diagnostic** establishing the misspecification is non-degenerate, not
a prediction of any arm's realised bias.

---

## 5. Regimes

| Regime | DGP(s) | Arms | Evidence for |
|---|---|---|---|
| **R1** | `A` | `doubletree`, `oracle_aipw`, `glm_main`, `glm_int`, `lasso_dml` | **C1** |
| **R2** | `A2` | `doubletree`, `doubletree_crossfit` | **C3** |
| **R3** | `B_eps0/05/10/15` | `doubletree`, `doubletree_crossfit`, `doubletree_anchor`, `glm_main`, `oracle_aipw` | **C2** |
| **R4** | `C` | `doubletree`, `glm_main`, `forest` | stress (a): overlap/boundary |
| **R5** | `A2`, `B_eps15` | `doubletree`, `..._crossfit_shallow`, `..._anchor_shallow`, `..._crossfit`, `..._anchor` | stress (b): anchor `ass:rate` violation |
| **R6** | `A2_L4/L8/L16` | `doubletree` | stress (c): `(2Lbar+1)/n` coverage boundary |

R4–R6 are the **required stress regimes** (constitution invariant: at least one
regime where the recommended method should struggle). Following Oracle's
redirection in spec §4, they target genuine failure modes in doubletree's
*inference* machinery — boundary weights, an anchor that does not satisfy its own
assumption, and the paper's own predicted finite-sample bias direction — rather
than a generic "hard DGP" aimed at the point estimate, which the anchor interval
is specifically designed to absorb.

**Deliberate additions to spec §4, flagged:**

- **R5 runs the shipped-default anchor alongside the underpowered one, on two
  DGPs.** Spec §4 asks only for the underpowered anchor and offers "A2 or B".
  Without a same-data reference anchor, an anchor-coverage shortfall cannot be
  attributed to the anchor's depth rather than to the DGP or the harness, so the
  default is included as a control. And the 2026-09-09 pilot showed **A2 alone
  cannot exhibit the failure mode at all** — see §10 R5.
- **R4 has no oracle arm.** Spec §4 lists flagship/GLM/RF. ST3's `e_0` reaches
  0.982, so an oracle RMSE ratio would be a ratio to a near-degenerate-weight
  benchmark — not the C1 statistic, and an invitation to over-read.
- **`ratio_ref` for R2 is `doubletree_crossfit`, not `oracle_aipw`.** C3 is the
  one comparison where a "win" is a paper claim, and its natural denominator is
  the fallback, not the oracle. (F1's remainder is exactly zero, so no GLM arm
  runs on it and no oracle-ratio question arises.)

### R4 / DGP-C calibration: measured, not pre-emptively escalated

Spec §8 decision 4 left ST3's exact recalibration open. Per the user's decision
(2026-09-09) it is **used as shipped**, and `verify_dgps.R` section 5 *measures*
whether it produces a **truly ill-defined** propensity weight (a leaf with zero
controls, so `m_0` is undefined and `e_hat` is driven to the clip) rather than
merely a large one. Measured 2026-09-09, 400 draws per `n`:

| `n` | E[cell size] | E[controls in `(X1,X2)=(1,1)`] | mean | min | **P(zero controls)** | P(≤2 controls) |
|---|---|---|---|---|---|---|
| 500 | 125 | 2.25 | 2.32 | **0** | **7.5 %** | 59.0 % |
| 1000 | 250 | 4.50 | 4.56 | **0** | 1.75 % | 18.25 % |
| 2000 | 500 | 8.99 | 8.96 | 2 | 0 % | 1.25 % |

**Verdict: adequate as-is, no escalation.** At `n = 500` the extreme cell is
control-*empty* in 7.5 % of draws and holds ≤2 controls in 59 % — genuinely
ill-defined, not merely extreme. If the full-scale R4 run shows no inference
degradation, the next lever is a *smaller* `n` or steeper propensity
coefficients, not a change to the shipped ST3 generator (which other studies
depend on).

---

## 6. Metrics (spec §5) — every number carries a Monte Carlo SE

- bias, relative bias, RMSE, and **RMSE ratio to the regime's reference arm**,
  computed **paired** (same replications) with a **paired-bootstrap MC SE** that
  resamples replications, preserving the strong positive correlation between two
  arms' errors on a shared draw. A delta-method or independent-samples SE would
  badly overstate the ratio's uncertainty.
- empirical coverage of the nominal 95 % CI, with MC SE. `analyze.R` reports the
  **empirical-proportion** SE `sqrt(p(1−p)/R)` *and* spec §6's sizing SE
  `sqrt(0.95·0.05/R)`; the former is the right thing for *reporting* a realised
  coverage (at `p = 0.62` the latter understates it ~2.2×), the latter is the
  right thing for *sizing* and stays visible next to it.
- CI width (mean and its own MC SE).
- `n_leaves_e`/`n_leaves_m0` (interpretability proxy, doubletree only) and
  `certified_frac`, plus **coverage conditional on `certified_*`** — per Oracle's
  Rashomon-multiplicity point: if several near-optimal partitions differ, the
  "auditable" claim is weaker than a single certified fit implies.
- clip fraction, `max_e_hat`, `min_e_hat` for **every** arm in **every** regime,
  so a coverage anomaly can be attributed to the clip biting rather than guessed
  at.
- wall-clock seconds per replication per arm (threads pinned — see §3.6).
- **`n_treated` in every row**: effective precision for the ATT scales with
  `n_treated`, not `n` (spec §5).
- `theta_0` is **exact and closed-form** for every DGP — none requires a
  Monte-Carlo "truth".

## 7. Replications

Spec §6: nsim ≈ **1000** per (regime, n) cell for the coverage-focused regimes
(R1, R3, R4, R5) and ≈ **2000** for R2 and R6, where the comparison is an RMSE
*ratio* between two methods on the same data. Encoded in `NSIM_FULL` and used by
`run_sweep.R`; `run_pilot.R` runs at ~30 and says so, and at 30 reps a coverage
MC SE is ≈ 0.04, so **a 0.95-vs-0.90 distinction is not resolvable from the
pilot** — directional sanity only.

## 8. Results format

Long: one row per `(regime, dgp, n, rep, arm)`. The arm *set* varies by regime
(5 arms in R1, 1 in R6), so a wide schema would be mostly `NA` and every analysis
would have to know which columns are meaningful where; long format makes "which
arms ran in this cell" a fact in the data. Arm-specific diagnostics that only one
or two arms produce live in a single `extras` string; the diagnostics every
analysis needs (leaves, certification, clip fraction, `n_treated`) are real
columns. Failure rows carry the condition **message** and the same columns as
success rows — never a shorter row that silently drops out of an `rbind`, and
never a plausible-looking `NA` estimate with no explanation.

## 9. Constitution compliance

- [x] **Stress regimes** — R4/R5/R6, targeted at doubletree's actual inference
  machinery.
- [x] **Design serves understanding, not impressive results** — §0 above is
  written *before* any replication runs and fixes in advance what conclusions are
  licensed regardless of the numbers.
- [x] **UQ propagates** — every regime reports CI coverage and width, not point
  estimates alone; every metric carries its own MC SE.
- [x] **Reproducibility** — one shared seed stream per replication across all
  arms; the whole sweep is a function of the single integer `SEED_MASTER`
  (`20260908`); `nsim` fixed in advance from a stated precision target.
- [x] **No quiet favoritism** — §0 and §5 state, before the run, where a
  competitor is *expected* to do as well as or better than doubletree (R3's
  GLM arm) and commit to reporting it that way.
- [x] **No silent fallbacks** — no `tryCatch`-to-`NULL`, no `%||%`, no
  `possibly()`/`safely()` in any estimator. `att_lasso_dml()` deliberately has
  **no** main-effects fallback (unlike `att_linear()`), because a fallback would
  report main-effects numbers under the `lasso_dml` label.

---

## 10. Pilot results (2026-09-09, nsim = 30, `results/pilot_*.rds`)

**Pilot scale. Coverage MC SE ≈ 0.04–0.09 throughout, so a 0.95-vs-0.90
distinction is _not_ resolvable here.** These numbers establish that every regime
executes end-to-end and that the directional predictions hold; they are not the
study's findings. Wall time: 13.3 min for 18 cells × 30 reps (+6.6 min for R5's
second DGP). Projected full sweep at spec §6's nsim over the full `n` grid,
sequential and single-threaded: **≈ 53 h**, dominated by R2 (22.8 h) and R5
(18.1 h) — both `estimate_att_crossfit()`-heavy, whose nested CV costs ~2–5 s/rep
against the flagship's ~0.02 s. Parallelise over **cells**.

### R1 — C1: no interpretability-for-efficiency cost. **Confirmed.**

RMSE ratio to `oracle_aipw`, paired, bootstrap MC SE:

| `n` | `doubletree` | `glm_main` | `glm_int` | `lasso_dml` |
|---|---|---|---|---|
| 500 | **1.030 ± 0.036** | 1.290 ± 0.114 | 0.989 ± 0.122 | 1.143 ± 0.070 |
| 1000 | **0.998 ± 0.029** | 1.364 ± 0.170 | 1.086 ± 0.052 | 1.068 ± 0.046 |
| 2000 | **1.039 ± 0.024** | 1.876 ± 0.229 | 1.079 ± 0.049 | 1.061 ± 0.025 |

Relative bias and coverage:

| `n` | `doubletree` | `oracle_aipw` | `glm_main` | `glm_int` | `lasso_dml` |
|---|---|---|---|---|---|
| 500 | 1.3 %, 1.000 | 1.9 %, 1.000 | **27.2 %, 0.967** | 4.8 %, 1.000 | 10.1 %, 1.000 |
| 1000 | −0.2 %, 0.967 | −0.8 %, 0.967 | **27.4 %, 0.833** | −0.7 %, 0.967 | 2.8 %, 0.967 |
| 2000 | 3.0 %, 0.967 | 2.9 %, 0.967 | **33.4 %, 0.767** | 4.2 %, 0.933 | 3.9 %, 0.933 |

- doubletree tracks the oracle to within MC error at every `n`, with a ratio CI
  containing 1 throughout — **the C1 finding**, in the form the theory licenses
  ("matched the oracle"), not as a win over anyone.
- `glm_main` is genuinely misspecified: relative bias 27–33 %, RMSE ratio rising
  1.29 → 1.88 in `n` (bias is fixed while the oracle's RMSE shrinks, so the ratio
  must diverge — the signature of a bias, not of noise), and coverage collapsing
  0.967 → 0.833 → 0.767 as the interval shrinks around the wrong value.
- Observed relative bias 27–33 % vs spec §3's predicted **45 %** asymptotic value.
  Expected and flagged in advance (§4): the 45 % follows from a
  *linear-probability* population projection, while `att_linear()` fits a *logit*
  main-effects model whose population limit differs; and the realised bias is
  still rising in `n` (0.0218 → 0.0220 → 0.0267 against the asymptotic 0.0361).
  Not a discrepancy to chase; the full sweep's larger `n` will resolve it.
- Adding the right basis fixes it: `glm_int` and `lasso_dml` both return to ~1.0.
  This is what separates "the GLM lacked the right basis" from "any
  linear-in-parameters method fails here" — precisely what `lasso_dml` was added
  to establish (spec §8 decision 3).

### R2 — C3: full sample vs its own cross-fit fallback. **Consistent, unresolved at this nsim.**

| `n` | RMSE ratio flagship/crossfit | CI width flagship | CI width crossfit | s/rep flagship | s/rep crossfit |
|---|---|---|---|---|---|
| 500 | 0.986 ± 0.063 | 0.1803 | 0.2016 | 0.019 | 3.08 |
| 2000 | 1.015 ± 0.023 | 0.0904 | 0.0917 | 0.031 | 5.21 |

"At least as precise" is consistent with both cells (ratio CIs contain 1) and the
flagship's mean CI width is narrower at both `n`, but at nsim = 30 the ratio's MC
SE (0.02–0.06) cannot resolve the small effect C3 predicts. This is exactly the
regime spec §6 assigns nsim = 2000. Both arms are unbiased and at/above nominal
coverage, as exact sparsity requires. The flagship is **~160× faster**.

### R3 — C2: anchor retains coverage as sparsity degrades. **Confirmed, and the honest GLM framing held.**

`n = 500`. Coverage (MC SE 0.046–0.086):

| ε | predicted bias | `doubletree` plain Wald | `doubletree_anchor` | `doubletree_crossfit` | `glm_main` | `oracle_aipw` |
|---|---|---|---|---|---|---|
| 0.00 | 0.0000 | 0.933 | 0.933 | 0.933 | 0.933 | 0.933 |
| 0.05 | 0.0123 | 0.933 | 1.000 | 0.967 | 0.967 | 0.933 |
| 0.10 | 0.0491 | 0.900 | 0.933 | 0.900 | 0.900 | 0.933 |
| 0.15 | 0.1105 | **0.667** | **1.000** | 0.967 | 0.967 | 1.000 |

- **C2 as stated:** the plain Wald interval degrades monotonically
  (0.933 → 0.667) while the anchor interval retains nominal coverage, paying for
  it in width (0.208 → 0.223 → 0.313 → **0.521**) — `cor:width`'s growth-then-plateau
  behaviour in the right direction.
- doubletree's realised bias at ε = 0.15 is 0.0913 against `dgps_delta_dial.R`'s
  closed-form pseudo-true prediction of **0.1105** — same order, right direction,
  and the gap is what one expects when not every replication's fit lands exactly
  at the pseudo-true partition at `n = 500`.
- **`glm_main` does not degrade** (0.933 / 0.967 / 0.900 / 0.967) and its RMSE
  ratio to the oracle stays near 1 (1.01 → 1.02 → 1.14 → 1.19) while doubletree's
  rises to 2.48. **This is the predicted, pre-committed outcome, not a
  doubletree failure**: DGP-B is additive on the probability scale at every ε
  (spec §3, verified), so a main-effects model is approximately correctly
  specified there. Recorded here in the terms §0 fixed before the run.

### R4 — stress (a), overlap/boundary. **Ran; the arm that breaks is the forest, not doubletree.**

| `n` | arm | rel. bias | RMSE | coverage | CI width |
|---|---|---|---|---|---|
| 500 | `doubletree` | 44.3 % | 0.1421 | 0.900 | 0.505 |
| 500 | `glm_main` | 31.1 % | 0.1823 | 0.900 | 0.593 |
| 500 | `forest` | **191.3 %** | **0.2794** | **0.767** | 0.497 |
| 1000 | `doubletree` | 10.8 % | 0.1073 | 0.900 | 0.406 |
| 1000 | `glm_main` | 13.9 % | 0.1196 | 0.900 | 0.385 |
| 1000 | `forest` | **140.9 %** | **0.2448** | **0.667** | 0.299 |

All three arms sit at or below nominal (`θ₀ = 0.10`, RMSE ≥ 0.10 — this is a hard
regime for everyone), but the **anticipated doubletree-specific inference failure
did not appear at pilot scale**: doubletree is the *least* biased and the
*best*-covering of the three. The forest is severely biased (141–191 %) with
coverage 0.667–0.767, because `ranger`'s smoothing cannot resolve the 0.982
propensity cell and its over-narrow CI then excludes `θ₀`. Reported as measured;
this regime was included to look for a doubletree failure and instead found a
competitor's, which §0 commits to reporting either way. Whether doubletree's
0.900 is a genuine shortfall needs nsim = 1000 (MC SE here is 0.055).

### R5 — stress (b), the anchor's own `ass:rate`. **Confirmed on `B_eps15`; `A2` cannot exhibit it.**

`n = 500`, coverage and CI width:

| DGP | `doubletree` | `..._anchor` (depth 4) | `..._anchor_shallow` (depth 1) | `..._crossfit` bias | `..._crossfit_shallow` bias |
|---|---|---|---|---|---|
| `A2` (flagship unbiased) | 1.000 / 0.183 | 1.000 / 0.223 | 1.000 / **0.231** | +0.0133 | +0.0149 |
| `B_eps15` (flagship biased) | 0.767 / 0.202 | **1.000** / 0.467 | **0.833** / **0.227** | −0.0185 | **+0.0753** |

**The pilot changed this regime's design.** On `A2` the flagship is unbiased
(exact sparsity), so degrading the anchor only makes
`|δ̂| = |θ̂_full − θ̂_anchor|` *larger*, which **widens** the anchor interval
(0.231 > 0.223) and *preserves* coverage. `thm:anchor`'s own conservatism absorbs
the `ass:rate` violation, so `A2` alone cannot show the failure mode this regime
exists to find. `B_eps15` was therefore added (spec §4 offers "A2 or B"), and it
shows the mechanism cleanly: the depth-1 anchor is **as biased as the flagship**
(+0.0753 vs +0.0762), so `δ̂ ≈ 0`, the anchor interval **collapses to a plain
Wald interval on a biased estimate** (width 0.227 ≈ the plain 0.202, against
0.467 for the properly-powered anchor), and coverage falls to **0.833**. This is
`thm:anchor` failing exactly where its own assumption fails, checked directly.
Both DGPs are kept: `A2` is the informative negative control.

### R6 — stress (c), the `(2L̄+1)/n` coverage boundary. **Not resolvable at nsim = 30.**

| `L̄` | `(2L̄+1)/n` at n=200 | coverage n=200 | coverage n=500 | CI width n=200 |
|---|---|---|---|---|
| 4 | 0.045 | 0.900 | 0.967 | 0.2755 |
| 8 | 0.085 | 0.933 | 0.967 | 0.2972 |
| 16 | 0.165 | 0.967 | 0.900 | 0.2794 |

No monotone pattern, which is expected: the MC SE is 0.033–0.055 and the effect
sought is a few percentage points. This is the regime spec §6 assigns nsim = 2000,
and it is cheap (0.23 projected hours) — run it at full scale before drawing any
conclusion. Bias is small and of both signs at every cell, so nothing suggests a
gross problem.

### Open items after the pilot

1. **`glm_main`'s realised bias (27–33 %) vs the predicted 45 %.** Expected
   (logit vs linear-probability projection, plus finite `n`); the full sweep's
   `n = 4000` cell will show whether it continues toward 0.0361. No action.
2. **R4's doubletree coverage of 0.900** — needs nsim = 1000 to distinguish from
   nominal. The DGP-C calibration is adequate (§5), so if 0.900 persists it is a
   real finding, and if it does not, the regime found a forest failure instead.
3. **R2 and R6 are the two regimes whose conclusions genuinely require their
   nsim = 2000 targets.** Neither is resolvable from this pilot, by design.
4. **Full sweep is ~53 h sequential**, 77 % of it in R2 and R5's
   `estimate_att_crossfit()` calls. `run_sweep.R` is resumable per cell; parallelise
   over cells, never within a fit.
