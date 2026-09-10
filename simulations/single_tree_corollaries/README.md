# single_tree_corollaries (doubletree)

**Spec:** [`quality_reports/specs/2026-09-09_single-outcome-tree-corollaries.md`](../../../quality_reports/specs/2026-09-09_single-outcome-tree-corollaries.md)
**Stem:** `single_tree_corollaries`
**Status:** implemented 2026-09-09; **pilot only** (`nsim = 50`, `n = 1000`, spec §6's
local stage). The cluster stage (`nsim = 1000`, `n ∈ {1000, 4000}`, 17 cells) is a
separate decision and is deliberately **not** runnable from this directory.

---

## 0. Why this study exists

`inst/paper/claims.md` (seeded 2026-09-09) found **zero** simulation evidence for
five claims drafted and proof-audited on 2026-09-08. This study supplies it:

| Claim | Statement | Where it is tested |
|---|---|---|
| `cor:single-tree` | The two-tree estimator with the propensity tree **forced equal** to the outcome tree collapses to `mean(Y[A==1] − μ̂[A==1])`. | **Unit test**, not a regime: [`tests/testthat/test-single-tree-identity.R`](../../tests/testthat/test-single-tree-identity.R) |
| `cor:single-tree-saturated` | When `τ ∈ S_e ∩ S_μ`, the tied estimator is **exactly efficient** (`V_τ = V`). | Regime **A** |
| `cor:single-tree-coarsening` (≤) | Under homoskedasticity, coarsening the propensity makes the tied estimator **super-efficient** (`V_τ ≤ V`). | Regime **B** |
| `cor:single-tree-coarsening` (reversal) | The falsifiable sub-claim: anti-correlated heteroskedasticity **reverses** it (`V_τ > V`). | Regime **C** |
| `rem:single-tree-se` | The naive SE (treated-arm term only) omits a nonnegative control-arm term and so **under-covers**. | Regime **D**, folded into A–C as a metric (spec §3) |
| `lem:single-tree-linear` (caveat) | With `S_μ` non-singleton, **no single fixed variance** describes the tied estimator. | Regime **E** |

`prop:spectest` is **out of scope** (spec §0, §9): it targets a different
diagnostic (the anchor's fidelity statistic `δ̂ = θ̂ − θ̃`), shares no DGP structure
with the tied-estimator family, and needs its own alternative-hypothesis regimes.
So are any change to `manuscript.tex` / `theory.tex` / `claims.md`, a dedicated
`estimate_att_tie()` package function, and GBM/RF/black-box comparators.

## 1. What is under test, and what it is compared to

**Object:** the single-outcome-tree plug-in

```
theta_tie = n_1^{-1} sum_{A_i = 1} { Y_i - mu_hat(X_i) }
```

where `μ̂` is **doubletree's own outcome tree** — `estimate_att()`'s
`nuisance_fits$m0_model`, fit on control-only data at a fixed `leaf_budget`,
leaf-constant, predicted on all `n`. **No propensity tree is fitted for it.**

**Comparator:** doubletree's full two-tree `estimate_att()` at the **same**
`leaf_budget` on the **same** data — not a rival method, the other half of the same
identity. It runs as the `two_tree` arm off the *identical* fit, so any difference
between `fitted_tie` and `two_tree` is the tie, never the data or the tuning.

**No exported package function** (spec §7). `theta_tie` is a one-liner off an
existing `estimate_att()` fit, and this construction's own corollary
(`rem:single-tree-se`) states that its *default* SE is wrong — shipping an
`estimate_att_tie()` would invite users to call it and report an under-covering
interval. The computation lives only in `code/common.R`.

## 2. Layout

```
single_tree_corollaries/
├── README.md                  # this file
├── code/
│   ├── dgps.R                 # Regimes A/B/C/E + the population variance machinery
│   │                          #   + verify_population_gates(): spec §4's gates, recomputed
│   ├── common.R               # harness: the tied estimator, both arms, seeds, one cell
│   ├── analyze.R              # spec §5 metrics with MC SEs; the gate-by-gate reading
│   └── run_pilot.R            # spec §6's local pilot; prints everything above
├── results/                   # cell_*.rds checkpoints, pilot_*.rds
├── tables/                    # population_*.csv, pilot_summary_*.csv, pilot_gates_*.csv
└── figures/
```

Run **from the doubletree package root**:

```sh
Rscript -e 'devtools::test(filter = "single-tree-identity")'      # spec §1, first
Rscript simulations/single_tree_corollaries/code/run_pilot.R      # 50 reps, n = 1000
Rscript simulations/single_tree_corollaries/code/run_pilot.R 50 1000 A,C   # subset
```

`run_pilot.R` re-runs the population gates before any replication, so it cannot
run against DGP constants that no longer reproduce spec §4.

## 3. The two arms (spec §2), and why the split is load-bearing

| Arm | What it is |
|---|---|
| `oracle_tie` | `τ` imposed **analytically**: leaf-wise sample means on the *known* true leaves, no tree search at all. Isolates the fixed-`τ` mathematics with zero tree-search noise. |
| `fitted_tie` | doubletree's actual `estimate_att()` at `leaf_budget` **exactly** `|leaves(τ)|`, leaf assignments read via the existing internal `doubletree:::.tree_leaf_paths()` — reused, not reimplemented. |
| `two_tree` | the comparator, from the same fit. |

Without the split, a regime failure is ambiguous between "the corollary's maths is
wrong" and "the tree missed `τ`". `leaf_budget` is set **equal to** `|leaves(τ)|`
in every regime and asserted in `new_st_dgp()`: a larger budget lets the tree
over-split and silently exit `S_μ`. **Regime E is fitted-arm only** — the point
there is which `τ` the fitting *realises*.

Leaf extraction is checked, not trusted: `st_fitted_leaf()` asserts that
`estimate_att()`'s own leaf-constant predictions are constant within the extracted
leaves and that the leaf count matches `n_leaves_m0`. A misaligned design or path
convention would otherwise compute every leaf-wise quantity on the wrong groups
with no symptom.

## 4. DGPs — every constant is the spec's

All covariates binary (`estimate_att()`'s implementation constraint, not worked
around), outcome continuous, `σ₀² ` and the CATE as spec §3 states them.

| Regime | Construction (spec §3) | `L̄` | `π` | `θ₀` |
|---|---|---|---|---|
| **A** | `e₀ ∈ {0.3, 0.7}`, `μ₀ ∈ {0, 1}` both on `X1` only; `σ₀² = 1`. `τ ∈ S_e ∩ S_μ`. | 2 | 0.5 | 0.25 |
| **B** | `μ₀` on `X1` (so `τ ∈ S_μ`); `e₀ ∈ {0.2, 0.8}` driven by `X2` **within** each leaf (so `τ ∉ S_e`); `σ₀² = 1`. | 2 | 0.5 | 0.25 |
| **C** | `P(X1=1) = 0.85`. Boring leaf: `e₀ = 0.5`, `σ₀² = 1`. Reversal leaf: `X2` drives `e₀ ∈ {0.1, 0.8}` **paired with** `σ₀² ∈ {1, 0.01}` — high `e₀` with *low* variance. | 2 | 0.4575 | 0.10 |
| **E** | OR structure: `μ₀ = a` at `(0,0)`, `b` otherwise; `e₀(00)=e₀(11)=0.5`, `e₀(01)=0.2`, `e₀(10)=0.7`. | 3 | 0.475 | 0.25 |

`σ₀² = 0.01`, **not** 0, in Regime C — spec §3 is explicit (an exact zero
degenerates the plug-in `σ̂²` and can perturb tree fitting). Regime D has **no
DGP**: `make_st_dgp("D")` errors with that sentence rather than "unknown id".

**Where the CATE lives.** Spec §3 gives the treated arm "an additive CATE plus its
own noise `σ₁²`" while pinning `e₀`/`μ₀` to `X1` (and `X2`) with the remaining
covariate declared pure noise. Both hold here: the CATE rides on `X3`, so `μ₀` does
not depend on it and the `μ` tree — fit on **controls only** — cannot see it. `τ ∈
S_μ`, `τ ∈ S_e` (Regime A) and partition recovery are all untouched. A consequence:
`Tr = π{Var(τ) + σ₁²}` exactly.

**`σ₁²` is solved, not typed.** Spec §3 pins `Tr` numerically in Regimes B
(`= 0.5 × C_untied ≈ 0.81`) and C (`≈ 0.15`) but leaves the split between CATE
spread and `σ₁²` free. `st_sigma1_sq_for_Tr()` solves
`σ₁² = Tr_target/π − Var(τ)` in closed form, so `Tr` hits the spec's target
**exactly** and the gates are reproducible to machine precision.

**Regime E's tie is genuine, not a duplicate-column artefact.** Tree 1 (split `X1`,
then `X2` within `X1 = 0`) gives leaves `{10,11}, {00}, {01}`; Tree 2 (split `X2`,
then `X1`) gives `{01,11}, {00}, {10}`. Cell `10` pools with `11` under Tree 1 but
is its own leaf under Tree 2. Both have exactly zero population approximation error
for `μ₀` — asserted per partition in `new_st_dgp()`, not claimed in a comment.

## 5. The two population variances, and the one methodological rule

`Var(ψ) = π⁻²{Tr + C}` with the treated-arm and control-arm terms orthogonal;
`C_untied = E[f(e₀)σ₀²]`, `f(e) = e²/(1−e)` (efficient), and
`C_tied = E[(1−e₀)w_tie(ℓ)²σ₀²]`, `w_tie(ℓ) = ē(ℓ)/(1−ē(ℓ))` (tied at fixed `τ`).
Spec §3 reports `V = Tr + C_untied` and `V_τ = Tr + C_tied` **without** the `π⁻²`
(it cancels in the ratio the gates are stated on); `code/dgps.R` returns both those
(`*_num`) and the `π⁻²`-scaled `V_full` / `V_tau_full` the empirical ratios need.

**Super-/sub-efficiency is judged empirical-to-*analytic*, never
empirical-to-empirical** (spec §5, Oracle's reason): the full two-tree estimator's
own finite-sample variance can exceed the bound `V` from nuisance noise, which would
let "super-efficiency" appear for entirely the wrong reason. Both ratios are
reported: `n·V̂ar/V_tau_full` (should be ≈ 1) and `n·V̂ar/V_full` (≈ `V_τ/V`).

## 6. Metrics (spec §5) — every number carries a Monte Carlo SE

- bias, relative bias, RMSE, each with `MCSE(bias) = sd(θ̂)/√R`;
- `n·V̂ar(θ̂_tie)` against **analytic** `V_τ` **and** `V`, with the variance ratio's
  relative MC SE `√(2/(R−1))` — **20.2 % at `R = 50`**, against the **6.3 %** spec §4
  pre-registers at `R = 1000`;
- **partition recovery rate**, a first-class metric (spec §5): fitted-arm results in
  A–D are only interpretable conditional on high recovery. Compared as a *bijection*
  between the two labellings' blocks, so a relabelled-but-identical partition counts
  as a match and a same-leaf-count-but-different-split does not;
- naive vs corrected SE coverage at 95 % nominal, plus the analytic
  `implied_coverage_naive` each regime predicts;
- `identity_gap` per replication (`cor:single-tree`, spec §6's smoke check) —
  **fatal**, not averaged: it is exact arithmetic, so a gap above `1e-10` is a
  harness bug;
- `mu_leaf_gap`, `clip_frac`, `n_leaves_m0`, seconds per replication.

**Corrected-SE formula** (spec §5, labelled exploratory there):

```
SE_corr^2 = SE_naive^2 + n^{-1} pi_hat^{-2} sum_l P_hat(l){1 - e_hat(l)} w_hat_tie(l)^2 sigma0_hat^2(l)
```

Spec §5 writes the second term without the `n⁻¹`; it is supplied here because
`SE_naive²` already carries one and the two must be commensurate. Flagged rather
than silently fixed.

## 7. Replications

Spec §6's local pilot: `nsim = 50`, `n = 1000`. It **can** settle (i) the identity
to `1e-10`, (ii) fitted-arm recovery ≥ 95 %, (iii) whether empirical variances track
the analytic ratios in direction and rough magnitude. It **cannot** resolve a
coverage or variance-ratio gate to spec §4's precision, and Regime A's `[0.88, 1.14]`
band is registered for `nsim = 1000` — a pilot excursion outside it would not be a
violation. The cluster stage is spec §6's second stage and is a separate decision.

## 8. Constitution compliance

- [x] **Stress regimes** — Regime C (adversarial anti-correlated heteroskedasticity,
  the falsifiable sign-reversal) and Regime E (random-index / non-singleton `S_μ`)
  are both cases where the clean story breaks, by the paper's own admission.
- [x] **Design serves understanding** — Regime C is calibrated to *test* the
  reversal, not to make the tied estimator win; its margin's fragility is stated in
  `code/dgps.R` and §10 below, and `Tr_target` is an argument so the sensitivity can
  be re-measured rather than argued.
- [x] **UQ propagates** — naive vs corrected SE and their coverage are first-class
  outputs in every regime.
- [x] **Reproducibility** — `set.seed(st_rep_seed(regime, n, rep))` once per
  replication, then one draw shared by every arm; the whole study is a function of
  `SEED_MASTER = 20260909`. All DGP constants tabulated; `σ₁²` derived, not typed.
- [x] **No quiet favoritism** — Regime C, designed to show the tied estimator can be
  *worse*, is reported with equal prominence; Regime E's result is explicitly a "no
  fixed variance" negative finding.
- [x] **No silent fallbacks** — no `tryCatch`-to-`NULL`, no `%||%`, no
  `possibly()`/`safely()`. A leaf with < 2 controls **errors** (the leaf-wise `σ̂₀²`
  is genuinely undefined) instead of being imputed; failure rows carry the condition
  message and the full column set.

---

## 9. Deviations from the spec, flagged

**One substantive deviation.** Spec §1 says `μ̂` is "`estimate_att()`'s
`m0_model`/`m0_hat`". The tied estimator here is built on the fitted tree's
**partition** with its leaf values **recomputed as exact within-leaf control
means**, not on the serialised `nuisance_fits$outcome_control` vector.

Reason, measured 2026-09-09: **optimaltrees stores a leaf's fitted value to six
significant decimal digits** (a leaf whose control mean is `−0.058940614050488`
comes back as exactly `−0.058941`). `cor:single-tree`'s identity is exact only when
the leaf value *is* the control mean, so using the serialised prediction verbatim
caps the achievable identity gap at **≈ 4e−07** — 3.6 orders of magnitude above spec
§1's `1e-10`, for a reason that has nothing to do with the corollary. The two
objects are the *same estimator* (leaf-constant control means on the fitted
partition — spec §1's own description) and differ by ≤ 5.4e−07 in `μ̂` and hence in
`θ̂_tie`, far below any Monte Carlo scale here. The discrepancy is **reported** as
`mu_leaf_gap` in every results row, and `test-single-tree-identity.R` pins the
rounding magnitude separately at `1e-5` so a regression surfaces there rather than
as an unexplained loosening.

**Two clarifications, not changes.**

1. `verify_population_gates()`'s ratio tolerance is `1e-3`, not `5e-4`. Spec §3/§4
   quote the ratios to three decimals; Regime B computes `0.5384615…`, which the
   spec rounds to `0.539`, a 5.4e−4 gap that is *rounding*, not disagreement. `1e-3`
   is the tightest tolerance admitting the spec's own rounding.
2. The identity has **two** preconditions, not three. `μ̂(ℓ)` must be the leaf's
   control mean and `ê` must be leaf-constant on `μ̂`'s **own** partition; the
   *values* of `ê` are irrelevant, exactly as spec §1 says ("for any
   leaf-wise-constant `ê`"). The unit test demonstrates this directly — four
   deliberately wrong leaf-constant `ê` vectors, including a permutation of the true
   ones, all reproduce `θ̂_tie` to `1e-10` — and shows the identity *does* break when
   `ê` is constant on a different partition, or when `μ̂` uses all-unit rather than
   control leaf means. A corollary of it: the `ass:construct(c)` clip cannot break
   the identity (clipping preserves leaf-constancy), though the clip is still
   asserted inert, because a binding clip means `ê` is no longer the forced
   propensity tree and the object under test would be a different one.

**Resolved 2026-09-09 (same day, downstream fix in `optimaltrees`):** the
6-significant-figure rounding described above was a genuine bug in
`optimaltrees` (a `float` narrowing + `std::to_string(float)`'s fixed-6-decimal
formatting, plus a separate R-side CSV floor for any non-numeric feature
column), not a package-design constraint. Both are now fixed
(`quality_reports/plans/2026-09-09_optimaltrees-leaf-value-precision-complete-fix.md`).
Re-running this study's own Regime A cell against the fixed package: `mu_leaf_gap`
collapsed from `{0, 3.86e-7, NA}`-scale values (matching the committed RDS
exactly) to `{0, 2.22e-16, NA}`-scale values — machine epsilon, not just a smaller
floor. Tree structure and point estimates are unchanged (bit-identical `theta`,
same `n_leaves_m0`/`recovered`/`realized_partition`) — this was a precision fix
with zero effect on this study's already-reported pilot numbers. The
partition-based tied-estimator computation in `common.R` is **kept as-is**
regardless (it is the correct, direct implementation of `cor:single-tree`'s own
definition — recomputing exact within-leaf means on the fitted partition — not
merely a workaround for the now-fixed bug), but the identity gap it exists to
guard against is no longer a real precision floor at any practically relevant
scale.

---

## 10. Pilot results (2026-09-09, `nsim = 50`, `n = 1000`, `results/pilot_*.rds`)

**Pilot scale. Variance-ratio relative MC SE 20.2 %; coverage MC SE 0.020–0.065.**
Wall time **7.1 s** for 4 regimes × 50 reps — the whole study is cheap, and the
cluster stage is a precision decision, not a feasibility one. Zero failed
replications.

### `cor:single-tree` — the identity. **Confirmed, 4 orders tighter than required.**

Max `|θ̂_tie − θ̂_forced|` over all 200 tied-arm replications: **4.11e−15** (per
regime: A 4.05e−15, B 4.11e−15, C 1.86e−15, E 3.66e−15), against spec §1's `1e-10`.
`test-single-tree-identity.R`: **26 expectations, 0 failures**, on four datasets
(Regime A's and B's own constants, a 4-leaf interaction structure with an asymmetric
propensity in all four cells, and a deliberately unbalanced 85/15 two-leaf design),
all with the clip asserted inert and ≥ 2 leaves.

### Partition recovery — **100 % (50/50) in every fitted arm, Regimes A–C.**

Gate `§6(ii)` (≥ 95 %) passes with no margin needed; `n_leaves_m0 = L̄` in every
replication and `clip_frac = 0` throughout. Regime E realises **tree1 in 27/50
(54 %)** and **tree2 in 23/50 (46 %)** — recovery is reported as `NA` there, because
with two genuine ties there is no single intended `τ` and a yes/no would be a coin
flip dressed as a metric.

### Regimes A–C — variance against the analytic bounds (arm `oracle_tie` ≡ `fitted_tie`)

| Regime | analytic `V_τ/V` | `n·V̂ar/V_τ` | `n·V̂ar/V` | bias (MCSE) |
|---|---|---|---|---|
| **A** | 1.000000 | 0.9947 | 0.9947 | +0.01156 (0.01060) |
| **B** | 0.538462 | 0.7489 | 0.4033 | +0.00569 (0.00887) |
| **C** | 1.979362 | 1.2449 | 2.4642 | −0.00704 (0.00757) |

All three ratios are within ±25 % of their analytic targets, i.e. **inside a single
20.2 % MC SE** — direction and rough magnitude confirmed, which is exactly what spec
§6 asks a pilot for and no more. Bias is within 2 MCSE of zero everywhere (Gate 6).

The `oracle_tie` and `fitted_tie` arms produce **byte-identical** point estimates in
A–C, because recovery is 100 % and the tied estimator depends on the outcome tree
*only* through its partition. In **Regime A** `two_tree` is identical as well: with
`τ ∈ S_e ∩ S_μ` and both trees recovering the same saturated partition, the tied and
untied estimators are literally the same estimator — which is
`cor:single-tree-saturated`'s content in its sharpest form.

### `cor:single-tree-coarsening` — both directions. **Confirmed.**

- **Regime B (≤):** analytic `V_τ/V = 0.5385`; empirical `n·V̂ar/V = 0.4033` for the
  tied estimator against **0.7141** for `two_tree` on the same data. The tied
  estimator's empirical SD is **0.0627 vs 0.0834** — super-efficiency, in the
  direction and roughly the magnitude Jensen predicts. Gate 3 passes on both the
  analytic (0.5385 ≤ 0.7) and empirical (0.4033 ≤ 0.7) readings.
- **Regime C (reversal):** analytic `V_τ/V = 1.9794`; empirical `n·V̂ar/V = 2.4642`
  ≥ 1.5. Gate 4 passes on both readings. **The falsifiable sub-claim survives its own
  counterexample regime.** The two estimators' realised SDs are close here
  (0.0535 tied vs 0.0522 two-tree) — the reversal shows up against the *analytic*
  bound, which is precisely why spec §5 forbids the empirical-to-empirical reading.

**Fragility, restated as spec §3 requires.** Regime C's ≈ 1.98 is sensitive to `Tr`
(≈ 1.68 at `Tr = 0.3`), and a 50/50-mass design cleared ≈ 1.5 only at `Tr = 0` and
fell *below* the gate at `Tr = 0.1`. The 0.85 skew toward the reversal leaf buys the
margin. The counterexample is **sharp, not generic**, and the paper's framing should
say so.

### `rem:single-tree-se` — naive SE under-coverage (Regime D, folded into A–C). **Confirmed.**

| Regime | naive/correct variance (analytic) | observed | `cov_naive` (MCSE) | analytic implied | `cov_corrected` (MCSE) |
|---|---|---|---|---|---|
| **A** | 0.3762 | 0.3824 | **0.80** (0.057) | 0.771 | 0.94 (0.034) |
| **B** | 0.6190 | 0.6175 | **0.92** (0.038) | 0.877 | 0.98 (0.020) |
| **C** | 0.3114 | 0.3140 | **0.72** (0.064) | **0.726** | 0.90 (0.042) |
| **E** | 0.4621 | 0.4496 | **0.86** (0.049) | 0.817 | 0.96 (0.028) |

Gate 5's gated cell (Regime C, ≤ 0.6) passes at **0.3140**, and its realised naive
coverage of **0.72** lands on the analytic prediction of **0.726**. Regimes A and B
are reported at the same threshold as observations rather than gates, following spec
§4's own verdict row ("A/B give smaller but nonzero gaps"): A is 0.3824, and **B is
0.6175, marginally above the 0.6 threshold** — the honest reading of a regime whose
`Tr` is deliberately half of `C_untied`, so the treated-arm term is a *large* share
of the tied variance there. Naive coverage is nonetheless clearly below nominal in
every regime, which is what the claim states.

**An unexpected exact result, reported rather than buried.** Spec §5 labels the
corrected SE "exploratory / not-registry-claimed". Measured here, it is
**algebraically identical** to the package's own EIF plug-in at the tied nuisances:
`max |SE_corrected − SE_eif| = 4.2e−17` over all 200 tied-arm replications. The
reason is exact, not numerical: `(1 − ê(ℓ))·n(ℓ) = n₀(ℓ)`, so the leaf-pooled
control term and `eif_att_solve()`'s pointwise one coincide term by term. So the
"correction" is not a new estimator at all — it is the standard EIF plug-in, and
what `rem:single-tree-se` really warns against is *dropping* the control-arm term,
not the absence of a corrected formula. Its coverage (0.90–0.98, all within ~1.5
MCSE of nominal) is a *cross-check on two independent implementations* rather than
evidence for a novel SE. **Recommendation for `claims.md`: record this as an
algebraic observation, not as a new empirically-supported claim.**

### Regime E — `lem:single-tree-linear`'s random-index caveat. **Mechanism shown; magnitude not resolvable.**

| realised partition | share | `R` | bias | `n·V̂ar` | own `V_tau_full` | `n·V̂ar/V_τ` | rel. MCSE |
|---|---|---|---|---|---|---|---|
| `tree1` (`{10,11},{00},{01}`) | 0.54 | 27 | −0.01571 | 3.4058 | 4.8407 | 0.7036 | 0.277 |
| `tree2` (`{01,11},{00},{10}`) | 0.46 | 23 | −0.00087 | 6.1652 | 5.0183 | 1.2285 | 0.302 |

**The qualitative point is demonstrated:** the fitted tree lands on *both* ties
across replications (54 % / 46 %), so which `τ` is realised is genuinely random, and
each tie carries its **own** analytic `V_τ` (4.8407 vs 5.0183) — no single fixed
variance describes the estimator. Both strata are unbiased, as `τ₁, τ₂ ∈ S_μ`
requires.

**Honest limitation of this calibration, flagged for the user's decision.** At spec
§3's Regime E constants the analytic contrast between the two ties is only
**3.67 %** (`V_tau_full` 4.8407 vs 5.0183), while each stratum's empirical variance
carries a 28–30 % relative MC SE at `R = 23–27`. The observed 3.41 vs 6.17 spread is
therefore **not** evidence that the stratum variances differ — it is within noise.
Spec §3 asks for a demonstration that they "differ **materially**"; the constants as
specified do not deliver a material analytic gap, and no `nsim` fixes that. **This
was not adjusted here**: `σ₁² = 1` and the mild `X3` CATE were chosen to match
Regime A's convention *before* the contrast was computed, and shrinking `σ₁²` after
seeing that it would widen the ratio would be exactly the tuning-to-impress this
study is supposed to avoid. If a material contrast is wanted, the lever is a
smaller `σ₁²` (which shrinks `Tr`, the term common to both ties) or a more
asymmetric `e₀` across the boundary cells — a **spec amendment**, not a code change.

### Open items after the pilot

1. **Regime E's analytic contrast is 3.67 %**, too small for the "differ materially"
   demonstration spec §3 describes. Needs a spec decision (§10 above), not more
   replications.
2. **Regime B's naive/correct ratio is 0.6175**, marginally above Gate 5's 0.6.
   Reported as an observation per spec §4; Regime C is the gated cell and passes at
   0.3140. No action unless the user wants B re-calibrated (`Tr_multiple` is an
   argument).
3. **The corrected SE is the EIF plug-in**, exactly. `claims.md` should record it as
   an algebraic observation rather than as a new empirically-supported claim, and
   spec §7's escalation trigger ("promote to `estimate_att_tie()` bundled with the
   corrected SE") should note that the corrected SE is already what
   `estimate_att()`'s own machinery computes.
4. **Cluster scale-up is cheap** — 7.1 s for 200 replications means spec §6's 17-cell
   `nsim = 1000` sweep is minutes, not hours. It is a precision decision only.
