# partition_recovery_clt (S1)

> ## ⚠ CORRECTION, 2026-09-01 (later the same day): read this before any number below
>
> **The F1/ST2 half of the first run was invalid.** `code/dgps.R` shipped with
> the two nuisances on **disjoint** covariate supports (`e_0` on `(X1, X2)`,
> `mu_0` on `(X3, X4)`), which makes `A` independent of `Y(0)`: the DGP had **no
> confounding**, a difference in means was already unbiased, and **no amount of
> partition-recovery failure in either nuisance could move `theta_hat`.** The
> defect was diagnosed in `quality_reports/2026-09-01_partition-recovery-clt-design-audit.md`
> §2 and pinned as item 1 of the spec's own §8 addendum, but the fix was **never
> applied to the code that actually ran**.
>
> This is fatal to precisely this study's reason for existing. Spec §4 reports
> recovery *alongside* coverage to separate "Condition (P) held and the CLT
> delivered" from "(P) failed but the bias happened to be small." On an
> unconfounded DGP the second branch is unreachable **by construction**, and the
> first run demonstrated exactly that pathology: at ST2 `n = 200` it reported
> `rec_both = 0.001` (recovery essentially never happened) with coverage
> `0.942` — indistinguishable from nominal.
>
> **Fix:** `mu_0` is now gated on `X1`, the coordinate that also gates `e_0`, so
> the supports share a covariate and the DGP is genuinely confounded
> (`E[Y(0)|A=1] = 0.55` vs `E[Y(0)|A=0] = 0.30` on F1; naive difference-in-means
> bias `+0.25`). Leaf values were recalibrated for the confounded geometry
> (`e = 0.20/0.45/0.75`, `mu = 0.10/0.50/0.90`; ST2 gap `0.20`), because sharing
> the gate depletes control mass where the propensity is high and thereby shrinks
> `Delta_mu`. `code/calibrate_margin.R` (new; a required deliverable of the spec's
> §8 addendum) recomputes every margin and its `lambda_n` crossing from the
> enumeration helper, so no constant here is a mystery number.
>
> **Status of the artifacts:**
>
> * **Authoritative results:** `tables/summary_20260901-185633.md` and the
>   `*_20260901-185633.csv` / `figures/*_20260901-185633.*` set. 18 cells x 1000
>   reps, 0 failures, all fitted partitions inside `T_Lbar`, 9.6 min of fit time.
> * The 12 `f1` / `f1_weak` checkpoints were deleted and re-run. **ST3's 6
>   checkpoints were never affected** — `st3_weak_overlap` uses
>   `generate_dgp_weak_overlap()`, whose two nuisances were always both on
>   `(X1, X2)` — and were kept, not re-run, so every ST3 number is unchanged.
> * **Everything below this notice, from "The question" through "Deviations", is
>   the write-up of the SUPERSEDED run.** Its ST3 numbers still hold; its F1,
>   ST1 and ST2 numbers do not, and neither do the conclusions drawn from them.
>   It is retained deliberately as the record of the defective run rather than
>   quietly overwritten. **The narrative has not been re-derived against the
>   corrected numbers** — that is interpretive work, not a mechanical patch, and
>   is left as the next task.
>
> **What actually changed in the headline numbers** (corrected run; recovery is
> the enumeration-based `tauhat_j in S_j` indicator, not a same-partition check):
>
> | regime | n | coverage | bias | rec_e | rec_mu | rec_both |
> |---|---|---|---|---|---|---|
> | F1  | 2000 | 0.946 |  0.0001 | 1.000 | 1.000 | 1.000 |
> | F1  | 4000 | 0.965 |  0.0009 | 1.000 | 1.000 | 1.000 |
> | F1  | 8000 | 0.944 |  0.0002 | 1.000 | 1.000 | 1.000 |
> | ST1 |  200 | 0.933 |  0.0002 | 0.376 | 0.587 | 0.228 |
> | ST1 |  500 | 0.955 |  0.0015 | 0.913 | 0.953 | 0.871 |
> | ST1 | 1000 | 0.952 | -0.0019 | 0.998 | 0.997 | 0.995 |
> | ST2 |  200 | 0.940 | -0.0042 | 0.088 | 0.097 | 0.011 |
> | ST2 |  500 | 0.942 | -0.0044 | 0.395 | 0.284 | 0.117 |
> | ST2 | 1000 | 0.954 |  0.0024 | 0.799 | 0.538 | 0.429 |
> | ST2 | 2000 | 0.960 |  0.0008 | 0.996 | 0.881 | 0.877 |
> | ST2 | 4000 | 0.945 | -0.0001 | 1.000 | 0.996 | 0.996 |
> | ST2 | 8000 | 0.952 |  0.0006 | 1.000 | 1.000 | 1.000 |
> | ST3 |  200 | 0.868 |  0.0282 | 0.381 | 0.000 | 0.000 |
> | ST3 |  500 | 0.886 |  0.0003 | 0.899 | 0.000 | 0.000 |
> | ST3 | 1000 | 0.920 |  0.0018 | 0.999 | 0.002 | 0.002 |
> | ST3 | 2000 | 0.935 |  0.0018 | 1.000 | 0.001 | 0.001 |
> | ST3 | 4000 | 0.944 |  0.0020 | 1.000 | 0.012 | 0.012 |
> | ST3 | 8000 | 0.943 | -0.0013 | 1.000 | 0.050 | 0.050 |
>
> Corrected exact population structure (`|S_e| = |S_mu| = 11` on both grid-sparse
> DGPs, `|T_4| = 356` at p = 5):
>
> | dgp | `Delta_kl(e)` | `Delta_kl(mu)` | `lambda_n < Delta` from (e / mu) |
> |---|---|---|---|
> | `f1`      | 0.011944 | 0.010599 | 1000 / 1000 |
> | `f1_weak` | 0.005087 | 0.002941 | 2000 / 4000 |
> | `st3_weak_overlap` (unchanged) | 0.021481 | 0.000191 | 500 / never |

**Study S1 of 3** for Paper B's Simulations section
(`quality_reports/specs/2026-09-01_paperb-simulations-decomposition.md`).
Spec: `quality_reports/specs/2026-09-01_partition-recovery-clt.md`.
Run 2026-09-01, 18 cells x 1000 replications = **18,000 `estimate_att()` fits**,
0 failures, 9.6 min of fit time sequential.

Covers §3.1-3.2 of `doubletree/inst/paper/manuscript.tex`: Condition (P),
Instantiation 1 (grid-exact sparsity). Claims exercised: `thm:main`,
`cor:variance`, `lem:selection`, `prop:selection-rate`, `lem:uniform`.

---

## The question, and why this study ran first

Every efficiency and coverage claim in §3 routes through Condition (P) holding.
S2 (the anchor-interval study) needs its exactly-sparse cell to be a valid null.
`propensity_loss_choice/tables/summary_20260821-113704.md` showed
`same_partition_rate_e = 0.53` at n = 2000 on an exactly-sparse DGP, with mean
leaf count still drifting upward -- suggesting (P1), the partition-recovery half
of Condition (P), might be too fragile at practical n for S2 to interpret
anything.

**That 0.53 was not a (P1) failure rate.** It is the rate at which two
*propensity-loss arms* selected the same partition as each other -- agreement
between arms, not membership in the sufficient class \(\mathcal S_j\). This study
measures the actual quantity, by exact enumeration.

---

## Answer, stated plainly

### Does Condition (P)'s selection-consistency claim (P1) hold cleanly at practical n?

**Yes -- when the population risk margin is real, and the leaf penalty is small
enough to permit the necessary split. Both conditions are checkable in advance,
and neither is automatic.**

| | |
|---|---|
| Wide margin (F1, \(\Delta^{\mathrm{KL}}_e = 0.0254\), \(\Delta^{\mathrm{KL}}_\mu = 0.0127\)) | \(\Pr(\hat\tau_j \in \mathcal S_j) = \mathbf{1.000}\) for **both** nuisances at n = 2000, 4000, 8000 (1000/1000 reps each); 0.991 already at n = 500 |
| Weak margin (ST2, same structure, \(\Delta^{\mathrm{KL}}\) ~10x smaller) | (P1) fails in **97% / 91% / 69% / 52%** of reps at n = 200 / 500 / 1000 / 2000; only reaches 1.000 at n = 8000 |
| Weak overlap (ST3) | \(\hat\tau_e\) recovers by n = 1000, but \(\hat\tau_\mu\) **never** does (0.0-5.0% across the whole grid) |

**The single quantity that decides it is \(\lambda_n\) versus \(\Delta_j\).**
`estimate_att()` minimises \(R_n(\tau) + \lambda_n|\tau|\) with
\(\lambda_n = \log n / n\). Buying the split that a sufficient partition requires
costs one leaf, i.e. \(\lambda_n\), and buys \(\Delta_j\). So whenever
\(\lambda_n \ge \Delta_j\), **the penalised population optimum lies outside
\(\mathcal S_j\)** and recovery is excluded structurally -- no sampling noise
required. Every cell's `penalised_opt_in_S` column records this, and the
empirical transitions line up with the crossings, not with the concentration
bound:

| regime | nuisance | crossing \(\lambda_n = \Delta_j\) | empirical failure rate around it |
|---|---|---|---|
| ST2 | e | between n = 2000 and 4000 | 0.280 -> **0.015** |
| ST2 | mu | between n = 4000 and 8000 | 0.039 -> **0.000** |
| ST3 | mu | never crosses (\(\Delta_\mu = 0.00019\)) | stays >= 0.95 at every n |

This is a mechanism `prop:selection-rate` does not model: its bound
\(|\mathcal T_{\bar L}|\exp(-c_1 n \Delta_j^2)\) contains no \(\lambda_n\) (the
proposition only assumes \(\lambda_n \to 0\)). At these sample sizes the
\(\lambda_n\)-vs-\(\Delta_j\) crossing predicts the empirical transition better
than the concentration bound does.

### What this means for S2 -- read this before running S2

1. **S2's exactly-sparse (\(\varepsilon = 0\)) cell *can* be a valid null.** On
   F1, (P1) holds at 1.000 with 1000 reps at every n >= 2000. The feared
   fragility is not a property of the estimator; it is a property of DGPs with
   small margins.
2. **But it is not a valid null for free.** S2 must compute \(\Delta_e\) and
   \(\Delta_\mu\) with `code/enumerate_sufficient_class.R` (the shared helper this
   study exists to build) and confirm \(\lambda_n < \Delta_j\) -- with margin to
   spare, a factor of ~3 or more -- at **every** n in its grid.
3. **Concretely: S2 must not reuse `generate_dgp_weak_overlap()` or
   `generate_dgp_simple()` as its \(\varepsilon = 0\) baseline.**
   `generate_dgp_weak_overlap()` has \(\Delta^{\mathrm{KL}}_\mu = 0.00019\), so
   \(\mu\)-recovery is structurally impossible at every n up to 8000 -- its
   "exactly sparse" cell is *not* a Condition-(P) cell in practice, even though
   \(\mathcal S_\mu \ne \emptyset\) formally. The additive-in-logit
   `generate_dgp_simple()` family has the same defect for the same reason (leaf
   values spanning only [0.38, 0.55]), which is what produced the 0.53 figure.
4. **S2 should report the recovery indicator per cell rather than infer it from
   coverage** -- because, as below, coverage was nominal even where (P1) failed
   in 97% of replications.

---

## The result that most changes how §3 should be read

**Coverage was at nominal 95% in every cell where propensities were not clipped
-- including cells where (P1) failed in 97% of replications.**

ST2 at n = 200: \(\Pr(\hat\tau_j \in \mathcal S_j)\) = 0.001 for the pair, yet
coverage = 0.942 (MC se 0.007) and bias = 0.0021. Same story at every ST2 cell.

So on these DGPs **(P1) is sufficient for the CLT but far from necessary.** The
reason is visible in the margin itself: an ST2 "failure" is a partition that
merges two leaf values 0.14 apart, so the resulting nuisance is *nearly* correct
(\(\delta_j\) small but nonzero), and the EIF's double robustness absorbs it.
Condition (P) is a convenient sufficient condition for the proof route, not a
sharp characterisation of when the interval works.

This cuts both ways for the paper, and both should be said:

* It is **good news for robustness** -- the interval is more forgiving than
  Condition (P) suggests.
* It is **bad news for using coverage as evidence** -- a coverage table alone
  cannot distinguish "Condition (P) held and the CLT delivered" from "(P) failed
  and the resulting bias happened to be negligible." Only the enumeration-based
  indicator separates them, which is exactly why spec §4 required both.

---

## ST1/ST3: does the required disentanglement actually distinguish the causes?

**Yes, decisively.** Spec §3 required ST3 to report, per replication, both the
recovery indicator and the propensity-clip fraction, so a coverage shortfall
could be attributed rather than left undifferentiated. It attributes cleanly:

| n | coverage (all) | clip fraction | coverage \| clipped | coverage \| unclipped | \(\Pr(\hat\tau_\mu\in\mathcal S_\mu)\) |
|---:|---:|---:|---:|---:|---:|
| 200 | 0.868 | 0.079 | **0.732** (310 reps) | 0.929 (690) | 0.000 |
| 500 | 0.886 | 0.081 | **0.833** (324) | 0.911 (676) | 0.000 |
| 1000 | 0.920 | 0.044 | 0.933 (179) | 0.919 (819) | 0.002 |
| 2000 | 0.935 | 0.019 | 0.909 (77) | 0.938 (922) | 0.001 |
| 4000 | 0.944 | 0.006 | 0.920 (25) | 0.956 (963) | 0.012 |
| 8000 | 0.943 | 0.0005 | 1.000 (2) | 0.966 (948) | 0.050 |

**The shortfall tracks clipping, not partition recovery.** Two independent
readings, and they agree:

* Coverage rises monotonically 0.868 -> 0.943 while the clip fraction falls
  0.079 -> 0.0005, and among *unclipped* replications coverage is already
  0.929-0.966 throughout.
* \(\Pr(\hat\tau_\mu \in \mathcal S_\mu)\) is ~0 at **every** n, including
  n = 8000 where coverage is 0.943. If \(\mu\)-recovery failure were the cause,
  coverage could not have improved by 7.5 points while that indicator stayed
  flat at zero.

This resolves the confound the spec flagged from `propensity_loss_choice`
(89.7-90.3% coverage with 36% of fits clipped at n = 500): **`ass:construct`'s
clipping bound biting, not (P1) failing.** Without the enumeration indicator both
studies would have had to report an undifferentiated "coverage is off."

**One caveat, stated because the numbers look alarming.** In ST3 the handful of
replications where \(\mu\) *was* recovered have coverage 0.00 (2, 1 and 12 reps at
n = 1000, 2000, 4000) with biases of 0.23, -0.16, -0.11. That is conditioning-on-a-
rare-event selection, not a finding: recovering the full 4-leaf \((X_1,X_2)\) grid on
the control subset under weak overlap requires an atypical control composition,
and those same samples give poor ATT estimates. Reported rather than dropped, but
it should not be read as "recovery hurts."

## ST1: small-n fragility did not materialise on a well-separated DGP

| n | \(\Pr(\text{both}\in\mathcal S)\) | coverage | bias | RMSE |
|---:|---:|---:|---:|---:|
| 200 | 0.723 | 0.954 | 0.0015 | 0.1099 |
| 500 | 0.991 | 0.940 | 0.0007 | 0.0666 |
| 1000 | 0.998 | 0.957 | -0.0005 | 0.0464 |
| 2000 (F1) | 1.000 | 0.953 | 0.0003 | 0.0319 |
| 4000 (F1) | 1.000 | 0.954 | 0.0003 | 0.0226 |
| 8000 (F1) | 1.000 | 0.959 | -0.0003 | 0.0158 |

RMSE falls by 0.71 and 0.70 per quadrupling of n (\(1/\sqrt 2 = 0.707\)) --
`thm:main`'s \(\sqrt n\) rate, on the nose. Coverage is within 1.96 MC se
(+-0.013) of 0.95 in every cell.

Note the n = 200 nuance: there \(\lambda_n = 0.0265\) exceeds *both* margins, so
the penalised population optimum was outside \(\mathcal S_j\) -- and recovery
still occurred in 72% of replications, because the empirical risk gap fluctuates
above \(\lambda_n\). The \(\lambda_n \ge \Delta_j\) condition predicts where
recovery is *not guaranteed*, not where it is impossible.

## `lem:uniform`'s signature, checked rather than assumed

Spec §7 asked whether coverage holds *while* \(\hat\tau_j\) varies across
replications within \(\mathcal S_j\). It does:

* In F1 at n = 2000/4000/8000, the \(\mu\) fit visits **11, 9 and 9 distinct
  sufficient partitions** across 1000 replications, while coverage is
  0.953/0.954/0.959.
* 1.4-1.8% of recovered \(\mu\) fits in F1 are strictly **refinements** -- more
  leaves than the 3-leaf minimum sufficient partition needs, i.e. a leaf spent on
  a coordinate the truth ignores. At n = 200 that reaches 12.0%, and in ST2 at
  n = 200, 15.8%.

Those refinement percentages are also the quantitative case for the
currency-correct indicator: **every one of those replications would have been
scored a (P1) FAILURE by a same-partition or same-leaf-count check**, which is why
spec §2 forbade one. `code/verify_enumeration.R` test B3 pins that
misclassification down as an assertion.

## `prop:selection-rate`: the bound holds, its shape is not confirmed

\(c_1\) is not identified by the proposition ("depending only on the outcome
bound, \(c\), and the positive leaf masses"), so per spec §4 it was **not fitted
post hoc**. The bound is reported at the stated convention \(c_1 = 1\), where it
is vacuous everywhere (\(n\Delta_j^2 \le 5.2 < \log|\mathcal T_{\bar L}| = 5.87\),
so \(|\mathcal T_{\bar L}|e^{-n\Delta^2} > 1\)). The checkable content is the
decay shape, via the \(c_1\) each observed rate implies:

| regime, nuisance | implied \(c_1\) across n |
|---|---|
| ST2, e | 4692, 1904, 991, 568, 400 (n = 200 -> 4000) |
| ST2, mu | 11280, 4580, 2385, 1323, 867 |
| ST1, mu | 231, 131, 75 (n = 200 -> 1000) |
| ST3, mu | 691619 ... 17466 |

**The bound is never violated, and by an enormous margin** (implied \(c_1 \gg 1\)
everywhere). But the implied constant **falls monotonically in n** in every
series, so \(\log(1/\text{rate})\) grows *slower* than linearly in
\(n\Delta_j^2\): the observed decay is sub-exponential on the bound's own scale.
The honest reading is that this study **confirms failure rates go to zero and the
bound is satisfied with room to spare, but does not confirm the bound's
exponential-in-\(n\Delta^2\) form** -- consistent with the \(\lambda_n\) mechanism
above, which the bound omits.

---

## Design

Binary covariates throughout (`estimate_att()` requires them; `ass:finite`).
\(\bar L = 4\), \(m_n = 1\), \(\lambda_n = \log n/n\), `propensity_loss =
"log_loss"`, `outcome_type = "binary"` -- i.e. the estimator exactly as it ships.
Y binary, hence bounded, as `prop:selection-rate` requires.

| regime | DGP | n | what it stresses |
|---|---|---|---|
| **F1** | `f1` | 2000, 4000, 8000 | nothing -- wide margin baseline |
| **ST1** | `f1` (same population DGP) | 200, 500, 1000 | selection not yet converged |
| **ST2** | `f1_weak` | 200 ... 8000 | \(\Delta_j\) near zero (near-tie) |
| **ST3** | `generate_dgp_weak_overlap()`, unmodified | 200 ... 8000 | positivity, with clip/recovery disentangled |

All four appear in the same tables and the same figure panels (spec §6): there is
no favourable-only headline.

**DGPs** (`code/dgps.R`) adapt the existing `simple`/`moderate` structure rather
than adding a fourth family, with two changes, both required:
leaf values are set **directly** rather than through a logit link (an
additive-in-logit truth needs the full 4-leaf grid, which makes
\(\mathcal S_j\) a singleton and the margin an accident of the coefficients), and
`e_gap`/`mu_gap` expose the **margin dial** ST2 needs. \(e_0\) is hierarchical on
\((X_1,X_2)\), \(\mu_0\) on \((X_3,X_4)\), \(X_5\) is pure noise; \(\theta_0\) is
computed in closed form (= \(\tau\) = 0.08 exactly, not a realised-sample
average). ST3 reuses `dgps_stress.R::generate_dgp_weak_overlap()` unmodified for
the draw; only its population counterpart is restated, and
`check_spec_vs_draw()` verifies the two agree to 1e-12.

**Exact population structure** (closed form, no simulation):

| dgp | nuisance | \(|\mathcal T_{\bar L}|\) | \(|\mathcal S_j|\) | min sufficient leaves | \(\Delta^{\mathrm{sq}}_j\) | \(\Delta^{\mathrm{KL}}_j\) |
|---|---|---:|---:|---:|---:|---:|
| f1 | e | 356 | 11 | 3 | 0.01000 | 0.02544 |
| f1 | mu | 356 | 11 | 3 | 0.00500 | 0.01272 |
| f1_weak | e | 356 | 11 | 3 | 0.00123 | 0.00251 |
| f1_weak | mu | 356 | 11 | 3 | 0.00079 | 0.00162 |
| st3_weak_overlap | e | 155 | 1 | 4 | 0.00338 | 0.02148 |
| st3_weak_overlap | mu | 155 | 1 | 4 | 0.00009 | 0.00019 |

\(|\mathcal S_j| = 11\) of 356 for F1/ST2 -- **not a singleton**, which is the
whole reason the recovery indicator must be the projection identity rather than a
partition comparison. ST3's classes *are* singletons, deliberately the opposite
corner.

**Two margins are reported on purpose.** `manuscript.tex` defines \(R^{(j)}\) as
a *squared-error* population risk (the display after `eq:risk-emp`), so
`prop:selection-rate`'s \(\Delta_j\) is \(\Delta^{\mathrm{sq}}\). But
`estimate_att()` *selects* with log-loss by default, so the margin that actually
governs these runs is \(\Delta^{\mathrm{KL}}\). \(\mathcal S_j\) itself is
loss-free, so only the margin's value depends on the choice. **This mismatch
between the stated theory and the shipped default is worth a sentence in §3.**

---

## The enumeration helper (shared infrastructure -- S2 imports this)

`code/enumerate_sufficient_class.R` enumerates every partition of \(\{0,1\}^p\)
realizable by a tree of \(\le \bar L\) leaves, deduplicated **as partitions**
(not topologies -- \(|\mathcal T_{\bar L}|\) in the bound counts partitions),
determines \(\mathcal S_j\) membership and \(\Delta_j\) exactly from closed-form
population probabilities, and provides the recovery indicator
\(\Pi_{\hat\tau_j}\gamma_{0,j} = \gamma_{0,j}\).

Three things it gets right that a shortcut would not:

1. **Recovery is evaluated from the definition, not by lookup** in the enumerated
   set, so it stays correct for a fit outside \(\mathcal T_{\bar L}\);
   `in_enumerated_class` is a separate diagnostic (it was TRUE for all 18,000
   fits).
2. **The check runs on the 2^p cell grid, not the sample.** At n = 200 with
   p = 5 several of the 32 cells are empty in a given sample; a sample-based
   check would skip them and could call a non-sufficient partition sufficient.
3. **The tightest competitor is found, not guessed.** For F1's \(\mu\) it is not
   the obvious "merge the two \(X_3 = 1\) cells" partition but one isolating
   \(\{X_1 = 1\}\), where the control weighting \(\nu_x = p_x\{1-e_0(x)\}\) makes
   the same merge **4x cheaper**. A two-point formula would have overstated
   \(\Delta_\mu\) by 4x, and F1 would have been mis-sized.

### It was verified before being trusted (spec §7)

`code/verify_enumeration.R`, **16 assertions, all passing**:

* **A** -- the constructive enumeration is compared against an **independent
  brute-force enumerator** (all set partitions via restricted-growth strings,
  filtered by a from-scratch tree-realizability predicate sharing no code) for
  p in {2,3} x \(\bar L\) in {2,3,4}. A hand-checked count for one configuration
  would not have caught an off-by-one in the budget split; this does.
* **B** -- all 8 members of \(\mathcal T_4\) at p = 2 are classified against a
  hand-worked table written out in the file, including the two
  redundant-refinement cases. **B3 asserts that the naive checks really do
  misclassify them**, so the reason for rejecting the naive check is a test, not
  a comment.
* **C** -- \(\Delta_{\mathrm{sq}} = 0.045\) and
  \(\Delta_{\mathrm{KL}} = (\log 2 - H(0.2))/2\) against closed-form hand
  arithmetic, plus a non-uniform-weight case (\(\Delta = 0.024\) via the
  two-point formula) because the \(\mu\) risk is control-weighted, never uniform.
* **D** -- the fitted-model path end to end against the real solver: TRUE for a
  recovered tree, FALSE for a collapsed one, TRUE for a redundant refinement, and
  a hard error on non-binary X.
* **E** -- bound/implied-\(c_1\) round trip, including both edges where the
  reported bound goes uninformative.
* **F** -- the DGP guards (below).

---

## Files

```
code/enumerate_sufficient_class.R  the shared helper: T_Lbar, S_j, Delta_j, recovery indicator
code/verify_enumeration.R          17 assertions incl. independent brute-force cross-check
code/dgps.R                        f1 / f1_weak (shared-X1 confounding + margin dial) + ST3 population side + guards
code/calibrate_margin.R            deterministic: every (Delta_sq, Delta_kl) + its lambda_n crossing
code/common.R                      design grid, seeds, one replication, one cell
code/run_pilot.R                   spec/draw checks, exact margins, timing projection
code/run_cell.R                    one cell per process:  Rscript ... <dgp> <n> <reps>
code/run_sweep.R                   process-per-cell driver, resumable
code/analyze.R                     metrics, attribution table, figures, tables
code/review_checks.R               review artifact: 7 independent re-derivations (not in the pipeline)
results/cell_<dgp>_n<N>_r1000.rds  18 per-cell checkpoints
tables/summary_<runid>.md          all tables in one document
tables/{population,primary,uniformity,selection_rate,attribution}_<runid>.csv
figures/coverage_bias_recovery_<runid>.{png,pdf}
figures/selection_failure_rate_<runid>.{png,pdf}
```

## Reproduce

```sh
cd doubletree
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
Rscript --vanilla simulations/partition_recovery_clt/code/verify_enumeration.R   # always first
Rscript --vanilla simulations/partition_recovery_clt/code/calibrate_margin.R     # margins + lambda_n crossings
Rscript --vanilla simulations/partition_recovery_clt/code/run_cell.R f1 8000 10  # timing gate at the LARGEST cell
Rscript --vanilla simulations/partition_recovery_clt/code/run_sweep.R            # 18 cells, ~10 min
Rscript --vanilla simulations/partition_recovery_clt/code/analyze.R
Rscript --vanilla simulations/partition_recovery_clt/code/review_checks.R        # optional: re-derivations
```

To re-run only some cells (as the correction above did for the 12 `f1`/`f1_weak`
cells, leaving ST3's 6 checkpoints untouched):
`PRC_CELLS="f1:200,...,f1_weak:8000" Rscript --vanilla .../run_sweep.R`.

Seeds are a pure function of one integer:
`rep_seed(dgp, n, rep) = xxhash64(20260901 | dgp | n | rep)`, the same scheme as
`propensity_loss_choice`. Re-running reproduces the numbers above exactly.

**Memory posture** (this project's convention;
`../docs/MEMORY_SAFE_SIMULATIONS.md`, `MEMORY.md` `[LEARN:rashomon-memory]`): no
Rashomon enumeration anywhere (single best tree only, `certified_* = TRUE` on
every fit), no parallelism (`worker_limit = 1`, BLAS/OpenMP pinned to 1), one
**process per cell**, `gc(full = TRUE)` on a heartbeat, checkpoint per cell.
Sequential total 11.7 min for 18,000 fits (~0.04 s/fit); a pilot measured this
before the full grid was committed to, per spec §5.

---

## Deviations from the spec, and why

1. **F1's n = 200 cell has \(\lambda_n \ge \Delta_j\) for both nuisances**
   (0.0265 vs 0.0254 and 0.0127), so at that one cell the penalised population
   optimum is outside \(\mathcal S_j\) even though the spec calls F1
   "comfortably" separated. This is **not fixable by choosing better leaf
   values**: at \(\bar L = 4\), p = 5 a competitor can always isolate a
   weight-1/4 subcube and merge inside it, capping \(\Delta_\mu\) near 0.0127
   given \(\mu_{\mathrm{hi}} + \tau \le 1\). Reported per cell
   (`penalised_opt_in_S_*`) rather than hidden, and it turned into a finding: the
   estimator's default \(\lambda_n\) rate, not the DGP, is what makes n = 200
   marginal.
2. **ST3 was swept over the full 6-point n grid**, not the {500, 2000} of its
   "own existing sweep" in `propensity_loss_choice`. Matched n across regimes is
   what makes the side-by-side comparison in §6 meaningful, and it is what let
   the clip-vs-recovery attribution be read as a trend rather than two points.
3. **Both \(\Delta^{\mathrm{sq}}\) and \(\Delta^{\mathrm{KL}}\) are reported**,
   resolving an ambiguity the spec did not anticipate: `manuscript.tex` states
   \(R^{(j)}\) as squared error while `estimate_att()` defaults to log-loss.
4. **`prop:selection-rate`'s \(c_1\) is not identified**, so the bound is shown
   at the stated convention \(c_1 = 1\) plus an implied-\(c_1\) column. Per spec
   §4, \(c_1\) was not fitted to make the curves agree.
5. No replication counts were reduced and no regime was skipped: 1000 reps in all
   18 cells, as specified.

## Review

`quality_reports/2026-09-01_partition-recovery-clt-r-review.md`. Score 92/100
(PR gate). Two Critical and one High finding, all fixed; the sweep was then re-run
from scratch with `PRC_FORCE=1` and `primary_*.csv` came out **bit-identical**,
which verifies both the seeding scheme and that the fixes were behaviour-preserving.

The review pass ran **inline**, not as an isolated subagent (the independent
subagent dispatch timed out), so it is a review of the code by its own author.
What was substituted for independence: `code/review_checks.R` re-derives seven
load-bearing numeric claims by routes sharing no code with the implementation --
the partition enumeration against a from-scratch realizability predicate
(two-sided, p = 4), \(\Delta_j\) against both an explicit per-leaf loop and hand
two-point arithmetic, \(\theta_0\) against a 400k Monte Carlo, \(\nu\) against
the empirical control distribution, and `fitted_cell_labels()` against an
independently written tree walker. All seven agree. Two of them failed on first
run and caught real defects (a vacuous random-probe design, and a wrong assumption
that the CI multiplier was 1.96 rather than `qnorm(0.975)`).

## Three defects found while building this

Recorded because all three are the kind that produce plausible numbers rather
than errors.

1. **The two nuisances were on disjoint covariate supports**, so the DGP had no
   confounding and `theta_hat` could not respond to nuisance misfit at all. This
   is the defect in the correction notice at the top of this file; it survived a
   design audit that diagnosed it correctly, because the audit's recommended fix
   was never applied to the code that ran. **A full 18-cell x 1000-rep sweep
   completed on it with 0 failures and every fitted partition inside
   \(\mathcal T_{\bar L}\)** — the only symptom was nominal-looking coverage in
   cells where recovery had essentially never happened. Now guarded:
   `verify_enumeration.R` F1b asserts both that the declared supports intersect
   and that \(E[Y(0)\mid A=1]\ne E[Y(0)\mid A=0]\) in closed form, and
   `calibrate_margin.R` aborts on a disjoint-support DGP.
2. **`mu_0` was keyed on \((X_1, X_3)\) while DECLARING \((X_3, X_4)\)** at an
   even earlier point, *consistently on both the population and the sample side*.
   Nothing failed: the two agreed with each other, every fit succeeded, every
   fitted partition landed inside \(\mathcal T_{\bar L}\), and the only symptom
   was \(\Delta_\mu\) coming out the wrong size. Note that
   \((X_1, X_3)\) is now the *intended* keying — sharing \(X_1\) is defect 1's
   fix — so `verify_enumeration.R` F2 was retargeted to probe \(X_2\) instead.
   `check_variable_roles()` runs on **every** `build_dgp()` call, asserting each
   nuisance depends on exactly the coordinates it declares;
   `check_spec_vs_draw()` covers the adjacent failure mode (population/sample
   disagreement) for every DGP, not just ST3.
3. **`same_partition_rate_e` in `propensity_loss_choice` is not a (P1) recovery
   rate**, and the 0.53 figure that motivated the S1-before-S2 ordering should
   not be cited as one. Worth correcting wherever it is referenced.

## Expected, benign warning

Every tree fit emits `High-dimensional data detected (estimated 155 binary
features after discretization). Setting model_limit=1000000`. The estimate is
made *before* discretization runs and so does not know that all-binary X takes
`discretize_features()`' no-op fast path; the `model_limit` it sets caps the
Rashomon *model count*, irrelevant on this single-tree path. It cannot be
silenced by passing `model_limit` explicitly (`optimaltrees()` then errors with
"formal argument model_limit matched by multiple actual arguments"), so warnings
are muffled and **recorded per replication** (`n_warnings`,
`warning_messages`), with distinct texts reported once per cell. (Two distinct
texts appear across the 18 cells, not one: `155` binary features at p = 5
(`f1`/`f1_weak`) and `124` at p = 4 (ST3).)

## Certification: `certified_* = TRUE` on 17,949 of 18,000 fits, not all of them

An earlier version of this README claimed `certified_e = certified_m0 = TRUE` on
**all** 18,000 fits. That was wrong, and it was wrong before the DGP correction
too: 30 of the 51 uncertified fits are in ST3's checkpoints, which the correction
did not touch. Audited directly over the checkpoints (corrected run):

| | count |
|---|---|
| `certified_e == FALSE` | 2 |
| `certified_m0 == FALSE` | 49 |
| `used_search_e == TRUE` | 19 |
| `used_search_m0 == TRUE` | 77 |
| by DGP | `f1` 15, `f1_weak` 6, `st3_weak_overlap` 30 (unchanged from the first run) |

All 51 are optimaltrees' documented **case (iii), "staircase skip"**
(`bisect_lambda_to_budget.R::classify_lambda_fit()`): the bisection landed at a
`lambda > lambda_n` whose fit uses strictly fewer than `leaf_budget` leaves, so
the certificate lemma's cases (i)/(ii) do not apply and a **computable
suboptimality slack** `gap = (lambda - lambda_n)(Lbar - L)` is returned instead.
`feasible` held on every one; the largest slack observed is `gap_m0 = 0.0237`,
and recovery/coverage among these fits (`rec_e = 0.63`, `cov = 0.94`) is
unremarkable. So this is a rare alternate path through the estimator's own
penalty search, not truncation and not a failure — but the blanket "all 18,000"
claim should not be repeated.

`certified_* = TRUE` together with `used_search_* = FALSE` remains
`prop:parsimony`'s operational signature on the other 99.7% of fits.
