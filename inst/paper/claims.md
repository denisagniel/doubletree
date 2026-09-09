# Claims Registry — doubletree manuscript

Seeded 2026-09-09. This registry did not exist before this date (a gap against this
project's own `paper-sequencing-gate.md` hard gate 2, which requires `claims.md` to be
updated in the same turn as any prose introducing a claim or numbered result — retroactively
backfilled here from a full read of `manuscript.tex`, not reconstructed from memory).

**Scope:** `doubletree/inst/paper/manuscript.tex` only. `theory.tex` is the companion proof
document for the same claims (Condition-(P), the general anchor argument) and is not a
second, independent claims surface. Claims belonging to the separate continuum-covariate /
Rashomon "companion paper" track (per `theory.tex`'s own crosswalk table) are noted where
relevant but not registered here — they are out of scope for THIS manuscript's evidentiary
needs.

---

## Table 1: Numbered claims and their simulation evidence

| Label | Type | Statement (one line) | Line | Sim evidence | Status |
|---|---|---|---|---|---|
| `ass:causal` | Assumption | Identification + **uniform** positivity bound $1-e_0(x)\ge c$ | 120 | S5 R4 (ST3 weak-overlap boundary) stresses the consequence of this bound nearly failing | Consequence tested |
| `ass:construct` | Assumption | Clipping keeps fitted $1-\hat e\ge c$ | 194 | Implicit in every estimator's clip convention (S1/S2/S5) | Definitional, not separately tested |
| `rem:two-trees-why` | Remark | No-splitting works because the candidate class is fixed/finite, not from the bilinear product-rate condition alone | 188 | S5 **C3** (flagship vs. its own cross-fit fallback, matched $n$) | **Covered** |
| `lem:selection` | Proposition | Selected partitions land in the sufficient class w.p. $\to1$ | 239 | S1 (`partition_recovery_clt`) | **Covered** |
| `prop:selection-rate` | Proposition | Exponential rate for that selection consistency | 257 | S1, if it reports the rate and not just the limit — **not independently re-verified this pass** | Covered (unverified rate check) |
| `cor:consistency` | Corollary | $\thetahat$ is consistent | 284 | S1 | **Covered** |
| `lem:uniform` | Proposition | Uniform asymptotic linearity over the whole sufficient class | 299 | S1 (underlies its CLT check) | **Covered** |
| `prop:P-instantiation1` | Proposition | **Headline CLT**: $\sqrt n(\thetahat-\theta_0)\rightsquigarrow N(0,V)$ under exact sparsity | 319 | S1 | **Covered** |
| `cor:variance` | Corollary | Plug-in $\Vhat$ consistent, nominal-coverage CIs; finite-sample downward bias $O((2\bar L+1)/n)$ | 347 | S1 (coverage), S5 **R6** ($\bar L\times n$ boundary sweep) | **Covered** |
| `cor:single-tree` | Corollary | Algebraic identity: single-outcome-tree plug-in **is** doubletree with tied trees | 396 | None (algebra, not data-dependent) | Not sim-testable as stated |
| `lem:single-tree-linear` | Lemma | Uniform asymptotic linearity of the **tied** estimator; random-index expansion | 415 | **None.** `single_tree_coverage` was believed to cover this but targets a different, excluded object (see Table 2) | **GAP** |
| `cor:single-tree-saturated` | Corollary | Under saturated leaves, tied estimator attains the same bound $V$ | 438 | **None** | **GAP** |
| `cor:single-tree-coarsening` | Corollary | Coarsening + homoskedasticity $\Rightarrow$ tied estimator super-efficient ($V_\tau\le V$); **heteroskedasticity can reverse this** (worked example: $\approx0.006$ vs $\approx0.30$, a $\sim$50x reversal) | 443 | **None.** No existing DGP varies $\sigma_0^2$ against $e_0$ to test the reversal | **GAP — includes a specific falsifiable sub-claim with zero coverage** |
| `rem:single-tree-se` | Remark | Naive SE for the tied estimator **under-covers** | 458 | **None** (see Table 2 — `single_tree_coverage` targets a different estimand) | **GAP** |
| `rem:single-tree-distinct` | Remark | Distinguishes this corollary family from a separate shared-structure-tree construction (LOO stability, working-model estimand) | 463 | N/A — a scoping remark, not itself a claim | N/A |
| `prop:bilinear` | Proposition | Exact bilinear-remainder bias decomposition | 474 | Population identity; numerically self-verifying wherever used (e.g. S5's DGP-A construction) but not separately validated for doubletree's own fitted nuisances | Used, not independently stress-tested |
| `lem:biasbound` | Lemma | $\lvert\text{bias}\rvert\le \pi^{-1}D_wD_\mu$ | 493 | Partial overlap with `remainder_terms_value_recovery`'s $D_w,D_\mu$ measurements, but that study tests a different (companion-track) empirical-process bound, not this population bound | Not directly tested |
| `ass:rate` | Assumption | Anchor's fold-wise nuisance rate condition | 519 | S2 regime (b), S5 **R5** (both needed a deliberately-underpowered anchor to actually exhibit violation) | **Covered** |
| `thm:anchor` | Theorem | **Headline validity claim**: anchor interval retains $\ge1-\alpha$ coverage regardless of sparsity | 544 | S2 (`honest_inference_sparsity_failure`) | **Covered** |
| `cor:width` | Corollary | Anchor width $O_p(\max\{n^{-1/2},D_wD_\mu\})$ | 579 | S2 | **Covered** |
| `prop:spectest` | Proposition | Fidelity diagnostic: size $\to0$, consistent under separated-margin alternative | 602 | **None found** | **GAP** |
| `thm:regular` | Theorem | $\thetahat$ attains the full nonparametric bound $V$ (regular, efficient) | 631 | S1 ($\Vhat\to_pV$), S5 **C1** (RMSE ratio to oracle-AIPW $\approx1$) | **Covered** |
| `prop:sparse-bound` | Proposition | The smaller submodel bound $V_\tau$ is **not** exploited; doubletree attains $V$, not $V_\tau$ | 646 | Not independently tested (would pair naturally with the `cor:single-tree-coarsening` gap above) | **GAP** (shares root cause with the single-tree gap) |
| `sec:lasso-comparison` (prose) | Discussion | doubletree vs.\ Lasso/GLM: same two-tier resolution shape, different mechanism; complementary representational trade-offs | 687–694 | S5 **C1** (GLM main-effects genuinely misspecified vs.\ doubletree/oracle) | **Covered** |
| `sec:discussion` (prose) | Discussion | Explicit limitation: no demonstration on real/simulated payment-triage-sentencing data; applied comparison "left to future work" | 699–701 | S5 partially discharges the simulated half; the sentence itself needs rescoping once S5 has results (tracked separately, not yet done) | Partially covered; manuscript edit still pending |

---

## Table 2: Simulation studies and what they actually validate

Built from a factual re-audit 2026-09-09 (session `ses_f7979be5dffeWBKqR31LGmBc9b`) —
**not** from study names, which were found to be actively misleading in three cases.

| Study | Claims targeted | Verdict |
|---|---|---|
| `partition_recovery_clt` (S1) | `lem:selection`, `prop:selection-rate`, `cor:consistency`, `lem:uniform`, `prop:P-instantiation1`, `cor:variance`, `thm:regular` | **Load-bearing, keep** |
| `honest_inference_sparsity_failure` (S2) | `ass:rate`, `thm:anchor`, `cor:width` | **Load-bearing, keep** |
| `head_to_head_comparison` (S5) | `rem:two-trees-why` (C3), `cor:variance` finite-sample note (R6), `ass:rate` violation (R5), `thm:regular`/`prop:sparse-bound` (C1), `sec:lasso-comparison` (C1) | **Load-bearing, keep** |
| `single_tree_coverage` | **None of the above.** Targets `notes/single-tree-inference.tex`'s working-model estimand $\theta^*_n$ via LOO-stability/Rashomon selection — the object `rem:single-tree-distinct` explicitly sets *aside* from this manuscript's corollary family. | **Companion-paper track, not this manuscript** |
| `single_tree_inference` | None — never implemented (all code files are TODO stubs, 2026-04-29). Superseded by `single_tree_coverage`. | **Dead, superseded** |
| `propensity_loss_choice` | None — a design/default-selection ablation (log-loss vs.\ squared-error for the propensity tree), explicitly self-described as "not adjudicating a theoretical validity question." | **This paper's estimator, not a claim validator** |
| `remainder_terms_value_recovery` | None directly. Tests the refining-grid/continuum-covariate "P5 (L2-closeness/value-recovery)" empirical-process bound, rooted in `refining-grid-standalone.tex`, which is explicitly the companion document to this fixed-grid manuscript. Coincidentally shares the control-weighted norm with `lem:biasbound` but tests a different bound. | **Companion-paper track, not this manuscript** |
| `six-approach-arbitration` | None (audited 2026-09-09, `quality_reports/2026-09-09_six-approach-arbitration-audit.md`) — all 7 arms either architecturally superseded (not `estimate_att()`) or Rashomon/msplit machinery on the companion track. | **Companion-paper track / superseded, pending harvest-then-archive decision** |
| `functional_consistency` | None — an msplit/averaged-tree architecture comparison from before the no-splitting flagship architecture existed. | **Orphaned, superseded** |
| `threshold_superconsistency` | None — tests data-adaptive threshold localization, which the manuscript's own setup (analyst-fixed, data-independent cutpoint grid, `manuscript.tex` L80–82) structurally excludes by construction. Self-labelled "not a production simulation study." | **Companion-paper track (`discretization-theory-section.tex`, `refining-grid-standalone.tex`), not this manuscript** |

---

## The gap this registry surfaces

**Zero existing simulations validate the `cor:single-tree` corollary family** (`lem:single-tree-linear`, `cor:single-tree-saturated`, `cor:single-tree-coarsening`, `rem:single-tree-se`) or `prop:sparse-bound`/`prop:spectest`. This is not a minor omission:

- `cor:single-tree-coarsening` makes a **specific, falsifiable, sign-reversing** claim (heteroskedasticity flips a super-efficiency result into a large loss, worked numerically to a ~50x swing in the manuscript's own prose) that no simulation currently probes in either direction.
- `rem:single-tree-se` makes a **direct coverage claim** (naive SE under-covers) — exactly the kind of claim this project's constitution requires stress-testing for, and exactly the kind `single_tree_coverage` was believed, before this audit, to already supply.
- `prop:spectest`'s diagnostic (used to decide, in practice, whether to trust the plain CI or fall back to the anchor) has never had its size or power checked empirically.

These corollaries were drafted and proof-audited 2026-09-08 (session notes,
`session_notes/doubletree-2026-09-08.md` §2.3) but a simulation study for them was never
built — the closest-named existing study (`single_tree_coverage`) targets a different,
excluded object entirely.

---

## Determining the "sims we are going to use" list — from this registry

**Definitely in active use (load-bearing for this manuscript, keep compliant):**
`partition_recovery_clt` (S1), `honest_inference_sparsity_failure` (S2),
`head_to_head_comparison` (S5).

**Not in use for this manuscript (companion-paper track or dead) — do not spend cluster-compliance effort here unless the companion paper is confirmed live:**
`single_tree_coverage`, `single_tree_inference`, `remainder_terms_value_recovery`,
`six-approach-arbitration`, `threshold_superconsistency`.

**Not a claim validator, low stakes either way:** `propensity_loss_choice`
(design-decision sim; small, non-urgent).

**Orphaned, candidate for archive regardless of companion-paper status:** `functional_consistency`.

**New gap, no existing study — needs a decision, not just a compliance fix:** a
single-outcome-tree simulation study covering `lem:single-tree-linear`,
`cor:single-tree-saturated`, `cor:single-tree-coarsening` (including the heteroskedasticity
reversal), `rem:single-tree-se`, and ideally `prop:spectest`. This does not exist anywhere in
the current `simulations/` tree under any name.
