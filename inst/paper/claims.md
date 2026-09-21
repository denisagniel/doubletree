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
| `ass:sparsity` | Assumption | **New, 2026-09-17, not yet in `manuscript.tex`.** Structural sparsity, $\cS_e\ne\emptyset$ and $\cS_\mu\ne\emptyset$ — formalizes what was previously six separate unlabeled prose restatements (`:226,:248,:322,:417,:445,:470`) into one `\ref`-able assumption, placed after `def:sufficient` (`:114`) per `outline.md`'s §2 row for `sec:trees`. The single most load-bearing condition in the paper; every efficiency claim below rests on it. | New (proposed at ~`:115`) | S1 (`partition_recovery_clt`)'s DGPs are constructed to satisfy this by design | **Proposed, accepted 2026-09-17 — not yet drafted** |
| `ass:pseudo` | Assumption | **New, drafted 2026-09-17, added to `manuscript.tex` §3.3 (head of section, before `prop:bilinear`).** Pseudo-true margin: for each $j\in\{e,\mu\}$, a unique $\tau_j^\dagger$ minimizes $R^{(j)}$ over $\cT_{\Lbar}$. Strictly weaker than structural sparsity. Promoted from `prop:spectest`'s previously-unlabeled inline hypothesis; adapted from `theory.tex`'s own `ass:pseudo` (theory.tex line 636), with $V_*$ deliberately dropped (see `prop:pseudo-consistency` row below). | New (~`:472–478`) | None — new, no simulation targets it yet. `theory.tex`'s own `ass:pseudo` is likewise untested empirically. | **Drafted, not yet audited against a live simulation** |
| `prop:pseudo-consistency` | Proposition | **New, drafted 2026-09-17.** Under `ass:pseudo`, selection lands on the pseudo-true partitions w.p.$\to1$, $\hat\eta\to_p\eta_*$, and $\thetahat-\theta_*=o_p(1)$ — consistency (not a CLT) for "the best effect fittable by a sparse tree" when structural sparsity fails but a unique best approximation exists. Audited by a `proof-auditor` pass (2026-09-17, session `ses_f4f2239a5ffeAooiC3YAPgU6Dj`): the mechanism (reusing `lem:selection`'s margin-comparison argument with $\Delta_j^\dagger$ standing in for $\Delta_j$) checks out; the audit caught and required fixing three defects in an earlier proof-block draft (a citation to a proof `lem:selection` does not have in this document; a false claim that $\psi$ is bounded under the wrong assumption pairing; an invalid same-data LLN interchange) — all three are why the final version states the mechanism in prose, matching §3.3's own zero-inline-proof convention, rather than as a `\begin{proof}` block. Confirmed: no CLT with a computable variance $V_*$ exists anywhere in `theory.tex` either — `ass:pseudo`'s own $V_*$ was deliberately not carried into the manuscript version. **Now also previewed in the Introduction** (folded into the contribution paragraph's non-sparse-setting discussion, drafted 2026-09-17, compressed into 4 sentences 2026-09-18 at the author's request — no longer a standalone paragraph), which required its own precision pass — a first draft conflated "the trees converge" with "the point estimate converges" and let an "exact tie" qualifier appear to scope over the anchor interval's unconditional validity; both fixed and preserved through the later compression (final sentence keeps the anchor's independence explicit: "requires neither sparsity nor a unique approximation to remain valid"). | New (~`:479–486`) | None — genuinely new result, no simulation targets the pseudo-true-consistency claim specifically. | **Drafted, audited, previewed in the Introduction, no simulation coverage — new gap, parallel to the single-tree family's** |
| `rem:saturated-general` | Remark | **New, drafted 2026-09-18.** At the saturated budget, the tree class represents exactly the grid-measurable functions, which strictly contains any generalized linear working model's fit whenever the truth departs from additivity on the link scale; the cost of that generality is leaves ($2^k$ for $k$ active covariates, robust to coefficient ties since leaves are axis-aligned unions of atoms, not arbitrary level sets). Corrects the manuscript's prior framing ("neither representation is uniformly preferable," `sec:lasso-comparison`/`sec:discussion`), which the author identified as understating trees' representational generality — general (unconstrained) trees strictly dominate GLMs; only doubletree's *chosen, legibility-motivated* leaf budget makes the comparison two-sided. Audited by `domain-reviewer` (2026-09-18): caught and fixed a critical overclaim (an earlier draft said a saturated tree represents "any function on the grid," contradicted by §2.2's own "$\cX$ need not be finite" — corrected to grid-measurable functions only), a link-scale subtlety (a log-link GLM can represent what looks like an interaction on the identity scale), and an unsupported empirical claim ("real nuisance functions rarely..." — downgraded to "an empirical question this paper does not settle"). | New (~`:103`, §2.2 `sec:trees`, after $\cT_{\bar L}$'s definition) | None — new representational claim, not itself the kind of thing a simulation study tests | **Drafted, audited, no simulation coverage** (not applicable — this is a deterministic/combinatorial claim, not a statistical one) |
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
| `lem:single-tree-linear` | Lemma | Uniform asymptotic linearity of the **tied** estimator; random-index expansion | 415 | **Updated 2026-09-17**: `single_tree_coverage` still targets a different, excluded object (unchanged finding — see Table 2). But `simulations/single_tree_corollaries/` (commit `8138ed2`, 2026-09-10 — postdates this row's original 2026-09-09 "None" finding) implements Regime E specifically for this lemma's stabilization caveat. **Pilot-scale only (nsim=50).** | **Pilot evidence, not scale-validated** (was: GAP) |
| `cor:single-tree-saturated` | Corollary | Under saturated leaves, tied estimator attains the same bound $V$ | 438 | **Updated 2026-09-17**: `single_tree_corollaries` Regime A tests this directly (`tables/pilot_summary_20260909-145931.csv`, regime A rows: `ratio_emp_V`≈0.995, consistent with $V_\tau=V$ under saturation). Pilot-scale only (nsim=50, n=1000). | **Pilot evidence, not scale-validated** (was: GAP) |
| `cor:single-tree-coarsening` | Corollary | Coarsening + homoskedasticity $\Rightarrow$ tied estimator super-efficient ($V_\tau\le V$); **heteroskedasticity can reverse this** (worked example: $\approx0.006$ vs $\approx0.30$, a $\sim$50x reversal) | 443 | **Updated 2026-09-17**: `single_tree_corollaries` Regime B tests the $\le$ direction (super-efficiency under homoskedasticity) and Regime C tests the falsifiable reversal sub-claim. Regime C's pilot `ratio_emp_V`≈2.46 (tied variance exceeds full $V$) — **runs in the direction the reversal claim predicts.** Pilot-scale only (nsim=50, n=1000); not yet promoted to the cluster stage (spec's nsim=1000, n∈{1000,4000}, 17 cells). | **Pilot evidence in the predicted direction, not scale-validated** (was: GAP — includes a specific falsifiable sub-claim with zero coverage) |
| `rem:single-tree-se` | Remark | Naive SE for the tied estimator **under-covers** | 458 | **Updated 2026-09-17**: `single_tree_corollaries`'s pilot table reports `cov_naive` vs. `cov_corrected`/`cov_eif` per regime (e.g. regime A: 0.80 vs. 0.94/0.94) — under-coverage in the predicted direction at pilot scale. Folded into regimes A–C as a metric, not a separate regime, per that study's own README §0. | **Pilot evidence, not scale-validated** (was: GAP) |
| `rem:single-tree-distinct` | Remark | Distinguishes this corollary family from a separate shared-structure-tree construction (LOO stability, working-model estimand) | 463 | N/A — a scoping remark, not itself a claim | N/A |
| `prop:bilinear` | Proposition | Exact bilinear-remainder bias decomposition | 474 | Population identity; numerically self-verifying wherever used (e.g. S5's DGP-A construction) but not separately validated for doubletree's own fitted nuisances | Used, not independently stress-tested |
| `lem:biasbound` | Lemma | $\lvert\text{bias}\rvert\le \pi^{-1}D_wD_\mu$ | 493 | Partial overlap with `remainder_terms_value_recovery`'s $D_w,D_\mu$ measurements, but that study tests a different (companion-track) empirical-process bound, not this population bound | Not directly tested |
| `ass:rate` | Assumption | Anchor's fold-wise nuisance rate condition | 519 | S2 regime (b), S5 **R5** (both needed a deliberately-underpowered anchor to actually exhibit violation) | **Covered** |
| `thm:anchor` | Theorem | **Headline validity claim**: anchor interval retains $\ge1-\alpha$ coverage regardless of sparsity | 544 | S2 (`honest_inference_sparsity_failure`) | **Covered** |
| `cor:width` | Corollary | Anchor width $O_p(\max\{n^{-1/2},D_wD_\mu\})$ | 579 | S2 | **Covered** |
| `prop:spectest` | Proposition | Fidelity diagnostic: size $\to0$, consistent under separated-margin alternative | 602 | **None found** | **GAP** — statement now cites `ass:pseudo`/`prop:pseudo-consistency` by label (2026-09-17 edit) instead of an inline, unlabeled hypothesis; evidentiary status unchanged, still no simulation. **2026-09-21:** the Discussion (`sec:discussion`, first paragraph) now explicitly discloses this GAP to the reader — states that the diagnostic's finite-sample size/power are unvalidated, notes its dependence on the anchor's own (unverifiable) rate condition `ass:rate`, and clarifies it does not establish coverage for an interval chosen on the diagnostic's outcome. Prose-level disclosure only; the simulation gap itself is unchanged. |
| `thm:regular` | Theorem | $\thetahat$ attains the full nonparametric bound $V$ (regular, efficient) | 631 | S1 ($\Vhat\to_pV$), S5 **C1** (RMSE ratio to oracle-AIPW $\approx1$) | **Covered** |
| `prop:sparse-bound` | Proposition | The smaller submodel bound $V_\tau$ is **not** exploited; doubletree attains $V$, not $V_\tau$ | 646 | Not independently tested (would pair naturally with the `cor:single-tree-coarsening` gap above) | **GAP** (shares root cause with the single-tree gap) |
| `sec:lasso-comparison` (prose) | Discussion | doubletree vs.\ Lasso/GLM: same two-tier resolution shape, different mechanism. **The representational-tradeoff content (trees vs.\ linear bases) executed 2026-09-18**: moved into the Discussion section (`sec:discussion`, a new paragraph after the audit-trail paragraph, before the "what this paper has not done" paragraph) — adapted, not a straight copy; not yet independently audited (flagged for a future pass, same as the other pending Discussion additions). `sec:lasso-comparison` itself ends on the rate-condition-comparison paragraph. | Body restructuring verified 2026-09-18 (43 unique labels, all `\ref`s resolve); Discussion paragraph itself not yet audited | S5 **C1** (GLM main-effects genuinely misspecified vs.\ doubletree/oracle) | **Covered** |
| *(ATT-vs-ATE scope, no label — prose)* | Scope decision | **Decided 2026-09-18; drafted 2026-09-21 (`manuscript.tex`, `sec:discussion`, third paragraph); updated same day once Appendix~\ref{app:ate} landed.** ATT stays the paper's primary estimand throughout. The bilinearity of the estimating equation's bias in the two nuisance errors (`prop:bilinear`'s structure) is not itself ATT-specific: `Appendix~\ref{app:ate}` (new, see its own row below) now sketches — does not prove — how the mechanism extends, splitting into one term per outcome regression. The Discussion prose was updated to point there and to state explicitly that the sketch "stops short of a selection-consistency result, a central limit theorem, or an anchor interval," so no claim of a proven extension is made. Also notes the ATE requires both outcome regressions (three trees, not two). **Explicitly NOT generalized, with a corrected attribution:** the single-tree corollary family (§3.5) is ATT-specific because it is stated for the ATT's efficient score and efficiency bound, drawing its *interest* (not its premises — `cor:single-tree` itself is a pure algebraic identity) from `hahn1998`'s asymmetry; only `rem:two-trees-why`'s *second* point (propensity-nuisance informativeness) rests on that asymmetry — its first point (propensity tree as an interpretability, not validity, choice) is not ATT-specific and is not claimed to be. | `manuscript.tex:750` | None — scope decision, no numbered result yet | **Drafted** (`domain-reviewer`-audited 2026-09-21; see `quality_reports/reviews/2026-09-21_discussion-completion.md` and `2026-09-21_ate-appendix.md`) |
| `eq:ate-id` | Equation (identification) | **New, drafted + audited 2026-09-21.** ATE identification, $\theta_0^{\mathrm{ATE}}=\E[m_1(\bX)-\mu_0(\bX)]$, under consistency, ignorability, and pointwise (not necessarily uniform) two-sided positivity. `domain-reviewer` independently re-derived and confirmed correct; caught and required fixing a false comparative claim in the first draft ("unlike the ATT's own identification result, this does not involve the propensity" — the ATT's `eq:att-id` is equally propensity-free; both identification results are, the EIF is where the propensity enters for both estimands). | New, `app:ate` | None — expository identification result, not a claim requiring simulation evidence | **Drafted, audited, correct** |
| `eq:score-ate` | Equation (efficient influence function) | **New, drafted + audited 2026-09-21.** Standard AIPW-form efficient influence function for the ATE in the nonparametric model, for candidate triple $\eta^{\mathrm{ATE}}=(e,m,\mu)$; cites `hahn1998`. `domain-reviewer` verified all four signs correct and confirmed mean-zero at the truth by direct substitution. | New, `app:ate` | None | **Drafted, audited, correct** |
| `rem:bilinear-ate` | Remark | **New, drafted + audited 2026-09-21.** Bilinear bias decomposition for the ATE's EIF, splitting `prop:bilinear`'s single (propensity-error)×(outcome-error) term into two, one per treatment arm. `domain-reviewer` independently re-derived the full identity from the tower property (both signs, including the control-arm term's sign, which does not flip) and found it exactly correct; also found — and this session fixed — a genuine missing hypothesis in the first draft (no stated bound on the candidate propensity $e$, needed for the two expectations to be finite and the tower-property exchange to be valid); added "$c'\le e(\bx)\le1-c'$" as an explicit premise. Also surfaced a non-obvious corroborating identity not previously stated anywhere: the ATE's control-arm bias term is exactly `prop:bilinear`'s own remainder, up to the ATT's $\pi^{-1}$ normalization — now noted in the remark's own text. Explicitly hedged, in the subjunctive, that a three-nuisance analogue of structural sparsity would zero both terms — not verified whether the fixed-finite-class selection argument, the CLT, or the anchor interval actually extend to three simultaneously selected partitions. | New, `app:ate` | None — new, no simulation targets a result not yet proven | **Drafted, audited, correct as far as it goes — explicitly not a proof of the further extension** |
| `sec:discussion` (prose) | Discussion | Explicit limitation: no demonstration on real/simulated payment-triage-sentencing data; applied comparison "left to future work" | 699–701 | S5 partially discharges the simulated half; the sentence itself needs rescoping once S5 has results (tracked separately, not yet done) | Partially covered; manuscript edit still pending |

---

## Table 2: Simulation studies and what they actually validate

Built from a factual re-audit 2026-09-09 (session `ses_f7979be5dffeWBKqR31LGmBc9b`) —
**not** from study names, which were found to be actively misleading in three cases.

| Study | Claims targeted | Verdict |
|---|---|---|
| `partition_recovery_clt` (S1) | `lem:selection`, `prop:selection-rate`, `cor:consistency`, `lem:uniform`, `prop:P-instantiation1`, `cor:variance`, `thm:regular` | **Load-bearing, keep.** Dated summary tables exist (`tables/summary_*.md`, 2026-09-01); no pilot caveat found on a quick 2026-09-17 audit — likely reportable in the new §4 Simulation evidence section, not independently reconfirmed line-by-line. |
| `honest_inference_sparsity_failure` (S2) | `ass:rate`, `thm:anchor`, `cor:width` | **Load-bearing, keep.** Result cells at `r300` scale (300 reps, n=500) plus a `sweep.log`, per a 2026-09-17 audit — more mature than a small pilot, completeness not independently confirmed. |
| `head_to_head_comparison` (S5) | `rem:two-trees-why` (C3), `cor:variance` finite-sample note (R6), `ass:rate` violation (R5), `thm:regular`/`prop:sparse-bound` (C1), `sec:lasso-comparison` (C1) | **Load-bearing, keep — but NOT yet reportable.** Own README states verbatim (2026-09-09): "Status: implemented 2026-09-09; **pilot only** (nsim ≈ 30). The full sweep at spec §6's nsim = 1000/2000 is a separate cluster/batch job." Found 2026-09-17 during the citations/simulations sequencing check (`outline.md` §5) — this is the long pole before the new §4 can report S5 at all. |
| `single_tree_corollaries` | `cor:single-tree` (unit test, not this study), `lem:single-tree-linear`, `cor:single-tree-saturated`, `cor:single-tree-coarsening` (incl. the heteroskedasticity reversal sub-claim), `rem:single-tree-se` | **NEW ROW, added 2026-09-17 — was missing from this table entirely, despite existing since commit `8138ed2` (2026-09-10).** Pilot-scale only (nsim=50, n=1000; spec's cluster stage is nsim=1000, n∈{1000,4000}, 17 cells, not yet run). Regime-C pilot results run in the direction the reversal claim predicts. Built specifically to fill the "zero simulation coverage" gap this registry flagged for the single-tree corollary family (2026-09-09) — see Table 1's updated rows above. `prop:spectest` is explicitly out of scope for this study (its own README §0, §9). |
| `single_tree_coverage` | **None of the above.** Targets `notes/single-tree-inference.tex`'s working-model estimand $\theta^*_n$ via LOO-stability/Rashomon selection — the object `rem:single-tree-distinct` explicitly sets *aside* from this manuscript's corollary family. | **Companion-paper track, not this manuscript** |
| `single_tree_inference` | None — never implemented (all code files are TODO stubs, 2026-04-29). Superseded by `single_tree_coverage`. | **Dead, superseded** |
| `propensity_loss_choice` | None — a design/default-selection ablation (log-loss vs.\ squared-error for the propensity tree), explicitly self-described as "not adjudicating a theoretical validity question." | **This paper's estimator, not a claim validator** |
| `remainder_terms_value_recovery` | None directly. Tests the refining-grid/continuum-covariate "P5 (L2-closeness/value-recovery)" empirical-process bound, rooted in `refining-grid-standalone.tex`, which is explicitly the companion document to this fixed-grid manuscript. Coincidentally shares the control-weighted norm with `lem:biasbound` but tests a different bound. | **Companion-paper track, not this manuscript** |
| `six-approach-arbitration` | None (audited 2026-09-09, `quality_reports/2026-09-09_six-approach-arbitration-audit.md`) — all 7 arms either architecturally superseded (not `estimate_att()`) or Rashomon/msplit machinery on the companion track. | **Companion-paper track / superseded, pending harvest-then-archive decision** |
| `functional_consistency` | None — an msplit/averaged-tree architecture comparison from before the no-splitting flagship architecture existed. | **Orphaned, superseded** |
| `threshold_superconsistency` | None — tests data-adaptive threshold localization, which the manuscript's own setup (analyst-fixed, data-independent cutpoint grid, `manuscript.tex` L80–82) structurally excludes by construction. Self-labelled "not a production simulation study." | **Companion-paper track (`discretization-theory-section.tex`, `refining-grid-standalone.tex`), not this manuscript** |

---

## The gap this registry surfaces

**Updated 2026-09-17 — narrowed, not closed.** A pilot study (`single_tree_corollaries`,
commit `8138ed2`, 2026-09-10) now exists for the `cor:single-tree` corollary family
(`lem:single-tree-linear`, `cor:single-tree-saturated`, `cor:single-tree-coarsening`,
`rem:single-tree-se`) and, in its regime-C results, runs in the direction the
heteroskedasticity-reversal sub-claim predicts. This is not the same as the claim being
validated:

- The pilot is nsim=50 at a single $n=1000$; the study's own spec calls for a cluster stage
  at nsim=1000, $n\in\{1000,4000\}$, 17 cells, not yet run. "Pilot-scale, direction-consistent"
  is the honest current status, not "covered."
- `cor:single-tree-coarsening`'s reversal sub-claim is still the one with the most riding on a
  single numerical worked example (`manuscript.tex:451`, "roughly a fiftyfold reversal") — the
  pilot's regime-C ratio (≈2.46, tied variance exceeding the full-model $V$) is qualitatively
  consistent but the magnitude has not been checked against the manuscript's own worked number.
- `prop:sparse-bound` and `prop:spectest` remain **fully GAP**, unaffected by this pilot
  (`single_tree_corollaries`'s own README explicitly scopes `prop:spectest` out).

Per `outline.md`'s 2026-09-17 restructuring, this corollary family is being relocated to its own
subsection (new §3.5, "A special case: the single-outcome-tree plug-in") specifically so that if
the cluster-stage run does not promote the pilot's direction, softening or hedging the claim is a
one-subsection edit rather than surgery on the headline CLT narrative it currently sits inside.

---

## Determining the "sims we are going to use" list — from this registry

**Definitely in active use (load-bearing for this manuscript, keep compliant):**
`partition_recovery_clt` (S1), `honest_inference_sparsity_failure` (S2),
`head_to_head_comparison` (S5 — **but see the sequencing note below: this one needs a
full-scale cluster run before it can be reported**), and, as of 2026-09-17,
`single_tree_corollaries` (pilot stage; cluster stage not yet run).

**Sequencing note, added 2026-09-17 (`outline.md` §5):** the author has decided the manuscript
needs a Simulation evidence section (new §4) and real citations before further drafting. Of the
three previously load-bearing studies, `head_to_head_comparison` is explicitly self-labeled
"pilot only" and cannot be reported until its cluster-stage sweep runs — this is the long pole,
not a writing task. `partition_recovery_clt` looks reportable now (dated full-run tables);
`honest_inference_sparsity_failure`'s `r300`-scale cells look plausibly reportable but
completeness has not been independently confirmed.

**Not in use for this manuscript (companion-paper track or dead) — do not spend cluster-compliance effort here unless the companion paper is confirmed live:**
`single_tree_coverage`, `single_tree_inference`, `remainder_terms_value_recovery`,
`six-approach-arbitration`, `threshold_superconsistency`.

**Not a claim validator, low stakes either way:** `propensity_loss_choice`
(design-decision sim; small, non-urgent).

**Orphaned, candidate for archive regardless of companion-paper status:** `functional_consistency`.

**Resolved 2026-09-17 (was: "New gap, no existing study"):** a single-outcome-tree simulation
study covering `lem:single-tree-linear`, `cor:single-tree-saturated`, `cor:single-tree-coarsening`
(including the heteroskedasticity reversal), and `rem:single-tree-se` now exists —
`single_tree_corollaries` (commit `8138ed2`, 2026-09-10), pilot stage. `prop:spectest` remains
explicitly out of scope for it (that study's own README §0, §9) and is still a fully open gap.

---

## Table 3: Finite-$X$, legibility, and coarsening bias — **accepted 2026-09-17**

Registered here per this project's hard gate 2 — **before** the corresponding prose exists,
not after — because the gate requires a claim to have an anchor before it is written, and
writing the anchor down first is what makes that checkable. All five rows below were **accepted
by the author 2026-09-17**, including `prop:coarsen-bound` (previously the one open decision).
See `outline.md`'s updated §2 for placement — note the "Proposed location" column below uses the
**pre-2026-09-17 section names**; `sec:finite-discussion` is now §3.4 "What the fixed candidate
class restricts" (content unchanged, heading and container changed) and `sec:discussion` is now
§5. None of these rows are in `manuscript.tex` yet.

| Label (tentative) | Type | Statement (one line) | Proposed location | Anchor it rests on | Status |
|---|---|---|---|---|---|
| *(regime ladder, no label — prose)* | Motivation, not a claim | Names three regimes (audit-complete / sparse / general covariate space) to motivate why legibility favors finite $X$, without asserting finite $X$ is required | Introduction | None needed — explicitly motivation, not a contribution; must not be written as a claim | **Accepted** |
| *(legibility-vs-audit-completeness remark, no label — prose)* | Remark | Distinguishes legibility (delivered by the fixed cutpoint grid regardless of $\cX$) from audit-completeness/fidelity (requires atom-measurable nuisances — restates `ass:sparsity`, does not add to it) | `sec:trees`, after `manuscript.tex:88` | `def:sufficient`; **now also `ass:sparsity` directly, by label** (new numbered assumption, accepted 2026-09-17 — see Table 1) | **Accepted** |
| `rem:coarsening` (tentative) | Remark | Discretizing a continuous covariate onto the grid is a special case of sparsity failure that `thm:anchor` already covers; explicitly disclaims "residual confounding" terminology (epidemiology usage: Becher; Brenner & Blettner) in favor of "coarsening bias" / leaf-constancy misspecification, because `ass:causal` conditions on the full $X$ and identification is unaffected — only the working-model error in `prop:bilinear`'s remainder is | `sec:honest-manuscript` (now §3.3), after `manuscript.tex:598` (`cor:width`) | `prop:bilinear`, `lem:biasbound`, `thm:anchor` — **no new theory required** | **Accepted** |
| `prop:coarsen-bound` (tentative) | Proposition | Under a Lipschitz sensitivity device on $e_0,\mu_0$ (constants $L_e,L_\mu$, local to this proposition only) and grid mesh $h$: $\delta_e\lesssim L_eh$, $\delta_\mu\lesssim L_\mu h$; combined with a finite-class ERM excess-risk step bridging $\delta$ to the realized $D_w,D_\mu$ that `cor:width` actually uses, gives anchor-interval width $=O_p(\max\{n^{-1/2},h^2\})$ | Same location as `rem:coarsening` (§3.3) | Builds on `lem:biasbound`, `cor:width`; the ERM-bridging step is the one piece of genuinely new technical work in this plan | **Accepted — needs its own simulation coverage; none exists yet (see `outline.md` §5)** |
| *(narrow, don't delete, `:666–671`)* | Amendment to existing prose | Points `sec:finite-discussion`'s (now §3.4's) existing "we do not minimize this cost in general" at `prop:coarsen-bound` | `sec:finite-discussion` (now §3.4, "What the fixed candidate class restricts") | `prop:coarsen-bound` | **Accepted** (no longer conditional — `prop:coarsen-bound` is accepted) |
| *(narrow, don't delete, `:701`)* | Amendment to existing prose | Narrows (not deletes) the "genuine misspecification... left to future work" clause in `sec:discussion` (now §5), since `prop:coarsen-bound` would be a first cut at exactly this | `sec:discussion` (now §5) | `prop:coarsen-bound` | **Accepted** (no longer conditional) |

**What is settled vs. not, as of this registry entry (updated 2026-09-17):**

- **Settled:** all five rows above, including `prop:coarsen-bound` — the only remaining open item
  is that it needs a dedicated simulation study (none exists) before it can be marked "Covered"
  rather than merely "drafted," per this table's own evidentiary standard.
- **Explicitly rejected, not merely deferred:** stating "finite $X$ is required" as a new
  assumption. This would contradict `manuscript.tex:79` and `theory.tex:594–618`'s existing,
  deliberate "escape hatch... does not survive" framing, and is not being carried forward in
  any of the rows above.
