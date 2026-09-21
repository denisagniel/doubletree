# Paper Outline — `doubletree`

**Instance path:** `inst/paper/outline.md`
**Governed by:** `paper-protocol.md` (via the `draft-paper-section` / `edit-paper-with-context` skills)
**Status:** draft, created 2026-09-16, restructured 2026-09-17. The 2026-09-17 pass is a full
renumbering of §3 plus a new §4, decided after a `structure-reviewer` critique
(session-local report, not filed under `quality_reports/reviews/` yet — see maintenance
contract) triggered by the author's observation that the old §3.4/§3.5 "were things we
discussed and then got stuck on the end." Nothing below reflects a change to `manuscript.tex`
itself yet — this file states the *target* structure; the moves are still pending.

> **Why this file exists.** It is a content-distribution plan, not a table of contents: which
> section owns which idea, so a claim is made once, in one place, at full strength.

**Backfill scope.** This is a full backfill of §2 for the existing manuscript (every section
and subsection in `manuscript.tex` gets a row), because the immediate task — adding a
legibility/finite-$X$/coarsening-bias discussion — needs to know what every existing section
already owns to avoid duplicating `manuscript.tex:79`'s or `manuscript.tex:594`'s existing
finite-$X$ material. It is *not* a backfill of `theory.tex` (the companion proof document,
4383 lines, out of scope for this pass) or of the appendix-track files under `inst/paper/notes/`
and the `_*-draft.tex` files, which are not currently `\input` into `manuscript.tex` at all.

---

## 1. Fixed before drafting

| Field | Value |
|---|---|
| **Venue class** | **unknown — reported as a finding, not guessed.** No `**Type:**` line in a governance file (no `CLAUDE.md`/`AGENTS.md` at the project root), no requirements spec under `quality_reports/plans/` naming a target outlet (`2026-09-15_real-data-pipeline-spec-shell.md` is a data-pipeline spec, not a venue declaration), no governance prose naming the paper's type. Resolution chain exhausted at step 4. This matters concretely for one thing: whether new assumptions (e.g. a possible Lipschitz sensitivity device, §2 below) get a formal numbered environment or prose treatment — methods-theory papers use formal environments 47.7% of the time, applied papers 9.0%. Recommend the author state this explicitly rather than have it inferred later. |
| **Target venue** | TBD — not stated anywhere found. |
| **Estimand, in one sentence** | The average treatment effect on the treated, $\theta_0=\E[Y(1)-Y(0)\mid A=1]$ (`manuscript.tex:58`, eq. `att`). Matches `notation.md` row 1. |
| **The one-sentence contribution** | *doubletree*: two size-budgeted, jointly-displayed decision trees (one propensity, one outcome), fit on the full sample with no sample splitting, that attain the semiparametric efficiency bound for the ATT under exact structural sparsity on the analyst's own pre-specified grid, and remain validly (if not efficiently) inferable via an anchor-interval fallback when that sparsity fails. |
| **What this paper does NOT do** | No demonstration on any real or simulated payment/triage/sentencing dataset (`manuscript.tex:701`); the partition-recovery guarantee is pointwise, not uniform, in the data-generating process (`manuscript.tex:44`, `:701`); no claim of matching an unrestricted black-box ML ensemble's flexibility — the comparison defended is against a GLM working model (`manuscript.tex:33`, `:44`). As of this pass, also does **not** claim finite $X$ is required (see §2 rows for `sec:trees`/Introduction below — this is the specific thing the pending edit must not accidentally assert). |

---

## 2. Section plan — TARGET structure (post-2026-09-17 restructuring)

**Status update, 2026-09-18: EXECUTED.** Everything below described the target structure as of 2026-09-17; as of 2026-09-18 the body moves are done (verified: 43 labels, all unique; every `\ref`/`\eqref` resolves; `grep` confirms the six-subsection order below). The Introduction's forward references (`\S\ref{sec:instantiation1}`, `\S\ref{sec:honest-manuscript}`, `\S\ref{sec:finite-discussion}`, `\S\ref{sec:lasso-comparison}`) now point at sections that actually contain what the Introduction says they do.

**All decisions below are settled** (author sign-off 2026-09-17): the Table 3 finite-$X$/
coarsening-bias additions (accepted, previously "proposed"); the §3.4/§3.5 restructuring
(`structure-reviewer`'s "Option B", accepted in full); promoting structural sparsity to a
numbered assumption `ass:sparsity` (accepted); relocating the single-tree corollary family
(accepted); and adding real citations plus a simulation-evidence section (accepted, scope
and sequencing tracked in §6 below — this is new work, not yet drafted, and larger than a
registry edit).

**Why the renumbering happened.** The old §3 ran 497 lines (≈70% of the body) and mixed
inference machinery with two discussion/positioning subsections (old 3.4 "is this basically
parametric?", old 3.5 the Lasso comparison) that a `structure-reviewer` critique found
misplaced for different reasons: old 3.5 had *no* dependency holding it at the end (it was
promised in the abstract and at `:40`, delivered in 8 lines 650 lines later); old 3.4 was
correctly *ordered* (its closing paragraphs need the anchor interval from old §3.3) but wrongly
*contained* — it buried `thm:regular`, arguably the paper's second-strongest result, under a
FAQ-style heading, and it topically answers §2.2 (`sec:trees`), not "bias and inference." The
critique also found the single-tree corollary family (old lines 396–466, inside old §3.2) is
the section's actual worst-placed block: it forward-cites `prop:bilinear` (old §3.3) and
`thm:regular` (old §3.4), a genuine ordering violation neither old 3.4 nor old 3.5 has.

| § | Working title (label) | Owns (target, post-restructuring) | Claims made (→ `claims.md`) | Results delivered | Depends on | Moved from |
|---|---|---|---|---|---|---|
| 1 | Introduction | The audit-trail problem; the gap; the contribution; roadmap. **Current state as of 2026-09-18** (five sequential edit passes, each audited by `domain-reviewer` + `proofreader`): ¶1 (`:38`) — citations for AIPW/DML/TMLE (`robins1994`, `bang2005`, `chernozhukov2018`, `vanderLaanTargetedLearningCausal2011`, `kennedy2016semiparametric`) and the empirical-process/Donsker completeness fix. ¶2 (`:40`) — rewritten by the author to argue cross-fitting itself breaks auditability even with an interpretable learner, sidestepping the citation-fidelity risk of characterizing specific competing papers (`athey2016recursive`/`ChenSyrgkanisAustern2022` were tried here, found to mischaracterize both papers per `domain-reviewer`, and removed — see §5 below for what happened to that content). ¶3 (contribution) — restored the ATT name (author decision 2026-09-18, see the new Appendix/Discussion row below), the pre-specified grid qualifier, a real quantified leaf-budget cost (`cor:width`), aligned the efficiency claim with the abstract, added a citation for "globally optimal decision tree" (`xuStatisticalOptimalityOptimal2026`, `mctavishFastSparseDecision2022`); also fixed a matching contradiction in the **abstract** ("no sample splitting or cross-fitting" → "...for either displayed tree"). **New, folded into ¶3's non-sparse-setting sentences (drafted 2026-09-17, compressed to 4 sentences 2026-09-18 at the author's request — no longer a standalone paragraph):** previews `prop:pseudo-consistency` (§3.3) — when sparsity fails but a unique best sparse-tree approximation exists, doubletree is consistent for it; audited twice before compression (once for fidelity against the body result, once for mechanics) to avoid conflating "the trees converge" with "the point estimate converges" and to stop an "exact tie" qualifier from appearing to scope over the anchor interval's unconditional validity — the compressed version preserves both fixes (ends on "requires neither sparsity nor a unique approximation to remain valid"). ¶4 (roadmap) — unchanged. **Still pending, not yet drafted:** the three-rung regime ladder (Table 3, accepted) and a roadmap pointer to the new §4 Simulation evidence. | Abstract-level claims (not separately tabled) | — | — | — |
| 2.1 | Data, estimand, and notation | Unchanged. | — | — | §1 | — |
| 2.2 | Covariate space, tree partitions, and optimal trees (`sec:trees`) | Unchanged existing content, **plus**: new numbered `ass:sparsity` (accepted) placed immediately after `def:sufficient` (`:114`), stating $\cS_e\ne\emptyset,\cS_\mu\ne\emptyset$ formally for the first time — collapses six unlabeled prose restatements currently at `:226,:248,:322,:417,:445,:470` into one `\ref`-able assumption the rest of the paper cites; the legibility-vs-audit-completeness remark (Table 3, accepted), now explicitly restating `ass:sparsity` by label, not just in spirit; **new forward cross-reference** at ~`:93` ("it is not merely a convention") pointing to §3.4 below, where the parametric-model objection is answered — closes the gap `structure-reviewer` found between where the objection is provoked and where it was previously answered (537 lines away). | `ass:sparsity` (new row, `claims.md` Table 1) | — | §2.1 | — |
| 2.3 | Identification and influence function | Unchanged. | `ass:causal` | — | §2.1, §2.2 | — |
| 2.4 | The doubletree algorithm | Unchanged. | `ass:construct`; `rem:two-trees-why` | — | §2.2, §2.3 | — |
| 3.1 | Partition recovery and the CLT (`sec:partition-recovery`) | Unchanged. **Flagged, not actioned this pass:** `structure-reviewer` found this subsection has zero numbered results and repeats its own scoping disclaimer five times (`:211–220,:230–234,:292–297,:313–314,:326–329`, the last *inside* the headline proposition's statement) — a real finding, deferred pending a separate decision, not part of the four items the author signed off on 2026-09-17. | — | — | §2.2–2.4 | — |
| 3.2 | Instantiation 1: fixed grid, exact sparsity (`sec:instantiation1`) | Exact sparsity setting (now citing `ass:sparsity` by label instead of restating $\cS_e,\cS_\mu\ne\emptyset$); `lem:selection`, `prop:selection-rate`, `cor:consistency`, `lem:uniform`, `prop:P-instantiation1` (headline CLT), `cor:variance`. **New here (hoisted):** `thm:regular` (regularity/efficiency in the full nonparametric model) and `prop:sparse-bound` (submodel bound not exploited), placed immediately after `cor:variance` (old `:394`) under a plain result-naming heading, not a question — fixes the one real ordering violation in the manuscript (old `cor:single-tree-coarsening` at `:455` cited `thm:regular` at `:631`, a forward reference; hoisting makes it backward). **Loses** the single-tree corollary family (→ new §3.5 below): it doesn't belong under this heading (it studies a *different estimator*, $\hat\theta_{\mathrm{tie}}$, for 70 lines with no signal in the section title) and it forward-cited both `prop:bilinear` (§3.3) and `thm:regular` (formerly §3.4) — the section's actual worst ordering problem, worse than the one the author originally flagged. | `lem:selection`, `prop:selection-rate`, `cor:consistency`, `lem:uniform`, `prop:P-instantiation1`, `cor:variance`, `thm:regular`, `prop:sparse-bound` | Prop. `P-instantiation1` (headline CLT); Cor. `variance`; Thm. `regular`; Prop. `sparse-bound` | §2.2 (`ass:sparsity`), §3.1 | `thm:regular`, `prop:sparse-bound` (from old §3.4) |
| 3.3 | When partition recovery fails: the anchor interval (`sec:honest-manuscript`) | **New, drafted 2026-09-17, at the head of the section (before `prop:bilinear`):** `ass:pseudo` (pseudo-true margin — unique population-risk minimizer per nuisance) and `prop:pseudo-consistency` (under `ass:pseudo`, $\hat\theta\to_p\theta_*$, the best-sparse-tree effect — consistency, not a CLT; audited by `proof-auditor`, session `ses_f4f2239a5ffeAooiC3YAPgU6Dj`, three defects in an earlier proof-block draft found and fixed). Then unchanged existing content: `prop:bilinear`, `lem:biasbound`, the anchor construction, `ass:rate`, `thm:anchor`, `cor:width`, `prop:spectest` — the latter's own hypothesis now cites `ass:pseudo`/`prop:pseudo-consistency` by label instead of an inline unlabeled restatement. **Plus** (Table 3, accepted): `rem:coarsening` after `cor:width` (old `:598`) — discretizing a continuous covariate is a special case of sparsity failure `thm:anchor` already covers, with the "coarsening bias, not residual confounding" terminology disclaimer; **and** `prop:coarsen-bound` — Lipschitz sensitivity device ($\delta_e\lesssim L_eh,\ \delta_\mu\lesssim L_\mu h$) plus the finite-class ERM bridging step connecting $\delta$ to the realized $D_w,D_\mu$ that `cor:width` uses, giving width $=O_p(\max\{n^{-1/2},h^2\})$ — **this is the one item with genuinely new technical work (the ERM bridging step) and needs its own simulation coverage; none exists yet (see §5).** | `prop:bilinear`, `lem:biasbound`, `ass:rate`, `thm:anchor`, `cor:width`, `prop:spectest`, `ass:pseudo`, `prop:pseudo-consistency`, `rem:coarsening`, `prop:coarsen-bound` | Thm. `anchor` (headline validity); Prop. `pseudo-consistency` (new, drafted); Prop. `coarsen-bound` (new, not yet drafted) | §2.2, §2.3, §3.2 (`lem:selection`'s margin-comparison mechanism, reused with $\Delta_j^\dagger$) | — |
| 3.4 | What the fixed candidate class restricts | Residual of old `sec:finite-discussion` after `thm:regular`/`prop:sparse-bound` move to §3.2: the framing paragraph (old `:629`), "declines the information rather than exploiting it" (old `:660–664`), "what does restrict the DGP is structural sparsity itself... we do not minimize this cost in general" (old `:666–671`) — **narrowed, not deleted, to point at `prop:coarsen-bound`** now that it's accepted — and the anchor-based closing answer to "what if the class is wrong" (old `:673–685`). Retitled away from the old rhetorical-question heading; this is the subsection §2.2's new cross-reference (row above) points to. | (prose; ties `ass:sparsity`'s cost to the DGP, not the model) | — | §3.2 (`thm:regular`), §3.3 (`thm:anchor`, `cor:width`, `prop:coarsen-bound`) | Renamed/shrunk from old `sec:finite-discussion` |
| 3.5 | A special case: the single-outcome-tree plug-in | Relocated single-tree corollary family: `cor:single-tree`, `lem:single-tree-linear`, `cor:single-tree-saturated`, `cor:single-tree-coarsening`, `rem:single-tree-se`, `rem:single-tree-distinct`. Now correctly ordered — sits after `prop:bilinear` (§3.3) and the hoisted `thm:regular` (§3.2), so both citations inside this family become backward references. Opens with a relocated/restated motivating sentence from `rem:two-trees-why` (old `:190`, "Corollary~\ref{cor:single-tree} and its consequences below make both points precise...") — that motivation currently sits 206 lines upstream of the material it motivates; this move puts it at the point of delivery. **Evidentiary note:** a pilot simulation (`simulations/single_tree_corollaries/`, commit `8138ed2`, 2026-09-10) exists and its regime-C results run in the direction `cor:single-tree-coarsening`'s heteroskedasticity-reversal claim predicts — but it postdates `claims.md`'s "GAP" labeling for this family (seeded 2026-09-09) and is pilot-scale only (nsim=50); see `claims.md` Table 1/2 updates. Isolating this family in its own subsection is what makes a future retreat (if the pilot doesn't promote) a one-subsection edit instead of surgery on §3.2's main narrative. | `cor:single-tree` … `rem:single-tree-distinct` (full list in `claims.md` Table 1) | — | §3.2 (`thm:regular`), §3.3 (`prop:bilinear`) | Relocated from old §3.2 (old `:396–466`) |
| 3.6 | Relation to penalized linear methods (`sec:lasso-comparison`) | **Shrunk, confirmed complete 2026-09-17** (not just planned): keeps the two-tier-resolution-shape-vs-mechanism paragraph and the rate-condition/margin comparison, now ending on "...not a free alternative to the compatibility conditions used for Lasso." **Lost** the representational-tradeoff paragraph → moved to §1 (Introduction), restored there with its concrete mechanism ("a leaf for every configuration of the active covariates") that an earlier draft of the move had dropped, per `domain-reviewer`'s audit. Removed from this location, not duplicated — verified no residual copy remains here. | (prose, `sec:lasso-comparison` row in `claims.md`) | — | §3.3 (`prop:selection-rate`), §3.4 | Lost old `:721` (→ §1) |
| 4 | Simulation evidence *(NEW — placeholder, content not yet drafted)* | Reports on the three load-bearing studies `claims.md` Table 2 already treats as evidence for claims made in §3: S1 `partition_recovery_clt` (selection consistency, headline CLT, variance/coverage — appears to have completed, dated full-run tables, not just a pilot); S2 `honest_inference_sparsity_failure` (anchor validity/width under sparsity failure — has `r300`-scale result cells, completeness not yet confirmed); S5 `head_to_head_comparison` (doubletree vs. forest-DML/linear-DML/oracle-AIPW, incl. the GLM-misspecification comparison that currently backs the §3.6 row above — **explicitly self-labeled "pilot only, nsim≈30" as of 2026-09-09; the full sweep is a separate, not-yet-run cluster job**). **This row cannot be filled in until the sequencing decision in §6 resolves** — specifically, whether S5 gets a full-scale run before this section is drafted. | (prose + tables, once drafted) | — | §3 (all subsections) | New section; previously this evidence existed only in `claims.md`, never in the manuscript |
| 5 | Discussion (`sec:discussion`) | Renumbered from old §4. Unchanged existing content (returns to payment/triage/sentencing; audit-trail framing; limitations; future work), **plus two limitations `structure-reviewer` found omitted, both drafted 2026-09-21:** (i) `prop:spectest`'s size/power is unvalidated (`claims.md`: GAP) — the first Discussion paragraph now states this explicitly, including the diagnostic's dependence on the anchor's own (unverifiable) rate condition, rather than presenting the diagnostic as an unqualified operational decision rule; (ii) the audit-trail claim is narrowed — the point estimate is auditable in every case, but running the diagnostic pulls the (possibly non-interpretable, cross-fit) anchor into the audit trail regardless of outcome, so the interval is fully anchor-free only when the analyst commits to structural sparsity outright, without running the diagnostic. **New, decided 2026-09-18, drafted 2026-09-21:** (iii) a new third paragraph on the ATT-vs-ATE scope decision, citing `hahn1998`, correctly splitting `rem:two-trees-why`'s two points (only its second, propensity-nuisance-informativeness point rests on Hahn's asymmetry) and noting the single-tree family (§3.5) is ATT-specific by construction, not because it depends on Hahn's result directly (it is a pure algebraic identity) — no promise of specific Appendix content, since none is drafted. Both (i)-(iii) audited by `domain-reviewer` + `proofreader` (2026-09-21); see `quality_reports/reviews/2026-09-21_discussion-completion.md`. **New, decided 2026-09-17, not yet drafted:** (iv) the trees-vs-linear-models representational trade-off + its concrete mechanism ("a leaf for every configuration of the active covariates") — this content was removed from both the Introduction and old §3.6 during sequential edits and currently lives nowhere in the manuscript; author has decided it belongs in the Discussion, not the Introduction. **Correction, 2026-09-21:** item (iv) is in fact already drafted — see the second Discussion paragraph (`rem:saturated-general` positioning argument, citing `rudin2019`), landed in the same 2026-09-18 commit that added the remark itself; this row previously listed it as pending in error. | (prose, `sec:discussion` row in `claims.md`) | — | §1, §3, §4 | Renumbered from old §4 |
| Appendix | — | **Drafted 2026-09-21** (`app:ate`, "Extension to the average treatment effect (sketch)"): states the ATE estimand, its identification via both outcome regressions (`eq:ate-id`), the efficient influence function (`eq:score-ate`, citing `hahn1998`), and a Remark (`rem:bilinear-ate`) deriving the bilinear bias structure — split across two outcome-regression terms instead of `prop:bilinear`'s one, verified correct by `domain-reviewer` (2026-09-21) via independent re-derivation. **Explicitly does not** derive a selection-consistency result, central limit theorem, or anchor interval for the resulting three-nuisance case — no numbered proposition/theorem/corollary appears, only equations and one remark, by design, matching the author's "minimal: estimand + identification + sketch" scope decision. Proof content for the ATT case continues to live in the separate companion document `theory.tex`, out of scope for this backfill. Decision unchanged from 2026-09-18: ATT stays the paper's primary estimand throughout; the single-tree corollary family (§3.5) and `rem:two-trees-why`'s justification stay **ATT-only**, since both lean on `hahn1998`'s asymmetry (knowing $e_0$ exactly lowers the efficiency bound for the ATT but not the ATE) — not re-derived for ATE, and the appendix explicitly leaves open (rather than answers) whether a propensity-free single-tree analogue exists for the ATE. | `eq:ate-id`, `eq:score-ate`, `rem:bilinear-ate` (`claims.md` rows, once added) | — | §2.1 (pointer sentence after `eq:att`), §5 Discussion (updated paragraph pointing here) | — |

**Resolved, no longer open:** the Table 3 Lipschitz-sensitivity proposition (`prop:coarsen-bound`)
is accepted, along with the rest of Table 3 — see `claims.md`'s updated Table 3 status column.

---

## 3. Dependency and ordering check

- [x] **No section depends on a result introduced later.** The restructuring in §2 above fixes
      the one violation `structure-reviewer` found (old `cor:single-tree-coarsening` at `:455`
      citing `thm:regular` at `:631`) by hoisting `thm:regular` into §3.2, ahead of its consumer's
      new position in §3.5. All other cross-section citations checked in the 2026-09-17 pass
      point backward: §3.4 cites §3.2/§3.3 only; §3.5 cites §3.2/§3.3 only; §3.6 cites §3.3 only;
      §4 (new) cites all of §3; §5 cites §1/§3/§4.
- [x] **Every assumption is introduced before the first result that uses it** — checked against
      `notation.md`'s "Used by" field this pass. The one gap this check surfaced is now fixed:
      structural sparsity was used (six times, unlabeled) before this pass ever gave it a
      `\ref`-able home; promoting it to `ass:sparsity` in §2.2 (row above) closes this.
- [ ] Range-citation and appendix-numbering checks (`Assumptions~\ref{first}-\ref{last}`,
      appendix-declared-but-body-used assumptions) — **still not checked in this pass.** Deferred;
      lower priority than the structural moves above.
- [x] The discussion (§5) already contains results from substantive sections (`prop:spectest`,
      the audit-trail claim tied to `thm:anchor`/`cor:width`), per the existing text at `:699`.
- [x] `\appendix` is used correctly (`:703`) and the appendix section is genuinely empty rather
      than competing with a body `\section{Appendix}`.

---

## 4. Explicitly not planned here

- Whether to include a roadmap paragraph in the introduction — already present (`:44`); not
  revisited by this edit.
- **Whether to add a dedicated `\section{Related work}` heading** — `structure-reviewer`
  recommends against this specifically (minority convention, not itself a defect) but found the
  underlying deficit real: the Introduction has zero citations while making the paper's
  positioning/novelty claims (`:38`, `:40`). The fix is citations at those two lines (§6 below),
  not a new heading.
- Whether every assumption carries a `[Name]` — mixed already (`ass:causal` named, `eq:score`
  unnamed as an assumption but is not one; `ass:sparsity`, new this pass, is named). Not otherwise
  revisited.
- A full backfill of `theory.tex`'s own section plan — flagged as future work, not blocking.
- **§3.1's five-times-repeated scoping disclaimer and its lack of numbered results**
  (`structure-reviewer` finding, `:211–220,:230–234,:292–297,:313–314,:326–329`) — flagged, not
  actioned. Not one of the four items the author signed off on 2026-09-17.
- **The proof-placement convention** (4 of ~13 numbered results have inline proofs, no stated
  rule for which) and the orphaned "Instantiation 1" label (no "Instantiation 2" anywhere in this
  file) — flagged by `structure-reviewer`, not actioned this pass.

---

## 5. Sequencing: citations and simulation evidence (decided 2026-09-17, not yet executed)

The author confirmed both are needed before further drafting, not deferred to a later revision.
This is new work, not a registry edit — recorded here so the outline stays the single place that
tracks what §1 and §4 (new) actually require before they can be drafted.

**Citations (§1, `:38` and `:40`).** Needed: real references for (i) the AIPW/DML/TMLE
robustness claim at `:38` (candidates not yet confirmed — e.g. Robins/Rotnitzky/Zhao-era AIPW,
van der Laan/Rubin TMLE; `chernozhukov2018` is already in the bib for a different point and may
cover part of this) and (ii) the tree-based-causal-methods gap claim at `:40` ("existing
constructions require sample splitting... or sacrifice interpretability or an asymptotic
guarantee" — currently asserted with zero citations against a literature the paper never names).
`ChenSyrgkanisAustern2022` is already in the bib (cited only at `:465`, 425 lines downstream) and
may belong at `:40` too. **Not yet done: an actual literature search.**

**Simulation evidence (new §4).** A quick audit of `simulations/` this pass found the three
studies backing `claims.md`'s "Covered" claims are at three different stages of readiness, which
directly determines whether drafting §4 is a write-up task or blocked on new compute:
- `partition_recovery_clt` (S1) — dated, sized summary tables exist (`tables/summary_*.md`, Sep
  1); no "pilot" caveat found. Likely reportable as-is, not independently reconfirmed.
- `honest_inference_sparsity_failure` (S2) — result cells at `r300` scale (300 reps, n=500) plus
  a `sweep.log`; more mature than a small pilot but completeness not independently confirmed.
- `head_to_head_comparison` (S5) — **its own README states, verbatim, "Status: implemented
  2026-09-09; pilot only (nsim ≈ 30). The full sweep... is a separate cluster/batch job."** This
  is the study backing the §3.6 GLM-misspecification comparison and part of `thm:regular`'s
  "Covered" status in `claims.md`. **§4 cannot responsibly report S5 until that cluster run
  happens** — this is compute, not writing, and is the long pole in the sequencing decision.

**Open, not yet decided:** whether §4 waits for S5's full-scale run, or drafts around S1/S2 first
and adds S5 once it exists, or reports S5 explicitly as pilot-scale with the caveat carried into
the text. Author has not yet chosen among these.

---

## 6. Maintenance contract

- **Structure changes mid-draft → stop drafting, update this file, then resume.**
- This file is current as of 2026-09-17. All four of the author's 2026-09-17 decisions (Option B
  restructuring; `ass:sparsity`; single-tree relocation; citations + simulations before further
  drafting) are reflected in §2 and §5 above. The next actual prose change should be either (a)
  the mechanical restructuring itself in `manuscript.tex` (moving/retitling existing text per the
  "Moved from" column in §2, no new content), or (b) the citation search / simulation-audit
  follow-through in §5, before the Introduction is drafted from scratch — drafting the roadmap
  paragraph before §5 resolves risks citing section labels or evidence that don't exist yet.
