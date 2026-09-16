# Paper Outline — `doubletree`

**Instance path:** `inst/paper/outline.md`
**Governed by:** `paper-protocol.md` (via the `draft-paper-section` / `edit-paper-with-context` skills)
**Status:** draft, created 2026-09-16. This file did not exist before today — created now, at the
author's request, to settle structure before the next manuscript edit rather than after it.

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

## 2. Section plan

One row per section/subsection currently in `manuscript.tex`. "Owns (pending)" marks what the
finite-$X$/legibility/coarsening-bias edit is proposed to add; nothing in that column is written
into the manuscript yet.

| § | Working title (label) | Owns (as of current manuscript) | Owns (pending, this edit) | Claims made (→ `claims.md`) | Results delivered | Depends on |
|---|---|---|---|---|---|---|
| 1 | Introduction | The audit-trail problem; the gap (existing methods can't jointly deliver flexible + interpretable + full-sample-valid nuisance estimation); the contribution; roadmap. | **The three-rung regime ladder** (audit-complete / sparse / general) motivating *why* the paper cares about $X$'s finiteness, stated as motivation, not as a new assumption. Cross-referenced from §`sec:trees`, not restated there. | Abstract-level claims (not separately tabled) | — | — |
| 2.1 | Data, estimand, and notation | $O=(X,A,Y)$, potential outcomes, $\theta_0$, nuisances $e_0,\mu_0$, weight $w_0$, $\Pn$, the two norms. | — (no change) | — | — | §1 |
| 2.2 | Covariate space, tree partitions, and optimal trees (`sec:trees`) | **Existing, must not be duplicated:** "$\cX$ need not be finite" (`:79`); the fixed, data-independent cutpoint grid; grid atoms $\cA$; leaf budget $\bar L$; candidate class $\cT_{\bar L}$; the sufficient class `def:sufficient`. | **New remark distinguishing legibility from audit-completeness/fidelity**: legibility (reviewer can trace a unit to a leaf and read its cutpoints) holds under any $\cX$ and is already delivered by the fixed grid; audit-completeness (the displayed leaf number *is* the unit's true nuisance value, nothing hidden) additionally requires atom-measurable nuisances — i.e., is exactly `ass:sparsity`/the sufficient-class condition restated in a motivating register, not a new assumption. One forward cross-reference to the regime ladder in §1. **Must not re-assert or contradict `:79`.** | — | — | §2.1 |
| 2.3 | Identification and influence function | **Existing, adjacent but distinct — do not conflate:** `ass:causal` (identification + *uniform* one-sided positivity $1-e_0\ge c$); `:128`'s own finite-vs-general-$\cX$ discussion — but that passage is about the uniform-vs-pointwise positivity bound's variance-decomposition term, a *different* finite/general-$\cX$ distinction than the sparsity/legibility one this edit adds. EIF `eq:score`. | None planned. If the pending §3.3 remark (below) needs to say "identification is unaffected by coarsening," it should cite `ass:causal` from here rather than restate it. | `ass:causal` | — | §2.1, §2.2 |
| 2.4 | The doubletree algorithm | Feasible sets `eq:feasible`; penalized-ERM selection `eq:select`; leaf refits `eq:refit`; no-splitting remark; `rem:two-trees-why`; `ass:construct` (clipping). | None planned. | `ass:construct`; `rem:two-trees-why` | — | §2.2, §2.3 |
| 3.1 | Partition recovery and the CLT (`sec:partition-recovery`) | Efficiency rests on partition recovery, not grid properties; Condition (P), deferred to the companion theory document; notes the analyst's grid is finite by construction so Condition (P) is not needed directly. | None planned. | — | — | §2.2–2.4 |
| 3.2 | Instantiation 1: fixed grid, exact sparsity (`sec:instantiation1`) | Exact structural sparsity setting; `lem:selection`, `prop:selection-rate`, `cor:consistency`, `lem:uniform`, `prop:P-instantiation1` (headline CLT), `cor:variance`; the single-tree corollary family (`cor:single-tree`, `lem:single-tree-linear`, `cor:single-tree-coarsening`, `rem:single-tree-se`, `rem:single-tree-distinct`). | None planned. **Note (pre-existing gap, not from this edit):** per `claims.md`, the single-tree corollary family has zero simulation coverage — flagging so the pending edit doesn't get blamed for it later. | `lem:selection` … `rem:single-tree-distinct` (full list in `claims.md` Table 1) | Prop. `P-instantiation1` (headline CLT); Cor. `variance` | §3.1 |
| 3.3 | When partition recovery fails: the anchor interval (`sec:honest-manuscript`) | Non-sparse setting; `prop:bilinear` (exact bilinear remainder); `lem:biasbound`; the anchor construction; `ass:rate`; `thm:anchor` (headline validity); `cor:width`; `prop:spectest` (fidelity diagnostic). | **New remark, placed after `cor:width` (~line 598):** discretizing a continuous covariate onto the grid is a *special case* of sparsity failing that `thm:anchor` already covers, with a terminology disclaimer (avoid "residual confounding"; see `notation.md` §3). **Optional, pending decision (see open question below):** a new proposition bounding the coarsening-specific piece of $\delta_e,\delta_\mu$ under a Lipschitz sensitivity device on $e_0,\mu_0$ and grid mesh $h$, combined with a finite-class ERM step connecting $\delta$ to the realized $D_w,D_\mu$ that `cor:width` actually uses — yielding width $=O_p(\max\{n^{-1/2},h^2\})$. | `prop:bilinear`, `lem:biasbound`, `ass:rate`, `thm:anchor`, `cor:width`, `prop:spectest` | Thm. `anchor` (headline validity); **new, if adopted:** Prop. `coarsen-bound` (tentative label) | §2.2, §2.3, §3.2 |
| 3.4 | Discussion: is a fixed candidate class basically a parametric model? (`sec:finite-discussion`) | `thm:regular` (regularity/efficiency in the full nonparametric model); `prop:sparse-bound` (submodel bound not exploited); argues `ass:sparsity`'s cost is about the DGP, not the model; explicitly states "we do not minimize this cost in general" (`:666–671`). | **If, and only if, the optional proposition in 3.3 is adopted:** point `:666–671`'s existing "we do not minimize this cost in general" sentence at that new proposition, rather than leaving it as an unaddressed promise. No change if the optional piece is not adopted. | `thm:regular`, `prop:sparse-bound` | Thm. `regular`; Prop. `sparse-bound` | §3.2, §3.3 |
| 3.5 | Relation to penalized linear methods (`sec:lasso-comparison`) | Two-tier resolution shape shared with Lasso/debiased-inference; different mechanism (no orthogonality property of the displayed estimator); representational trade-offs (hierarchical/region-varying vs. additive structure). | None planned. | (prose, `sec:lasso-comparison` row in `claims.md`) | — | §3.3, §3.4 |
| 4 | Discussion (`sec:discussion`) | Returns to payment/triage/sentencing; states the audit trail is available regardless of sparsity; explicit limitations (no real/simulated demonstration; pointwise not uniform; not matching black-box flexibility); future work (finite-sample behavior, "genuine misspecification," applied comparison, `:701`). | **If the optional proposition in 3.3 is adopted:** narrow — not delete — the "genuine misspecification... left to future work" clause at `:701`, since that proposition would be a first cut at exactly this. If not adopted, no change; the promise stands as-is. | (prose, `sec:discussion` row in `claims.md`) | — | §1, §3.3, §3.4 |
| Appendix | — | Empty in `manuscript.tex` itself (`\appendix` at `:703` precedes only the bibliography). Proof content lives in the separate companion document `theory.tex`, out of scope for this backfill. | None planned. | — | — | — |

**Open question this table surfaces, not yet resolved:** whether to adopt the optional
Lipschitz-sensitivity proposition in row 3.3 (and its downstream touches in rows 3.4 and 4).
This is the one item in the plan with real technical work (a finite-class ERM excess-risk step);
everything else is prose. See `claims.md`'s new Table 3 for the claims-registry consequence of
each choice.

---

## 3. Dependency and ordering check

- [x] No section depends on a result introduced later, as far as this backfill can tell —
      the pending edit's new material (row 1, 2.2, 3.3) only depends on rows already earlier in
      the table (§2.2, §2.3, `thm:anchor`, `prop:bilinear`), and its one internal forward
      dependency (3.4/4 pointing at the optional 3.3 proposition) is recorded conditionally above,
      not silently assumed.
- [ ] Every assumption is introduced before the first result that uses it — **not yet checked
      against `notation.md`'s "Used by" field**, since `notation.md` did not exist before this
      pass either; see that file's own scope note.
- [ ] Range-citation and appendix-numbering checks (`Assumptions~\ref{first}-\ref{last}`,
      appendix-declared-but-body-used assumptions) — **not checked in this pass.** `manuscript.tex`
      has an empty appendix, so the appendix-crossing risk is likely moot, but this was not
      verified line-by-line.
- [x] The discussion (§4) already contains results from substantive sections (`prop:spectest`,
      the audit-trail claim tied to `thm:anchor`/`cor:width`), per the existing text at `:699`.
- [x] `\appendix` is used correctly (`:703`) and the appendix section is genuinely empty rather
      than competing with a body `\section{Appendix}`.

---

## 4. Explicitly not planned here

- Whether to include a roadmap paragraph in the introduction — already present (`:44`); not
  revisited by this edit.
- Whether to have a related-work section — the paper does not have one; not revisited.
- Whether every assumption carries a `[Name]` — mixed already (`ass:causal` named, `eq:score`
  unnamed as an assumption but is not one); not revisited except for the new sensitivity device
  in 3.3, if adopted, which per Oracle's review must **not** be given the weight of a standing
  named assumption.
- A full backfill of `theory.tex`'s own section plan — flagged as future work, not blocking.

---

## 5. Maintenance contract

- **Structure changes mid-draft → stop drafting, update this file, then resume.**
- This file is now current as of 2026-09-16. The next person to draft against it should update
  the "Owns (pending)" column to "Owns (as of current manuscript)" once the corresponding prose
  actually lands, and delete the "pending" framing at that point.
