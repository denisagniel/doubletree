# Review log — drafting `ass:sparsity` (§2.2), 2026-09-21

**Scope:** `inst/paper/manuscript.tex` §2.2 (`sec:trees`) — new `\begin{assumption}[Structural
sparsity]\label{ass:sparsity}` placed after `def:sufficient`, plus a new Remark distinguishing
legibility from a complete audit trail, plus six existing sites elsewhere in the manuscript that
restate or invoke $\cS_e\ne\emptyset,\cS_\mu\ne\emptyset$ in prose. This closes the gap the
2026-09-21 09:05 session found: `outline.md` had marked this restructuring row "accepted" and
even "EXECUTED" (elsewhere in the file), but `ass:sparsity` did not exist as a label anywhere in
`manuscript.tex` before this pass.

Audit commissioned **after** first drafting (not before, since the content did not yet exist to
audit): `domain-reviewer` (session `ses_f3ab34850ffeYoRX4KDn8O1oxk`) and `proofreader` (session
`ses_f3ab34745ffegDbFjez3SdWyLl`), run in parallel against the first draft. Fixes below were
applied after both reports returned, then re-verified by a fresh 4-pass XeLaTeX+bibtex compile
(zero undefined refs/citations, zero multiply-defined labels, stable 27 pages).

## What was drafted

1. `ass:sparsity`: $\cS_e\ne\emptyset$ and $\cS_\mu\ne\emptyset$, immediately after
   `def:sufficient`, matching `theory.tex`'s own `ass:sparsity` (the formal companion document
   already had this exact condition named; the manuscript did not).
2. A Remark ("Legibility does not by itself make the audit trail complete") distinguishing
   legibility (delivered by the fixed grid/leaf budget, independent of whether the assumption
   holds) from a complete audit trail (requires the assumption) — the "legibility-vs-audit-
   completeness remark" `claims.md` Table 3 had marked accepted but, like the assumption itself,
   never actually drafted.
3. Four of the six prose restatements the outline's plan targeted were converted to cite
   `Assumption~\ref{ass:sparsity}` by label: the opening of `sec:instantiation1`, `lem:selection`'s
   hypothesis, `prop:P-instantiation1`'s hypothesis, and `prop:pseudo-consistency`'s exclusion
   clause.
4. **The other two were deliberately left as direct prose, not converted** — `lem:single-tree-
   linear` (needs only $\cS_\mu\ne\emptyset$) and `cor:single-tree-coarsening` (adds
   $\cS_e\ne\emptyset$ individually to that). Converting either to the joint-assumption label would
   have overstated what those specific results actually assume. This judgment call is exactly what
   the audit was commissioned to stress-test — see below.

## What the audit found, and what was fixed

**BLOCKING, all fixed:**

- **`lem:selection`'s hypothesis was stated jointly but a later lemma (`lem:single-tree-linear`)
  legitimately invokes it one-sidedly.** `domain-reviewer` traced the actual optimization
  (`eq:select`/`eq:feasible`) and confirmed selection is separable per nuisance — disjoint
  objectives, disjoint feasible sets, no coupling term — so the joint hypothesis was stronger than
  the mechanism requires. Independently re-verified against `eq:select`/`eq:feasible` and
  `prop:selection-rate`'s own per-$j$ $\Delta_j$ before accepting. **Fixed:** restated `lem:selection`
  per $j$ ("if $\cS_j\ne\emptyset$, then $\Prob(\tauhat_j\in\cS_j)\to1$"), with a joint corollary
  sentence for consumers that need both nuisances at once. This is what actually makes leaving
  `lem:single-tree-linear` unconverted (item 4 above) a *correct* citation rather than a citation of
  a proposition whose stated hypothesis it does not satisfy.
- **Pre-existing error, newly exposed: "`ass:pseudo` is strictly weaker than structural sparsity"
  is false.** Independently re-derived against `eq:risk-pop`: for every $\tau\in\cS_\mu$,
  $R^{(\mu)}(\tau)=\E[(Y-\mu_0(\bX))^2\mid A=0]$ — the *same* value — so a non-singleton $\cS_\mu$
  (which the manuscript repeatedly says is the generic case, e.g. lines 330, 365) is exactly an
  exact tie, and `ass:pseudo`'s required unique minimizer fails to exist. Sparsity does not imply
  `ass:pseudo`; they are non-nested. **Fixed:** rewrote the paragraph to state non-nestedness and
  the tie mechanism directly, dropping the now-unsupported "finest partition in $\cS_j$" claim.
- **New Remark overclaimed twice.** "Legibility is therefore free" (unqualified) is contradicted
  verbatim by §3.4's own "It can cost something in representation, however" and by
  `rem:saturated-general`'s $2^k$-leaf cost. "The best tree-representable approximation" presupposes
  `ass:pseudo` (existence *and* uniqueness), which the Remark never invokes, and which
  `prop:pseudo-consistency`'s own surrounding prose says can fail ("the selected partitions can
  oscillate with no fixed limit at all"). **Fixed:** scoped both claims — legibility is free *of
  Assumption~\ref{ass:sparsity}*, not free of representational cost; the displayed trees are *a*
  tree-representable approximation, population-risk-optimal only in the limit and only when
  `ass:pseudo` supplies a unique one.

**MINOR, all fixed:**

- `cor:single-tree-coarsening`'s hypothesis ($\cS_\mu\ne\emptyset$ plus $\cS_e\ne\emptyset$
  individually) is **not** actually one-sided — combined, it *is* the joint condition, and
  `prop:P-instantiation1` (cited two lines later) needs exactly that. Fixed the parenthetical from
  "structural sparsity holds for the propensity too" (ill-formed once "structural sparsity" is a
  named joint condition) to "jointly, Assumption~\ref{ass:sparsity}." **This is the asymmetry in
  the original (d) judgment call:** `lem:single-tree-linear` (item 4, first case) is genuinely
  one-sided and was right to leave alone; `cor:single-tree-coarsening` (item 4, second case) is not
  one-sided, and only its wording — not its citation status — needed fixing.
- Three remaining locators for "structural sparsity" still pointed at `\S\ref{sec:instantiation1}`
  (its old, unlabeled home) rather than the new `Assumption~\ref{ass:sparsity}`; fixed all three
  (§3.4's "what does restrict the DGP" sentence, and both mentions in the new ATE appendix).
  `prop:P-instantiation1`'s hypothesis line named the assumption twice (once transitively via
  "the conditions of Proposition~\ref{lem:uniform}," once explicitly); reworded as an
  "in particular" clause rather than deleting the emphasis, since `lem:uniform`'s own conclusion is
  vacuous rather than false when sparsity fails.
- `proofreader` caught: a stray "topology recovery" (the manuscript's own term, used everywhere
  else, is "partition recovery" — fixed); a missing `\label` on the new Remark, the only titled
  remark in the document without one (added `rem:legibility-not-complete`); a singular/plural drift
  ("the displayed tree ... the true nuisance" mid-paragraph against "the displayed trees ... they"
  one sentence later — fixed to plural throughout, which also resolved a number-agreement slip
  domain-reviewer independently flagged in the same sentence); a self-referencing `\S\ref{sec:trees}`
  from inside `sec:trees` (removed); a stray `\noindent`, the only occurrence in the 858-line file
  (removed); a duplicated forward pointer to `sec:honest-manuscript` twelve lines apart (trimmed to
  one).
- **Terminology drift, resolved by checking actual usage rather than either reviewer's
  first-pass guess.** `domain-reviewer` suggested standardizing on "audit trail" by changing the
  *abstract*; `proofreader` suggested conforming the new Remark to "auditable trail." Checked
  actual usage first: "audit trail" appears 3 times pre-existing in §5 Discussion versus one
  "auditable trail" in the abstract — the new Remark now uses "audit trail," matching the dominant
  existing form, and the abstract was left untouched (out of scope for this pass; already carries
  both forms in the same paragraph, a pre-existing minor inconsistency, not one introduced today).

## Deliberately not fixed (out of scope, flagged for a future pass)

- `cor:single-tree-coarsening`'s heteroskedasticity worked example (line ~747, "$\approx0.30$")
  does not reproduce under either natural reading of $\sigma^2_{0,\ell}$ that `domain-reviewer`
  tried (got $\approx0.18$ or $\approx0.37$, not $0.30$) — pre-existing, not part of today's
  content, needs its own numeric check.
- `\bar L` vs. `\Lbar` and `\that` vs. `\tauhat` are both used interchangeably throughout the
  document (they render identically via `common-defs.tex`'s macro, so this is invisible to
  readers, only to source-level grep/redefinition) — pre-existing, paper-wide, not introduced by
  this pass.
- `proofreader` also flagged a missing `~` before `\ref` at (pre-existing) line ~171 and a
  double-comma at the newly-fixed `lem:selection` hypothesis line — both noted as low-value,
  pre-existing/cosmetic in their own report; left for a future proofreading pass rather than
  folded into a substance-fix commit.

## Registries updated

`claims.md` (Table 1 row for `ass:sparsity`, Table 3 row for the legibility remark),
`outline.md` (§2.2 row), `notation.md` (moved the "Structural sparsity" entry from §2b
"implicit" to a new §2a "formal" entry; updated the `$\cS_j$` symbol-table row) — all now say
**Drafted**, with the corrected count of converted vs. deliberately-unconverted sites, superseding
this morning's "six locations, one collapse" framing.
