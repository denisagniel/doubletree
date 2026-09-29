MODE: independent subagent audit (proof-auditor) — context isolated, drafted by a separate
oracle-consultation pass earlier in the same session.

# Proof Audit: `inst/paper/theory.tex` — new `sec:onesided` ("One-sided sufficiency")
**Date:** 2026-09-28
**Reviewer:** proof-auditor agent (model: Claude Opus 5)
**Target:** new material added this session — macros `\etatil,\wtil,\mtil,\cF`; subsection
`sec:onesided` with `ass:pseudo-onesided`, `ass:clip-margin`, `lem:onesided`,
`cor:onesided-mu`, `cor:onesided-e`, and eight remarks; the corresponding proofs appended to
`app:main`.

## Coverage Manifest

Read directly: preamble macros (1–120); setup/score/assumptions A1–A8 (143–402);
`ass:budget`…`prop:bilinear`, `lem:tangent`/`eq:vdecomp` (556–875); `lem:margin`,
`lem:selection` statements (946–1039); `sec:main` in full (`lem:skeleton`, `lem:U`,
`lem:cross`, `lem:uniform`, `thm:main`, `cor:variance`, `eq:vhat`, `lem:findim`, 1274–1463);
the new `sec:onesided` in full (1464–1827); `prop:sparse-bound` (1880–1939); `lem:clip`,
`lem:normequiv` and proofs (2961–3035); `lem:margin`/`lem:selection` proofs (3405–3593);
`lem:skeleton`/`lem:findim`/`lem:U`/`lem:cross` proofs (3941–4219); the new proofs of
`lem:onesided` and both corollaries in full (4220–4552).

Not read: `manuscript.tex`; `lem:invariance`'s proof; `sec:condP`/`sec:crossfit`/`sec:honest`
beyond the specific lines the new material cites (`thm:anchor`, `prop:spectest`).

## Verdict: CONDITIONAL PASS

**All mathematics in the new material is correct.** Independent re-derivation confirmed both
exact variance formulas (Corollaries `cor:onesided-mu`/`cor:onesided-e`), the `V_*-V`
difference formula, the `rem:onesided-bound` minimization and its equality condition, and every
`O_p(\cdot)` rate claim in the proof of `lem:onesided` (nine order claims recomputed from
scratch, not checked for plausibility). No sign error, no circular reasoning, no silent
weakening of any claim.

Findings: 0 critical, 5 major, 17 minor, 6 suggestions. The five major findings are
presentation/consistency issues, not mathematical errors:

- **M1** — `V_*` is reused for a different quantity than the one `ass:pseudo` already defines
  it to mean (asymptotic variance at the purged nuisance `\etatil`, vs. `ass:pseudo`'s variance
  at the pseudo-true `\eta_*`); the two coincide only in the collapse case, which masks the
  clash. Fix: rename the new quantity.
- **M2** — `V_*>0` is used by both corollaries' coverage claims but never stated as a
  hypothesis or established. Fix: add as a hypothesis (`cor:onesided-mu`) or prove it from
  standing assumptions (`cor:onesided-e`, where a two-line argument suffices).
- **M3** — the mean-zero property in case (i) is attributed to `prop:bilinear`, but that
  proposition's own hypothesis (`1-e\ge c`) need not hold for the purged `\wtil`; the mean-zero
  fact is in fact established directly two lines later and does not need the citation. Fix:
  rewrite the citation to point at the direct verification.
- **M4** — two false leaf-constancy claims inside the proof (`b` and `a` asserted leaf-constant
  on the pseudo-true side's own partition, when by construction they are not); not load-bearing
  anywhere downstream. Fix: delete the two false clauses.
- **M5** — the two new assumptions are declared inside `sec:onesided` rather than `sec:setup`,
  conflicting with the document's own stated convention that assumptions live only in
  `sec:setup`. Pre-existing violation (`ass:rate` already breaks this), aggravated here. Fix:
  relocate, or amend the convention statement to acknowledge locally-declared assumptions.

Full findings (including all 17 minor items and 6 suggestions, each with exact line numbers,
independent recomputation, and a specific recommended fix) are in the proof-auditor's full
report, retained in this session's transcript rather than duplicated here.

## Disposition

All five major findings and the cheap minor one-liners (m1, m4, m8, m9, m13, m14, m15, m16,
m17) were applied to `theory.tex` in this pass, in the auditor's suggested cheapest-first order.
M5 was resolved by amending the `sec:setup` convention statement to acknowledge the two
pre-existing exceptions (`ass:rate`, the rectangular-covariate assumption) and the two new ones,
rather than relocating four assumptions that are used only locally. The remaining minor items
(m2, m3, m5, m6, m7, m10, m11, m12) and suggestions (s1–s6) are lower-value notation/cross-
reference polish; not applied in this pass, left for a future proofreading sweep.

## Registries updated

`notation.md`, `claims.md`, `outline.md` — see this session's other changes to those files.

## Follow-up audit (same day, second pass)

A second `proof-auditor` pass verified the fix pass above, specifically because a bulk
line-range `sed` substitution (the `V_*` rename) and several hand-edits are exactly the kind of
change that can fix one defect while introducing another. Verdict: **CONDITIONAL PASS** — 0
critical, 4 major, 5 minor findings, again **zero new mathematical errors**. Two of the four
majors were genuinely new defects introduced by the first pass's own fixes:

- **M3 was only partially applied.** The statement-side prose was corrected, but the proof of
  `lem:onesided` still invoked `prop:bilinear` outside its stated hypothesis in case (i) (the
  defect had been relocated, not removed). Fixed: the proof's citation is now restricted to
  case (ii), where it is legitimate, with an explicit note on why it does not extend to case (i).
- **M5's amendment introduced a false cross-reference and an unmet convention.** The rewritten
  scope statement pointed to "the rectangular-covariate assumption of Section~\ref{sec:condP}"
  when that assumption (`ass:bridge-mesh`) is actually in `sec:discussion`; and it claimed every
  locally-declared assumption is "flagged as such at the point of definition" when only one of
  the four actually was. Fixed: corrected the section reference and cited `ass:bridge-mesh` by
  label; added the missing locality sentence (matching `ass:bridge-mesh`'s own wording) to
  `ass:pseudo-onesided`, `ass:clip-margin`, and `ass:rate`; softened the scope statement's claim
  to match what is now actually true of all four; and updated the standing-assumptions taxonomy
  paragraph, which had gone stale in the same way (an outdated `\ref` range, and a now-false
  claim that `ass:pseudo` is "never in force elsewhere" given `rem:onesided-random-index`
  invokes it as a sufficient condition for part (c)).

Also fixed: `cor:onesided-e`'s $\Vtil>0$ derivation (mathematically correct in the first pass,
but referenced an undisplayed "fourth term"/"second term" and skipped a required step) — now
displays the four-term variance expansion explicitly and walks through all four terms. Three
cheap minors (an imprecise "preserve the sup-norm bounds on $w_*$" that should reference $\ez$;
a one-sided `$e_*\le1-c$` left inconsistent with a two-sided bound introduced elsewhere; an
indexing-base wording issue) were also fixed. Recompiled clean after each pass (0 undefined
references, 0 multiply-defined labels, 64 pages, same 11 pre-existing overfull hboxes).

Not re-litigated, per the second pass's own instruction: m2, m3, m5, m6, m7, m10, m11, m12, m17,
s1–s6 from the first pass. The second pass additionally proposed `S1` (weakening
`prop:bilinear`'s hypothesis to arbitrary bounded $(w,\mu)$, which its own proof already
supports and would make the case-(i)/case-(ii) asymmetry moot) as the highest-leverage remaining
edit — not applied this pass; left for a future revision since it touches a proposition with
other consumers (`lem:biasbound`, `cor:convert`, `thm:crossfit`) that would need re-checking.
</content>
