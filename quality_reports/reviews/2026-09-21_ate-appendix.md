# Review log — ATE appendix (sketch) added, 2026-09-21

**Scope:** `inst/paper/manuscript.tex` — new `\section{Extension to the average treatment effect
(sketch)}\label{app:ate}` (inside `\appendix`, previously empty except for the bibliography); a
new pointer sentence in §2.1 right after `eq:att`; an update to the Discussion's existing
ATT-vs-ATE scope paragraph to reference the new appendix. Requested by the author explicitly:
"add ATE theory to the appendix, and state after introducing the estimand as the ATT that ATE
theory is available in the appendix." Scope was clarified with the author before drafting — three
options were presented (minimal sketch / full proved parallel treatment / bare identification
only); the author chose **minimal: estimand + identification + sketch**, explicitly not a full
proof. Commissioned `domain-reviewer` (session `ses_f3b7ba358ffeIY7ZRtqC77Ai4c`) and `proofreader`
(session `ses_f3b7ba2a8ffeGDnfLLyjvMTGBR`) audits against the first draft before finalizing.

## Why this needed care beyond ordinary prose editing

This is new derived mathematical content, not an editorial change — the appendix states the ATE's
identification result, its efficient influence function, and a remark deriving a bilinear bias
decomposition by direct calculation. Per the project's own escalation logic, new technical content
gets checked for correctness regardless of how "minimal" its scope is; a wrong formula in a
one-paragraph sketch is exactly as wrong as one in a ten-page proof.

## What the audit found — and the headline result

**No error in the mathematics.** `domain-reviewer` independently re-derived every displayed
result from scratch:
- The identification result $\theta_0^{\mathrm{ATE}}=\E[m_1(\bX)-\mu_0(\bX)]$ — confirmed correct
  (standard g-formula), and confirmed it genuinely does not require the propensity, only
  consistency + ignorability + arm-wise positivity.
- The efficient influence function `eq:score-ate` — confirmed to be the standard AIPW-form EIF for
  the ATE in the nonparametric model, with all four signs correct, and confirmed mean-zero at the
  truth by direct substitution.
- **The bilinear bias identity in `rem:bilinear-ate`** (the result this session's own drafting was
  least certain about, and the one whose sign errors would be easiest to make and hardest to
  catch by eye) — independently re-derived from the tower property, in full, by the reviewer, and
  confirmed **exactly correct**, including the control-arm term's sign, which does not flip
  (a genuine risk point: the naive "fix" of inserting a minus sign there would have introduced an
  error). The reviewer additionally found, and the manuscript now states, a non-obvious
  corroborating fact not previously known anywhere in this project: the ATE's control-arm bias
  term is *exactly* `prop:bilinear`'s own ATT remainder, up to the ATT's $\pi^{-1}$ normalization
  — an independent internal consistency check between brand-new content and a pre-existing,
  already-audited result.

## What was fixed (three MAJOR, several MINOR)

**MAJOR — a false comparative claim.** The first draft said "unlike the ATT's own identification
result \eqref{eq:att-id}, [the ATE's] does not involve the propensity at all" — but the ATT's own
`eq:att-id` is equally propensity-free (this is exactly `thm:regular`'s own point, that
differentiating it "needs only the control regression $\mu_0$"). Fixed to state the true contrast:
both identification results are propensity-free; the genuine difference is that the ATE needs a
*second* outcome regression ($m_1$, in addition to $\mu_0$).

**MAJOR — conflated identification-level and efficiency-level positivity.** The first draft
attributed the two-sided *uniform* positivity bound to *identification*, justified by an
*estimation*-side reason (denominators in the EIF). Identification of the ATE needs only pointwise
two-sided positivity ($0<e_0(\bx)<1$ a.s., so $m_1$ is defined); the uniform bound is what efficient
*estimation* needs, for exactly the denominator reason originally given. This distinction is one
the main text is already scrupulous about for the ATT's own one-sided bound (`ass:causal`'s
surrounding prose, `:137`) — the appendix's first draft fell below the paper's own established
standard on precisely the point the paper is careful about elsewhere. Fixed by splitting into two
sentences, one per level, each with its correct justification.

**MAJOR — a missing hypothesis on a displayed identity.** `rem:bilinear-ate`'s bias decomposition
has the candidate propensity $e$ in both denominators, but the first draft stated no bound on it
— without one, the two expectations need not be finite and the tower-property exchange used to
derive the identity is unjustified. `prop:bilinear` states the analogous hypothesis explicitly;
this appendix's version did not. Fixed by opening the remark with the explicit premise
$c'\le e(\bx)\le1-c'$.

**MINOR, all fixed:** "structural sparsity" mis-attributed to the wrong section (`sec:trees`,
which supplies the *sufficient class* definition, rather than `sec:instantiation1`, which is where
the rest of the paper points for the *condition* "structural sparsity" itself); no citation on the
EIF display (added `\citep{hahn1998}`, matching the main text's own practice for `eq:score`); the
candidate outcome-regression symbol $m$ silently overloading the pre-existing leaf-size floor
$m_n$ (added a disambiguating registry note; left the manuscript prose as-is since every use is
subscript-disambiguated in context); the constant $c$ reused for a strictly different range
without a new symbol (renamed the ATE's uniform bound constant to $c'$ throughout); the phrase
"therefore three trees" overclaiming that the tree count *follows from* the two displayed
equations, when it actually follows from the paper's established one-tree-per-nuisance convention
and leaves open (rather than forecloses) whether a single-tree ATE analogue could exist — softened
and connected explicitly to the open question the Discussion paragraph already raises; the new
label `eq:att-id-ate` (near-collision with `eq:att-id` in a document where the ATT/ATE distinction
is load-bearing) renamed to `eq:ate-id`; a display using $\eta^{\mathrm{ATE}}$ before defining it
(moved the definition earlier); minor prose polish (an awkward "does not complete them for a
different estimand," a mid-sentence parenthetical, three uses of "sketch"/"sketches" in one
Discussion paragraph reduced to two, a clunky possessive, "zeroes" → "zeroes out").

## Scope discipline — confirmed respected

`domain-reviewer` checked specifically whether the appendix, anywhere, implies more than was
authorized (a proven three-nuisance extension). Verdict: no. The appendix uses no numbered
proposition, theorem, lemma, or corollary — only two labeled equations and one remark, and every
forward-looking sentence about the unproven parts (selection consistency, the CLT, the anchor
interval for three nuisances) is explicitly hedged ("is not verified here," "if every nuisance
*were*... both terms *would* vanish" in the subjunctive). The updated Discussion paragraph
describes the appendix's actual contents accurately and, if anything, *under*-states them (it does
not mention the appendix also supplies the EIF) — the reviewer noted this is the safe direction for
any discrepancy to run.

## Registry deltas (same turn)

- `outline.md`: Appendix row updated from "empty... not yet drafted" to reflect the drafted
  content, its audit status, and what it explicitly does not claim.
- `claims.md`: ATT-vs-ATE scope row updated to point at the appendix instead of pure "future
  work"; three new rows added (`eq:ate-id`, `eq:score-ate`, `rem:bilinear-ate`), each recording
  what was verified and how.
- `notation.md`: new symbol-registry block for the appendix's ATE-local symbols
  ($\theta_0^{\mathrm{ATE}}$, $\eta^{\mathrm{ATE}}$, $\psi^{\mathrm{ATE}}$, $c'$, $\cS_m$),
  explicitly noting each is Appendix-local and does not overwrite the corresponding main-text
  symbol's meaning — verified by `domain-reviewer` via exhaustive grep for any bare (unsuperscripted)
  occurrence of the ATT-reserved symbols inside the new appendix text; none found.

## Verification

Label/`\ref`/`\eqref` integrity re-checked after every edit in this batch: 49 unique labels, 38
unique refs, zero missing, zero duplicates. No LaTeX compile run this session.

## Explicitly not done / not claimed

- No selection-consistency result, central limit theorem, or anchor interval for the ATE's
  three-nuisance case — by design, per the author's own scope decision, and stated as such in
  the manuscript itself, not just in this log.
- Whether a propensity-free, single-outcome-tree analogue exists for the ATE (paralleling
  `cor:single-tree`'s ATT-only result) is explicitly left open, not answered.
- LaTeX compilation itself was not verified this session, for any of today's three edit passes
  (Discussion completion, §2.1 audit, this ATE appendix) — needs a `compile-latex`/`verifier` run
  before any is treated as fully final.
- The `ass:sparsity` gap from the earlier session-status review, and the broader lowercase
  $x$/$\bx$ notation drift flagged during the §2.1 audit, remain open and unrelated to this log.
