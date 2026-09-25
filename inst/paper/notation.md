# Notation and Assumption Registry — `doubletree`

**Instance path:** `inst/paper/notation.md`
**Governed by:** `paper-protocol.md`
**Status:** full backfill of `manuscript.tex`, current as of 2026-09-16.

> **Scope.** Full backfill of `manuscript.tex` only, matching `claims.md`'s own precedent
> ("`theory.tex` is the companion proof document for the same claims... not a second,
> independent claims surface"). `theory.tex` (4383 lines; Condition (P), `ass:margin-floor`,
> `ass:pseudo`, `ass:solver`, `ass:finite`, `ass:budget`, `ass:global`, `ass:moments`,
> `ass:iid`, `ass:sparsity`) is **out of scope**, confirmed with the author 2026-09-16. Where
> `manuscript.tex` cites a `theory.tex` object without restating it (e.g. Condition (P),
> $\delta_e,\delta_\mu$), that is recorded below as "deferred to companion document," not
> backfilled from `theory.tex` itself.

---

## 1. Symbol registry

One row per symbol appearing in a formal statement in `manuscript.tex`, in the order it is
first introduced.

| Symbol | Meaning | Defined at | Notes / forbidden variants |
|---|---|---|---|
| $\theta_0$ | Estimand — the ATT, $\E[Y(1)-Y(0)\mid A=1]$ | `:58` (eq. `att`) | Never write bare $\theta$ for the estimand; bare $\theta$ is the generic candidate value inside $\psi(O;\theta,\eta)$. **2026-09-21:** the identified (observed-data) form is now labeled `eq:att-id` (`:136`) rather than left un-`\ref`-able; `thm:regular` (`:412`) cites it by label instead of restating it. |
| $\bO=(\bX,A,Y)\sim\Prob$ | One observation; $\bO_1,\dots,\bO_n$ iid realizations | `:53` | iid sampling itself is never pulled into a numbered assumption in this document — see §2b. **2026-09-21:** normalized to bold ($\bO$, not $O$) to match the manuscript's dominant convention (`domain-reviewer` finding, `quality_reports/reviews/2026-09-21_section2.1-audit.md`); `eq:score` (`:141`) and `ass:causal`'s ignorability clause (`:127`) were fixed to match in the same pass. |
| $\bX\in\cX\subset\R^p$ | Covariate vector | `:53` | |
| $A\in\{0,1\}$ | Binary treatment | `:53` | |
| $Y\in\R$ | Outcome | `:53` | **2026-09-21:** the standalone symbol $\cY$ (previously "$Y\in\cY\subset\R$") was dropped — confirmed zero downstream uses anywhere in the manuscript, so it was dead weight rather than a symbol worth defining. Boundedness of $Y$ remains a separate, bundled condition stated locally where needed (`lem:selection`'s moment condition, `prop:selection-rate`'s stronger boundedness), not a property of an outcome-space symbol. |
| $Y(a)$ | Potential outcome under $A=a$ | `:53` | |
| $e_0(\bx)=\Prob(A=1\mid\bX=\bx)$ | Propensity nuisance | `:53` | |
| $\mu_0(\bx)=\E[Y\mid A=0,\bX=\bx]$ | Control-outcome nuisance | `:53` | |
| $\pi=\Prob(A=1)$ | Treatment prevalence | `:53` | $\pi>0$ is part of `ass:causal`. |
| $m_1(\bx)=\E[Y\mid A=1,\bX=\bx]$ | Treated-outcome regression | `:53` | **2026-09-21, corrected:** the prior note ("used only in prose") was itself wrong — a full-document grep found zero uses anywhere after its definition. Fixed by actually using the symbol at its natural site, `thm:regular` (`:423`, "never $m_1(\bx)$, which would be undefined wherever $e_0(\bx)=0$"), which previously wrote $\E[Y\mid A=1,\bX=\bx]$ longhand. Still not used in any displayed formula among the manuscript's main results, which work with $\mu_0$ alone. |
| $\sigma_a^2(\bx)=\Var(Y\mid A=a,\bX=\bx)$ | Conditional outcome variance | `:53` | $\sigma_0^2(\bx)$ (the $a=0$ case) recurs in the single-tree corollary family (`:451`, `:455`) under the name $\sigma_{0,\ell}^2$ for its leaf-average — see below; do not conflate the pointwise and leaf-average versions. |
| $w_0(\bx)=e_0(\bx)/\{1-e_0(\bx)\}$ | Odds-scale propensity, true | `:61` | |
| $w(\bx)=e(\bx)/\{1-e(\bx)\}$ | Odds-scale propensity, candidate $e$ | `:61` | |
| $\eta_0=(e_0,\mu_0)$, $\eta=(e,\mu)$ | True / candidate nuisance pair | `:62` | |
| $j\in\{e,\mu\}$; $\gamma_{0,e}=e_0,\ \gamma_{0,\mu}=\mu_0$ | Index over the two nuisances | `:62–63` | |
| $\Prob_n$, $\Prob_nf=n^{-1}\sum_{i\le n}f(\bO_i)$ | Empirical measure | `:65` | |
| $\hat\pi=\Prob_nA$ | Empirical treatment prevalence | `:65` | The event $\{\hat\pi>0\}$, whose probability $\to1$, is the qualifier under which `thm:anchor` and `cor:variance` hold (`:369`, `:573`) — not a separate standing assumption, a probability-one event. |
| $\twonorm{f}^2=\E[f(\bX)^2]$ | Unweighted $L^2(P_X)$ norm | `:67–68` | |
| $\|f\|_{2,w}^2=\E[\{1-e_0(\bX)\}f(\bX)^2]$ | Control-weighted norm | `:67–70` | The "$w$" subscript here means *control-weighted*, unrelated to the propensity weight $w_0,w$ above — same letter, two roles. Flagging per the template's own instruction to record symbol collisions rather than let them stay latent. |
| $\rightsquigarrow$ | Convergence in distribution | `:72` | |
| $\to_p$ | Convergence in probability | `:73` (added 2026-09-21) | Used 6× downstream (`:368,:387,:390,:459,:547,:611`) — more often than $\rightsquigarrow$'s 5×, yet was previously undefined anywhere in the manuscript despite $\rightsquigarrow$ being defined in the same sentence. Gap found by `domain-reviewer`, fixed by adding three words to §2.1's existing sentence rather than a new sentence. |
| $z_{1-\alpha/2}$ | Standard normal quantile | `:73` | |
| $\cX$ | Covariate space. **Need not be finite.** | `:53`, elaborated `:82` | Do not write "$\cX$ is finite" unqualified anywhere — contradicts `:82` directly. |
| $\cA$ | Finite family of **grid atoms** from the analyst's pre-specified cutpoints; $M:=\card\cA<\infty$ | `:82–90` | The thing that is actually required finite — not $\cX$. Do not conflate. |
| $\tau$ | A tree partition of $\cX$, built from pre-specified cutpoints only | `:93–96` | |
| $\card\tau$ | Leaf count of $\tau$ | `:95` | |
| $\ell_\tau(\bx)$ | The leaf of $\tau$ containing $\bx$ | `:96` | **Corrected 2026-09-25** — was $\ell_\tau(x)$ (unbolded); a document-wide pass bolded every bare covariate symbol to $\bX$/$\bx$ since $\bX\in\cX\subset\R^p$ is a vector. Zero bare instances remain anywhere in `manuscript.tex` (verified by an independent `domain-reviewer` audit, `quality_reports/reviews/2026-09-25_sec2.2-saturated-general-broadening.md`). |
| $|\ell|,\ |\ell|_0$ | Leaf sample size; control ($A=0$) sample size within leaf $\ell$ | `:96` | $|\ell|_1:=|\ell|-|\ell|_0$ (treated count) introduced later. |
| $\Pi_\tau^\nu$, $\Pi_\tau:=\Pi_\tau^{\Prob}$ | $\nu$-weighted projection onto leaf-constant functions; unweighted case is the ordinary projection under $\bX$'s own law $\Prob$ | `:96` | **Corrected 2026-09-25**, replacing the row below. Previously written $\Pi_\tau^{p}$ with an atom-mass factor $p_x$ that was never given a reference measure — an independent `domain-reviewer` audit found this could make `:232`'s "these projections are exactly the squared-error minimizers" false under the literal (pointwise-weight) reading. Aligned with `theory.tex:206–207`'s own $\Pi_\tau=\Pi_\tau^P$ convention: the author confirmed $\nu$ is a **measure**, not a pointwise weight combined with an atom-mass density, and that this manuscript should match `theory.tex` throughout. |
| $d\nu=\{1-e_0(\bX)\}\,d\Prob$ | Control-weighted measure used for the outcome nuisance's projection | `:96` | **Corrected 2026-09-25**, replacing the row below (was $\nu_x=p_x\{1-e_0(x)\}$, a pointwise weight with the same unresolved-reference-measure problem). Now matches `theory.tex:207`'s $d\nu=\{1-\ez\}dP$ exactly, no atom-mass factor. |
| ~~$P(a)$, written $p_x$~~ | ~~Probability mass of an atom / covariate value~~ | — | **Removed 2026-09-25.** The symbol $p_x$ no longer appears anywhere in `manuscript.tex` — it was defined solely to serve the now-removed pointwise-weight reading of $\nu$ (see the two corrected rows above) and had no other use in the document. `theory.tex:186–187` reserves this exact shorthand for the special case where $\cX$ itself is finite (atoms coincide with points); the manuscript had inherited it into the general-$\cX$ setting, which is what an independent `domain-reviewer` audit flagged. Row kept, struck through, rather than deleted outright, so a reader of this registry's history can see the symbol existed and why it was removed — do not reintroduce $p_x$ without re-deriving which reference measure it would need. |
| $\bar L$ | Leaf budget (analyst-chosen, before data) | `:96` | Saturated budget $\bar L=M$. |
| $\cT_{\bar L}$ | Candidate class of tree partitions, $\card\tau\le\bar L$ | `:96–101` | Fixed and finite regardless of whether $\cX$ itself is finite — depends only on $\cA,\bar L$. |
| $\cS_j$ | **Sufficient class** — $\tau\in\cT_{\bar L}$ on whose leaves $\gamma_{0,j}$ is $P$-a.s. constant | `:112–121` (`def:sufficient`) | The manuscript-level stand-in for `theory.tex`'s formal `ass:sparsity`; here it is a *definition*, and "structural sparsity" is now the numbered `\label{ass:sparsity}` (`manuscript.tex:123`, drafted 2026-09-21) built from it — see §2a. |
| $c\in(0,1)$ | Uniform one-sided positivity constant | `:120–126` (`ass:causal`) | |
| $\psi(O;\theta,\eta)$ | Efficient influence function (EIF) for the ATT | `:132–139` (eq. `score`) | $\psi_0:=\psi(\cdot;\theta_0,\eta_0)$. |
| $\cT_n^{(\mu)},\ \cT_n^{(e)}$ | Feasible sets — partitions whose every leaf carries $\ge m_n$ (control) observations | `:181–185` (eq. `feasible`) | Data-dependent (through the realized leaf counts); each nonempty w.p.$\to1$ — added 2026-09-25 after a `domain-reviewer` finding that non-emptiness/randomness were previously unaddressed. |
| $m_n$ | Minimum per-leaf sample-size threshold | `:181–185` | **Corrected 2026-09-25** — now directly conditioned in `manuscript.tex` ($1\le m_n\to\infty$, $m_n=o(n)$), not just via the imported $\Lbar m_n\lesssim n$ (see the `theory.tex:3167` rate-mismatch note below); the parent-leaf no-controls fallback (`theory.tex`'s `ass:construct`(d)) is also now stated here. |
| $\lambda_n$ | Penalty on leaf count in the selection criterion | `:188–192` (eq. `select`) | Governed by $\lambda_n\to0$ — see §2c. |
| $\tauhat_j$ | Selected partition for nuisance $j$ | `:188–192` | **Corrected 2026-09-25** — `eq:select`'s own display previously used `\that_j` (renders as $\hat t_j$, hat over Latin *t*, per `common-defs.tex:291`), a broken symbol chain against every downstream use of $\tauhat_j$/`\tauhat`; fixed at `:189` and at the one other stray occurrence (`:428`, the degrees-of-freedom divisor). |
| $R_n^{(j)}(\tau)$ | Empirical risk of $\tau$ for nuisance $j$ | `:193–198` (eq. `risk-emp`) | **Corrected 2026-09-25** — $R_n^{(\mu)}$ was on the control-*conditional* scale ($\div n_0$, $\mid A{=}0$), which `theory.tex:3116–3121` explicitly rules out by name (silently rescales `prop:selection-rate`'s $\Delta_\mu$ by $(1-\pi)^{-1}$ relative to `theory.tex`'s Step 4). Now $R_n^{(\mu)}(\tau)=\Prob_n[(1-A)\{Y-\hat\mu_\tau(\bX)\}^2]$, unnormalised and control-*weighted*, matching `eq:risk-pop`'s population form and §2.1's $\|\cdot\|_{2,w}$ directly. |
| $\hat e_\tau(\bx),\ \hat\mu_\tau(\bx)$ | Leaf-wise empirical plug-ins (empirical proportion treated / control-outcome mean); both unclipped as a family in $\tau$ | `:199–204` (eq. `refit`) | `manuscript.tex` **overloads** $\hat e_{\tauhat_e}$ (only at the *selected* partition) to denote the clipped value from `ass:construct` — see that assumption's own registry entry above for the 2026-09-25 fix. |
| $\hat\tau_e,\hat\tau_\mu$ | Selected propensity / outcome partitions | throughout §2.4 | |
| $\hat\theta=\hat\theta(\hat\tau_e,\hat\tau_\mu)$ | doubletree point estimate | `:168–175` | |
| $\hat w_{\hat\tau_e}=\hat e_{\hat\tau_e}/(1-\hat e_{\hat\tau_e})$ | Fitted odds-scale propensity | `:175` | |
| $R^{(j)}(\tau)$ | Population risk (population counterpart of $R_n^{(j)}$) | `:176–182` (eq. `risk-pop`) | |
| Condition (P) | General partition-recovery condition | referenced `:211–220`, **not stated in this document** | Deferred entirely to the companion theory document; `manuscript.tex` uses it only to say the direct sparsity argument does not need to pass through it. Out of scope for this backfill (theory.tex excluded). |
| $\Delta_j$ | Population risk margin for nuisance $j$ — gap between the best-in-class and best-outside-$\cS_j$ risk | `:257–270` (`prop:selection-rate`) | |
| $c_1>0$ | Constant in the exponential selection-consistency bound | `:271–279` | Local to `prop:selection-rate`; depends on the outcome bound, $c$, and the sufficient partitions' leaf masses. Not reused elsewhere — included per the template's literal rule ("one row per symbol in a formal statement") but flagged as non-recurring, unlike every other row in this table. |
| $\hat V$ | Plug-in variance estimator | `:347–362` (`cor:variance`) | $\hat V\to_pV$ where $V=\E[\psi_0^2]$ (defined via `prop:P-instantiation1`'s hypothesis, `:322`). |
| $n_1=\sum_{i\le n}A_i$ | Treated-arm sample size | `:404` | |
| $\hat\theta_{\mathrm{tie}}$ | Single-outcome-tree "tied" plug-in estimator (no propensity tree) | `:396–406` (`cor:single-tree`) | |
| $|\ell|_1=|\ell|-|\ell|_0$ | Treated count within leaf $\ell$ | `:406` | |
| $\hat w_{\mathrm{tie}}(\bx)=|\ell_{\hat\tau_\mu}(\bx)|_1/|\ell_{\hat\tau_\mu}(\bx)|_0$ | Empirical treated-to-control odds within the outcome tree's leaf | `:406` | |
| $e_\tau(\bx)=\Prob(A=1\mid\bX\in\ell_\tau(\bx))$ | Leaf-averaged (not leaf-refit) true propensity | `:413` | Distinct from $\hat e_\tau$ (the empirical leaf plug-in): $e_\tau$ is a population object, the true propensity's leaf average, used only in the single-tree corollary family's population-side analysis. |
| $w_\tau=e_\tau/(1-e_\tau)$; $\mu_\tau=\Pi_\tau^\nu\mu_0$ | Leaf-average propensity odds; control-weighted leaf projection of $\mu_0$ | `:413` | |
| $\tau_j^\dagger$ ($j\in\{e,\mu\}$) | The unique population-risk minimizer over $\cT_{\Lbar}$ for nuisance $j$, when it exists | New, `ass:pseudo` (~`:472`) | Existence/uniqueness is the content of `ass:pseudo`, not automatic. Under structural sparsity, $\tau_j^\dagger$ is the finest partition in $\cS_j$. |
| $\Delta_j^\dagger$ | Margin between $\tau_j^\dagger$'s population risk and every other candidate's | New, `ass:pseudo` (~`:472`) | Redundant given uniqueness + finiteness of $\cT_{\Lbar}$ (a minimum over finitely many strictly positive numbers is strictly positive) — kept in the assumption's statement for readability/parallelism with `prop:selection-rate`'s $\Delta_j$, not as an independent requirement. |
| $\eta_*=(e_*,\mu_*)$, $\theta_*$ | Pseudo-true nuisance pair ($e_*=e_{\tau_e^\dagger}$, $\mu_*=\Pi^\nu_{\tau_\mu^\dagger}\mu_0$, two generally distinct partitions) and the effect it implies (solves $\E[\psi(\bO;\theta,\eta_*)]=0$) | New, `ass:pseudo` (~`:472`) | Do not confuse with $\eta_\tau,\theta_\tau$ (`:413`), which are indexed by one shared $\tau$; $\eta_*$'s two components carry their own, generally different, partitions. Under structural sparsity, $\eta_*=\eta_0$. |
| $\eta_\tau=(e_\tau,\mu_\tau)$; $\theta_\tau$ solves $\E[\psi(\bO;\theta,\eta_\tau)]=0$; $\psi_\tau=\psi(\bO;\theta_\tau,\eta_\tau)$ | Population objects at a fixed partition $\tau$ | `:413` | |
| $\eta_{\mathrm{tie}}=\eta_{\hat\tau_\mu}$, $\theta_{\mathrm{tie}}=\theta_{\hat\tau_\mu}$, $\psi_{\mathrm{tie}}=\psi_{\hat\tau_\mu}$ | The above, evaluated at the *realized* (random) partition | `:413` | Random-index objects — `lem:single-tree-linear`'s whole point is that these need not be a single fixed number across realizations of $\hat\tau_\mu$ unless $\cS_\mu$ is a singleton or $\hat\tau_\mu$ stabilizes. |
| $\bar Y_1(\ell)$ | Treated-arm sample mean within leaf $\ell$ | `:435` (inside `lem:single-tree-linear`'s proof) | |
| $V_\tau=\Var(\psi_\tau)$ | Fixed-$\tau$ asymptotic variance | `:447` (`cor:single-tree-coarsening`) | $V_\tau\le V$ under within-leaf homoskedasticity; the ordering can reverse under heteroskedasticity (`:451`) — this is the corollary's own worked numerical example, not a separate symbol. |
| $\bar e_\ell=\E[e_0(\bX)\mid\bX\in\ell]$ | Leaf-average true propensity | `:451` | |
| $f(e)=e^2/(1-e)$ | Strictly convex function used in the Jensen's-inequality variance comparison | `:451` | |
| $\sigma_{0,\ell}^2$ | Leaf-average control-outcome variance (leaf-average of $\sigma_0^2(\bx)$) | `:451` | See the flag on $\sigma_a^2(\bx)$ above — do not conflate the pointwise and leaf-average versions. |
| $\tau_\star$ | Deterministic partition $\hat\tau_\mu$ is assumed to stabilize to | `:453` | Local to `cor:single-tree-coarsening`'s stabilization branch. |
| $D_w=\|\hat w-w_0\|_{2,w}$, $D_\mu=\|\hat\mu-\mu_0\|_{2,w}$ | **Realized** nuisance errors (random, properties of the actual fit) | `:497–500` (`lem:biasbound`) | **Do not conflate with $\delta_e,\delta_\mu$** — `:515` explicitly warns bias is controlled by $D_w,D_\mu$, not by the population approximation errors. **2026-09-21, corrected:** this row previously used a `\Lw{\cdot}` macro borrowed from `theory.tex:55`, which does not exist in `common-defs.tex` — pasting it into `manuscript.tex` would fail to compile. Fixed to the manuscript's actual notation, `\|\cdot\|_{2,w}$. Also flagged: at `$D_w=\|\hat w-w_0\|_{2,w}$` the letter $w$ carries three roles at once (propensity odds $\hat w,w_0$; the $_{2,w}$ control-weight subscript; $D_w$'s own subscript) — §2.1 (`:70–71`) now has a clarifying half-sentence disambiguating the $_{2,w}$ subscript from the propensity-odds $w_0,w$, added in the same 2026-09-21 pass. |
| $\delta_e,\delta_\mu$ | **Population approximation** errors — best achievable leaf-constant error at budget $\bar L$, a deterministic functional of $P,\cT_{\bar L}$ | referenced `:515`, formally defined only in `theory.tex:270–276` | Out of scope to backfill from `theory.tex`, but flagged: `manuscript.tex` uses the term without giving it its own display equation in this document — a reader relying on `manuscript.tex` alone must take the distinction from $D_w,D_\mu$ on faith at `:515` rather than seeing $\delta_e,\delta_\mu$ defined. Worth a forward-reference to the companion document if this is ever revised. |
| $\widetilde\theta$ | The **anchor** estimator — any $\sqrt n$-consistent, asymptotically normal estimator | `:517` | Need not be interpretable or tree-based. |
| $\widetilde e^{(-k)},\widetilde\mu^{(-k)}$ | Fold-wise anchor nuisance fits, $K$ folds fixed | `:519–530` (`ass:rate`) | |
| $\widetilde V$ | Anchor's asymptotic variance | `:544–551` (`thm:anchor`) | |
| $\widetilde\sigma^2$ | Anchor's variance estimator, $\widetilde\sigma^2\to_p\widetilde V$ | `:544–553` | |
| $\widehat\delta=\hat\theta-\widetilde\theta$ | Anchor-interval bias proxy | `:555` | |
| $C_n$ | Anchor confidence interval | `:557–568` (`thm:anchor`) | |
| $\eta_*$ | Probability limit of $\hat\eta$ under a sparsity-violating alternative | `:616–620` (`prop:spectest`) | |
| $b:=\E[\psi(\bO;\theta_0,\eta_*)]$ | Population bias under that alternative | `:620` | `prop:spectest`'s consistency claim is exactly that the diagnostic detects $b\ne0$. |
| $\widehat{\mathrm{SE}}\asymp n^{-1/2}$ | Generic standard-error rate used by the fidelity-diagnostic test | `:609–614` | Not a specific estimator — any SE at the parametric rate qualifies. |
| $\cP_\tau$ | Submodel in which $e_0$ is exactly constant on the leaves of $\tau$ | `:646–649` (`prop:sparse-bound`) | |
| **New, only if the optional coarsening-bound proposition (`outline.md` row 3.3) is adopted:** | | | |
| $L_e,L_\mu$ | Lipschitz constants for $e_0,\mu_0$ on $\cX$ | not yet in the manuscript — proposed | **Sensitivity-device parameters, not estimable from data.** Must be labeled local to this one proposition, not folded into the standing assumption list — see §2c. |
| $h$ | Maximum grid-atom diameter (mesh) | not yet in the manuscript — proposed | Analyst-chosen, same status as $\bar L$ and the cutpoints — known exactly, not estimated. |
| **New, 2026-09-21: Appendix~\ref{app:ate}'s ATE-extension symbols. All are Appendix-local by design and do not overwrite the main-text meanings of $\theta_0$, $\eta$, $\psi$, $c$, or $m$ — verified by `domain-reviewer` (2026-09-21) via exhaustive grep of the appendix for any bare (unsuperscripted) occurrence of these symbols; none found.** | | | |
| $\theta_0^{\mathrm{ATE}}=\E[Y(1)-Y(0)]$ | The ATE — distinct from the paper's primary estimand $\theta_0$ (the ATT) | `app:ate` (new) | Never write bare $\theta_0$ for this — that symbol is reserved for the ATT everywhere else in the document. |
| $\eta^{\mathrm{ATE}}=(e,m,\mu)$ | Candidate nuisance **triple** for the ATE — distinct arity from the main text's pair $\eta=(e,\mu)$ | `app:ate` (new) | $m$ here is a candidate for $m_1$ (the treated-outcome regression, already defined at `:53`) — not to be confused with the leaf-size floor $m_n$ (eq. `feasible`), which is unrelated and coexists in the same document. |
| $\psi^{\mathrm{ATE}}(\bO;\theta,\eta^{\mathrm{ATE}})$ | Efficient influence function for the ATE in the nonparametric model, standard AIPW form, `\citep{hahn1998}` | `app:ate`, eq. `score-ate` (new) | `domain-reviewer` (2026-09-21) verified all four signs against the standard result and confirmed mean-zero at the truth by direct substitution. No $\pi^{-1}$ normalization (that factor is ATT-specific, from conditioning on $A=1$; absent here since the ATE's estimating equation is not conditioned on either arm). |
| $c'\in(0,1/2)$ | Uniform two-sided positivity constant for the ATE, $c'\le e_0(\bx)\le1-c'$ | `app:ate` (new) | **Deliberately a distinct symbol from `ass:causal`'s $c\in(0,1)$** — the two-sided bound is a strictly stronger, generally smaller constant; `domain-reviewer`/`proofreader` (2026-09-21) both flagged that reusing the bare letter $c$ for both would read as claiming they are the same constant. |
| $\cS_m$ | Hypothetical third sufficient class, for $m_1$, in a three-nuisance analogue of structural sparsity | `app:ate`, `rem:bilinear-ate` (new) | **Not defined as a formal object** — invoked only inside a subjunctive ("if every nuisance were exactly tree-representable..."), alongside the already-defined $\cS_e,\cS_\mu$. Whether it should use the treated-weighted measure ($d\nu=e_0(\bX)\,d\Prob$, the natural analogue of $\cS_\mu$'s control-weighted $d\nu=\{1-e_0(\bX)\}\,d\Prob$ — **updated 2026-09-25** to match the corrected, measure-based `\Pi^\nu_\tau` convention above; no longer $p_x$-based) is not settled here. |

---

## 2. Assumption registry

**Environment name(s) in use:** `assumption` — every formal instance in `manuscript.tex`
(`ass:causal`, `ass:construct`, `ass:rate`) uses `\begin{assumption}[...]`. No other
environment name (`assum`, `cond`, `hyp`, etc.) occurs in this document.

### 2a. Formal — explicitly numbered in the paper

#### Structural sparsity — `\label{ass:sparsity}` (drafted 2026-09-21)

- **Stated at:** `:123–128`
- **Statement:** $\cS_e\ne\emptyset$ and $\cS_\mu\ne\emptyset$ — some tree in $\cT_{\Lbar}$
  represents the propensity exactly, and some tree in $\cT_{\Lbar}$ represents the
  control-outcome regression exactly, each in the sense of `def:sufficient`.
- **Discharges:** `lem:selection` (selection lands in the sufficient class w.p.$\to1$),
  `prop:P-instantiation1` (the headline CLT), and everything in §3.2 built on them; also the
  new legibility-vs-audit-completeness remark immediately following it (`:139–150`), which
  restates it by label rather than in spirit.
- **Role:** The single most load-bearing condition in the manuscript — delivers full
  semiparametric efficiency at the parametric rate. Its failure is the entire subject of
  `sec:honest-manuscript` (the anchor interval).
- **When reasonable:** The two nuisances' additive structure is concentrated enough, relative
  to the leaf budget $\Lbar$, that a tree of at most $\Lbar$ leaves represents each exactly on
  the analyst's pre-specified grid — cheap when dependence is hierarchical (the active
  coordinate set varies by region), expensive when it is additive across many active
  coordinates.
- **When unreasonable:** The truth's additive structure is spread widely enough across active
  covariates that no tree within the leaf budget represents it exactly. `thm:anchor` remains
  valid regardless; only the plain interval's efficiency claim is lost, absorbed into the
  anchor interval's width (`cor:width`).
- **Verifiable from data?** No — a property of the unobservable population nuisances, like
  `ass:rate`. `prop:spectest`'s diagnostic offers indirect evidence, not verification.
- **Used by:** `lem:selection`, `prop:P-instantiation1`. **Not** by
  `lem:single-tree-linear`/`cor:single-tree-coarsening` in the single-tree corollary family —
  those need only $\cS_\mu\ne\emptyset$, or add $\cS_e\ne\emptyset$ as a separate hypothesis,
  rather than the joint condition, so citing `ass:sparsity` there would overstate what they
  actually assume; left as direct $\cS_e/\cS_\mu$ prose (`:711,:739`) on purpose.
- **Cited as a range anywhere?** No.
- **Declared in:** body (§2.2).

#### Causal identification — `\label{ass:causal}`

- **Stated at:** `:120–126`
- **Statement:** Consistency $Y=AY(1)+(1-A)Y(0)$; ignorability $\{Y(0),Y(1)\}\perp A\mid X$;
  $\pi>0$; $\exists\,c\in(0,1)$ with $1-e_0(x)\ge c$ for every $x\in\cX$.
- **Discharges:** Identification of $\theta_0=\E\{Y-\E(Y\mid X,A=0)\mid A=1\}$; the *uniform*
  bound $c$ (as opposed to pointwise $e_0(x)<1$) specifically keeps the variance-decomposition
  term $\E[e_0^2\sigma_0^2/(1-e_0)]$ finite (`:128`).
- **Role:** Conditions on the **full** $X$ — the assumption a coarsening/discretization
  discussion must be checked against, since it is unaffected by how the estimator later
  represents $X$ for tree-fitting purposes. See §3 below.
- **When reasonable:** Standard causal-inference setting where all confounders are captured in
  the measured $X$, of whatever dimension or type.
- **When unreasonable:** Genuine unmeasured confounding — a variable outside $X$ affecting both
  $A$ and $Y$. This is the failure mode "confounding" language properly refers to, and it is
  categorically different from coarsening a measured $X$ (§3).
- **Verifiable from data?** No — see §2e.
- **Used by:** `eq:score`; every downstream result in §§2.3–4 (essentially the whole paper).
- **Cited as a range anywhere?** No.
- **Declared in:** body (§2.3).

#### Clipping — `\label{ass:construct}`

- **Stated at:** `:206–208`
- **Statement:** **Corrected 2026-09-25** — was an unscoped "fitted propensities are clipped for
  every $\bx$," which literally applied to `eq:refit`'s leaf refit and hence to `eq:select`'s
  argmin and $R^{(e)}_n$ itself, contradicting `theory.tex:234–235`'s explicit "structure
  selection uses the unclipped leaf refit" and leaving $\thetahat$ (displayed 20 lines earlier)
  technically undefined wherever a leaf is all-treated. Now scoped: `manuscript.tex` overloads
  $\hat e_{\tauhat_e}$ (evaluated at the *selected* partition only) to mean the clipped value,
  $1-\hat e_{\tauhat_e}(\bx)\ge c$ for every $\bx$; the general family $\hat e_\tau$ used inside
  the argmin and $R^{(e)}_n$ stays unclipped throughout. Matches `theory.tex:235`'s own
  "$\ehat=\ehat_{\that_e}$" overloading move exactly, rather than inventing new notation.
- **Discharges:** Keeps the fitted nuisances inside `ass:causal`'s positivity region, so
  downstream bias/variance bounds (`lem:biasbound`, `thm:anchor`, `cor:variance`'s proof)
  apply to the estimator actually reported, not an unclipped one.
- **Role:** Reconciles estimation with `ass:causal`'s positivity bound.
- **Not (yet) restored:** `theory.tex`'s `ass:construct` is a four-part bundle (Loss, Leaf refit,
  Clipping, Leaf mass); the manuscript's version under this same label is Clipping only. Coverage
  of the other three parts is now distributed elsewhere in §2.4 (Loss: `eq:risk-emp`'s own
  squared-error statement; Refit: `eq:refit`; Leaf mass: the new $m_n$ condition text after
  `eq:feasible`, `:185`) rather than bundled under `ass:construct` itself — flagged by a
  `domain-reviewer` audit (2026-09-25,
  `quality_reports/reviews/2026-09-25_sec2.4-algorithm-audit.md`) as MINOR/not-yet-actioned: a
  reader tracing "the conditions of Proposition~\ref{lem:selection}" via `ass:construct` alone
  cannot recover the leaf-mass floor from that citation, only from the surrounding prose.
- **When reasonable:** Whenever a fixed clip constant consistent with domain knowledge about
  achievable overlap is chosen.
- **When unreasonable:** If $c$ is misspecified too aggressively, clipping could itself
  introduce bias — not separately analyzed in this manuscript.
- **Verifiable from data?** Partially — enforced by construction, not verified.
- **Used by:** `lem:selection`, `lem:biasbound`, `thm:anchor`, `cor:variance`.
- **Cited as a range anywhere?** No.
- **Declared in:** body (§2.4).

#### Anchor nuisance rate — `\label{ass:rate}`

- **Stated at:** `:519–534`
- **Statement:** $\max_{k\le K}\twonorm{\widetilde e^{(-k)}-e_0}\,
  \|\widetilde\mu^{(-k)}-\mu_0\|_{2,w}=o_p(n^{-1/2})$ for fold-wise fits, $K$ fixed.
- **Discharges:** The product-rate condition making the anchor itself $\sqrt n$-efficient,
  feeding `prop:spectest`'s diagnostic and `cor:width`'s "efficient anchor" branch.
- **Role:** Lets the anchor be built by *any* method (not necessarily trees) while still
  supporting the diagnostic and efficient-width results.
- **When reasonable:** Flexible black-box nuisance estimation (e.g. cross-fit boosting) is
  available even though it would not itself be displayed.
- **When unreasonable:** No flexible nuisance estimator satisfying the product-rate condition
  is available — then `thm:anchor` alone (needing only $\sqrt n$-consistency and asymptotic
  normality of the anchor) still gives validity, just not efficient width or the diagnostic.
- **Verifiable from data?** No — see §2e.
- **Used by:** `prop:spectest`; `cor:width`'s efficient-anchor branch.
- **Cited as a range anywhere?** No.
- **Declared in:** body (§3.3).

#### Pseudo-true margin — `\label{ass:pseudo}`

- **Stated at:** ~`:472–478` (new, inserted 2026-09-17 at the head of §3.3, before `prop:bilinear`)
- **Statement:** For each $j\in\{e,\mu\}$, a unique $\tau_j^\dagger\in\cT_{\Lbar}$ minimizes $R^{(j)}$
  over $\cT_{\Lbar}$ (the margin $\Delta_j^\dagger>0$ is then automatic — see the symbol table's
  note).
- **Discharges:** Gives `prop:pseudo-consistency` a fixed target $(\tau_e^\dagger,\tau_\mu^\dagger)$
  to land selection on when structural sparsity fails but a well-defined best approximation still
  exists; supplies `prop:spectest`'s consistency-direction hypothesis as a derived consequence
  ($\hat\eta\to_p\eta_*$) rather than an assumed one.
- **Role:** Strictly weaker than structural sparsity — requires uniqueness of the best
  tree-representable approximation, not exactness. Promoted from a condition previously bundled,
  unlabeled, into `prop:spectest`'s own hypothesis line (see §2c below, now superseded).
- **When reasonable:** Whenever the population risk surface over the fixed finite class
  $\cT_{\Lbar}$ has no exact ties among its best candidates — a strictly weaker requirement than
  structural sparsity, and one violated only by an exact tie (a near-tie leaves the assumption
  intact and only affects the rate via `prop:selection-rate`'s exponent).
- **When unreasonable:** An exact tie between two or more risk-minimizing partitions — then
  $\tau_j^\dagger$ is undefined, $\eta_*,\theta_*$ do not exist, and the selected partitions can
  oscillate with no fixed limit. `thm:anchor` remains valid regardless; only the pseudo-true
  target and `prop:spectest`'s consistency direction are lost.
- **Verifiable from data?** No — a property of the population risk surface, like `ass:rate`.
- **Used by:** `prop:pseudo-consistency`; `prop:spectest`'s consistency direction.
- **Cited as a range anywhere?** No.
- **Declared in:** body (§3.3, at the head of the section, before `prop:bilinear`).

### 2b. Implicit — carried in constraints or prose, not numbered

| Assumption | Where it hides | Role | Should it be promoted to formal? |
|---|---|---|---|
| **Structural sparsity** — **MOVED to §2a, 2026-09-21.** Promoted to `\label{ass:sparsity}` (`manuscript.tex:123`); this row's own earlier finding ("Yes, arguably" it should be promoted) is now resolved. See §2a's new entry for the full record. | — | — | **Resolved 2026-09-21** |
| **iid sampling** | Bundled into `lem:selection`'s hypothesis list, `:241` ("the observations are iid") | Standard sampling assumption underlying every asymptotic result in the paper. | Optional — 9.5% of methods-theory papers name every assumption (per `paper-protocol.md`'s measurement); bundling a sampling assumption into a proposition's hypothesis rather than a separate numbered environment is common practice, not a defect. |
| **Uniformly bounded conditional second moments of $Y$** | Bundled into `lem:selection`'s hypothesis (`:242`), `cor:variance`'s proof (`:376–377`) | Needed for the LLN/CLT arguments underlying `lem:selection`, `cor:variance`. | Optional, same reasoning as above. |
| **$Y$ bounded** (stronger than the second-moment condition) | Bundled into `prop:selection-rate`'s hypothesis, `:259` | Needed for the exponential (rather than merely polynomial) selection-consistency rate. | Optional. |
| **Exact global minimization (solver fidelity)** | Bundled into `lem:selection`'s hypothesis, `:246–247` ("$\hat\tau_j$ is an exact global minimizer... with $\lambda_n\to0$") | The margin-based selection-consistency argument (`prop:selection-rate`) needs the *actual* global minimizer, not merely a good local one — this is `theory.tex`'s `ass:solver`, restated here only inline. | Optional to promote in this document; `theory.tex` already promotes it. |

### 2c. Regularity conditions — technical requirements

| Condition | Stated at | Role | Relaxable? |
|---|---|---|---|
| $\bar L$ fixed (not growing with $n$) | `:96–102`, restated in `lem:selection`'s hypothesis `:245` | Keeps $\cT_{\bar L}$ a fixed finite set, the property the whole no-splitting argument leans on. | Not within this framework — a growing $\bar L$ is exactly the relaxation `theory.tex`'s Condition (P) / general-anchor track handles instead, out of scope here. |
| $\bar Lm_n\lesssim n$ | `lem:selection`'s hypothesis, `:246` | Keeps the minimum-leaf-mass feasibility constraint from binding asymptotically. | Yes — a rate condition, standard. |
| $\lambda_n\to0$ | `lem:selection`'s hypothesis, `:247`; `lem:uniform` notes exact global optimization and $\lambda_n$'s rate are *not* required for its own (weaker) uniform statement, `:313–314` | Penalty vanishing fast enough not to distort selection asymptotically. | Yes — a tuning-parameter rate condition. |
| $V:=\E[\psi_0^2]>0$ | `prop:P-instantiation1`'s hypothesis, `:322` | Non-degeneracy — rules out a degenerate limiting variance. | Not really "relaxable"; a genuine non-degeneracy requirement for the CLT to be non-trivial. |
| **Superseded 2026-09-17** — unique, margin-separated population risk minimizer under the alternative | `prop:spectest`'s hypothesis, `:616–618` (pre-edit line numbers) | Promoted to formal `ass:pseudo` (§2a above), at the head of §3.3. `prop:spectest`'s own text now cites `ass:pseudo` by label instead of restating the condition inline. | Done — see §2a. |
| **Lipschitz nuisances on grid atoms** (proposed, not yet in the manuscript) | Would be new, in §3.3 | Would bound the coarsening-specific component of $\delta_e,\delta_\mu$ via mesh $h$: $\delta_e\lesssim L_eh$, $\delta_\mu\lesssim L_\mu h$. | **Yes, and must be presented as relaxable/local** — a sensitivity-analysis device attached to one proposition, not a standing assumption that would re-impose the smoothness restriction the fixed-grid framework was built to avoid (Oracle's explicit warning, carried over from the prior pass). |

### 2d. Summary — what is assumed and why

| Assumption(s) | Role / why needed |
|---|---|
| `ass:causal` | Identifies $\theta_0$ from observed data; the uniform positivity bound keeps the efficiency bound's variance decomposition finite. |
| `ass:construct` | Keeps fitted nuisances inside the positivity region `ass:causal` requires. |
| Structural sparsity ($\cS_e,\cS_\mu\ne\emptyset$) | Delivers full parametric-rate efficiency with no splitting; its absence is what §3.3–3.4 are about. |
| `ass:rate` | Lets the anchor be $\sqrt n$-efficient, powering the fidelity diagnostic and the efficient-width branch of `cor:width`. |
| Sampling/moment/solver conditions (§2b) | Standard technical scaffolding for the LLN/CLT machinery; bundled into individual propositions rather than separately numbered. |
| $\bar L$-fixed / $\lambda_n$-rate / margin conditions (§2c) | Keep the no-splitting, fixed-finite-class argument mechanically valid. |

### 2e. Blanket statements

| Statement | Covers (assumption labels) | Stated at | Reason |
|---|---|---|---|
| Not empirically verifiable | `ass:causal` (ignorability clause), `ass:rate` | Not currently stated together as a single sentence anywhere in `manuscript.tex` — recorded here as a registry-level observation about the manuscript's content, not yet reflected in its prose. | Ignorability given $X$ can never be checked from observed data; the anchor's rate condition concerns unobservable population quantities ($e_0,\mu_0$ themselves). |

---

## 3. Identification status

- **Point-identified from observed data alone, given `ass:causal`:** $\theta_0$, via
  $\theta_0=\E\{Y-\E(Y\mid X,A=0)\mid A=1\}$ (`:130`) — using the **full**, possibly continuous
  $X$, not the discretized grid-atom version.
- **Partially identified, given `<what>`:** none in the manuscript.
- **Not identified without further assumptions:** $\theta_0$, if `ass:causal`'s ignorability
  clause fails (genuine unmeasured confounding — a variable outside $X$ affecting both $A$ and
  $Y$).

**If an assumption fails, what happens?**

| Assumption | If it fails | Consequence for the estimand |
|---|---|---|
| `ass:causal` (ignorability given full $X$) | Genuine unmeasured confounding | $\theta_0$ is **not identified**; no estimator, anchor or otherwise, repairs this. The failure mode "confounding" language properly refers to. |
| Structural sparsity ($\cS_e,\cS_\mu\ne\emptyset$) — includes, as a special case, coarsening a continuous $X$ so $e_0,\mu_0$ are not exactly atom-constant | Sparsity fails | $\theta_0$ **remains identified and point-estimable**; only the doubletree point estimate's efficiency and the plain CI's validity are lost. `thm:anchor`'s interval remains valid with no rate/smoothness condition on the doubletree nuisances. Estimation-side working-model misspecification, not an identification failure. |
| `ass:construct` (clipping) | Clip constant misspecified | Not separately analyzed in the manuscript as of this pass — recorded as an open question, not a resolved consequence. |
| `ass:rate` (anchor nuisance rate) | Anchor's fold-wise fits do not meet the product-rate condition | `thm:anchor`'s coverage guarantee is **unaffected** (it needs only $\sqrt n$-consistency/normality of the anchor); `cor:width`'s efficient-width conclusion and `prop:spectest`'s power guarantee are what is lost. |

---

## 4. Maintenance contract

- A drafting or editing skill that introduces a symbol, assumption, or numbered result updates
  this file in the same turn.
- Divergence between manuscript and registry is a finding, not a cue to silently sync.
- This file is now a full backfill of `manuscript.tex` (confirmed scope, 2026-09-16); if
  `theory.tex` is ever brought into scope, that is a new, separate backfill, not an incremental
  extension of this file's existing rows.
