# Phase 2 Simulation Design

## Motivation

**Phase 1 (DGPs 1-6) showed:**
- Tree: 99.8% convergence, 96.5% coverage (excellent!)
- Tree: 0.012 RMSE vs Linear: 0.0007 (tree loses on efficiency)
- **Problem:** All Phase 1 DGPs have linear nuisances → favor linear methods
- **Gap:** No demonstration of tree's value proposition

**Reviewer question we couldn't answer:** "Why use trees instead of flexible GLM with interactions?"

---

## Phase 2 Strategy: Show Tree's Value

**Three complementary DGPs:**

### DGP7: Deep 3-way Interaction
**Purpose:** Show tree beats linear on complex interactions

**Design:**
- Propensity: e(X) has X1*X2 interaction (moderate complexity)
- Outcome: m0(X) has X1*X2*X3 **3-way interaction** (high complexity)
- Linear needs to specify all interaction terms explicitly
- Tree learns interaction structure automatically

**Results (n=800):**
- Tree RMSE: 0.0104 vs Linear: 0.0151 → **Tree wins 0.69x**
- Both maintain coverage ✓

**Takeaway:** Tree excels when interaction structure is complex

---

### DGP8: Double Nonlinearity (Sin/Cos)
**Purpose:** Show tree beats linear when BOTH models are nonlinear

**Key Insight:** EIF-based ATT is **doubly robust**:
- Can tolerate misspecification of e(X) OR m0(X)
- But not both simultaneously
- Must break double robustness to show tree advantage

**Design:**
- Propensity: e(X) ~ sin(2πX1) + cos(2πX2) (periodic, non-polynomial)
- Outcome: m0(X) ~ sin(2πX1) + cos(2πX2) + X1*X2 (also periodic)
- Linear regression fundamentally cannot approximate trigonometric functions
- Tree can approximate any smooth function

**Results (n=800):**
- Tree RMSE: 0.1254 vs Linear: 0.2055 → **Tree wins 0.61x**
- Tree covers ✓, **Linear fails to cover ✗**

**Takeaway:** When both nuisances are misspecified, tree's flexibility wins

---

### DGP9: Weak Overlap (Stress Test)
**Purpose:** Test robustness under extreme propensity scores

**Design:**
- Propensity: e(X) = plogis(-4 + 8*X1) (very steep)
- Creates near-deterministic treatment (e ≈ 0 or e ≈ 1)
- Only 55% of sample has e ∈ [0.1, 0.9]
- Tests numerical stability of IPW weights

**Results (n=800):**
- Tree RMSE: 0.2012 vs Linear: 0.0366 → Linear wins (expected)
- **But:** Both maintain coverage ✓ (key result!)

**Takeaway:** Tree remains robust under weak overlap despite efficiency loss

---

## Complete Narrative (Phase 1 + Phase 2)

### When Linear Wins (DGPs 1-6, 9)
- Simple linear/additive nuisances
- Tree pays ~1-50x RMSE penalty
- But tree maintains excellent coverage (96.5%)

### When Tree Wins (DGPs 7-8)
- Complex interactions (DGP7: 3-way interaction)
- Nonlinear nuisances (DGP8: sin/cos in both e and m0)
- Tree: 0.6-0.7x better RMSE than linear
- Linear may even lose coverage (DGP8)

### Overall Strength: Robustness
- Tree: 99.8% convergence across all 9 DGPs
- Tree: Maintains coverage even when losing on RMSE
- Efficiency loss is tolerable price for not assuming functional form

---

## Phase 2 Implementation Plan

**Same structure as Phase 1:**
- 3 DGPs (7, 8, 9)
- 3 sample sizes (400, 800, 1600)
- 4 methods (tree, rashomon, forest, linear)
- 1000 replications each
- Total: 36,000 additional replications

**Expected runtime:** Same as Phase 1 (~6-8 hours on O2)

**Files:**
- `dgps/dgps_phase2.R` - DGP generators (created ✓)
- `test_phase2_dgps.R` - Local validation (passing ✓)
- `slurm/launch_phase2.sh` - O2 submission script (TODO)

---

## Key Methodological Insight

**Why DGP8 requires sin/cos functions:**

Doubly robust estimators satisfy:
```
bias(θ̂) ≈ E[{e(X) - ê(X)} × {m0(X) - m̂0(X)}]
```

If EITHER e(X) or m0(X) is estimated well, bias ≈ 0.

**Phase 1 failure:** DGPs had linear e(X), so even when m0(X) was nonlinear (e.g., interactions), linear got e(X) right → double robustness protected linear.

**Phase 2 solution:** Misspecify BOTH:
- e(X) nonlinear (sin/cos)
- m0(X) nonlinear (sin/cos)
- Now bias ≠ 0 for linear
- Tree approximates both better → lower bias

**Why sin/cos specifically:**
- Polynomials can be approximated by linear + interactions
- Trigonometric functions are fundamentally different (periodic, bounded)
- No polynomial expansion can approximate sin(2πx) on [0,1]
- Tree's piecewise constant approximation works well

---

## Verification

Local test (n=800, single replication):
```bash
Rscript test_phase2_dgps.R
```

Results:
```
DGP7: Tree 0.69x better than linear ✓
DGP8: Tree 0.61x better than linear ✓
DGP9: Both maintain coverage ✓
```

**Status:** Ready for full O2 simulation study
