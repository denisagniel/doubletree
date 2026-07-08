# Lambda Sensitivity Test - Preliminary Findings

**Date:** 2026-05-22
**Status:** In Progress (50% complete for theory lambda)

---

## Test Configuration

**DGP:** Complex (dgp=4)
**n:** 2000
**Replications:** 100 per lambda
**Approaches:** 1-3 (full_sample, crossfit, doubletree)

**Lambda values tested:**
- λ = 0.1 (current baseline) ✓ Complete
- λ = 0.05 (intermediate) ✓ Complete
- λ = log(2000)/2000 ≈ 0.0038 (theory) ⏳ 50/100 reps

---

## Critical Finding: Theory Lambda Computationally Infeasible

### Model Limit Warnings

**Total warnings:** 54 (and counting)
**Distribution by lambda:**
- λ = 0.1: **0 warnings** ✓
- λ = 0.05: **0 warnings** ✓
- λ = 0.0038 (theory): **54 warnings** (36% of fits so far)

**Pattern:** All computational issues are with theory lambda.

### What "Model limit exceeded" Means

**From TreeFARMS/GOSDT solver:**
- Solver explores tree space up to a model limit (default: 10,000 models)
- When limit hit: returns best tree found so far (may be suboptimal)
- Indicates search space is too large to explore exhaustively

**Why theory lambda hits limit:**
- λ ≈ 0.0038 is **very weak regularization**
- Allows extremely complex trees (potentially 100+ leaves)
- Combinatorial explosion in tree space
- Solver cannot explore all possibilities within computational budget

### Practical vs Theoretical Trade-off

**Theory says:** λ ~ log(n)/n minimizes worst-case risk (minimax optimal)

**Practice shows:** λ = 0.0038 at n=2000 exceeds computational limits

**Implications:**
1. **Theory assumes unlimited computation** - oracle can find optimal tree in T_{M_n}
2. **Practice has constraints** - solver has finite time/memory budget
3. **Gap widens with n** - at n=5000, theory would prescribe λ ≈ 0.0017 (even worse)

---

## Preliminary Observations

### Computational Feasibility

**λ = 0.1 (current):**
- Fast: ~7 minutes for 300 fits (100 reps × 3 approaches)
- No solver warnings
- Trees fit reliably

**λ = 0.05 (half current):**
- Fast: ~7 minutes for 300 fits
- No solver warnings
- Trees fit reliably

**λ = 0.0038 (theory):**
- Slow: ~50+ minutes expected for 300 fits
- **36% of fits hit model limit** (54/150 so far)
- Some trees may be suboptimal due to incomplete search
- Practical feasibility questionable

### Runtime Comparison

| Lambda | Time per 300 fits | Warnings | Feasibility |
|--------|-------------------|----------|-------------|
| 0.1 | 7 min | 0 | ✓ Excellent |
| 0.05 | 7 min | 0 | ✓ Excellent |
| 0.0038 | ~60 min | 54+ | ⚠️ Problematic |

**Ratio:** Theory lambda is **8-9× slower** than current lambda.

---

## Revised Hypothesis

### Original Hypothesis

λ = 0.1 is too large → over-regularization → under-coverage

### Updated Hypothesis (Pending Full Results)

**Two competing effects:**

1. **Statistical:** Larger λ → simpler trees → worse approximation → SE underestimation
2. **Computational:** Smaller λ → search space explosion → solver fails → unpredictable quality

**Possibility:** λ = 0.05 might be the **practical optimum**
- 2× weaker than current (more flexibility)
- Still computationally feasible (no solver warnings)
- Fast enough for production use

### Questions to Answer

**When full results available:**

1. **Does λ = 0.05 improve coverage over λ = 0.1?**
   - If yes: confirms over-regularization hypothesis
   - If no: problem is elsewhere

2. **Does λ = 0.0038 improve coverage further?**
   - If yes and no warnings: theory is right, just slow
   - If yes but many warnings: unstable, not practical
   - If no: theory lambda doesn't help (unexpected!)

3. **Do solver warnings correlate with poor coverage?**
   - Check if fits with "Model limit exceeded" have worse coverage
   - If yes: incomplete search degrades inference

---

## Practical Implications

### For Future Simulations

**Don't use theory lambda directly at large n:**
- λ ~ log(n)/n becomes computationally infeasible
- 36%+ failure rate at n=2000
- Would be worse at n=5000, n=10000

**Consider modified theory lambda:**
- Add computational feasibility constraint
- Example: λ = max(log(n)/n, λ_min) where λ_min ensures tractability
- Or: λ = c · log(n)/n with c > 1 (e.g., c=5-10)

### For Paper

**Need to address theory-practice gap:**

**Current claim (line 183):**
> "regularization parameter λ_n is chosen by 5-fold cross-validation"

**Reality:**
- Simulations: λ = 0.1 (fixed)
- Theory: λ ~ log(n)/n
- Practice: λ ~ log(n)/n hits solver limits

**Options:**

1. **Acknowledge constraint:**
   - "Theory prescribes λ ~ log(n)/n for minimax optimality"
   - "In practice, computational constraints may require λ > log(n)/n"
   - "We use λ = 0.1 as computationally feasible approximation"

2. **Propose modified theory:**
   - "Theory assumes oracle access to optimal tree in T_{M_n}"
   - "With approximate solvers, add constraint: λ ≥ λ_min(solver)"
   - "Investigate rate implications of λ = c·log(n)/n with c > 1"

3. **Use CV (as claimed):**
   - Actually implement CV selection
   - But constrain grid to computationally feasible range
   - E.g., λ ∈ [0.01, 0.2] instead of [(log n)/n, ...]

---

## Next Steps

### Immediate (When Simulation Completes)

1. **Analyze full results**
   - Coverage by lambda
   - Tree complexity by lambda
   - Correlation: solver warnings vs coverage

2. **Determine practical optimum**
   - Is λ = 0.05 better than λ = 0.1?
   - Is λ = 0.0038 better than λ = 0.05 (despite warnings)?

3. **Make recommendation**
   - If λ = 0.05 improves coverage: update simulations
   - If λ = 0.0038 needed: investigate solver tuning
   - If no improvement: look elsewhere for under-coverage cause

### Future Work

1. **Solver tuning for small lambda**
   - Increase model_limit (default 10,000 → 50,000?)
   - Add time_limit parameter
   - Test if longer search helps

2. **Modified theory lambda**
   - Test λ = 5·log(n)/n, λ = 10·log(n)/n
   - Find practical constant that balances theory and computation

3. **Sensitivity analysis**
   - Coverage as function of λ
   - Identify acceptable λ range

---

## Current Status

**Runtime:** 35+ minutes
**Progress:**
- λ = 0.1: Complete (100 reps, ~7 min)
- λ = 0.05: Complete (100 reps, ~7 min)
- λ = 0.0038: 50/100 reps (~25+ min so far, ~30 min remaining)

**Estimated completion:** ~15-20 more minutes

**Key finding already clear:** Theory lambda is computationally problematic (54+ solver warnings).
