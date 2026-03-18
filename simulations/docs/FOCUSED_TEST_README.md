# Focused Test: Getting Theory to Work

**Goal:** Validate theory under favorable conditions before expanding

**Strategy:** Large n + smooth DGP + theory-guided epsilon

---

## Design

### DGP
- **DGP3 only:** Smooth, high signal-to-noise ratio
- Most favorable for approximation theory
- Binary features X1, X2 (median-split from continuous)

### Sample Sizes
- n ∈ {800, 1600, 3200}
- Large enough for theory to kick in
- β > d/2 condition more likely to hold

### Epsilon (Theory-Guided)
Formula: **ε_n = c√(log n/n)**

From /tmp/epsilon_guidance_summary.md:
- Satisfies ε_n = o(n^{-1/2}) deterministically
- Conservative relative to empirical process noise
- c ∈ {1, 2, 3} provides reasonable range

**Computed values:**
- n = 800, c = 1: ε = 0.0914
- n = 800, c = 2: ε = 0.183
- n = 800, c = 3: ε = 0.274
- n = 1600, c = 1: ε = 0.068
- n = 1600, c = 2: ε = 0.136
- n = 1600, c = 3: ε = 0.204
- n = 3200, c = 1: ε = 0.048
- n = 3200, c = 2: ε = 0.096
- n = 3200, c = 3: ε = 0.144

Note: Epsilon decreases with n (as theory predicts)

### Grid
- 1 DGP × 3 n × 3 c × 100 reps = **900 runs**
- Estimated time: ~45 minutes

---

## What This Tests

1. **Coverage:** Does 95% CI achieve ~95%?
   - Oracle should work (validates DML framework)
   - Fold-specific should be close
   - Rashomon: target 90-95%

2. **Intersection Existence:**
   - With smooth DGP, expect high success (85-95%)
   - Should improve with larger n

3. **Epsilon Sensitivity:**
   - c = 1 (tight): Lower bias, may have empty intersections
   - c = 2 (moderate): Balance bias vs intersection success
   - c = 3 (loose): Higher bias, more successful intersections

4. **Sample Size Effects:**
   - Coverage should improve with n
   - Epsilon automatically tightens with n (by design)
   - RMSE should decrease

---

## Success Criteria

**Minimum (proceed to full grid):**
- Oracle: 93-97% coverage
- Rashomon: ≥90% coverage (at least for c=1 or c=2)
- Intersection: ≥70% success rate

**Ideal:**
- Rashomon: 93-97% coverage
- Intersection: ≥85% success rate
- Low bias relative to fold-specific

---

## Execution

```r
# Run simulations (~45 min)
source("simulations/run_focused_test.R")

# Analyze results (~1 min)
source("simulations/analyze_focused_results.R")
```

---

## If It Works

✓ Theory validated under favorable conditions
✓ Ready to expand to full grid:
  - Add DGP 1-2 (binary features, varying confounding)
  - Maybe add DGP 4 (rough) as stress test
  - Keep large n (800, 1600)
  - Use theory-guided epsilon

---

## If It Doesn't Work

Debug in order:
1. Check Rashomon metadata (pct_nonempty still NA?)
2. Tune tree complexity (lambda parameter)
3. Try even larger n (4000, 8000)
4. Check DGP3 implementation (is it really smooth?)
5. Verify cross-fitting implementation

---

## Advantages of This Approach

1. **Fast iteration:** 45 min vs 4-6 hours for full grid
2. **Clear signal:** Smooth DGP should work if theory is sound
3. **Theory-aligned:** Epsilon choice follows manuscript guidance
4. **Focused debugging:** One DGP = easier to diagnose issues
5. **Expandable:** Can add DGPs/conditions once working
