# Strategic Subset Simulation (Option B)

**Purpose:** Fast iteration to validate core claims before full grid

**Scope:** 2 DGPs × 3 n × 2 ε × 100 reps = **1,200 runs** (~1 hour)

---

## What This Tests

### DGPs (2 total - extremes of spectrum)
- **DGP 3:** Smooth, high SNR (ideal case for Rashomon)
- **DGP 4:** Non-smooth, low SNR (stress test)

### Sample Sizes (3 total)
- n ∈ {200, 400, 800}
- Skip n = 1,600 for now (can add in full grid)

### Tolerances (2 total)
- ε ∈ {0.05, 0.1}
- Skip extreme values (0.01, 0.2) for now

### Methods (3 total)
- Fold-specific optimal trees (baseline)
- Rashomon intersection (interpretable)
- Oracle (true nuisances - performance ceiling)

---

## Execution

### Run simulations
```r
source("simulations/run_strategic_subset.R")
```

**Time:** ~1 hour
**Output:** `simulations/results_strategic/result_*.rds`

### Analyze results
```r
source("simulations/analyze_strategic_results.R")
```

**Time:** ~1 minute
**Output:**
- `simulations/figures_strategic/` - 3 PDF figures
- `simulations/results_strategic/summary_strategic.csv` - detailed results

---

## What We'll Learn

### Figure 1: Coverage
- Does 95% CI achieve ~95% coverage?
- Validates DML theory (o_p(n^{-1/4}) rate)

### Figure 2: Intersection Existence
- What % of replications have non-empty intersection?
- Expected: DGP 3 (smooth) > DGP 4 (rough)
- Target: 80-90% success rate

### Figure 3: Bias Comparison
- How much overhead does Rashomon incur?
- Is interpretability-validity tradeoff acceptable?

---

## Decision Points After Results

### If results look good:
→ Proceed to full grid (Option A)
  - Add DGP 1-2 (binary features)
  - Add n = 1600
  - Add extreme tolerances (0.01, 0.2)
  - Total: 6,400 runs (~4-6 hours)

### If issues found:
→ Iterate on design
  - Adjust DGPs if needed
  - Tune optimaltrees parameters
  - Revise Rashomon tolerance range
  - Re-run strategic subset

### If coverage issues:
→ Check:
  - Nuisance estimation quality
  - Cross-fitting implementation
  - Variance estimation

### If intersection rate too low:
→ Consider:
  - Wider tolerances
  - Different DGP designs
  - Looser complexity constraints

---

## Notes

- Strategic subset deliberately tests extremes (smooth vs rough)
- If both work well → full grid should be fine
- If DGP 4 (rough) struggles → may need design refinement
- Iteration is fast (~1 hr per round)
