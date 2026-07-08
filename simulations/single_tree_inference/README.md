# Single-Tree Inference Simulation Studies

**Created:** 2026-04-29
**Purpose:** Determine if single full-sample tree can be used for interpretation with valid inference
**Status:** Setup complete, ready to run

---

## Quick Start

**To run pilot studies (quick test):**
```r
source("code/study1_similarity/run_study1_pilot.R")
source("code/study2_coverage/run_study2_pilot.R")
```

**To run full studies:**
```r
source("code/study1_similarity/run_study1.R")    # ~3-4 hours
source("code/study2_coverage/run_study2.R")      # ~4-5 hours
```

**To analyze results:**
```r
source("code/study1_similarity/analyze_study1.R")
source("code/study2_coverage/analyze_study2.R")
```

**Results:** Check `results/latest/` symlinks for most recent runs

---

## Overview

We need to decide: **Can we use a single full-sample tree for interpretation while maintaining valid inference?**

Two approaches being tested:

### Approach A: Show Full-Sample Tree, Report Cross-Fitted Estimates
- **Idea:** Fit one tree on all data for visualization, but use cross-fitted estimates for inference
- **Requires:** Full-sample tree must be similar to cross-fitted trees (RMSE < 5%)
- **Tested by:** Study 1 (Similarity)

### Approach B: Report Full-Sample Estimate with Bias-Adjusted CIs
- **Idea:** Report estimate from full-sample tree, but inflate CIs to account for overfitting bias
- **Requires:** Bias-adjusted CIs must achieve 95% coverage without being too wide
- **Tested by:** Study 2 (Coverage)

---

## Directory Structure

```
single_tree_inference/
├── README.md                 # This file
├── SPEC.md -> ../../quality_reports/specs/2026-04-29_single-tree-inference-simulation.md
│
├── code/
│   ├── dgps.R                # Data-generating processes (simple, moderate, complex)
│   ├── estimators.R          # Full-sample and cross-fitted ATT estimators
│   ├── metrics.R             # Compute similarity and coverage metrics
│   ├── utils.R               # Plotting and summary helpers
│   │
│   ├── study1_similarity/    # Study 1: Are full-sample trees similar to cross-fitted?
│   │   ├── run_study1.R      # Full study (500 reps × 9 settings)
│   │   ├── run_study1_pilot.R # Quick test (50 reps × 2 settings)
│   │   └── analyze_study1.R   # Generate tables and plots
│   │
│   └── study2_coverage/      # Study 2: Do bias-adjusted CIs maintain coverage?
│       ├── run_study2.R      # Full study (500 reps × 9 settings)
│       ├── run_study2_pilot.R # Quick test (50 reps × 2 settings)
│       └── analyze_study2.R   # Generate tables and plots
│
├── results/
│   ├── study1_similarity/
│   │   ├── YYYYMMDD_HHMM/    # Timestamped run
│   │   │   ├── results.rds   # Full results object
│   │   │   ├── summary.csv   # Summary table
│   │   │   ├── plots/        # Diagnostic plots
│   │   │   └── run_info.txt  # Timestamp, system info, elapsed time
│   │   └── latest -> YYYYMMDD_HHMM/  # Symlink to most recent run
│   │
│   ├── study2_coverage/
│   │   ├── YYYYMMDD_HHMM/
│   │   │   ├── results.rds
│   │   │   ├── summary.csv
│   │   │   ├── plots/
│   │   │   └── run_info.txt
│   │   └── latest -> YYYYMMDD_HHMM/
│   │
│   └── DECISION_SINGLE_TREE_INFERENCE.md  # Final decision (created after both studies)
│
└── archive/                  # Old runs (when re-running)
```

---

## Study 1: Similarity Test

**Research Question:** Are full-sample trees similar enough to cross-fitted trees that we can use them for interpretation?

**Design:**
- 3 DGPs (simple, moderate, complex)
- 3 sample sizes (500, 1000, 2000)
- 500 replications per setting
- Compare structure and leaf values

**Metrics:**
- Structure match rate (% exact structure equality)
- Leaf RMSE: √mean[(μ_full - μ_crossfit)²]
- Max absolute difference in leaf values

**Success Criterion:**
- RMSE < 5% in ≥80% of replications (for n ≥ 1000)
- Structure match rate ≥ 70% (for simple DGP)

**If successful:** Approach A is viable

---

## Study 2: Coverage Test

**Research Question:** Can we report full-sample tree estimates with bias-adjusted CIs that maintain 95% coverage?

**Design:**
- 3 DGPs (simple, moderate, complex)
- 3 sample sizes (500, 1000, 2000)
- 500 replications per setting
- Test 3 CI adjustment methods

**CI Methods Tested:**
1. **Standard (naive):** θ̂_full ± 1.96 × SE
2. **Additive:** θ̂_full ± (1.96 × SE + |B̂|)
3. **Conservative:** θ̂_full ± (1.96 × SE + 2 × |B̂|)

Where B̂ = θ̂_full - θ̂_crossfit (estimated bias)

**Metrics:**
- Empirical coverage (should be ≥ 95%)
- Average CI width
- Width ratio vs standard cross-fitted CI
- Bias magnitude: |B̂|/SE

**Success Criterion:**
- At least one method achieves coverage ≥ 95%
- Successful method has width ratio < 2.0
- Bias |B̂| < 1.5 × SE in ≥80% of replications

**If successful:** Approach B is viable

---

## DGP Specifications

All DGPs have:
- Binary covariates: X₁, ..., X₅ ~ Bernoulli(0.5)
- Binary treatment: A ~ Bernoulli(e(X))
- Binary outcome: Y ~ Bernoulli(μ(X, A))
- True ATT = 0.15

**Simple DGP (Few interactions, ~2-3 split trees optimal):**
```
e(X) = expit(-0.5 + 0.3×X₁ + 0.3×X₂)
μ₀(X) = 0.2 + 0.15×X₁ + 0.15×X₃
```

**Moderate DGP (Some interactions, ~4-5 split trees):**
```
e(X) = expit(-0.5 + 0.3×X₁ + 0.2×X₂ + 0.3×X₁×X₂)
μ₀(X) = 0.2 + 0.2×X₃ + 0.15×X₄ + 0.2×X₃×X₄
```

**Complex DGP (Many interactions, ~6-8 split trees):**
```
e(X) = expit(-0.5 + 0.2×(X₁+X₂+X₃) + 0.3×X₁×X₂ + 0.2×X₂×X₃)
μ₀(X) = 0.2 + 0.15×(X₃+X₄+X₅) + 0.2×X₃×X₄ + 0.15×X₄×X₅
```

---

## Expected Timeline

**Phase 1: Setup and Code (4-5 hours)**
- ✓ Create directory structure (done)
- ✓ Write specification (done)
- Write DGP functions (1 hour)
- Write estimator functions (1 hour)
- Write metrics and utils (1 hour)
- Write run scripts (1 hour)
- Test with pilots (1 hour)

**Phase 2: Pilot Runs (1 hour)**
- Run study1_pilot.R (10 min)
- Run study2_pilot.R (15 min)
- Review results, adjust if needed
- Debug any issues

**Phase 3: Full Runs (7-9 hours, can run overnight)**
- Study 1: ~3-4 hours (9 settings × 500 reps)
- Study 2: ~4-5 hours (9 settings × 500 reps)
- Can run in parallel if multiple cores available

**Phase 4: Analysis (3 hours)**
- Generate summary tables (30 min)
- Create diagnostic plots (1 hour)
- Interpret results (1 hour)
- Write decision document (30 min)

**Total:** 15-18 hours over 2-3 days

---

## Outputs

After running both studies, you will have:

### Study 1 Outputs
1. **results.rds:** Full simulation results (all reps, all settings)
2. **summary.csv:** RMSE mean/median/p90, structure match rates by DGP and n
3. **Plots:**
   - `rmse_by_dgp_n.png`: Distribution of RMSE across settings
   - `structure_match.png`: Structure match rates by DGP
   - `scatter_example.png`: Full vs cross-fitted predictions (example)

### Study 2 Outputs
1. **results.rds:** Full simulation results
2. **summary.csv:** Coverage, width ratio, bias/SE by CI method, DGP, n
3. **Plots:**
   - `coverage_by_method_dgp.png`: Coverage heatmap
   - `width_ratio.png`: Width ratio by method and setting
   - `bias_distribution.png`: Distribution of |B̂|/SE

### Decision Document
`DECISION_SINGLE_TREE_INFERENCE.md` synthesizing both studies with clear recommendation:
- Which approach(es) work
- Under what conditions
- Recommended strategy for doubletree package
- Manuscript framing implications

---

## Interpretation Guide

### If Study 1 Succeeds (RMSE < 5%)
→ **Approach A viable:** Show full-sample tree, report cross-fitted estimates
→ **Manuscript:** "Full-sample tree shown for visualization; cross-fitted estimates used for inference (trees are similar, RMSE < 5%)"
→ **Package:** Add function to visualize full-sample tree

### If Study 2 Succeeds (Coverage ≥ 95%, Width < 2×)
→ **Approach B viable:** Report full-sample estimate with bias-adjusted CIs
→ **Manuscript:** "Single-tree estimate with bias-corrected confidence intervals"
→ **Package:** Add option for full-sample estimation with CI adjustment

### If Both Succeed
→ **Choose based on user preference:**
- Approach A: Simpler (show different tree than used)
- Approach B: More honest (report what you used)
→ **Recommend:** Probably Approach A (simpler communication)

### If Neither Succeeds
→ **Stick with current approach:** "One structure, cross-fitted leaf values"
→ **Manuscript:** Honest framing about K refits being necessary
→ **Value:** Validated that alternatives don't work

---

## How to Re-Run

If you need to re-run studies (e.g., after code changes):

```r
# Move old results to archive
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
system(paste0("mv results/study1_similarity/latest archive/", timestamp, "_study1"))
system(paste0("mv results/study2_coverage/latest archive/", timestamp, "_study2"))

# Re-run
source("code/study1_similarity/run_study1.R")
source("code/study2_coverage/run_study2.R")

# Re-analyze
source("code/study1_similarity/analyze_study1.R")
source("code/study2_coverage/analyze_study2.R")
```

Results will be in new timestamped directories with updated `latest` symlinks.

---

## Troubleshooting

**If pilot runs fail:**
- Check that optimaltrees and doubletree packages are installed
- Verify paths are correct (run from project root)
- Check for package function changes

**If full runs are too slow:**
- Reduce replications from 500 to 300 (still sufficient)
- Run overnight or on cluster
- Use parallel::mclapply for multi-core

**If results are inconclusive:**
- Try larger sample sizes (n = 5000, 10000)
- Try different tree regularization (λ)
- Check if specific DGPs drive results

---

## Contact

**Questions or issues:** See SPEC.md or session notes for rationale and design decisions

**Related documents:**
- Full specification: `quality_reports/specs/2026-04-29_single-tree-inference-simulation.md`
- Decision framework: `quality_reports/2026-04-29_averaging-msplit-decision-framework.md`
- Session notes: `session_notes/2026-04-29.md`

---

## Status

- [✓] Directory structure created
- [✓] Specification written and approved
- [ ] DGP functions implemented
- [ ] Estimator functions implemented
- [ ] Metrics functions implemented
- [ ] Pilot runs completed
- [ ] Full Study 1 run completed
- [ ] Full Study 2 run completed
- [ ] Results analyzed
- [ ] Decision document written

**Next:** Implement code functions and run pilots
