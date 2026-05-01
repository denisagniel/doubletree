# Six-Approach Comparison Study

**Created:** 2026-05-01
**Purpose:** Systematic comparison of 6 tree-based causal inference approaches
**Status:** Ready to run on cluster

---

## Quick Start

### On Cluster

```bash
# 1. Navigate to directory
cd doubletree/simulations/six_approach_comparison

# 2. Launch all 90 jobs
bash slurm/launch_all.sh

# 3. Monitor progress
bash slurm/check_progress.sh

# 4. When complete, combine results
Rscript code/combine_results.R

# 5. Analyze
Rscript code/analyze_results.R
```

---

## The Six Approaches

| # | Name | Structure From | Leaf Values From | Interpretation |
|---|------|---------------|------------------|----------------|
| **(i)** | Full-sample | All n obs | Single fit (all n) | "One tree" (simplest) |
| **(ii)** | Cross-fit separate | Per fold (no intersection) | Cross-fitted (K fits) | "K trees" (standard) |
| **(iii)** | Doubletree | Rashomon intersection (K folds) | Cross-fitted (K refits) | "One structure, K values" (current) |
| **(iv)** | Doubletree + single fit | Rashomon intersection | Single fit (all n) | "One tree + stable structure" |
| **(v)** | M-split | Modal across M splits | M-split averaged (M×K fits) | "One structure + stability" |
| **(vi)** | M-split + single fit | Modal across M splits | Single fit (all n) | "One tree + stability evidence" |

---

## Study Design

**DGPs:** 3 (simple, moderate, complex) + 1 continuous covariate DGP
- Simple: Linear functions (~2-3 optimal splits)
- Moderate: 2-way interactions (~4-5 splits)
- Complex: 3-way interactions (~6-8 splits)
- Continuous: Mixed binary + continuous covariates

**Sample sizes:** n ∈ {500, 1000, 2000}

**Replications:** 500 per (approach × DGP × n)

**Total:** 6 approaches × 4 DGPs × 3 n × 500 reps = 36,000 estimations

---

## Cluster Job Structure

**90 jobs split across 3 arrays:**

### Array 1: Fast Approaches (i, iv, vi) - 27 jobs
- 3 approaches × 3 DGPs × 3 n = 27 jobs
- 500 reps per job
- ~0.7-1.3 hours per job

### Array 2: Medium Approaches (ii, iii) - 18 jobs
- 2 approaches × 3 DGPs × 3 n = 18 jobs
- 500 reps per job
- ~1.7-2.1 hours per job

### Array 3: M-Split Approach (v) - 45 jobs
- 3 DGPs × 3 n × 5 batches = 45 jobs
- 100 reps per job (5 batches of 100 = 500 total)
- ~1.7-2.0 hours per job

**Wall time:**
- 90 cores: ~2 hours
- 45 cores: ~4 hours
- 30 cores: ~6 hours

---

## Directory Structure

```
six_approach_comparison/
├── README.md                       # This file
├── SPEC.md -> quality_reports/specs/...  # Full specification
│
├── code/
│   ├── dgps.R                      # 4 DGPs
│   ├── estimators.R                # 6 estimator functions
│   ├── metrics.R                   # Coverage, similarity metrics
│   ├── run_single_replication.R    # Main simulation worker
│   ├── combine_results.R           # Merge 90 job outputs
│   └── analyze_results.R           # TODO: Generate plots and tables
│
├── slurm/
│   ├── run_fast_approaches.sh      # Array 1
│   ├── run_medium_approaches.sh    # Array 2
│   ├── run_msplit_approach.sh      # Array 3
│   ├── launch_all.sh               # Submit all arrays
│   └── check_progress.sh           # Monitor completion
│
├── results/
│   ├── raw/                        # 90 individual job outputs
│   ├── combined/                   # Merged results
│   │   ├── all_results.rds
│   │   ├── summary_inference.csv
│   │   └── summary_timing.csv
│   └── plots/                      # Diagnostic plots
│
└── logs/                           # SLURM stdout/stderr
```

---

## Key Research Questions

### 1. Inference Validity
- Which approaches maintain valid inference (bias ≈ 0, coverage ≈ 95%)?
- Do single-fit approaches (i, iv, vi) have substantial bias?
- Can bias be corrected with adjusted CIs?

### 2. Interpretability
- How similar are structures across approaches?
- Can we show approach X but report approach Y?
- Is "one tree" feasible with valid inference?

### 3. Efficiency
- Computational cost tradeoffs?
- CI width comparison?
- Is M-split worth the extra computation?

### 4. Stability
- How stable are structures (for v, vi)?
- Does M-split reduce variance vs doubletree?
- When does Rashomon effect dominate?

---

## Expected Outcomes

### Scenario 1: Approach (i) has negligible bias
→ Simplest approach works! Recommend (i).

### Scenario 2: Approaches (i, iv, vi) have bias, but similar structures
→ Show simple tree, report valid estimates from (ii) or (iii).

### Scenario 3: Only cross-fitted approaches (ii, iii, v) are valid
→ Recommend (iii) or (v) based on stability/cost tradeoff.

### Scenario 4: M-split (v) provides clear stability benefits
→ Recommend (v) as enhanced doubletree.

---

## Monitoring

### Check progress
```bash
bash slurm/check_progress.sh
```

Shows:
- Completed jobs (out of 90)
- Error count
- Running jobs
- Estimated completion

### Check specific job
```bash
# View output
cat logs/fast_1.out

# View errors
cat logs/fast_1.err

# Check all jobs for errors
grep -i error logs/*.err
```

### Resubmit failed jobs
If some jobs fail, resubmit specific array indices:
```bash
# Example: Resubmit fast jobs 5-8
sbatch --array=5-8 slurm/run_fast_approaches.sh
```

---

## Analysis Pipeline

After all jobs complete:

1. **Combine results** (5 min)
```bash
Rscript code/combine_results.R
```

2. **Analyze** (TBD - need to implement)
```bash
Rscript code/analyze_results.R
```

3. **Review outputs**
- `results/combined/summary_inference.csv` - Bias, RMSE, coverage by setting
- `results/combined/summary_timing.csv` - Computational costs
- `results/plots/` - Diagnostic plots

---

## Troubleshooting

### Jobs failing immediately
- Check module availability: `module avail R`
- Check paths in SLURM scripts
- Test locally first: `Rscript code/run_single_replication.R --approach 1 --dgp 1 --n 500 --reps 1 --output test.rds`

### Out of memory errors
- Increase --mem in SLURM scripts
- M-split jobs may need more (currently 12G)

### Long queue times
- Check partition availability
- Consider different partition (change --partition)
- Split into smaller batches

### Missing results files
- Check logs for errors
- Resubmit failed jobs manually
- Verify output directory exists

---

## Timeline

**Day 1 (Complete):** Implementation
- ✓ DGP functions
- ✓ 6 estimator functions
- ✓ SLURM scripts
- ✓ Combine script

**Day 2 (Next):** Local testing
- Test all 6 approaches with 1 rep
- Fix any bugs
- Verify outputs

**Day 3:** Cluster launch
- Submit 90 jobs
- Monitor progress
- ~2-6 hours depending on cores

**Day 4:** Analysis
- Combine results
- Generate plots
- Write decision document

---

## Files Status

- [✓] DGPs implemented (4 DGPs)
- [✓] Estimators implemented (6 approaches)
- [✓] Metrics implemented
- [✓] Main worker script
- [✓] SLURM scripts (3 arrays)
- [✓] Launch and monitoring scripts
- [✓] Combine results script
- [ ] Analysis script (TODO)
- [ ] Local testing
- [ ] Cluster run
- [ ] Results analysis

---

## Contact

**Questions:** See SPEC.md for detailed design rationale
**Issues:** Check logs/ directory and session notes
**Related:** See quality_reports/specs/2026-05-01_six-approach-comparison.md

---

**Next step:** Local testing before cluster launch
```bash
# Test one rep of each approach
for i in {1..6}; do
  Rscript code/run_single_replication.R --approach $i --dgp 1 --n 500 --reps 1 --output test_approach_${i}.rds
done
```
