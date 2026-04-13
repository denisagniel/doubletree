# Deployment Checklist

**Date:** 2026-04-13
**Project:** Functional Consistency Simulation Study
**Total replications:** 67,500 (135 configs × 500 reps)

---

## Pre-Deployment: Local Testing ✓

- [x] All packages load correctly
- [x] Standard M-split method works
- [x] Averaged tree method works
- [x] Pattern aggregation method works
- [x] Functional consistency metrics computed correctly
- [x] All three methods tested with same seed (results differ as expected)

**Test Results (n=200, K=2, M=10):**
```
              method att_est coverage max_diff_e max_diff_m0
     standard_msplit  0.6838     TRUE   1.51e-01    8.56e-02
       averaged_tree  0.6778     TRUE   0.00e+00    0.00e+00
 pattern_aggregation  0.6789     TRUE   0.00e+00    0.00e+00
```

✓ Standard M-split shows expected FC gap
✓ Both averaged tree and pattern aggregation achieve perfect FC
✓ All methods provide reasonable ATT estimates

---

## Deployment Steps

### Step 1: Commit and Push

```bash
cd ~/RAND/rprojects/global-scholars
git add doubletree/simulations/functional_consistency/
git commit -m "Add functional consistency simulation with three methods"
git push
```

### Step 2: Connect to O2

```bash
ssh username@o2.hms.harvard.edu
```

### Step 3: Pull Repository on O2

```bash
cd ~/global-scholars
git pull
```

### Step 4: Install Packages (if first time)

```bash
module load gcc/9.2.0 R/4.2.1
R CMD INSTALL optimaltrees
R CMD INSTALL doubletree
```

### Step 5: Test on O2

```bash
cd doubletree/simulations/functional_consistency/slurm
bash quick_test.sh
```

Expected: 3 successful replications with summary statistics

### Step 6: Launch Full Simulation

```bash
bash launch_all_simulations.sh
```

This will submit 135 SLURM jobs.

---

## Monitoring

### Regular Progress Checks

```bash
cd ~/global-scholars/doubletree/simulations/functional_consistency/slurm
bash check_progress.sh
```

Run every 30-60 minutes.

### Expected Timeline

- **Launch:** ~5 minutes (submit all jobs)
- **Queue time:** Variable (depends on O2 load)
- **Execution:** ~8-12 hours total
  - 135 jobs × ~30 min each
  - Run in parallel (depends on available nodes)
- **Completion:** Check when all 67,500 files present

---

## After Completion

### Step 1: Combine Results

```bash
cd ~/global-scholars/doubletree/simulations/functional_consistency/slurm
Rscript combine_results.R
```

Expected output:
- `results/combined_fc_simulations.rds`
- `results/combined_fc_simulations.csv`
- Summary statistics printed to console

### Step 2: Transfer Results

**Option A: Git (if small enough)**
```bash
git add results/combined_fc_simulations.csv
git commit -m "Add FC simulation results"
git push
```

**Option B: SCP**
```bash
# On local machine
scp username@o2.hms.harvard.edu:~/global-scholars/doubletree/simulations/functional_consistency/results/combined_fc_simulations.rds ./
```

---

## Critical Questions to Answer

### Q1: Perfect Functional Consistency

**Test:** max_diff ≈ 0 for averaged tree and pattern aggregation across all n, K, DGP?

**Expected:** Yes (both achieve perfect FC by construction)

### Q2: Coverage Degradation

**Test:** Does coverage stay at 95% as n increases from 200 → 3200?

**Three scenarios:**

**A. Coverage stays valid (surprising!):**
- Both averaged tree and pattern aggregation maintain 95% coverage at n=3200
- Implications: Bias vanishes fast enough, or partial cross-fitting sufficient
- Conclusion: **Perfect FC comes for free!**

**B. Coverage degrades for both:**
- Coverage drops below 90% at n=3200 for both methods
- Confirms bias = O(1/n_ℓ) or O(1/n_x) dominates
- Conclusion: **Trade-off exists: perfect FC vs valid inference**

**C. One succeeds, one fails:**
- Averaged tree maintains coverage, pattern aggregation doesn't (or vice versa)
- Suggests mechanism (direct vs indirect bias) or group size matters
- Conclusion: **Use the successful method**

### Q3: Method Comparison

**Test:** Compare averaged tree vs pattern aggregation coverage at n=3200

**If pattern aggregation < averaged tree:**
- Confirms n_x < n_ℓ matters
- Group size tuning helps
- Use averaged tree with appropriate regularization

**If pattern aggregation ≈ averaged tree:**
- Group size doesn't matter much
- Methods effectively equivalent
- Use simpler approach (pattern aggregation)

### Q4: Effect of K

**Test:** Within averaged tree/pattern aggregation, does K matter?

**Expected:**
- Larger K → more cross-fitting → better coverage
- Effect should be stronger for pattern aggregation (relies on cross-fit predictions)

### Q5: DGP Robustness

**Test:** Do results hold across simple/complex/sparse DGPs?

**Expected:**
- Sparse DGP (small n_x) hardest for pattern aggregation
- Should see largest gap between methods in sparse DGP

---

## Analysis Plan

After results are combined, key analyses:

### 1. Coverage by Method and n

```r
library(ggplot2)
coverage_summary <- aggregate(coverage ~ method + n, data = results, mean)

ggplot(coverage_summary, aes(x = n, y = coverage, color = method)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  labs(title = "Coverage vs Sample Size",
       y = "Coverage", x = "Sample Size") +
  theme_minimal()
```

### 2. Functional Consistency

```r
fc_summary <- aggregate(cbind(max_diff_e, max_diff_m0) ~ method,
                        data = results,
                        FUN = function(x) c(mean=mean(x), max=max(x)))
```

### 3. Standardized Bias

```r
ggplot(results, aes(x = n, y = abs(standardized_bias), color = method)) +
  geom_point(alpha = 0.1) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Standardized Bias vs Sample Size",
       y = "|√n × bias|", x = "Sample Size (log scale)")
```

### 4. Method Comparison at n=3200

```r
large_n <- subset(results, n == 3200)
aggregate(cbind(coverage, bias, standardized_bias) ~ method + K + dgp,
          data = large_n, mean)
```

---

## Success Criteria

**Minimum success:**
- All 67,500 replications complete
- No systematic errors in any configuration
- Can answer Q1 (perfect FC achieved?)

**Full success:**
- Can definitively answer Q2-Q5
- Clear recommendation on which method to use when
- Understanding of bias-FC trade-off

---

## Troubleshooting

### Jobs fail immediately

**Check:**
1. Packages installed? `Rscript -e "library(optimaltrees); library(doubletree)"`
2. Paths correct? `ls ~/global-scholars/doubletree/simulations/functional_consistency/`
3. SLURM logs? `cat slurm/logs/*.err`

### Jobs stuck pending

**Check:**
- Queue status: `squeue -u $USER`
- Partition availability: `sinfo`
- Consider switching to medium/long partition if short is full

### Incomplete results

**Check:**
- Which configs missing? `Rscript slurm/combine_results.R` shows incomplete
- Re-run specific configs by modifying `launch_all_simulations.sh`

---

## Timeline Estimate

- **Monday:** Deploy to O2 (morning)
- **Tuesday:** Jobs complete, combine results (morning)
- **Tuesday-Wednesday:** Analysis and interpretation
- **Thursday:** Write up findings, update documentation

---

## Contact

For O2 issues: https://rc.hms.harvard.edu/

For simulation questions: See `METHODS_COMPARISON.md` and `SIMULATION_SUMMARY.md`
