# Deployment Instructions: Six-Approach Comparison on O2 Cluster

**Date:** 2026-06-09
**Status:** Ready to deploy — packages updated, parameter audit complete

---

## Changes Since Last Run (critical — must reinstall before running)

### doubletree (HEAD: latest — see git log)

⚠️ **CI formula bug fixed in approaches 4, 5, 6** — old results INVALID (22x narrow CIs)
- `estimate_att_averaged.R`: approaches 4 and 6 now use `att_ci(theta, sigma)` (was `sigma/sqrt(n)`)
- `estimate_att_msplit.R`: approach 5 same fix
- `dgps.R`: DGP_complex intercept corrected (0.2 → 0.05) to prevent mu0+ATT clipping

### simulation scripts (2026-06-09)

⚠️ **Parameter audit — all 6 approaches now use uniform adaptive CV**
- `code/estimators.R`: approaches 1 and 2 now use `cv_regularization_adaptive` with `max_lambda = 15*log(n)/n` (was fixed-grid `cv_regularization`)
- `code/estimators.R`: approach 3 now uses `auto_tune_intersecting = TRUE`; hard `stop()` on failed intersection (was fallback to fold-specific trees)
- `code/estimators.R`: production parameter constants block added (SIM_K, SIM_M, SIM_K_MSPLIT, EPS_N_C)
- SLURM: full redesign for 1000 reps; 6 per-approach scripts (`run_approach{1-6}.sh`) replacing 3-script grouping
- SLURM: total jobs increased from 120 to 6624; total reps from 6000 to 72,000

### optimaltrees (HEAD: fe357e1)

✅ **Package install fixed** (R CMD INSTALL now succeeds)
- NAMESPACE stale exports removed (export("for"), export(Plot), export(optimaltrees_model))
- `auto_tune_regularization_for_intersection`: lambda_min=NULL (auto-computes 0.5*log(n)/n)

---

## Pre-Deployment Steps (On Cluster)

### 1. SSH to O2 Cluster

```bash
ssh yourusername@o2.hms.harvard.edu
```

### 2. Load Modules

```bash
module load gcc/14.2.0 R/4.4.2
```

**Verify:**
```bash
module list
# Should show: gcc/14.2.0, R/4.4.2
```

### 3. Navigate to Project Directory

```bash
cd /path/to/your/global-scholars
```

### 4. Pull Latest optimaltrees Changes

```bash
cd optimaltrees
git pull origin main
```

**Expected output:**
```
From github.com:denisagniel/treefarmr
 * branch            main       -> FETCH_HEAD
   961b4cc..94e79dc  main       -> origin/main
Updating 961b4cc..94e79dc
Fast-forward
 14 files changed, 635 insertions(+), 102 deletions(-)
```

**Verify commit:**
```bash
git log -1 --oneline
# Should show: 94e79dc Remove obsolete man file
```

### 5. Reinstall optimaltrees

```bash
# From optimaltrees directory
R CMD INSTALL . --preclean
```

**Expected:** Package compiles successfully (C++ compilation may take 2-3 minutes)

**Verify installation:**
```bash
R --quiet -e "cat('optimaltrees version:', as.character(packageVersion('optimaltrees')), '\n')"
# Should show: optimaltrees version: 0.4.0
```

### 6. Reinstall doubletree

```bash
cd ../doubletree
R CMD INSTALL . --preclean
```

**Verify installation:**
```bash
R --quiet -e "library(doubletree); library(optimaltrees); cat('Both packages loaded successfully\n')"
# Should load without errors
```

### 7. Test Partition-Based Comparison

```bash
R --quiet --vanilla <<'EOF'
library(optimaltrees)

# Quick test: fit two trees, check structure comparison
set.seed(123)
X <- data.frame(x1 = rbinom(100, 1, 0.5), x2 = rbinom(100, 1, 0.5))
y <- rbinom(100, 1, 0.5)

m1 <- fit_tree(X, y, loss_function = "log_loss", seed = 1)
m2 <- fit_tree(X, y, loss_function = "log_loss", seed = 2)

s1 <- extract_tree_structure(m1)
s2 <- extract_tree_structure(m2)

# Check that partition_hash property exists
if (is.null(s1@partition_hash) || is.null(s2@partition_hash)) {
  stop("ERROR: partition_hash not found")
} else {
  cat("✓ Partition-based comparison active\n")
  cat("  Hash 1:", s1@partition_hash, "\n")
  cat("  Hash 2:", s2@partition_hash, "\n")
  cat("  Same structure:", compare_structures(s1, s2), "\n")
}
EOF
```

**Expected output:**
```
✓ Partition-based comparison active
  Hash 1: <16-char hash>
  Hash 2: <16-char hash>
  Same structure: TRUE or FALSE
```

---

## Deployment

### 8. Navigate to Simulation Directory

```bash
cd simulations/six_approach_comparison
```

### 9. Check Directory Structure

```bash
ls -l
# Should see: code/, slurm/, results/, logs/
```

### 10. Create Output Directories

```bash
mkdir -p results/raw logs results/combined results/plots
```

### 11. Test Single Job Locally (Recommended)

```bash
Rscript code/run_single_replication.R \
  --approach 1 \
  --dgp 1 \
  --n 500 \
  --reps 10 \
  --output results/test_local.rds
```

**Expected:** Completes successfully, creates `results/test_local.rds`

**Check results:**
```bash
R --quiet -e "readRDS('results/test_local.rds') |> str(max.level=1)"
# Should show list of 10 replication results
```

### 12. Launch All Jobs

```bash
bash slurm/launch_all.sh
```

**Expected output:**
```
==============================================
Six-Approach Comparison Study
==============================================

Total jobs: 120
  Array 1 (Fast):   36 jobs (approaches i, iv, vi)
  Array 2 (Medium): 24 jobs (approaches ii, iii)
  Array 3 (M-split): 60 jobs (approach v)

Submitting Array 1: Fast approaches (i, iv, vi) - 36 jobs...
  Job ID: 12345678

Submitting Array 2: Medium approaches (ii, iii) - 24 jobs...
  Job ID: 12345679

Submitting Array 3: M-split approach (v) - 60 jobs...
  Job ID: 12345680

==============================================
All jobs submitted successfully!
==============================================

Monitor progress with:
  bash slurm/check_progress.sh
```

---

## Monitoring

### Check Job Status

```bash
# Quick check
squeue -u $USER

# Detailed progress
bash slurm/check_progress.sh
```

### Check Logs

```bash
# Recent errors (per approach)
tail logs/approach1_*.err
tail logs/approach3_*.err

# Recent outputs
tail logs/approach1_*.out
```

### Check Results

```bash
# Count completed jobs
ls results/raw/*.rds | wc -l
# Should reach 6624 when all complete
```

---

## Expected Runtime

Per-rep timing estimates (measured smoke test + cluster scaling):

| Approach | DGP | n | s/rep local | s/rep cluster est. | Reps/batch | Max batch time |
|----------|-----|---|-------------|-------------------|-----------|---------------|
| 1 | all | 500 | 0.8 | 0.4 | 500 | 3 min |
| 1 | all | 2000 | 6 | 3 | 500 | 25 min |
| 2 | all | 2000 | 76 | 38 | 10 | 6 min |
| 3 | 1-3 | 2000 | ~480 | ~240 | 10 | 40 min |
| 3 | 4 | 500 | ~8000 (fail) | ~4000 | 1 | 1.1 h |
| 3 | 4 | 2000 | ~32000 (fail) | ~16000 | 1 | 4.4 h |
| 4 | 1-3 | 2000 | ~110 | ~55 | 10 | 9 min |
| 4 | 4 | 500 | ~6500 (fail) | ~3250 | 1 | 0.9 h |
| 4 | 4 | 2000 | ~26000 (fail) | ~13000 | 1 | 3.6 h |
| 5 | all | 2000 | ~1320 | ~660 | 5 | 55 min |
| 6 | all | 2000 | ~208 | ~104 | 20 | 35 min |

**Key:** Approaches 3/4 on DGP4 exhaust all auto-tune tiers (~50 GOSDT calls × 2 nuisances) before hitting the hard `stop()` on failed intersection. Times marked `(fail)` are for the failure path. DGPs 1-3 usually succeed at Tier 1 and are much faster.

**Wall time budget:**
- DGPs 1-3, 6h: all approaches comfortably within budget
- DGP4, approaches 3/4: 1 rep/batch with 8h wall covers worst case (n=2000, ~4-5h cluster)

---

## Troubleshooting

### Package Loading Errors

```bash
# Check package versions
R -e "cat('optimaltrees:', as.character(packageVersion('optimaltrees')), '\n')"
R -e "cat('doubletree:', as.character(packageVersion('doubletree')), '\n')"

# Reinstall if versions don't match
cd /path/to/optimaltrees && R CMD INSTALL . --preclean
cd /path/to/doubletree && R CMD INSTALL . --preclean
```

### Job Failures

```bash
# Check error logs
grep -i "error\|fatal" logs/*.err

# Resubmit failed jobs
# (Manual resubmission: adjust array indices in slurm scripts)
```

### Missing partition_hash

If you see `Can't find property @partition_hash`:
- optimaltrees not properly installed
- Using old cached version
- Reinstall: `R CMD INSTALL optimaltrees --preclean`

---

## Post-Processing

### After All Jobs Complete

```bash
# Combine results
Rscript slurm/combine_results.R

# Check combined file
ls -lh results/combined/all_results.rds
```

### Download Results

```bash
# From your local machine
scp -r yourusername@o2.hms.harvard.edu:/path/to/six_approach_comparison/results ./
```

---

## Verification Checklist

Before launching:
- [ ] SSH to O2 cluster
- [ ] Modules loaded (gcc/14.2.0, R/4.4.2)
- [ ] optimaltrees updated and reinstalled
- [ ] doubletree updated and reinstalled
- [ ] Partition-based comparison tested and working
- [ ] All 6 approaches use `cv_regularization_adaptive` (verify: `grep -n "cv_regularization[^_]" code/estimators.R` returns zero matches in function bodies)
- [ ] Approach 3 uses `auto_tune_intersecting = TRUE` with no fallback (verify: `grep -n "auto_tune_intersecting = FALSE" code/estimators.R` returns zero matches)
- [ ] Test job runs successfully (`Rscript code/run_single_replication.R --approach 1 --dgp 1 --n 500 --reps 5 --output results/test.rds`)
- [ ] Output directories created
- [ ] Ready to launch

After launching:
- [ ] All 12024 jobs submitted successfully (15 sbatch calls)
- [ ] Monitor logs for errors
- [ ] Approaches 3/4 DGP4: may have elevated failure rate; errors are logged and informative
- [ ] Wait for completion
- [ ] Combine results (note: approach{3,4}_dgp4_*.rds files may all be errors)
- [ ] Download and analyze

---

## Questions?

If issues arise:
1. Check logs: `logs/*.err`
2. Test locally: `Rscript code/run_single_replication.R ...`
3. Verify packages: `R -e "library(doubletree); library(optimaltrees)"`
4. Check module versions: `module list`

---

## Expected Improvements

With partition-based comparison:

**Approach 4 (doubletree_averaged):**
- More robust modal structure selection
- Fewer failures due to split order differences
- Should see improved coverage in some settings

**Approach 6 (msplit_averaged):**
- Similar improvements
- Better consistency across replications

**Approaches 1-3:**
- Unaffected (don't use structure comparison)

---

**Status:** Ready to deploy
**Date:** 2026-05-21
**Commits:** 94e79dc (optimaltrees), partition-based comparison
