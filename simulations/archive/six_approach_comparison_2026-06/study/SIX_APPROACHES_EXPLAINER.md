# Six ATT Estimators: What Each Approach Does

**Date:** 2026-06-10
**Context:** Six-approach comparison simulation study (ATT estimation with optimal decision trees)

---

## Background

All six approaches estimate the **Average Treatment Effect on the Treated (ATT)**
using the Efficient Influence Function (EIF):

$$\hat\theta = \frac{1}{n} \sum_{i=1}^n \psi_i, \qquad \psi_i = \frac{A_i}{\hat\pi}\bigl(Y_i - \hat m_0(X_i) - \hat\theta\bigr) - \frac{(1-A_i)\,\hat e(X_i)}{\hat\pi\,(1-\hat e(X_i))}\bigl(Y_i - \hat m_0(X_i)\bigr)$$

Every approach needs two nuisance functions:
- **e(X)** — propensity score P(A=1|X), fit on all n observations
- **m0(X)** — control outcome E[Y|A=0,X], fit on control units only

The approaches differ in **how** those nuisance functions are estimated, which affects:
- Whether predictions are cross-fitted (determines bias)
- How many trees are fit (determines stability and interpretability)
- Whether a single presentable tree is returned

---

## The Four Design Dimensions

| Dimension | Question | Values in use |
|-----------|----------|---------------|
| Structure selection | How do we pick which tree shape to use? | CV-best (per fold), Rashomon intersection, modal vote |
| Leaf values | How do we fill in the leaf predictions? | Single fit, K cross-fits, M×K average |
| Cross-fitting | Are predictions out-of-sample for each obs? | Yes (valid inference) / No (may have bias) |
| Interpretability | How many trees to show? | 1 averaged tree, K fold-trees, predictions only |

---

## The Six Approaches

### Approach 1 — Full-Sample Tree

**In one sentence:** Fit one tree on all the data; predict all observations in-sample.

**Algorithm:**
1. CV-select λ on all n observations; fit propensity tree e(X) on all n.
2. CV-select λ on control units; fit outcome tree m0(X) on controls.
3. Predict e(X_i) and m0(X_i) for every i using these same trees.
4. Compute ATT via EIF.

**Cross-fitting:** None — each observation was in the training set.

**Output:** 2 trees (e and m0), one each.

**Expected behavior:**
- Fast (2 CV runs).
- Positive bias because trees overfit their training data (predictions are too
  close to actual outcomes for in-sample observations). The bias is structural
  and doesn't vanish as n → ∞.

---

### Approach 2 — Standard Cross-Fit

**In one sentence:** Fit K independent trees (one per fold) and predict each observation out-of-sample.

**Algorithm:**
1. Partition n observations into K = 5 folds.
2. For each fold k:
   - Train e_k and m0_k on the other K−1 folds.
   - Predict on fold k (out-of-sample).
3. Pool the K sets of predictions into one length-n vector.
4. Compute ATT via EIF.

**Cross-fitting:** Yes — every observation is predicted by a tree that never
saw it during training.

**Output:** K=5 propensity trees + K=5 outcome trees (each potentially with a
different structure); pooled predictions.

**Expected behavior:**
- Asymptotically unbiased (cross-fitting eliminates first-order bias).
- Coverage valid at large n.
- No single interpretable tree; each fold can select a different structure.

---

### Approach 3 — Doubletree (Rashomon Intersection + Cross-Fit)

**In one sentence:** Find the one tree structure that appears near-optimal in every fold,
then use K cross-fitted refits of that structure for inference.

**Algorithm:**
1. Partition into K = 5 folds.
2. For each fold k: fit trees with CV-selected λ, collect the Rashomon set
   (all trees within ε_n = 2√(log n/n) of optimal). At n=500, ε_n ≈ 0.22.
3. Intersect the K Rashomon sets → common structures that are near-optimal
   in every fold. Auto-tune ε_n upward if intersection is empty.
4. For each fold k: refit the selected common structure on the K−1 training
   folds → e_k(X), m0_k(X).
5. Predict on fold k (out-of-sample).
6. Compute ATT via EIF.

**Cross-fitting:** Yes — same structure, different leaf values per fold;
predictions are out-of-sample.

**Output:** K=5 propensity trees + K=5 outcome trees, **all with the same
structure** (different leaf values). The structure is genuinely interpretable:
it's the one that works well across all folds.

**Expected behavior:**
- Asymptotically valid inference (same guarantees as Approach 2).
- Returns a single tree structure that can be shown to an audience.
- Slower (K Rashomon set fits + K refits for each nuisance).

---

### Approach 4 — Doubletree Averaged

**In one sentence:** Use the Rashomon intersection structure from Approach 3, but
average the K sets of leaf values into one tree and predict all observations with it.

**Algorithm:**
1. Run Approach 3 Steps 1–4 to obtain K trees e_1,...,e_K (common structure,
   different leaves).
2. **Average leaf values** across e_1,...,e_K weighted by fold sample sizes
   → one averaged propensity tree; same for m0.
3. Predict **all n observations** using the averaged tree (no cross-fitting).
4. Compute ATT via EIF.

**Cross-fitting:** None — all observations are predicted by a tree that
averaged over K fits, most of which used observation i in training.

**Output:** 1 averaged propensity tree + 1 averaged outcome tree.

**Expected behavior:**
- Returns a single maximally-stable tree (leaf values are averages, less noisy
  than any single fold).
- Has in-sample contamination: for K=5, each leaf average was trained on the
  observation's own fold 4/5 of the time. This introduces structural positive
  bias ≈ (K−1)/K × (leaf overfitting), which is **smaller** than Approach 1's
  bias because the Rashomon intersection constrains the structure.
- Coverage may be below nominal at finite n.

---

### Approach 5 — M-Split (Modal Structure)

**In one sentence:** Vote on the tree structure across M independent splits, then
estimate with proper cross-fitting on every split.

**Algorithm:**

**Stage 1 — Structure selection:**
1. Repeat M = 10 times (each with a fresh random fold assignment):
   - CV-select λ on fold 1's training set; fit e and m0; record their structures.
2. Find the **modal structure** (most common across M candidates) separately
   for e and m0, excluding stumps (constant-prediction trees).

**Stage 2 — Refit M×K times:**
3. For each split m = 1,...,M (same fold assignments as Stage 1):
   - For each fold k = 1,...,K:
     - Refit the modal structure on the K−1 training folds.
     - **Predict on fold k** (out-of-sample) → store in column m of an n×M matrix.

**Stage 3 — Average predictions:**
4. For each observation i, average its M out-of-sample predictions (one per split).
5. Compute ATT via EIF using these averaged predictions.

**Cross-fitting:** Yes — every observation's prediction in each split is from
a tree trained without it.

**Output:** Averaged predictions (one number per observation); modal structure
can be shown for interpretation.

**Expected behavior:**
- Asymptotically valid inference (M independent cross-fits, averaged → lower variance).
- Structure more stable than Approach 2 (modal vote reduces Rashomon ambiguity).
- Computationally expensive (M×K = 50 tree refits per nuisance).

---

### Approach 6 — M-Split Averaged

**In one sentence:** Use the M-split modal structure but, instead of averaging
predictions, average the M×K leaf values into one tree and predict all
observations with it.

**Algorithm:**

**Stage 1 — Structure selection:** (identical to Approach 5)
1. Find modal structure across M splits.

**Stage 2 — Refit M×K times:**
2. For each split m = 1,...,M and each fold k = 1,...,K:
   - Refit modal structure on K−1 training folds.
   - Store the refitted tree (not predictions).
3. Also store out-of-sample predictions during each refit in an n×M matrix
   (available via `predictions_all_splits` for comparison).

**Stage 3 — Average leaf values:**
4. Weighted-average leaf values across all M×K trees → one averaged propensity
   tree; same for m0. (Weights = training fold sample size per leaf.)

**Stage 4 — Predict and compute ATT:**
5. Apply discretization; predict all n observations with the averaged tree.
6. Compute ATT via EIF.

**Cross-fitting:** None — the averaged tree's leaf values were estimated from
M×K refits, M(K−1) of which used each observation in training.

**Output:** 1 averaged propensity tree + 1 averaged outcome tree (maximally
stable — M×K refits averaged vs. K for Approach 4).

**Expected behavior:**
- Most stable/interpretable single tree (M×K = 50 refits averaged).
- In-sample contamination: for K=5, (K−1)/K = 80% of each leaf estimate came
  from refits that included observation i. This contamination is constant
  across M (M cancels in the average), so it does not vanish with more splits.
  The resulting bias = structural (does not decrease with n or M at fixed K).
- The `predictions_all_splits` return field stores properly cross-fitted
  predictions (as in Approach 5) to enable empirical comparison of both
  inference paths.

---

## Side-by-Side Summary

| | 1 Full | 2 Crossfit | 3 Doubletree | 4 DT-Avg | 5 M-split | 6 M-split-Avg |
|---|:---:|:---:|:---:|:---:|:---:|:---:|
| Structure selection | CV (all n) | CV (per fold) | Rashomon ∩ | Rashomon ∩ | Modal vote | Modal vote |
| Trees fit | 2 | 2K | K+Rashomon | K+Rashomon | M×K | M×K |
| Cross-fitted? | No | **Yes** | **Yes** | No | **Yes** | No |
| Single tree output | Yes | No | No (struct only) | **Yes** | No (struct only) | **Yes** |
| Bias expected? | Yes | No | No | Some | No | **Yes** |
| Inference valid? | No | **Yes** | **Yes** | Partial | **Yes** | No |

---

## Data-Generating Processes (Simulation Study)

All DGPs have n ∈ {500, 1000, 2000}, true ATT = 0.15, binary outcome Y.

| DGP | Covariates | Propensity e(X) | Outcome m0(X) | Complexity |
|-----|-----------|----------------|--------------|-----------|
| **Simple** | x1,x2,x3 binary | expit(−0.5 + 0.3x1 + 0.3x2) | 0.2 + 0.15x1 + 0.15x3 | ~2–3 splits |
| **Moderate** | x1–x4 binary | expit(−0.5 + 0.3x1 + 0.2x2 + 0.3x1x2) | 0.2 + 0.2x3 + 0.15x4 + 0.2x3x4 | ~4–5 splits |
| **Complex** | x1–x5 binary | expit(linear + 0.3x1x2 + 0.2x2x3) | 0.05 + linear + 0.2x3x4 + 0.15x4x5 | ~6–8 splits |
| **Continuous** | x1,x2 binary; x3 ∈ [−1,1]; x4 ∼ N(0,1) | expit(−0.5 + 0.3x1 + 0.4x3 + 0.2x4 + 0.2x1x3) | 0.2 + 0.15x2 + 0.2x3 + 0.15x4²/2 + 0.1x2x3 | ~4–6 splits |

The Continuous DGP is the hardest for tree methods: x4² requires multiple thresholds
to capture the parabolic shape. With `discretize_bins = 2` (median split), only one
threshold is placed on x4, collapsing the quadratic to a linear step — this is the
root cause of the DGP4 coverage failure fixed in this session.

---

## Production Parameters (Cluster Runs)

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| K (cross-fit folds) | 5 | Standard; balances bias-variance |
| M (modal splits) | 10 | Enough for stable modal vote |
| K per split (approaches 5–6) | 5 | Same as K above |
| λ via | `cv_regularization_adaptive` | Theory-aligned adaptive search |
| λ cap | 20 × log(n)/n | Prevents over-regularization to stumps |
| λ floor | log(n)/n | Prevents under-regularization (inside package) |
| ε_n (Rashomon bound) | 2√(log n/n) | Theory rate: ε_n = o(n^{−1/2}) |
| `auto_tune_intersecting` | TRUE (approaches 3–4) | Relaxes ε_n if intersection empty |
| `discretize_bins` | "adaptive" | ceil(log(n)/3) thresholds per continuous feature |
| `discretize_method` | "quantiles" | Quantile-based placement |
| `max_depth` | 4 | Caps Rashomon set size; prevents overflow |
| n values | 500, 1000, 2000 | Three sample sizes for rate verification |
| Reps per config | 500–1000 | Depends on approach cost |
