# Six Approaches: Detailed Specification

**Date:** 2026-05-26
**Status:** Authoritative specification

---

## Overview

Six tree-based causal inference approaches for estimating ATT, differing in:
1. How tree **structure** is selected
2. How **leaf values** are determined
3. Whether cross-fitting is used

---

## The Six Approaches

| # | Name | Structure Selection | Leaf Values | Cross-Fitting | Result |
|---|------|-------------------|-------------|---------------|--------|
| **1** | Fullsample | All n obs | Single fit (all n) | No | 2 trees (e, m0) |
| **2** | Crossfit | Per fold (independent) | K cross-fits | Yes | 10 trees (5e, 5m0) |
| **3** | Doubletree | Rashomon intersection | K cross-fits | Yes | K trees with common structure |
| **4** | Doubletree Averaged | Rashomon intersection | Average K trees | No | 1 averaged tree |
| **5** | M-split | Modal (vote M×K) | M cross-fits | Yes | Averaged predictions across M |
| **6** | M-split Averaged | Modal (vote M×K) | Average M×K trees | No | 1 averaged tree |

---

## Detailed Algorithms

### Approach 1: Fullsample

**Goal:** Simplest possible - one tree per nuisance

**Algorithm:**
1. Fit propensity tree e(X) on all n observations with CV-selected λ
2. Fit outcome tree m0(X) on control units with CV-selected λ
3. Predict e(X) and m0(X) for all observations
4. Compute ATT via EIF: ψ_i = (A_i/π̂) × [Y_i - m0(X_i)] + ...
5. θ̂ = mean(ψ_i)

**Output:**
- 2 trees: e(X), m0(X)
- Single cross-section predictions (no cross-fitting)
- **Inference:** May have bias (overfitting)

---

### Approach 2: Crossfit (Standard)

**Goal:** Valid inference via cross-fitting with independent trees per fold

**Algorithm:**
1. Split data into K=5 folds
2. For each fold k:
   - Train e_k(X) and m0_k(X) on other 4 folds with CV-selected λ
   - Predict on fold k (test set)
3. Pool cross-fitted predictions
4. Compute ATT via EIF with pooled predictions

**Output:**
- 10 trees (5 propensity + 5 outcome), each with potentially different structure
- Cross-fitted predictions
- **Inference:** Valid (unbiased)

---

### Approach 3: Doubletree (Rashomon Intersection)

**Goal:** Single interpretable structure + valid inference

**Algorithm:**
1. For each fold k = 1, ..., K:
   - Fit trees with CV-selected λ (floored at sqrt(log(n)/n))
   - Collect Rashomon set using theory-justified ε_n = 2*sqrt(log(n)/n)
     (≈ 0.22 at n=500, ≈ 0.17 at n=1000, ≈ 0.12 at n=2000)
2. Find COMMON structure via intersection of K Rashomon sets
3. Refit this common structure K times (once per fold) with cross-fitting
4. Each fold gets predictions from its refitted tree
5. Compute ATT via EIF with cross-fitted predictions

**Output:**
- K trees, all with **same structure**, different leaf values
- Cross-fitted predictions
- **Inference:** Valid (unbiased)
- **Interpretability:** Can show single structure

---

### Approach 4: Doubletree Averaged

**Goal:** Single averaged tree (maximum stability)

**Algorithm:**
1. Obtain K trees with common structure from Approach 3
2. **Average leaf values** across K trees → 1 averaged tree
3. Predict for ALL observations using single averaged tree (no cross-fitting)
4. Compute ATT via EIF

**Output:**
- 1 averaged tree per nuisance
- Non-cross-fitted predictions
- **Inference:** May have bias (uses all data)
- **Interpretability:** Single tree, averaged values

**Note:** Trades valid inference for maximum interpretability/stability

---

### Approach 5: M-split (Modal Structure)

**Goal:** Stable structure + valid inference via repeated splitting

**Algorithm:**
1. **Stage 1 - Structure Selection:**
   - Repeat M=10 times:
     - Split data into K=5 folds independently
     - Fit tree on fold 1's training set with CV-selected λ
     - Extract structure
   - Vote: Find **modal** (most common) structure across M trees

2. **Stage 2 - Refit Modal Structure:**
   - For each split m = 1, ..., M:
     - Split data into K folds (independent fold assignments per split)
     - For each fold k:
       - Refit modal structure on training folds
       - Predict on test fold k
   - Result: n×M matrix of predictions (M predictions per observation)

3. **Stage 3 - Average Predictions:**
   - For each observation i:
     - Average predictions across M splits (only test-set predictions)
   - Compute ATT via EIF with averaged predictions

**Output:**
- Modal structure selected from M×K = 50 trees
- Cross-fitted predictions, averaged across M splits
- **Inference:** Valid (unbiased, cross-fitted)
- **Stability:** Structure from voting reduces Rashomon effect

**Key:** Each observation gets M predictions (from M different splits), averaged.

---

### Approach 6: M-split Averaged

**Goal:** Single highly-averaged tree (maximum averaging)

**Algorithm:**
1. **Stage 1 - Structure Selection:** (Same as Approach 5)
   - Find modal structure across M trees

2. **Stage 2 - Refit Modal Structure:** (Same as Approach 5)
   - Refit modal structure M×K times
   - But **store trees**, not predictions

3. **Stage 3 - Average Leaf Values:**
   - Average leaf values across ALL M×K trees → **1 averaged tree**

4. **Stage 4 - Predict:**
   - Predict for ALL observations using single averaged tree (no cross-fitting)
   - Compute ATT via EIF

**Output:**
- 1 highly-averaged tree per nuisance (average of M×K trees)
- Non-cross-fitted predictions
- **Inference:** May have bias (uses all data)
- **Interpretability:** Single tree with maximal averaging
- **Stability:** Most stable leaf values (averaged across M×K)

**Note:** Trades valid inference for maximum stability/interpretability

---

## Key Differences Summary

### Structure Selection:
- **1, 2:** Each tree independent
- **3, 4:** Rashomon intersection
- **5, 6:** Modal voting

### Leaf Values:
- **1:** Single fit (all n)
- **2, 3, 5:** Cross-fitted (K or M refits)
- **4:** Average K trees
- **6:** Average M×K trees

### Cross-Fitting:
- **Yes (valid inference):** 2, 3, 5
- **No (may have bias):** 1, 4, 6

### Number of Trees:
- **2 trees:** 1, 4, 6
- **K trees:** 3 (same structure)
- **10 trees:** 2 (independent structures)
- **Predictions from M×K:** 5 (averaged)

---

## Interpretation vs Inference Tradeoff

**Valid Inference (Cross-Fitted):**
- Approaches 2, 3, 5
- Unbiased estimates
- Correct coverage
- But: Multiple trees (2, 3) or only predictions (5)

**Single Interpretable Tree:**
- Approaches 1, 4, 6
- Can show one tree to audience
- But: May have bias from overfitting

**Best of Both Worlds?**
- **Approach 3:** K trees with same structure → show structure, report valid estimate
- **Approach 5:** Valid estimate from averaged predictions, modal structure for interpretation

---

## Computational Cost

| Approach | Trees Fitted | CV Calls | Relative Cost |
|----------|-------------|----------|---------------|
| 1 | 2 | 2 | 1× (baseline) |
| 2 | 10 | 10 | 5× |
| 3 | K + Rashomon | K×5 (for Rashomon) + K | 8× |
| 4 | K + Rashomon | K×5 + K | 8× (same as 3) |
| 5 | M×K | M×K | 25× (M=10, K=5) |
| 6 | M×K | M×K | 25× (same as 5) |

---

## Recommended Use Cases

**Approach 1:** Quick exploration, hypothesis generation
**Approach 2:** Standard valid inference, don't care about interpretability
**Approach 3:** Want single interpretable structure + valid inference
**Approach 4:** Want single averaged tree, okay with potential bias
**Approach 5:** Want most stable valid inference
**Approach 6:** Want single maximally-averaged tree, okay with potential bias

---

## Implementation Status

- [x] Approach 1: Implemented and tested
- [x] Approach 2: Implemented and tested
- [x] Approach 3: Implemented (via `estimate_att` with `use_rashomon=TRUE`)
- [x] Approach 4: Implemented (`estimate_att_doubletree_averaged`)
      Uses Rashomon intersection to find K trees with common structure, then averages leaf values
- [x] Approach 5: Implemented (`estimate_att_msplit`)
- [x] Approach 6: Implemented (`estimate_att_msplit_averaged`)
      Refits modal structure M×K times, averages ALL leaf values across M×K trees

**Last updated:** 2026-06-09

---

## Production Parameters

| Parameter | Value | Used In |
|-----------|-------|---------|
| lambda selection | `cv_regularization_adaptive`, max_lambda = 15·log(n)/n | **All 6 approaches** |
| K (cross-fitting folds) | 5 | All approaches |
| M (modal structure splits) | 10 | Approaches 5-6 |
| K per split | 5 | Approaches 5-6 |
| epsilon_n (Rashomon bound) | `2*sqrt(log(n)/n)` | Approaches 3-4 |
| auto_tune_intersecting | TRUE | Approaches 3-4 |
| Fallback on empty intersection | None — hard stop() | Approaches 3-4 |
| lambda floor (after CV) | `sqrt(log(n)/n)` | Approaches 3-6 (inside package) |
| max_depth | 4 | All approaches (inside package) |
| n values | 500, 1000, 2000 | All |
| Reps per config | 1000 | All |

epsilon_n approx 0.22 (n=500), 0.17 (n=1000), 0.12 (n=2000).
Theory: satisfies o(n^{-1/2}) rate for valid EIF-ATT inference (Appendix A.5).
Auto-tuning starts at epsilon_n and increases only if needed; hard failure (stop()) if no intersection
found after exhausting all tuning attempts. Reps that fail are logged as errors, not silently skipped.
