# Theoretical Gap: Discretization and Rate Guarantees

**Date**: 2026-03-03
**Issue**: Implementation-theory mismatch regarding continuous features
**Status**: Identified and resolved in implementation; manuscript needs clarification

---

## The Problem

### What the Manuscript Says

The manuscript (lines 87-89) states:

> **Approximation:** The support of $X$ is contained in $[0,1]^d$ (or a compact set, after transformation). The true nuisances $\eta_0 = (e_0, m_{0,0}, m_{1,0})$ lie in a Hölder class $H^\beta([0,1]^d)$ for some $\beta > 0$.

And (line 89):

> **Complexity:** The tree complexity (number of leaves) satisfies $s_n \asymp n^{d/(2\beta+d)}$.

**Interpretation**: The theory assumes continuous $X \in [0,1]^d$ and that trees can achieve $s_n \asymp n^{d/(2\beta+d)}$ leaves.

### What the Implementation Does

TreeFARMS (via GOSDT) requires **binary features**: each feature must be in $\{0, 1\}$.

**Without discretization**: Users must manually binarize continuous features. No guidance on how many bins.

**With fixed discretization** (initial implementation):
- Median: 1 threshold → 1 binary indicator per feature
- Quantiles with `n_bins=k`: `k-1` thresholds → `k-1` binary indicators per feature
- With $d$ features and $b$ binary indicators per feature → tree has at most $2^{d \cdot b}$ leaves
- This is **constant in $n$**, violating $s_n \asymp n^{d/(2\beta+d)} \to \infty$

### The Gap

**Theory requires**: $s_n \to \infty$ as $n \to \infty$
**Fixed bins give**: $s_n \le 2^{d \cdot b} = O(1)$ (constant in $n$)

This is a **fundamental mismatch** between theory and implementation.

---

## Why It Matters

### Approximation Error

With $b$ binary indicators per feature (from $b-1$ thresholds), the tree can partition $[0,1]^d$ into at most $2^{d \cdot b}$ regions.

**Approximation error**: For a Hölder-$\beta$ function, partitioning into $M$ axis-aligned rectangles gives approximation error:
$$
\inf_{f \in \mathcal{T}_M} \|f - \eta_0\|_{L^2}^2 \lesssim M^{-2\beta/d}
$$

With fixed $b$ (thus fixed $M = 2^{d \cdot b}$), this error is **constant**, independent of $n$.

### Oracle Inequality

The excess risk bound is:
$$
R(\hat{\eta}) - R(\eta_0) \lesssim \inf_{f \in \mathcal{T}_{s_n}} [R(f) - R(\eta_0)] + \frac{s_n \log n}{n}
$$

- **Approximation term**: $\inf_{f} [R(f) - R(\eta_0)] \gtrsim s_n^{-2\beta/d}$ (assuming $s_n$ controls partition size)
- **Estimation term**: $\frac{s_n \log n}{n}$

**Optimal balance**: Set $s_n^{-2\beta/d} \asymp \frac{s_n \log n}{n}$, which gives $s_n \asymp n^{d/(2\beta+d)}$.

**With fixed $s_n$**: The approximation error dominates for large $n$, and the rate cannot improve beyond the discretization error.

### DML Rate Condition

DML requires $\|\hat{\eta} - \eta_0\|_{L^2} = o_p(n^{-1/4})$.

From the oracle inequality with optimal $s_n$:
$$
\|\hat{\eta} - \eta_0\|_{L^2} = O_p(n^{-\beta/(2\beta+d)})
$$

This is $o_p(n^{-1/4})$ if and only if $\beta/(2\beta+d) > 1/4$, i.e., $\beta > d/2$.

**But this requires $s_n \to \infty$!** With fixed discretization, the rate is bounded by the discretization error, which may not be $o_p(n^{-1/4})$.

---

## The Solution: Adaptive Discretization

### Mathematical Requirement

To allow $s_n \asymp n^{d/(2\beta+d)}$ leaves, we need enough binary indicators.

With $b_n$ binary indicators per feature (from $b_n - 1$ thresholds), the tree can create up to $2^{d \cdot b_n}$ leaves.

**Requirement**: $2^{d \cdot b_n} \gtrsim n^{d/(2\beta+d)}$

Taking logarithms:
$$
d \cdot b_n \gtrsim \frac{d}{2\beta+d} \log n
$$

Therefore:
$$
b_n \gtrsim \frac{\log n}{2\beta+d}
$$

**Practical choice**: $b_n = \max(2, \lceil \alpha \log n \rceil)$ for some constant $\alpha > 0$.

In implementation: $b_n = \max(2, \lceil \log n / 3 \rceil)$ gives $\alpha \approx 1/3$.

### Why This Works

**As $n \to \infty$**:
- $b_n \sim \log n \to \infty$
- Maximum possible leaves: $2^{d \cdot b_n} = n^{d \cdot \alpha}$ (assuming $\alpha$ is the constant)
- Tree can achieve complexity $s_n \asymp n^{d/(2\beta+d)}$ if $n^{d/(2\beta+d)} \le n^{d \cdot \alpha}$
- This holds if $\alpha \ge 1/(2\beta+d)$, i.e., $\alpha(2\beta+d) \ge 1$

**With $\alpha = 1/3$**: This works for $\beta \ge (1 - d/3)/2$. For example:
- $d=3$, $\beta > 3/2$: Need $\alpha(2(3/2)+3) = \alpha \cdot 6 \ge 1$, so $\alpha \ge 1/6$ ✓
- $d=5$, $\beta > 5/2$: Need $\alpha(2(5/2)+5) = \alpha \cdot 10 \ge 1$, so $\alpha \ge 1/10$ ✓

The choice $\alpha = 1/3$ is conservative and works for most practical $(\beta, d)$ pairs satisfying $\beta > d/2$.

### Approximation Error with Adaptive Bins

With $b_n \sim \log n$ binary indicators per feature:
- Partition into $\approx n^{d \cdot \alpha}$ regions (in expectation)
- Approximation error: $\approx n^{-2\beta \alpha}$
- For $\alpha$ large enough, this can be made $o(n^{-1/2})$, which is faster than needed ($o(n^{-1/4})$ suffices)

**Key point**: The discretization error does not bottleneck the convergence rate.

---

## Implementation

### Code

```r
# Adaptive discretization (theory-aligned)
model <- treefarms(
  X, y,
  discretize_method = "quantiles",
  discretize_bins = "adaptive"  # max(2, ceiling(log(n)/3))
)
```

### Growth Pattern

| Sample size $n$ | Bins $b_n$ | Max leaves $2^{d \cdot b_n}$ (d=3) |
|-----------------|------------|-------------------------------------|
| 50              | 2          | 64                                   |
| 200             | 2          | 64                                   |
| 1000            | 3          | 512                                  |
| 10,000          | 4          | 4,096                                |

As $n$ grows, the tree's representational capacity grows, allowing it to approximate increasingly complex functions.

### Backward Compatibility

**Default**: `discretize_bins = 2` (fixed, for simplicity in small samples)

**For theory**: Use `discretize_bins = "adaptive"` to satisfy rate conditions

---

## Manuscript Implications

### What Needs to Change

The manuscript should clarify:

1. **TreeFARMS operates on binary features**, not continuous $X$ directly
2. **Discretization is required** for continuous covariates
3. **Adaptive discretization** (bins $\sim \log n$) is necessary for theoretical rates

### Suggested Additions

#### In "Optimal decision trees as nuisance learners" (Section 2.3, after line 75)

Add:

> **Discretization of continuous features.** TreeFARMS, like GOSDT, operates on binary features $X \in \{0,1\}^p$. When the original covariates $X \in [0,1]^d$ are continuous, we discretize each feature using quantile-based thresholds. For feature $j$, we compute $b_n - 1$ thresholds (the $(k/(b_n))$-quantiles for $k=1,\ldots,b_n-1$), creating $b_n$ binary indicators. The resulting $p = d \cdot b_n$ binary features are passed to TreeFARMS.
>
> To ensure the tree can achieve complexity $s_n \asymp n^{d/(2\beta+d)}$, we require $b_n \to \infty$ as $n \to \infty$. We use $b_n = \max(2, \lceil \log n / c \rceil)$ for a constant $c > 0$ (we use $c=3$ in practice). This ensures the tree can partition $[0,1]^d$ finely enough that the discretization error is $o(n^{-1/2})$, which does not affect the DML rate. The choice $b_n \sim \log n$ is sharp: slower growth would bottleneck the approximation, while faster growth would be unnecessarily costly.

#### In "Main result" (Section 2.4, before line 89)

Modify the **Complexity** condition:

> - **Complexity:** The tree complexity (number of leaves) satisfies $s_n \asymp n^{d/(2\beta+d)}$. When continuous features are discretized with $b_n \sim \log n$ binary indicators per feature, the tree has sufficient representational capacity to achieve this complexity.

#### In Appendix (after Lemma on approximation, ~line 190)

Add:

> **Remark (Discretization).** The approximation bound assumes axis-aligned trees on $[0,1]^d$. When $X$ is continuous, TreeFARMS operates on discretized binary features. With $b_n$ binary indicators per feature, the tree partitions each coordinate into at most $b_n$ intervals, giving at most $(b_n)^d$ possible rectangular regions. For a Hölder-$\beta$ function, this yields approximation error $\lesssim b_n^{-2\beta}$ (in each coordinate). With $b_n = \Theta(\log n)$, the discretization error is $(\log n)^{-2\beta} = o(n^{-1/2})$ for any $\beta > 0$, hence negligible relative to the statistical error $n^{-\beta/(2\beta+d)}$ required for DML.

---

## Simulation Hypotheses

### What to Expect

**Fixed bins** (e.g., $b=2$ or $b=4$):
- Small $n$: Comparable performance to adaptive (both bottlenecked by estimation error)
- Large $n$: Performance plateaus as approximation error dominates
- Test MSE: Stops improving beyond a threshold
- Tree leaves: Constant (bounded by $2^{d \cdot b}$)

**Adaptive bins** ($b_n \sim \log n$):
- Small $n$: Similar to fixed (both use $b=2$)
- Large $n$: Continues improving as $n$ grows (approximation error shrinks)
- Test MSE: Decreases as $n^{-2\beta/(2\beta+d)}$ (minimax optimal)
- Tree leaves: Grows with $n$, approaching $n^{d/(2\beta+d)}$

### Key Comparison

**Critical metric**: Excess risk (test MSE - oracle MSE) vs. sample size on log-log scale.

- **Theory predicts**: Slope should be $-2\beta/(2\beta+d)$ for adaptive, but flatten for fixed bins
- **With $\beta = 2$, $d=3$**: Slope should be $-2(2)/(2(2)+3) = -4/7 \approx -0.57$

**Visual test**: On log-log plot of excess risk vs $n$:
- Adaptive: Straight line with slope $\approx -0.57$
- Fixed-2: Starts parallel, then flattens (becomes horizontal)
- Fixed-4: Better than Fixed-2 but still flattens

---

## Conclusion

**Bottom line**: Adaptive discretization is **required** for the theoretical guarantees to hold. Fixed bins create a glass ceiling on performance that violates the theory's assumptions.

**Implementation**: Now supports `discretize_bins = "adaptive"` ✓

**Manuscript**: Needs clarification that discretization is part of the method and must be adaptive ✓ (draft above)

**Simulation**: Running to empirically validate the theoretical predictions ⏳
