# 🎯 Non-Parametric Tests

> Status: ✅ Complete — [Open the notebook →](07_non_parametric_tests.ipynb)

Topic 07 of the Statistical Inference & Hypothesis Testing repo. Topics 04 and 06 built the t-test and ANOVA, both of which assume the data comes from a Normal (or at least symmetric, finite-variance) population. Non-parametric tests make no such assumption — they work on **ranks** instead of raw values. This notebook builds the three most common rank-based tests from scratch, validates them against `scipy.stats`, and directly measures the actual tradeoff: non-parametric tests are less powerful under clean Normal data, but far more robust once real-world contamination (outliers, skew) enters the picture.

## 📑 Table of Contents

1. [Concept & Intuition](#-concept--intuition)
2. [Why This Topic Matters](#-why-this-topic-matters)
3. [Mathematical Explanation](#-mathematical-explanation)
4. [Formula Reference](#-formula-reference)
5. [Common Pitfalls & Gotchas](#️-common-pitfalls--gotchas)
6. [Function Reference](#-function-reference-used-in-the-notebook)
7. [Self-Test Exercises](#-self-test-exercises)
8. [Notebook](#-notebook)
9. [Further Reading](#-further-reading)

---

## 📖 Concept & Intuition

```
Parametric test          Non-parametric alternative        What changes
-------------------------------------------------------------------------------
t-test (2 groups)     -> Mann-Whitney U                  raw values -> ranks
paired t-test         -> Wilcoxon signed-rank             raw diffs -> ranked diffs
one-way ANOVA (k>=3)  -> Kruskal-Wallis                   raw values -> ranks

                No Normality assumption needed --
                but typically somewhat LESS power under clean Normal data,
                and MORE robust once outliers/heavy tails/skew show up.
                This is a genuine data-dependent tradeoff, not a fixed ranking.
```

---

## 🎯 Why This Topic Matters

- **All three tests were built from scratch and validated against `scipy.stats`** — §1's Mann-Whitney U (269.0 manual vs 269.0 scipy), §2's Wilcoxon signed-rank (W=3.0 both), and §3's Kruskal-Wallis (H=15.9405 manual vs scipy, exact match) all agree closely or exactly with the library implementations.
- **A single-sample outlier-robustness demonstration produced a genuinely surprising, honestly-reported result** — §4 injected one extreme outlier (value=500) into a clean dataset and found, on that *specific* sample, Mann-Whitney's p-value actually shifted *more* (0.0841) than the t-test's (0.0168) — the opposite of the textbook expectation.
- **Repeating the same experiment 2,000 times revealed the real, reliable pattern** — §4b averaged the outlier-injection experiment across 2,000 independent samples and found the t-test's average p-value shift was **0.3921** versus Mann-Whitney's **0.0699** — a nearly 6x difference confirming rank-based tests genuinely are far more outlier-robust *on average*, even though §4's single anecdote happened to land on the atypical side of that distribution — a direct, concrete lesson that any single illustrative example (no matter how intuitive it looks) can mislead, and only repeated trials reveal the real pattern.
- **The power tradeoff was measured in both directions, not just one** — §5 found the t-test held a power edge under clean Normal data (33.27% vs 32.77%), matching the well-known asymptotic relative efficiency theory, then found Mann-Whitney regaining the advantage under heavy-tailed, Cauchy-contaminated data (26.93% vs 14.33%) — the tradeoff genuinely flips based on the data's actual shape.
- **Tie-handling was demonstrated concretely, not just mentioned as a caveat** — §6 showed continuous data with 60/60 unique values (no ties) versus heavily rounded data with only 7/60 unique values (heavy tying), confirming `scipy`'s automatic tie correction keeps the test well-behaved (p=0.9639, a sensible near-1 result for two overlapping rounded distributions).
- **A working Bonferroni-corrected post-hoc procedure was built for Kruskal-Wallis** — §7 found only **1 of 3** pairwise comparisons remained significant after correction (A vs C), directly previewing Topic 11's full treatment of multiple-testing correction.
- **Everything was re-applied to the real, meaningfully skewed `tips` dataset** — §8 found tip percentage's skewness was **3.33** (a genuinely strong candidate for rank-based testing), and both Kruskal-Wallis vs. one-way ANOVA (day-of-week comparison) and Mann-Whitney vs. Welch's t-test (lunch vs. dinner) reached the *same* conclusions at α=0.05 despite testing fundamentally different things (rank distributions vs. means).

---

## 🧮 Mathematical Explanation

### 1. Mann-Whitney U

$$U_1 = R_1 - \frac{n_1(n_1+1)}{2}, \qquad U = \min(U_1, U_2)$$

Ranks all observations from both groups together; a large asymmetry in rank sums between groups indicates one group tends to have higher values than the other, without assuming any particular distribution shape.

### 2. Wilcoxon signed-rank

$$W = \min(W_+, W_-), \qquad W_+ = \sum_{d_i>0} R_i, \qquad W_- = \sum_{d_i<0} R_i$$

Ranks the *absolute* paired differences, then checks whether positive or negative differences dominate — the rank-based analogue of the paired t-test.

### 3. Kruskal-Wallis

$$H = \frac{12}{N(N+1)}\sum_{i=1}^{k}\frac{R_i^2}{n_i} - 3(N+1), \qquad df = k-1$$

Ranks all observations across all $k$ groups together, then compares each group's average rank — the rank-based analogue of one-way ANOVA's F-statistic.

### 4. Tie correction

When observations tie, they share the average of the ranks they would otherwise occupy; the variance formula used to compute a rank-sum test's p-value must be adjusted downward slightly to remain valid — §6's manual `tie_correction` term (in §1's from-scratch implementation) mirrors what `scipy` applies automatically.

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| Mann-Whitney U | $R_1 - n_1(n_1+1)/2$ | §1 |
| Wilcoxon W | $\min(W_+, W_-)$ | §2 |
| Kruskal-Wallis H | $\frac{12}{N(N+1)}\sum R_i^2/n_i - 3(N+1)$ | §3 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Trusting a single illustrative example over a repeated simulation** — §4's single-sample result (Mann-Whitney moved *more* than the t-test) directly contradicted §4b's 2,000-sample average (t-test moved nearly 6x more) — this notebook's own honest reporting is the pitfall demonstration.
2. **Assuming non-parametric tests are always "safer" or "more powerful"** — §5 showed the t-test genuinely holding a power edge under clean Normal data; non-parametric tests trade some power for robustness, they don't dominate unconditionally.
3. **Ignoring ties in rank-based tests on rounded or discrete-like data** — §6's heavily rounded data (7 unique values from 60 observations) needs the tie correction to stay valid; skipping it would drift the Type I error rate away from nominal.
4. **Treating a significant Kruskal-Wallis result as identifying which groups differ** — exactly like ANOVA's F-test (Topic 06), it only says *at least one* group's rank distribution differs; §7's Bonferroni-corrected post-hoc found only 1 of 3 pairs actually significant.
5. **Assuming non-parametric tests answer the exact same question as their parametric counterparts** — Mann-Whitney tests whether one distribution tends to produce larger values than another (roughly, a difference in medians/location under symmetry), not literally a difference in means; the two can occasionally disagree on genuinely skewed data even when they usually agree, as they happened to in §8.
6. **Reaching for a non-parametric test only because data "looks a little skewed"** — always check the actual robustness/power tradeoff for the situation (§5's dual-direction result) rather than defaulting to "non-parametric is always the safe choice."

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `scipy.stats.mannwhitneyu` | Mann-Whitney U test reference |
| `scipy.stats.wilcoxon` | Wilcoxon signed-rank test reference |
| `scipy.stats.kruskal` | Kruskal-Wallis test reference |
| `scipy.stats.rankdata` | Manual rank computation (with tie handling) |
| `scipy.stats.ttest_ind` / `.f_oneway` | Parametric counterparts used for direct power/robustness comparisons |
| `np.random.default_rng().standard_cauchy` | Heavy-tailed contamination for the power-tradeoff simulation |

---

## 📝 Self-Test Exercises

1. Section 4 found a single-sample result (Mann-Whitney moved more) that contradicted the 2,000-sample average (t-test moved far more). Using basic sampling variability reasoning, explain why a single random sample is not reliable evidence for "which test is more robust," even though the illustration felt concrete and specific.
2. Using §3 of the Mathematical Explanation, explain why the Wilcoxon signed-rank test operates on ranks of *absolute* differences rather than ranks of the raw before/after values directly.
3. Section 5 found the t-test's power edge under clean Normal data (33.27% vs 32.77%) was fairly small, while Mann-Whitney's edge under heavy-tailed data (26.93% vs 14.33%) was much larger. What does this asymmetry suggest about the relative "cost" of picking the wrong test in each direction?
4. Section 6 found only 7 unique values among 60 rounded observations. Using the tie-correction formula referenced in §4 of the Mathematical Explanation, explain in your own words why more ties should shrink the variance estimate used in the test statistic's standardization.
5. Section 8 found Kruskal-Wallis and one-way ANOVA reached the same conclusion on the real `tips` day-of-week comparison. Propose one scenario (in words, not code) where you'd expect these two tests to disagree on a real dataset, based on what you learned about their different underlying assumptions.

---

## 📓 Notebook

32 executed code cells: Mann-Whitney U, Wilcoxon signed-rank, and Kruskal-Wallis tests each built from scratch and validated against `scipy.stats`, a single-sample outlier-injection demonstration that (honestly) initially contradicted the expected direction, followed by a 2,000-repeat averaged simulation revealing the true, reliable pattern (t-test average shift 0.392 vs Mann-Whitney's 0.070), a bidirectional power-tradeoff simulation (t-test ahead under clean Normal data, Mann-Whitney ahead under heavy-tailed data), a concrete demonstration of tie-correction behavior on heavily rounded data, a Bonferroni-corrected pairwise post-hoc procedure following a significant Kruskal-Wallis result, and a full application of Kruskal-Wallis, Mann-Whitney, ANOVA, and Welch's t-test to the real, meaningfully skewed (skewness=3.33) `tips` dataset:

➡️ **[07_non_parametric_tests.ipynb](07_non_parametric_tests.ipynb)**

---

## 📚 Further Reading

- [scipy.stats.mannwhitneyu / wilcoxon / kruskal documentation](https://docs.scipy.org/doc/scipy/reference/stats.html)
- [Wasserman, *All of Statistics*, Ch. 10: Hypothesis Testing (non-parametric section)](https://link.springer.com/book/10.1007/978-0-387-21736-9)
- [Mann & Whitney (1947): On a Test of Whether One of Two Random Variables is Stochastically Larger](https://www.jstor.org/stable/2236101)
- [Kruskal & Wallis (1952): Use of Ranks in One-Criterion Variance Analysis](https://www.jstor.org/stable/2280779)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.07_Non_Parametric_Tests&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
