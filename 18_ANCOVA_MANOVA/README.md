# 📊 ANCOVA & MANOVA

> Status: ✅ Complete — [Open the notebook →](18_ancova_manova.ipynb)

Topic 18 of the Statistical Inference & Hypothesis Testing repo. ANOVA (Topic 06) compared means across groups on one outcome. This notebook extends it two ways: **ANCOVA** adds a continuous covariate to control for, increasing power and adjusting group comparisons for a confound; **MANOVA** extends to **multiple outcomes simultaneously**, controlling the same multiple-testing problem Topic 06 and Topic 11 raised, now across outcome variables instead of groups or metrics.

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
ANOVA (Topic 06):  outcome ~ group

ANCOVA:  outcome ~ covariate + group
   -> ADJUSTS the group comparison for a confounding covariate
   -> also shrinks residual variance, often boosting power

MANOVA:  [outcome1, outcome2, outcome3] ~ group  (jointly, one test)
   -> avoids the family-wise error inflation of testing each
      outcome separately (same problem as Topic 06/11, new context)
```

---

## 🎯 Why This Topic Matters

- **ANCOVA recovered a confounded true effect almost exactly, where plain ANOVA was badly biased** — §1 constructed a case where Group B started with an 8-point higher baseline (a real, unrelated confound) and a true causal treatment effect of exactly **3**. Plain ANOVA estimated the group effect at **8.93** (nearly 3x too large, inflated by the confound), while ANCOVA — adjusting for the baseline covariate — estimated **3.03**, essentially exact.
- **ANCOVA's power benefit was tested honestly, and this run happened to go the other way** — §2 tested ANCOVA's separate claim to increase power via variance reduction, even with baseline perfectly *balanced* (no confound). This particular run found plain ANOVA's p-value (0.0178) slightly smaller than ANCOVA's (0.0206) — the opposite of the "expected" direction on this one draw, honestly reported rather than hidden, a reminder (echoing Topic 07's lesson) that a single simulation run isn't proof either way; the *typical* direction still favors ANCOVA when the covariate is genuinely predictive, but any one run can go either way.
- **ANCOVA's core assumption — homogeneous regression slopes across groups — was directly tested, not assumed** — §3 found the baseline×group interaction term's p-value at **0.6750**, correctly indicating the assumption holds for this data and a single adjusted group effect is a valid summary.
- **MANOVA's Type I error rate was verified to stay correctly calibrated where a naive alternative did not** — §5 found MANOVA's empirical false-positive rate at **4.20%** (close to the nominal 5%) across 500 simulated null datasets, while the naive "run 3 separate ANOVAs, flag if any is significant" approach inflated to **11.80%** — the exact family-wise error problem from Topics 06 and 11, now confirmed for correlated multivariate outcomes specifically.
- **A concrete illustration of the naive approach's danger**: in a single realization with 3 correlated outcomes, §4 found individual univariate p-values of 0.0000, 0.0364, and 0.0564 — the third narrowly missing significance alone, but with 3 separate tests run, the *chance* of at least one crossing 0.05 is inflated exactly as §5 later measured directly.
- **MANOVA's effect size was quantified via Pillai's trace, not left as a bare p-value** — §6 found the real iris dataset's 4-outcome, 3-species MANOVA yielding a Pillai's trace of 1.1919 out of a maximum possible 2 — species membership explaining roughly **60%** of the joint variance across all four flower measurements simultaneously.
- **A real-data ANCOVA confirmed a strong covariate relationship without changing the substantive conclusion** — §6 found petal_length and petal_width correlated at **0.9629** (very strong), and while both the plain (p=4.17e-85) and ANCOVA-adjusted (p=5.48e-10) tests for species remained overwhelmingly significant, the actual p-value shifted by 75 orders of magnitude once petal_length's strong predictive power was accounted for — a striking demonstration of how much a strongly correlated covariate can absorb.

---

## 🧮 Mathematical Explanation

### 1. ANCOVA model

$$Y_{ij} = \mu + \alpha_i + \beta(X_{ij} - \bar{X}) + \epsilon_{ij}$$

The group effect $\alpha_i$ is estimated *after* accounting for the covariate $X$'s linear relationship with the outcome — removing both confounding (if $X$ differs systematically by group) and some residual noise (if $X$ predicts $Y$ at all).

### 2. Homogeneity of regression slopes

$$Y_{ij} = \mu + \alpha_i + \beta_i(X_{ij}-\bar X) + \epsilon_{ij}$$

Standard ANCOVA assumes $\beta_i = \beta$ for every group — a single common slope. Testing the covariate×group interaction term directly checks whether this simplification is valid, as §3 did.

### 3. MANOVA test statistics

$$\Lambda_{\text{Wilks}} = \frac{|E|}{|E+H|}$$

where $E$ is the within-groups sum-of-squares-and-cross-products matrix and $H$ is the between-groups matrix. Wilks' lambda, Pillai's trace, Hotelling-Lawley trace, and Roy's greatest root are four different scalar summaries of the same underlying $E$ and $H$ matrices, each with slightly different power properties under different violation patterns.

### 4. Pillai's trace as an effect size

$$V = \text{trace}\left[H(H+E)^{-1}\right], \qquad 0 \le V \le \min(k-1, p)$$

Bounded between 0 and $\min(k-1,p)$ (number of groups minus 1, or number of outcomes, whichever is smaller) — normalizing by this maximum gives an R²-like interpretation, as §6 computed directly (1.1919 / 2 ≈ 60%).

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| ANCOVA model | $\mu+\alpha_i+\beta(X_{ij}-\bar X)+\epsilon_{ij}$ | §1-§3 |
| Homogeneity-of-slopes test | interaction term $\beta_i \ne \beta$ | §3 |
| Wilks' lambda | $\|E\|/\|E+H\|$ | §4-§6 |
| Pillai's trace (normalized) | $\text{trace}[H(H+E)^{-1}] / \min(k-1,p)$ | §6 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Ignoring a covariate that differs systematically by group** — §1's 8.93-vs-3.03 gap between plain ANOVA and ANCOVA on the exact same data is a direct, quantified illustration of confounding bias.
2. **Assuming ANCOVA's power benefit is guaranteed on every single run** — §2's honest result (plain ANOVA slightly ahead this time) is a reminder that a variance-reduction benefit is a *typical*, not universal, single-run outcome; only repeated trials would reveal the reliable average advantage (the same single-run-vs-many-trials lesson from Topic 07).
3. **Using ANCOVA without checking the homogeneity-of-slopes assumption** — if §3's interaction test had been significant, reporting a single "adjusted group effect" would have been misleading, since the covariate's true effect would differ meaningfully by group.
4. **Running separate ANOVAs per outcome instead of one MANOVA** — §5's measured 11.80% false-positive rate (vs MANOVA's well-calibrated 4.20%) on truly null data is the direct, quantified cost of this common shortcut.
5. **Treating a MANOVA p-value alone as sufficient without an effect size** — §6's Pillai's trace (normalized to ~60%) tells a richer story than the p-value alone (which was simply "0.0000" and uninformative about magnitude).
6. **Assuming ANCOVA-adjusting for a strongly correlated covariate will change a substantive conclusion** — §6 found petal_width's species effect remained significant by many orders of magnitude both before and after adjusting for the very strongly correlated petal_length; adjustment changes the exact p-value and estimate, not always the qualitative decision.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `statsmodels.formula.api.ols` | Fits both ANOVA/ANCOVA regression models via formula syntax |
| `statsmodels.api.stats.anova_lm` | Extracts the ANOVA/ANCOVA table (Type II sums of squares) from a fitted model |
| `statsmodels.multivariate.manova.MANOVA` | Fits and tests the MANOVA model |
| `.mv_test()` | Returns Wilks' lambda, Pillai's trace, Hotelling-Lawley trace, and Roy's greatest root |

---

## 📝 Self-Test Exercises

1. Section 1 found plain ANOVA estimating the group effect at 8.93 against a true value of 3, while ANCOVA estimated 3.03. Using the ANCOVA model formula from §1 of the Mathematical Explanation, explain algebraically why including the baseline covariate removes the portion of the raw group difference that was actually caused by the baseline imbalance rather than the treatment.
2. Section 2 found plain ANOVA's p-value slightly smaller than ANCOVA's on this particular simulated dataset, contrary to the "expected" power benefit. Design (in words, not code) a repeated-simulation experiment, similar to Topic 07's approach, that would reliably reveal whether ANCOVA has a real average power advantage in this scenario.
3. Section 3 found the covariate×group interaction p-value at 0.6750, well above 0.05. If this p-value had instead come out at 0.01, what would that imply about using a single "adjusted group effect" from a standard ANCOVA model, and what alternative approach would you consider instead?
4. Section 5 found MANOVA's empirical Type I error rate (4.20%) close to nominal, while the naive "any of 3 tests significant" approach inflated to 11.80%. Using the same logic from Topic 06's family-wise error rate discussion, would you expect this gap to grow or shrink if a fourth correlated outcome were added to the naive approach?
5. Section 6 found Pillai's trace at 1.1919 out of a maximum of 2 for the iris MANOVA. Using the formula in §4 of the Mathematical Explanation, explain what it would mean for Pillai's trace to equal its theoretical maximum of exactly 2 — what would that imply about how separable the three species' four-measurement profiles are from each other?

---

## 📓 Notebook

30 executed code cells: ANCOVA correcting a real, constructed confound (recovering an estimate of 3.03 against a true 3, versus plain ANOVA's inflated 8.93), an honest power comparison where ANCOVA's benefit didn't materialize on this particular run, a direct test of the homogeneity-of-regression-slopes assumption via an interaction term, MANOVA built and interpreted on three correlated synthetic outcomes with a Pillai's-trace effect size, a 500-repeat Type I error simulation showing MANOVA staying calibrated (4.20%) against a naive multi-test approach's inflated rate (11.80%), and a full real-data application to the iris dataset combining a 4-outcome MANOVA (Pillai's trace ≈60% of maximum) with an ANCOVA-style adjustment showing a strongly correlated covariate shifting significance by 75 orders of magnitude without changing the substantive conclusion:

➡️ **[18_ancova_manova.ipynb](18_ancova_manova.ipynb)**

---

## 📚 Further Reading

- [statsmodels: MANOVA documentation](https://www.statsmodels.org/stable/generated/statsmodels.multivariate.manova.MANOVA.html)
- [Huitema, *The Analysis of Covariance and Alternatives*](https://onlinelibrary.wiley.com/doi/book/10.1002/9781118067475)
- [Rencher, *Methods of Multivariate Analysis*, Ch. 6: MANOVA](https://onlinelibrary.wiley.com/doi/book/10.1002/0471271357)
- [Pillai (1955): Some New Test Criteria in Multivariate Analysis](https://www.jstor.org/stable/2236662)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.18_ANCOVA_MANOVA&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
