# 🗂️ Categorical Data Analysis

> Status: ✅ Complete — [Open the notebook →](14_categorical_data_analysis.ipynb)

Topic 14 of the Statistical Inference & Hypothesis Testing repo. Topic 05 already built chi-square goodness-of-fit, independence, Cramer's V, and Fisher's Exact Test. This notebook extends the categorical data toolkit into territory Topic 05 didn't cover: odds ratios and relative risk with proper confidence intervals, McNemar's test for **paired** (not independent) categorical data, Simpson's Paradox and the Cochran-Mantel-Haenszel test that resolves it, and Cohen's Kappa for inter-rater agreement.

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
Odds Ratio & Relative Risk  -> two DIFFERENT measures of a 2x2 association,
                                 nearly equal only when the outcome is RARE

McNemar's Test              -> chi-square independence assumes independent
                                 samples; McNemar handles PAIRED categorical data

Simpson's Paradox           -> pooling across a hidden confounder can REVERSE
                                 what every subgroup individually shows

Cochran-Mantel-Haenszel     -> the fix: analyze WITHIN strata, don't naively pool

Cohen's Kappa                -> raw agreement overstates real agreement;
                                 kappa subtracts out what chance alone predicts
```

---

## 🎯 Why This Topic Matters

- **Odds ratio and relative risk were shown numerically diverging on common outcomes and converging on rare ones — the "rare disease assumption" demonstrated directly, not just cited** — §1 found OR=2.6667 vs RR=2.0000 (a gap of 0.6667) on a table with 20% baseline risk, versus OR=2.0081 vs RR=2.0000 (a gap of just 0.0081) on a table with the same 2x relative effect but a much lower ~0.6% baseline risk.
- **A from-scratch delta-method confidence interval was built for both OR and RR** — §1 found the rare-outcome table's OR 95% CI ([0.60, 6.69]) wide enough to include 1 despite a point estimate of 2.0 — a small-sample reminder that a "doubled" effect isn't automatically statistically significant.
- **McNemar's test was built from scratch and matched `statsmodels` exactly** — §2 found a manual chi-square statistic (13.3333) and p-value (0.000261) matching `statsmodels.stats.contingency_tables.mcnemar` to machine precision, correctly using only the 25-vs-5 discordant pairs and ignoring the 70 concordant ones.
- **The chi-square approximation and the exact binomial version of McNemar's test were cross-checked directly** — §2 found them agreeing closely (0.000261 vs 0.000325, a gap of just 0.000064) at a discordant count of 30 — mirroring Topic 05's approximation-vs-exact theme, now for paired data.
- **Simpson's Paradox was constructed with real, checkable numbers, not just described** — §3 built two hospitals where Treatment A won at **both** individually (93.10% vs 86.67% at Hospital 1; 73.00% vs 68.75% at Hospital 2), yet the naively pooled data showed Treatment B winning (82.57% vs 78.00%) — a genuine, verified reversal.
- **The Cochran-Mantel-Haenszel test correctly resolved the paradox by controlling for the hospital confounder** — §4 found the CMH-adjusted pooled odds ratio (1.4468, favoring Treatment A) correctly agreeing with what both individual hospitals showed — directly reversing the naive pooled conclusion from §3.
- **Cohen's Kappa was built from scratch and matched `sklearn` exactly, then shown correctly separating real skill from illusory agreement** — §5 found two moderately skilled raters (85% individual accuracy) reaching 79.00% raw agreement with kappa=0.6469 ("substantial" agreement), while §6 found two **purely chance-guessing** raters reaching a deceptively similar 76.00% raw agreement but a kappa of just **0.0964** — correctly revealing that most of that raw agreement was just both raters favoring the same imbalanced majority category, not real concordance.
- **Everything closed with a striking, real historical effect** — §7 found the Titanic dataset's female-vs-male survival odds ratio at **12.35** (95% CI [8.90, 17.14], firmly excluding 1) — the "women and children first" evacuation protocol, quantified directly from actual passenger records.

---

## 🧮 Mathematical Explanation

### 1. Odds ratio and relative risk

$$OR = \frac{ad}{bc}, \qquad RR = \frac{a/(a+b)}{c/(c+d)}$$

$$SE_{\log OR} = \sqrt{\frac{1}{a}+\frac{1}{b}+\frac{1}{c}+\frac{1}{d}}, \qquad
SE_{\log RR} = \sqrt{\frac{1}{a}-\frac{1}{a+b}+\frac{1}{c}-\frac{1}{c+d}}$$

Both confidence intervals are built on the log scale (via the delta method) then exponentiated back, since odds ratios and relative risks are naturally right-skewed and bounded below by 0.

### 2. McNemar's test

$$\chi^2 = \frac{(b-c)^2}{b+c}, \qquad df=1$$

Uses only the two *discordant* cells $b$ (changed one direction) and $c$ (changed the other) — exactly analogous to the paired t-test (Topic 04) using only within-subject differences.

### 3. Cochran-Mantel-Haenszel test

$$OR_{CMH} = \frac{\sum_k a_k d_k / n_k}{\sum_k b_k c_k / n_k}$$

A weighted combination of the odds ratio *within each stratum* $k$ — properly controlling for a confounder that a naive pooled 2x2 table ignores entirely, which is exactly what causes Simpson's Paradox.

### 4. Cohen's Kappa

$$\kappa = \frac{p_o - p_e}{1 - p_e}$$

$p_o$ is raw observed agreement; $p_e$ is the agreement expected if both raters' category choices were independent, computed from the marginal distributions. Kappa is 0 when agreement is no better than chance, 1 for perfect agreement.

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| Odds ratio | $ad/bc$ | §1 |
| Relative risk | $\frac{a/(a+b)}{c/(c+d)}$ | §1 |
| McNemar's chi2 | $(b-c)^2/(b+c)$ | §2 |
| CMH pooled OR | $\sum a_kd_k/n_k \big/ \sum b_kc_k/n_k$ | §4 |
| Cohen's Kappa | $(p_o-p_e)/(1-p_e)$ | §5-§6 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Treating odds ratio and relative risk as interchangeable** — §1's 0.6667 gap on a common outcome (20% baseline risk) is a real, practically significant difference; only use OR as an RR approximation when the outcome is genuinely rare.
2. **Applying chi-square independence (Topic 05) to paired/matched data** — this violates the independence assumption outright; §2's McNemar's test is the correct tool whenever the same subjects are measured twice.
3. **Trusting a naively pooled 2x2 table without checking for a lurking confounder** — §3's Simpson's Paradox is a real, constructed reversal, not a contrived edge case; always consider whether a natural stratifying variable (site, time period, subgroup) could be hiding the true relationship.
4. **Reporting raw percent agreement as evidence of real inter-rater reliability** — §6's 76.00% raw agreement from two *completely independent, chance-guessing* raters is the clearest possible warning; always report Cohen's Kappa (or a similar chance-corrected measure) instead.
5. **Ignoring how category imbalance inflates expected chance agreement** — §5-§6 both used an imbalanced 60-85% majority category, and in both cases $p_e$ was substantial (0.4053 and 0.7344 respectively); the more imbalanced the categories, the more raw agreement can mislead.
6. **Assuming a "significant" odds ratio automatically means a large practical effect** — §1's rare-outcome table had OR=2.0081 but a 95% CI of [0.60, 6.69], entirely consistent with no real effect at all; always report the CI alongside the point estimate.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `statsmodels.stats.contingency_tables.mcnemar` | McNemar's test reference (both chi-square and exact variants) |
| `statsmodels.stats.contingency_tables.StratifiedTable` | Cochran-Mantel-Haenszel test and pooled odds ratio |
| `sklearn.metrics.cohen_kappa_score` | Cohen's Kappa reference |
| `scipy.stats.norm.ppf` | Delta-method confidence intervals for OR and RR |
| `pandas.crosstab` | Building the real-data Titanic 2x2 table |

---

## 📝 Self-Test Exercises

1. Section 1 found the OR-RR gap shrinking from 0.6667 (20% baseline risk) to 0.0081 (~0.6% baseline risk) for the same underlying 2x relative effect. Using the OR and RR formulas from §1 of the Mathematical Explanation, explain algebraically why $a+b \approx b$ (i.e., $a$ is small relative to $b$) is what makes OR and RR converge as the outcome becomes rare.
2. Section 3's Simpson's Paradox showed Treatment A winning at both hospitals individually but Treatment B winning pooled. Using the raw counts in §3, identify which hospital contributed disproportionately more patients to Treatment B in the pooled total, and explain how that imbalance drives the reversal.
3. Section 4 found the CMH-adjusted pooled odds ratio (1.4468) favoring Treatment A, correctly matching both individual hospitals. Using the CMH formula from §3 of the Mathematical Explanation, explain in your own words why weighting by within-stratum information (rather than just adding raw counts across strata) prevents the confounding that caused the paradox in the first place.
4. Section 6 found two chance-only raters reaching kappa=0.0964, not exactly 0. Using the kappa formula, explain why kappa for two truly independent, random raters is expected to be *close to* zero but not necessarily *exactly* zero on any single simulated dataset.
5. Section 7 found the Titanic odds ratio for female survival at 12.35 with a 95% CI of [8.90, 17.14]. Using what you learned in §1 about the OR-RR relationship, would you expect the corresponding relative risk (3.93, computed in the same cell) to be smaller or larger than the odds ratio here, and why does that direction make sense given how common survival was in this dataset (not a rare outcome)?

---

## 📓 Notebook

30 executed code cells: odds ratio and relative risk built from scratch with delta-method confidence intervals, a direct common-vs-rare-outcome comparison confirming the "rare disease assumption," McNemar's test built from scratch and matched exactly to `statsmodels` (plus an exact-vs-chi-square-approximation cross-check), a fully constructed Simpson's Paradox with two hospitals each favoring Treatment A individually but Treatment B winning pooled, a Cochran-Mantel-Haenszel test correctly reversing that naive pooled conclusion, Cohen's Kappa built from scratch and matched exactly to `sklearn.metrics.cohen_kappa_score`, a direct demonstration of kappa correctly flagging two chance-only raters' deceptively high raw agreement as near-zero real agreement, and a real historical odds-ratio analysis of Titanic survival by sex:

➡️ **[14_categorical_data_analysis.ipynb](14_categorical_data_analysis.ipynb)**

---

## 📚 Further Reading

- [statsmodels.stats.contingency_tables documentation](https://www.statsmodels.org/stable/contingency_tables.html)
- [Simpson (1951): The Interpretation of Interaction in Contingency Tables](https://www.jstor.org/stable/2984065)
- [Mantel & Haenszel (1959): Statistical Aspects of the Analysis of Data from Retrospective Studies of Disease](https://academic.oup.com/jnci/article-abstract/22/4/719/912994)
- [Landis & Koch (1977): The Measurement of Observer Agreement for Categorical Data](https://www.jstor.org/stable/2529310)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.14_Categorical_Data_Analysis&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
