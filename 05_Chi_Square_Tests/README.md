# 🔢 Chi-Square Tests

> Status: ✅ Complete — [Open the notebook →](05_chi_square_tests.ipynb)

Topic 05 of the Statistical Inference & Hypothesis Testing repo. Topics 01-04 all dealt with continuous measurements (means, variances). Chi-square tests instead work on **counts** of categorical outcomes — does an observed distribution of categories match an expected one (goodness-of-fit), or are two categorical variables associated (independence)? Both are built from scratch here, validated against `scipy.stats`, and the notebook then does what a typical intro treatment skips: directly demonstrates the chi-square approximation's behavior at small expected cell counts, and shows Fisher's Exact Test as the fix.

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
Two chi-square questions, same core statistic:

GOODNESS-OF-FIT: does one categorical variable match an expected distribution?
   observed counts vs expected counts  ->  chi2 = sum((O-E)^2 / E)

INDEPENDENCE: are two categorical variables associated?
   observed contingency table vs the table expected IF independent  ->  same chi2 formula

                          v
        chi2 is only an APPROXIMATION to the true discrete sampling distribution --
        reliable when expected counts are large, unreliable when they are small
                          v
        Fisher's Exact Test: no approximation, exact for 2x2 tables, always valid
```

---

## 🎯 Why This Topic Matters

- **Both tests were built from scratch and matched `scipy.stats` exactly** — §1 and §3 show manual chi-square statistics for goodness-of-fit and independence both agreeing with `scipy.stats.chisquare` / `.chi2_contingency` to `0.0e+00` difference.
- **Type I error calibration was verified empirically, the same discipline applied to the t-test in Topic 04** — §2 found an empirical false-positive rate of **5.31%** across 8,000 simulations of a truly fair die, against a 5% nominal target.
- **A genuinely loaded die was tested, and the result was reported honestly even when the test missed it** — §2 built a die where face 6 has true probability 0.25 (not 1/6) and got p=**0.106** on this particular sample — a real, correctly-reported Type II error (the bias existed but wasn't detected on this draw), a direct reminder that a single non-significant test never proves the null.
- **Yates' correction's actual size was measured, not just described as "a small adjustment"** — §4 found it shrank the chi-square statistic by **9.75%** at n=100, but only **0.50%** at n=2,000 — the correction genuinely matters most exactly where intro courses warn it does: small tables.
- **Cramer's V was shown to be sample-size-independent while raw chi2 is not, mirroring Cohen's d from Topic 04** — §5 scaled the same table by 1x, 5x, 20x, and 100x and found chi2 growing from 16.67 to 1666.67 (linear in n) while Cramer's V stayed essentially fixed at **0.408** across every scale.
- **The chi-square approximation's real behavior at small expected counts was measured against a genuine simulated exact p-value, not assumed** — §6 built a table with 2 of 4 expected counts below 5 and found the chi-square p-value (0.00165) sat noticeably below a 20,000-simulation exact-permutation p-value (0.00490) — a real, non-trivial shift even though both values were small in absolute terms.
- **Fisher's Exact Test was validated against chi-square at both small and large expected counts** — §7 found a **3.3x** proportional gap between chi-square (p=0.00165) and Fisher's Exact (p=0.00548) at small counts, versus near-perfect agreement (0.000045 vs 0.000083) once expected counts were comfortably above 5.
- **Power was simulated directly across sample sizes for a fixed true association** — §8 found power climbing from **44.03%** (n=20) to **85.00%** (n=50) to essentially **100%** (n=300) for the same underlying moderate association.
- **Everything was re-applied to real data** — §9 found smoking status significantly associated with day of week in the real `tips` dataset (χ²=25.79, p<0.0001, Cramer's V=0.325 — a medium-or-larger effect), while smoking status and lunch/dinner timing showed no significant association (p=0.391), with chi-square and Fisher's Exact agreeing closely since all expected counts were comfortably above 5.

---

## 🧮 Mathematical Explanation

### 1. Goodness-of-fit

$$\chi^2 = \sum_{i=1}^{k} \frac{(O_i-E_i)^2}{E_i}, \qquad df = k-1$$

Tests whether one categorical variable's observed counts match an expected distribution.

### 2. Test of independence

$$\chi^2 = \sum_{i,j} \frac{(O_{ij}-E_{ij})^2}{E_{ij}}, \qquad E_{ij} = \frac{(\text{row}_i)(\text{col}_j)}{n}, \qquad df=(r-1)(c-1)$$

Both formulas are structurally the same statistic — only how $E$ is computed differs.

### 3. Yates' continuity correction

$$\chi^2_{\text{Yates}} = \sum_{i,j} \frac{(|O_{ij}-E_{ij}|-0.5)^2}{E_{ij}}$$

Compensates for approximating a discrete count distribution with a continuous chi-square curve, specifically for 2x2 tables. §4 measured this correction's real, shrinking-with-n magnitude directly.

### 4. Cramer's V

$$V = \sqrt{\frac{\chi^2/n}{\min(r-1,c-1)}}$$

A sample-size-independent standardized association strength, ranging 0 (none) to 1 (perfect) — the categorical-data analogue of Cohen's d.

### 5. Fisher's Exact Test

Computes the exact probability of a 2x2 table (or one more extreme) directly from the hypergeometric distribution, conditional on the observed row/column margins — no large-sample approximation at all, valid at any sample size.

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| Goodness-of-fit chi2 | $\sum(O_i-E_i)^2/E_i$, $df=k-1$ | §1 |
| Independence chi2 | $\sum(O_{ij}-E_{ij})^2/E_{ij}$, $df=(r-1)(c-1)$ | §3 |
| Yates' correction | $\sum(\|O_{ij}-E_{ij}\|-0.5)^2/E_{ij}$ | §4 |
| Cramer's V | $\sqrt{(\chi^2/n)/\min(r-1,c-1)}$ | §5 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Treating a single non-significant result as proof of no effect** — §2's loaded-die test (a real, constructed bias) still produced p=0.106 on one draw; a single test missing an effect is expected some fraction of the time, exactly what power (§8) quantifies.
2. **Ignoring the expected-count-≥5 rule of thumb** — §6 showed a genuine, measurable divergence between the chi-square approximation and the true exact p-value when 2 of 4 expected counts were below 5.
3. **Trusting chi-square blindly on small 2x2 tables instead of using Fisher's Exact** — §7's 3.3x proportional p-value gap at small counts is real, not a theoretical curiosity; Fisher's Exact costs nothing to compute for 2x2 tables and is always valid.
4. **Comparing raw chi-square statistics across tables of different sample sizes** — §5 showed chi2 growing 100x from a 100x sample-size increase alone, with the true association strength completely unchanged; always report Cramer's V (or another normalized effect size) alongside chi2.
5. **Forgetting Yates' correction exists (or over-relying on it)** — §4 showed it matters most at small n (9.75% shrinkage at n=100) and becomes nearly irrelevant at larger n (0.50% at n=2,000); neither "always apply it" nor "always ignore it" is correct.
6. **Assuming goodness-of-fit and independence are different statistics** — both use the identical $\sum(O-E)^2/E$ formula; only the source of the expected counts differs (a hypothesized distribution vs. row/column margins).

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `scipy.stats.chisquare` | Goodness-of-fit test reference |
| `scipy.stats.chi2_contingency` | Independence test reference (with/without Yates' correction) |
| `scipy.stats.fisher_exact` | Exact 2x2 test, no large-sample approximation |
| `scipy.stats.chi2.cdf` | Manual p-value computation from the chi-square statistic |
| `np.random.default_rng().multinomial` | Simulated categorical count data |
| `np.random.default_rng().multivariate_hypergeometric` | Exact-simulation null distribution under fixed margins |
| `pandas.crosstab` | Building a real-data contingency table |

---

## 📝 Self-Test Exercises

1. Section 2 found the loaded die's test produced p=0.106, failing to detect a real bias on that particular sample. Using §8's power-simulation logic, what would you expect to happen to the chance of correctly detecting this same bias if the number of rolls were increased from 300 to 3,000?
2. Using the Yates' correction formula from §3 of the Mathematical Explanation, explain why the correction's effect (the $-0.5$ term) matters proportionally less as the table's total count grows.
3. Section 5 found Cramer's V stayed at ~0.408 across a 100x range of sample sizes while raw chi2 grew 100x. Using the Cramer's V formula, explain algebraically why dividing chi2 by $n$ cancels out this sample-size effect.
4. Section 6 found the chi-square approximation's p-value (0.00165) was noticeably below the simulated exact p-value (0.00490) at small expected counts. Based on this direction, does the chi-square approximation tend to overstate or understate significance when expected counts are small — and why would this matter for someone deciding whether to reject $H_0$ at $\alpha=0.05$?
5. Section 9 found smoking status significantly associated with day of week (p<0.0001) but not with lunch/dinner timing (p=0.391). Propose one plausible real-world explanation for why day-of-week might correlate with smoking status while time-of-day does not, using the actual contingency table counts from §9 as evidence.

---

## 📓 Notebook

30 executed code cells: a goodness-of-fit test built from scratch and matched exactly to `scipy.stats.chisquare` on a simulated fair-die roll, a Type I error rate simulation (5.31% empirical vs 5% nominal over 8,000 runs) plus an honestly-reported single-run miss on a genuinely loaded die, a test of independence built from scratch and matched exactly to `scipy.stats.chi2_contingency`, a direct measurement of Yates' correction's shrinking-with-n effect, a Cramer's V implementation shown staying constant (~0.408) across a 100x sample-size range while raw chi2 scaled linearly, a small-expected-count table compared against a 20,000-simulation exact null distribution, a Fisher's Exact Test comparison at both small (3.3x p-value gap) and large (near-identical) expected counts, a chi-square independence power simulation climbing from 44% to 100% power across four sample sizes, and a full application of goodness-of-fit, independence, Cramer's V, and Fisher's Exact methods to the real 244-row `tips` dataset:

➡️ **[05_chi_square_tests.ipynb](05_chi_square_tests.ipynb)**

---

## 📚 Further Reading

- [scipy.stats.chisquare / chi2_contingency / fisher_exact documentation](https://docs.scipy.org/doc/scipy/reference/stats.html)
- [Wasserman, *All of Statistics*, Ch. 10: Hypothesis Testing (chi-square section)](https://link.springer.com/book/10.1007/978-0-387-21736-9)
- [Yates (1934): Contingency Tables Involving Small Numbers and the χ² Test](https://www.jstor.org/stable/2983604)
- [Fisher (1922): On the Interpretation of χ² from Contingency Tables](https://www.jstor.org/stable/2340521)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.05_Chi_Square_Tests&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
