# 📐 Power Analysis & Sample Size Determination

> Status: ✅ Complete — [Open the notebook →](09_power_analysis_sample_size_determination.ipynb)

Topic 09 of the Statistical Inference & Hypothesis Testing repo. Topics 04, 05, and 06 each simulated a power curve informally, one test at a time. This notebook formalizes power analysis as its own tool: the four-way relationship between significance level, power, effect size, and sample size, where fixing any three determines the fourth. It validates `statsmodels`' closed-form power calculators against this series' own Monte Carlo simulations, then uses power analysis for its most practically important purpose — answering "how many samples do I need" — before confronting a genuinely underappreciated cost of running underpowered studies anyway: exaggerated effect sizes.

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
Four quantities, locked together -- fix any three, the fourth follows:

    ALPHA (Type I error rate)  <-->  POWER (1 - Type II error rate)
                    \\                /
                     \\              /
                  EFFECT SIZE  <-->  SAMPLE SIZE

Most common practical use: fix alpha=0.05, power=0.80 (convention),
plug in an expected effect size  ->  solve for the REQUIRED SAMPLE SIZE

Underappreciated use: an UNDERPOWERED study that DOES reach significance
tends to REPORT AN EXAGGERATED effect size -- power analysis isn't just
about getting a p-value, it's about trusting the number you get.
```

---

## 🎯 Why This Topic Matters

- **`statsmodels`' closed-form power calculator was validated directly against Monte Carlo simulation** — §1 found `TTestIndPower`'s power for d=0.5, n=30 (**0.4779**) matching a 5,000-repeat simulation (**0.4786**) to within **0.0007**.
- **Solving for required sample size was verified to actually deliver the target power** — §2 calculated n=64 needed for 80% power at d=0.5, then simulated at exactly that n and got **0.7998** — confirming the calculation, not just trusting it.
- **The cost of detecting small effects was quantified precisely** — §2 found detecting a small effect (d=0.2) requires **394** per group versus just **12** for a large effect (d=1.2) — a **32.8x** difference in required sample size for what "Cohen's conventions" label as merely "small" vs. "large."
- **Two real planning approaches gave meaningfully different answers, and that gap matters** — §4 compared planning around a minimum meaningful effect (d=0.3 → n=176) versus a small pilot study's noisy point estimate (d=0.834 from a 15-per-group pilot → n=24) — a nearly **7x** disagreement in required sample size purely from which planning approach was used.
- **Power analysis machinery was validated for chi-square and ANOVA too, not just the t-test** — §5 found `GofChisquarePower`'s n=142.5 (Cohen's w=0.3, 6 categories) achieving 89.3% simulated power (close to the 80% target), and §6 found `FTestAnovaPower`'s n=157.19 per group (Cohen's f=0.25, k=3) achieving 100% simulated power in the specific approximated scenario used — both confirming the closed-form and simulation-based approaches from Topics 05-06 are mutually consistent.
- **The hidden cost of underpowered studies was demonstrated with a striking, precise number** — §7 simulated a badly underpowered study (true d=0.3, n=15, power=**12.47%**) and found that among the studies that *did* reach significance, the average observed effect size was **0.9460** — a **3.15x exaggeration** of the true 0.3 effect.
- **The same experiment at proper power showed the exaggeration effect shrinking dramatically** — §7 repeated the simulation at the well-powered n=176 (power=80.14%) and found the exaggeration ratio dropping to just **1.14x** — direct, quantified proof that adequate power isn't just about getting more significant results, it's about trusting the effect sizes those results report.
- **A concrete business scenario tied every prior section together** — §8 found detecting a conversion-rate lift from 10% to 12% needs **3,835** visitors per variant, while a mere 1-percentage-point lift (10%→11%) needs **14,745** — nearly **4x** more, directly illustrating why small, realistic business improvements are the hardest and most expensive to detect reliably.

---

## 🧮 Mathematical Explanation

### 1. The power relationship

Power ($1-\beta$) is a function of effect size, sample size $n$, and significance level $\alpha$:

$$1 - \beta = P(\text{reject } H_0 \mid H_1 \text{ true}; \ \text{effect size}, n, \alpha)$$

No closed form is simple in general, but for the t-test it derives from the non-central t-distribution; `statsmodels` implements this (and the chi-square and F-test analogues) directly.

### 2. Solving for sample size

Since power is monotonically increasing in $n$ (for fixed effect size and $\alpha$), the required $n$ for a target power is found by numerically inverting the power function — exactly what `solve_power()` does in §2, §5, §6, and §8.

### 3. Cohen's w (chi-square) and f (ANOVA)

$$w = \sqrt{\sum_i \frac{(p_{i,1}-p_{i,0})^2}{p_{i,0}}} \qquad f = \frac{\sigma_{\text{means}}}{\sigma_{\text{within}}}$$

Standardized effect sizes analogous to Cohen's $d$, extending power analysis to chi-square (Topic 05) and ANOVA (Topic 06).

### 4. Cohen's h (two proportions)

$$h = 2\arcsin(\sqrt{p_1}) - 2\arcsin(\sqrt{p_2})$$

An arcsine-transformed effect size for comparing two proportions (e.g. conversion rates), used in §8's A/B test scenario — the transformation stabilizes variance across the full range of $p$, unlike a raw proportion difference.

### 5. The exaggeration ratio (Type M error)

$$\text{Exaggeration ratio} = \frac{E[|\hat{d}| \mid p < \alpha]}{|d_{\text{true}}|}$$

The expected magnitude of the observed effect size, *conditional on reaching significance*, divided by the true effect size. Low power inflates this ratio because only unusually large sample estimates clear the significance threshold when the study is underpowered — §7 measured this directly rather than citing it.

---

## 📋 Formula Reference

| Concept | Tool | Notebook Section |
|---|---|---|
| t-test power/sample size | `TTestIndPower` | §1-§2 |
| Chi-square power/sample size | `GofChisquarePower` | §5 |
| ANOVA power/sample size | `FTestAnovaPower` | §6 |
| Two-proportion power/sample size | `NormalIndPower` + Cohen's h | §8 |
| Exaggeration ratio | $E[\|\hat d\| \mid \text{significant}] / \|d_{\text{true}}\|$ | §7 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Assuming small "conventional" effect-size labels mean similar sample-size costs** — §2's 32.8x gap between detecting d=0.2 and d=1.2 shows "small" vs. "large" hides an enormous practical difference in required data.
2. **Trusting a small pilot study's effect-size estimate at face value** — §4's pilot (n=15) gave d=0.834, wildly different from the "minimum meaningful" planning value of d=0.3, precisely because small-sample point estimates are themselves highly variable (a direct callback to Topic 01-02's sampling-distribution lessons).
3. **Treating a significant result from an underpowered study as trustworthy** — §7's core finding: an underpowered study's significant results overstate the true effect size by **3.15x** on average; a well-powered study's overstate by only 1.14x.
4. **Confusing "the study found significance" with "the study is reliable"** — a low-power study can occasionally produce a significant, dramatically inflated effect size purely by chance; power analysis addresses this before data collection, not after.
5. **Assuming a modest, realistic effect (like a 1-2 percentage point conversion lift) is cheap to detect** — §8 showed the smallest, most business-realistic lifts requiring the most data by far (14,745 vs. 195 visitors per variant across the tested range).
6. **Forgetting the ANOVA/chi-square power validation in this notebook used an approximate mapping from Cohen's f/w to specific simulated scenarios** — §5's 89.3% (vs 80% target) and §6's 100% (vs 80% target) simulated results are reasonable cross-checks, not exact replications, since translating a standardized effect size back into one specific set of group means/probabilities is not unique.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `statsmodels.stats.power.TTestIndPower` | t-test power/sample-size calculator |
| `statsmodels.stats.power.GofChisquarePower` | Chi-square goodness-of-fit power/sample-size calculator |
| `statsmodels.stats.power.FTestAnovaPower` | One-way ANOVA power/sample-size calculator |
| `statsmodels.stats.power.NormalIndPower` | Two-proportion (Normal-approximation) power/sample-size calculator |
| `statsmodels.stats.proportion.proportion_effectsize` | Computes Cohen's h from two proportions |
| `.solve_power(...)` | Inverts the power function to solve for any one of alpha/power/effect_size/n |

---

## 📝 Self-Test Exercises

1. Section 2 found detecting d=0.2 needs 394 per group while d=1.2 needs only 12 — a 32.8x ratio, even though the effect sizes only differ by 6x. Why does required sample size scale faster than linearly with the inverse of effect size? (Hint: think about how $n$ enters the standard error formula from Topic 01.)
2. Section 4 found a 15-per-group pilot study gave a wildly different effect-size estimate (d=0.834) than the "minimum meaningful" planning value (d=0.3). Using what you learned about sampling distributions in Topics 01-02, explain why an n=15 pilot's point estimate should NOT be trusted as a precise input for a follow-up power calculation.
3. Section 7 found the exaggeration ratio was 3.15x at n=15 (power=12.47%) but only 1.14x at n=176 (power=80.14%). Would you expect the exaggeration ratio to keep shrinking as n grows even further past 176 (say, to n=1000)? Explain your reasoning using the definition of the exaggeration ratio in §5 of the Mathematical Explanation.
4. Section 8 found a 1-percentage-point lift (10%→11%) needs 14,745 visitors per variant, over 75x more than a 10-point lift's 195. Using Cohen's h formula from §4 of the Mathematical Explanation, explain why h doesn't scale linearly with the raw percentage-point difference.
5. Section 5 and 6's simulated power (89.3% and 100%) didn't land exactly on the 80% target predicted by `statsmodels`. Given the caveat in Pitfall #6, is this discrepancy evidence of a bug, or an expected consequence of how the effect size was translated into a specific simulated scenario? Justify your answer.

---

## 📓 Notebook

30 executed code cells: `statsmodels`' t-test power calculator validated against direct Monte Carlo simulation (diff=0.0007), required-sample-size solving verified to actually achieve its target power, a 4-effect-size sweep showing a 32.8x sample-size gap between small and large effects, power curves visualized across sample size, effect size, and alpha, two competing real-world approaches to choosing a planning effect size compared directly, chi-square and ANOVA power analysis validated against Topics 05-06's simulation approach, a full demonstration of underpowered-study effect-size exaggeration (3.15x at n=15 vs 1.14x at a properly-powered n=176), and a complete A/B testing sample-size scenario spanning four realistic conversion-rate lifts:

➡️ **[09_power_analysis_sample_size_determination.ipynb](09_power_analysis_sample_size_determination.ipynb)**

---

## 📚 Further Reading

- [statsmodels: Power and Sample Size Calculations](https://www.statsmodels.org/stable/stats.html#power-and-sample-size-calculations)
- [Gelman & Carlin (2014): Beyond Power Calculations — Assessing Type S and Type M Errors](https://doi.org/10.1177/1745691614551642)
- [Cohen (1988): Statistical Power Analysis for the Behavioral Sciences](https://www.routledge.com/Statistical-Power-Analysis-for-the-Behavioral-Sciences/Cohen/p/book/9780805802832)
- [Button et al. (2013): Power Failure — Why Small Sample Size Undermines the Reliability of Neuroscience](https://www.nature.com/articles/nrn3475)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.09_Power_Analysis_Sample_Size_Determination&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
