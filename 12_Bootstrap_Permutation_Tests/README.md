# 🔄 Bootstrap & Permutation Tests

> Status: ✅ Complete — [Open the notebook →](12_bootstrap_permutation_tests.ipynb)

Topic 12 of the Statistical Inference & Hypothesis Testing repo. Topic 03 introduced the bootstrap percentile confidence interval briefly. This notebook goes deeper on two related but distinct resampling philosophies: the **bootstrap** (resample *with replacement* to estimate a statistic's uncertainty) and **permutation tests** (shuffle labels *without replacement* to test a null hypothesis directly, with no distributional assumptions at all). It builds a more accurate bootstrap CI method (BCa) from scratch, builds permutation tests from scratch and validates them against Topic 04's t-test, and — honestly — shows a case where the bootstrap itself breaks down.

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
BOOTSTRAP: resample the data WITH REPLACEMENT
   -> estimates a statistic's SAMPLING DISTRIBUTION / uncertainty
   -> answers: "how much would this statistic vary if I redid the study?"

PERMUTATION TEST: shuffle group LABELS WITHOUT REPLACEMENT
   -> simulates the NULL HYPOTHESIS directly
   -> answers: "if there were truly no effect, how extreme would this
                statistic look just from random labeling?"

Both need NO distributional assumption -- but they answer DIFFERENT questions,
and the bootstrap has a real, known failure mode for extreme-value statistics.
```

---

## 🎯 Why This Topic Matters

- **The bootstrap standard error was validated against the classical formula where one exists** — §1 found the bootstrap SE (1.7732) matching the textbook $s/\sqrt{n}$ formula (1.8078) to within **0.0346** — reassuring agreement establishing the bootstrap works where it can be checked.
- **BCa was built from scratch and shown converging to the same answer as the simpler percentile method when the statistic is well-behaved** — §2 found BCa and percentile 95% CIs for a sample mean differing by only **0.16** on the lower bound and **0.02** on the upper — nearly identical, as expected for a symmetric statistic.
- **BCa's corrections were shown to actually matter on a genuinely skewed statistic** — §3 found BCa and percentile diverging meaningfully on a sample variance from skewed data (skewness=1.99): lower bounds differed by **0.47** and upper bounds by **0.63** — a real, non-trivial difference exactly where BCa's bias/acceleration corrections are designed to kick in.
- **The bootstrap was applied to a statistic with no simple closed-form uncertainty formula, and cross-validated against an independent method** — §4 found the bootstrap 95% CI for a correlation coefficient ([0.370, 0.759]) closely matching Topic 08's Fisher's z-transformation CI ([0.333, 0.758]) — two completely different derivations agreeing closely.
- **The permutation test was built from scratch and matched the t-test almost exactly** — §5 found a permutation test's p-value (0.0041) nearly identical to the reference two-sample t-test's (0.0039), a difference of just **0.0002**, despite making zero Normality assumption.
- **Permutation test calibration was verified directly, not assumed** — §6 found an empirical Type I error rate of **7.20%** against a 5% nominal target (with the gap consistent with expected noise from a modest 500-simulation, 500-permutation budget), and found the permutation test's power under skewed data (**36.60%**) closely matching the t-test's (**35.80%**) — genuinely competitive even without the Normality assumption the t-test relies on.
- **The permutation test extended naturally to correlation, again agreeing closely with a parametric alternative** — §7 found a permutation-based correlation test's p-value (0.0000, i.e. no permuted correlation among 10,000 was as extreme) consistent with the parametric significance test's p-value (0.000076).
- **The bootstrap's real, known failure mode was demonstrated directly, not just described** — §8 showed that for the sample maximum, **every single one** of 5,000 bootstrap resamples' maxima fell at or below the observed sample max (95.98) — none could ever exceed it, since resampling with replacement can't invent a value not already in the data — producing a 95% "CI" of [93.68, 95.98] that structurally cannot capture the true population bound of 100.
- **Everything was re-applied to real data** — §9 built a BCa CI for the median tip percentage (a statistic with no closed form) and found weekend-vs-weekday tipping showing no significant difference by either the permutation test (p=0.667) or Welch's t-test (p=0.597) — full agreement between the assumption-free and parametric approaches.

---

## 🧮 Mathematical Explanation

### 1. Bootstrap standard error

$$\widehat{SE}_{\text{boot}} = \sqrt{\frac{1}{B-1}\sum_{b=1}^{B}\left(\hat\theta^{*}_b - \bar{\hat\theta}^{*}\right)^2}$$

The standard deviation of the statistic computed across $B$ bootstrap resamples — no formula for the statistic's true sampling variance required.

### 2. BCa confidence interval

$$\hat{z}_0 = \Phi^{-1}\left(\frac{\#\{\hat\theta^*_b < \hat\theta\}}{B}\right), \qquad
\hat{a} = \frac{\sum_i(\bar\theta_{(\cdot)} - \theta_{(i)})^3}{6\left[\sum_i(\bar\theta_{(\cdot)}-\theta_{(i)})^2\right]^{3/2}}$$

$\hat{z}_0$ (bias-correction, from how often the bootstrap distribution sits below the original estimate) and $\hat{a}$ (acceleration, from jackknife resampling — leave-one-out estimates $\theta_{(i)}$) together adjust which percentiles of the bootstrap distribution are used as CI endpoints, correcting for skewness the plain percentile method ignores.

### 3. Permutation test (two-sample)

$$p = \frac{1}{P}\sum_{p=1}^{P} \mathbb{1}\left[|\bar{X}^{*}_{p,A} - \bar{X}^{*}_{p,B}| \ge |\bar{X}_A - \bar{X}_B|\right]$$

Shuffle the combined data's group labels $P$ times, recompute the test statistic each time, and see how often the *label-shuffled* statistic is at least as extreme as the *actually observed* one — a direct simulation of the null hypothesis of exchangeability.

### 4. Why the bootstrap fails for extreme values

For the sample maximum $X_{(n)}$, every bootstrap resample is drawn *from the observed data itself*, so no resample's maximum can exceed $X_{(n)}$. The bootstrap distribution of the maximum is therefore degenerate at its upper end by construction — a structural limitation, not a sampling artifact that more resamples would fix.

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| Bootstrap SE | $\text{std}(\hat\theta^*_1,\dots,\hat\theta^*_B)$ | §1 |
| BCa bias-correction | $\Phi^{-1}(\#\{\hat\theta^*<\hat\theta\}/B)$ | §2 |
| BCa acceleration | jackknife-based $\hat{a}$ | §2 |
| Permutation p-value | fraction of shuffles $\ge$ observed | §5, §7 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Using the plain percentile bootstrap CI on a skewed statistic without checking BCa** — §3's 0.47-0.63 shift in CI bounds is a real, measurable consequence, not a theoretical nuance.
2. **Confusing what bootstrap and permutation tests actually estimate** — bootstrap quantifies a statistic's *uncertainty*; permutation simulates the *null hypothesis* directly. They're not interchangeable tools for the same job, even though both resample the data.
3. **Applying the standard bootstrap to extreme-value statistics** — §8's sample-maximum demonstration is a genuine, known failure mode; the bootstrap CI is guaranteed too narrow and biased low for statistics like min/max, which need specialized extreme-value theory instead.
4. **Assuming a permutation test needs Normality to be trustworthy** — §6's power comparison under skewed data (36.60% vs the t-test's 35.80%) shows it staying competitive precisely where the t-test's assumption is violated.
5. **Trusting a single Type I error simulation without accounting for its own Monte Carlo noise** — §6's 7.20% (vs 5% nominal) empirical rate came from only 500 repeats; a wider simulation budget would narrow the noise band, consistent with every earlier topic's calibration checks in this series.
6. **Forgetting the bootstrap's inherent conservativism for skewed statistics extends beyond variance** — any statistic whose sampling distribution is asymmetric (ratios, extreme quantiles, etc.) is a candidate for BCa over plain percentile, not just the variance example shown here.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `np.random.default_rng().choice(..., replace=True)` | Bootstrap resampling |
| `np.random.default_rng().permutation` | Label shuffling for permutation tests |
| `scipy.stats.norm.ppf` / `.cdf` | BCa's bias-correction and acceleration adjustments |
| `np.percentile` | Percentile and BCa-adjusted-percentile CI endpoints |
| `scipy.stats.ttest_ind` / `.pearsonr` | Parametric reference tests for cross-validation |

---

## 📝 Self-Test Exercises

1. Section 3 found BCa and percentile CIs diverging by 0.47-0.63 on a skewed sample variance, while §2 found them nearly identical (0.02-0.16) on a roughly-symmetric sample mean. Using the BCa formula from §2 of the Mathematical Explanation, explain what property of the bootstrap distribution's shape causes $\hat{a}$ (the acceleration term) to matter more in the skewed case.
2. Section 5 found the permutation test's p-value (0.0041) nearly identical to the t-test's (0.0039) on Normally-distributed data. Would you expect this close agreement to hold if the underlying data were instead drawn from a heavily skewed distribution? Use what you learned about the permutation test's assumption-free construction to justify your answer.
3. Section 8 showed the bootstrap's 95% "CI" for the sample maximum ([93.68, 95.98]) entirely missing the true population bound of 100. Using the explanation in §4 of the Mathematical Explanation, describe what would need to be true about the underlying resampling scheme for a bootstrap-like method to correctly capture the true maximum's uncertainty.
4. Section 4 found the bootstrap CI for a correlation coefficient closely matching Topic 08's Fisher's z-transformation CI. Given that these two methods make completely different assumptions (Fisher's z assumes approximate Normality of a transformed statistic; the bootstrap assumes nothing), what does their close agreement suggest about how well the Normality assumption holds for this particular sample size and correlation strength?
5. Section 9 found the real weekend-vs-weekday tipping comparison not significant by either the permutation test (p=0.667) or the t-test (p=0.597) — a fairly large p-value gap between the two methods despite both reaching the same qualitative conclusion. Propose a reason why the two methods' exact p-values might differ somewhat even when the underlying decision is the same.

---

## 📓 Notebook

30 executed code cells: a bootstrap standard error validated against the classical formula, a BCa confidence interval built from scratch (including the jackknife-based acceleration term) and shown converging to the percentile method on a symmetric statistic but diverging meaningfully on a skewed one, a bootstrap CI for a correlation coefficient cross-validated against Topic 08's Fisher's z method, a two-sample permutation test built from scratch and matched almost exactly to the reference t-test, a Type I error and power calibration check for the permutation test (including a direct power comparison against the t-test under skewed data), a permutation test for correlation cross-validated against the parametric significance test, an honest demonstration of the bootstrap's real failure mode for the sample maximum, and a full application of BCa and permutation methods to the real 244-row `tips` dataset:

➡️ **[12_bootstrap_permutation_tests.ipynb](12_bootstrap_permutation_tests.ipynb)**

---

## 📚 Further Reading

- [Efron & Tibshirani, *An Introduction to the Bootstrap*](https://www.taylorfrancis.com/books/mono/10.1201/9780429246593)
- [Efron (1987): Better Bootstrap Confidence Intervals (the original BCa paper)](https://www.jstor.org/stable/2289144)
- [Good, *Permutation, Parametric, and Bootstrap Tests of Hypotheses*](https://link.springer.com/book/10.1007/b138696)
- Fisher (1935): *The Design of Experiments* — the original source of the permutation test idea

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.12_Bootstrap_Permutation_Tests&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
