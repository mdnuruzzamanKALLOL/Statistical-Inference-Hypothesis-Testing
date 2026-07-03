# 📏 Confidence Intervals

> Status: ✅ Complete — [Open the notebook →](03_confidence_intervals.ipynb)

Topic 03 of the Statistical Inference & Hypothesis Testing repo. Topics 01-02 built the sampling distribution and confirmed $SE = \sigma/\sqrt{n}$ empirically. A confidence interval turns that standard error into a *range* of plausible parameter values — this notebook builds one from scratch, then does the one thing most explanations skip: directly simulates thousands of intervals to check whether a "95% confidence interval" actually contains the true parameter 95% of the time, and confronts the single most common misinterpretation of what a CI means head-on.

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
Point estimate (X_bar)  +/-  margin of error  =  Confidence Interval
                                  |
                    margin = critical_value * SE
                                  |
   critical value from:  z (sigma known)  or  t (sigma estimated, wider)
                                  |
   THE KEY QUESTION: if you repeated this whole procedure 10,000 times,
   what fraction of the resulting intervals would contain the true parameter?
                                  |
        That fraction is what "95% confidence" actually means --
        NOT "95% probability the true value is in THIS ONE interval"
```

---

## 🎯 Why This Topic Matters

- **The 95% coverage claim was tested directly, not assumed** — §2 built 10,000 independent 95% CIs from repeated sampling and found **95.05%** of them actually contained the true mean — a difference of just 0.05 percentage points from the nominal target.
- **A visual sample of that simulation made the 5% failure rate concrete** — §2 displayed 100 of those intervals color-coded by whether they hit or missed; **6 out of 100** missed the true mean, close to the expected 5%.
- **The single most common CI misinterpretation was addressed directly, not just warned about** — §3 took one already-realized interval and showed it is now simply a **fact** (contains the true mean: `True`), pointing out the "95% probability" language only ever applies to the *procedure* measured across many samples (§2's 95.05%), never to one already-computed interval.
- **The t-distribution's extra width was shown to matter most exactly where intuition says it should** — §4 found the ratio of t-critical to z-critical values shrinking from **1.42x at n=5** down to **1.01x at n=100** — the correction for estimating $\sigma$ is large at small samples and nearly vanishes by n=100.
- **A concrete, practical "how many samples do I need" calculation was included** — §5 computed that halving a target margin of error from 10 to 5 requires the sample size to roughly **quadruple** (from 9 to 35, both using the real $\sigma=15$ from this notebook's population) — the $1/\sqrt{n}$ relationship made directly actionable.
- **Bootstrap intervals were shown to agree closely with the parametric method where both apply, and to extend naturally where the parametric method has no formula** — §6 found the bootstrap 95% CI for the mean (width 9.70) matched the t-interval (width 10.24) to within 0.54, then built a CI for the **median** — a statistic with no simple closed form — directly from the same resampling machinery.
- **The Wald proportion interval's failure was measured, not cited from a textbook** — §7 found the Wald interval's *empirical coverage* was only **77.74%** against a nominal 95% target (n=30, true p=0.05), while the Wilson score interval reached **93.79%** — a 16-point real coverage gap from a formula still taught as the default in many intro courses.
- **An even starker Wald failure was constructed** — §7 found that with 0 successes out of 10 trials, the Wald interval collapses to the single point **[0.0000, 0.0000]** — a zero-width "confidence interval" — while the Wilson interval correctly stayed wide at [0.0000, 0.2775], acknowledging that observing zero events doesn't prove the true rate is exactly zero.
- **Every method was re-validated on real, messy data** — §9 applied the t-interval, bootstrap interval, and Wald/Wilson proportion intervals to seaborn's `tips` dataset (244 real dining parties), finding the t and bootstrap CIs for mean tip percentage agreed to within **0.0071** of each other, and Wald/Wilson nearly coincided at this larger, non-extreme real sample size — directly confirming §7's finding that Wald's failure mode is specific to small n and extreme p.

---

## 🧮 Mathematical Explanation

### 1. The basic confidence interval

$$\bar{X} \pm z_{\alpha/2} \cdot \frac{\sigma}{\sqrt{n}} \qquad \text{(} \sigma \text{ known)}$$

### 2. The correct frequentist interpretation

A $(1-\alpha)$ CI's guarantee is about the **construction procedure**, evaluated across hypothetical repeated sampling:

$$P(\hat{\theta}_{\text{lower}} \le \theta \le \hat{\theta}_{\text{upper}}) = 1 - \alpha \quad \text{(before the sample is drawn)}$$

Once a specific sample produces a specific numeric interval, $\theta$ either is or isn't inside it — there is no more randomness left to attach a probability to. §2's 95.05% empirical coverage across 10,000 *repeated* constructions is the correct object the "95%" describes; §3's single realized interval is not.

### 3. The t-interval when $\sigma$ is unknown

$$\bar{X} \pm t_{\alpha/2, \, n-1} \cdot \frac{s}{\sqrt{n}}$$

Substituting the sample standard deviation $s$ for the unknown $\sigma$ introduces extra uncertainty; the t-distribution's heavier tails (relative to Normal) widen the interval to compensate. §4 showed this correction shrinking from a 42% inflation at n=5 to a 1.2% inflation at n=100.

### 4. One-sided intervals

$$\bar{X} - z_{\alpha} \cdot \frac{\sigma}{\sqrt{n}} \quad \text{(lower bound only, full } \alpha \text{ in one tail)}$$

Putting the entire error budget in one tail makes a one-sided bound tighter than the corresponding two-sided bound's edge — §4's addendum confirmed the one-sided 95% lower bound (92.997) sits strictly above the two-sided lower bound (92.250) for the same data.

### 5. Minimum sample size for a target margin of error

$$n = \left(\frac{z_{\alpha/2} \cdot \sigma}{E}\right)^2$$

Solved directly for $n$ rather than for the interval — §5 used this to show precision requirements grow quadratically as the target margin shrinks.

### 6. The bootstrap percentile interval

Resample the observed data with replacement $B$ times, compute the statistic $\hat{\theta}^{*}_b$ each time, and take

$$\left[\hat{\theta}^{*}_{(\alpha/2)}, \ \hat{\theta}^{*}_{(1-\alpha/2)}\right]$$

the empirical $\alpha/2$ and $1-\alpha/2$ percentiles of the bootstrap distribution — no Normality assumption required. §6 applied this to both the mean (validating against the t-interval) and the median (where no simple closed form exists).

### 7. Wald vs. Wilson score interval for a proportion

$$\text{Wald: } \hat{p} \pm z_{\alpha/2}\sqrt{\frac{\hat{p}(1-\hat{p})}{n}} \qquad
\text{Wilson: } \frac{\hat{p} + \frac{z^2}{2n}}{1+\frac{z^2}{n}} \pm \frac{z}{1+\frac{z^2}{n}}\sqrt{\frac{\hat{p}(1-\hat{p})}{n} + \frac{z^2}{4n^2}}$$

Wald's margin term vanishes entirely when $\hat{p}=0$ or $\hat{p}=1$ (since $\sqrt{\hat{p}(1-\hat{p})}=0$), producing the zero-width failure §7 demonstrated directly. Wilson's more complex form avoids this degenerate case by design.

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| z-interval (σ known) | $\bar{X} \pm z_{\alpha/2}\sigma/\sqrt{n}$ | §1 |
| t-interval (σ estimated) | $\bar{X} \pm t_{\alpha/2,n-1}s/\sqrt{n}$ | §4 |
| Minimum n for margin E | $(z_{\alpha/2}\sigma/E)^2$ | §5 |
| Bootstrap percentile CI | $[\hat\theta^*_{(\alpha/2)}, \hat\theta^*_{(1-\alpha/2)}]$ | §6 |
| Wald proportion CI | $\hat{p}\pm z_{\alpha/2}\sqrt{\hat{p}(1-\hat{p})/n}$ | §7 |
| Wilson score CI | see Mathematical Explanation §7 | §7 |
| Difference of means (Welch) | $(\bar{X}_1-\bar{X}_2)\pm t_{\alpha/2,df}\cdot SE_{\text{diff}}$ | §8 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Interpreting a realized CI as "95% probability the true value is in here"** — §3 addressed this directly: the probability statement applies to the *procedure* (§2's 95.05% coverage across repeated sampling), not to any single already-computed interval, which either does or doesn't contain the truth.
2. **Using the z-interval when $\sigma$ is estimated, not known** — §4 showed this understates the interval width, especially badly at small n (42% too narrow at n=5).
3. **Trusting the Wald interval for proportions at small n or extreme p** — §7's 77.74% actual coverage (vs. 95% nominal) and the zero-successes collapse to `[0,0]` are both real, measured failures, not edge-case trivia.
4. **Assuming a smaller margin of error just needs "a bit more data"** — §5 showed halving the margin requires roughly quadrupling n, not doubling it.
5. **Forgetting a one-sided interval needs a different critical value** — §4's addendum used the full $\alpha$ in one tail ($z_{0.95}$, not $z_{0.975}$); reusing the two-sided critical value for a one-sided claim overstates the confidence level.
6. **Assuming Wald's failure always applies** — §9 showed Wald and Wilson nearly coinciding on the real `tips` proportion (n=244, p̂≈0.16) — the failure mode is specific to small n combined with extreme p, not universal.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `scipy.stats.norm.ppf` / `.t.ppf` | Critical values for z- and t-based intervals |
| `scipy.stats.t.pdf` | Visualizing the t-distribution's heavier tails vs Normal |
| `np.percentile` | Bootstrap percentile interval endpoints |
| `np.random.default_rng().choice(..., replace=True)` | Bootstrap resampling |
| `np.random.default_rng().binomial` | Simulating repeated proportion samples for the Wald/Wilson coverage test |
| `seaborn.load_dataset("tips")` | Real dataset used for the Section 9 validation |

---

## 📝 Self-Test Exercises

1. Section 2 found 95.05% empirical coverage from 10,000 simulated 95% CIs. If you instead simulated 99% CIs the same way, would you expect the empirical coverage to be higher, lower, or the same — and would the *width* of each interval be wider or narrower than the 95% version? Use §5's confidence-level-vs-width chart to check your answer.
2. Using §4's t/z ratio table, explain in your own words why the ratio approaches but never quite reaches 1.0 as n grows, even at n=100.
3. Section 7 found the Wald interval's actual coverage was 77.74% against a 95% target. Using the Wald formula from the Mathematical Explanation, explain why $\hat{p}=0$ (as in the zero-successes case) makes the margin of error term vanish completely, regardless of how large or small n is.
4. Section 6 built a bootstrap CI for the median because "there is no simple closed-form formula" for one. Using the same resampling logic, describe (in words, no code needed) how you would build a bootstrap CI for a statistic like the *interquartile range* — another statistic without a simple closed form.
5. Section 8 built a 95% CI for the difference of two means that contained 0 (interval: [-6.94, 2.85]). What does "contains 0" tell you, informally, about whether the two groups' true means might be equal — and why does the notebook describe this as "informal evidence" rather than a formal conclusion?

---

## 📓 Notebook

30 executed code cells: a first z-based 95% CI construction, a 10,000-repeat Monte Carlo coverage simulation confirming 95.05% empirical coverage against a 95% nominal target, a 100-interval visual coverage display, a direct confrontation of the "95% probability" misinterpretation using one realized interval, a t-vs-z critical value comparison across four sample sizes (ratio shrinking from 1.42x to 1.01x), a one-sided confidence interval construction, a confidence-level-vs-width and sample-size-vs-width comparison, a minimum-sample-size-for-target-margin calculator, a bootstrap percentile CI for both the mean (validated against the t-interval) and the median (no closed form), a Wald-vs-Wilson proportion interval comparison including a zero-successes collapse-to-a-point failure case and a 10,000-repeat coverage simulation (Wald 77.74% vs Wilson 93.79% against a 95% target), a Welch's-t difference-of-two-means interval, and a full re-validation of the t/bootstrap/Wald/Wilson methods on the real 244-row `tips` dataset:

➡️ **[03_confidence_intervals.ipynb](03_confidence_intervals.ipynb)**

---

## 📚 Further Reading

- [scipy.stats.t / scipy.stats.norm documentation](https://docs.scipy.org/doc/scipy/reference/stats.html)
- [Wasserman, *All of Statistics*, Ch. 6: Confidence Intervals](https://link.springer.com/book/10.1007/978-0-387-21736-9)
- [Newcombe (1998): Two-sided confidence intervals for the single proportion (Wilson score)](https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1097-0258(19980430)17:8%3C857::AID-SIM777%3E3.0.CO;2-E)
- [Efron & Tibshirani, *An Introduction to the Bootstrap*](https://www.taylorfrancis.com/books/mono/10.1201/9780429246593)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.03_Confidence_Intervals&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
