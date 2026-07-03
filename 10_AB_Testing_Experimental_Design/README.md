# 🧬 A/B Testing & Experimental Design

> Status: ✅ Complete — [Open the notebook →](10_ab_testing_experimental_design.ipynb)

Topic 10 of the Statistical Inference & Hypothesis Testing repo. Topic 09 answered "how many samples do I need" for an A/B test. This notebook covers what happens around that number: designing a proper randomized experiment, running the two-proportion test correctly, and — the real substance of this notebook — directly measuring two of the most common real-world ways A/B tests silently produce wrong conclusions: peeking at results before the planned sample size is reached, and undetected sample ratio mismatches in the randomization itself.

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
DESIGN:  randomize -> balances everything, known and unknown, between groups
                |
RUN:     collect data to the PRE-PLANNED sample size (Topic 09), don't stop early
                |
VERIFY:  Sample Ratio Mismatch check -- is the actual split what it should be?
                |
ANALYZE: two-proportion z-test (or the appropriate test for the metric type)
                |
DECIDE:  one look, at the planned n -- not a moving target updated every day

    Silent failure modes this notebook measures directly:
    - PEEKING: checking repeatedly and stopping early inflates false positives
    - SRM: a broken randomizer invalidates EVERYTHING before you even look at the outcome
    - NOVELTY EFFECTS: an overall average can hide a fading (or growing) effect over time
    - MULTIPLE METRICS: testing many outcomes multiplies the chance of a false "winner"
```

---

## 🎯 Why This Topic Matters

- **Randomization was verified to actually balance covariates, not just assumed to** — §1 randomly split 4,000 synthetic users and found all three tracked covariates (prior spend, tenure, mobile usage) balanced between groups (all balance-check p-values > 0.59).
- **The core two-proportion z-test was built from scratch and matched `statsmodels` exactly** — §2 found a genuine 10%→12% simulated lift reaching significance (p=0.0078), with the manual z-statistic matching `statsmodels.stats.proportion.proportions_ztest` to machine precision.
- **Peeking's real cost was measured directly, and it's substantial** — §3 found checking a null (no true effect) experiment just once produced a **4.30%** false-positive rate (correctly near nominal), but checking it **5 times** and stopping at the first p<0.05 inflated that to **13.77%** — nearly **3x** the nominal rate, with **20 checks** pushing it to **23.95%**.
- **A working fix for peeking was built and validated** — §4 applied a Bonferroni-style correction (α/5 per check) and brought the 5-check false-positive rate back down from **13.8%** to **3.65%** — close to the nominal 5% target.
- **Sample Ratio Mismatch detection was demonstrated on both a healthy and a broken randomizer** — §5 found a genuinely healthy 50/50-ish split (19,850 vs 20,150) correctly passing (p=0.134), while a subtly skewed split (19,400 vs 20,600 — just a 3% relative skew) was caught decisively (p=1.97e-09) — a real, practical bug detector that runs *before* trusting any outcome metric.
- **SRM sensitivity at scale was quantified precisely** — §5 found that at n=100,000 total, even a **0.1 percentage point** skew (49,900/50,100) was borderline undetectable (p=0.527) while a **1.0 percentage point** skew was caught with overwhelming certainty (p=2.5e-10) — large experiments can catch very subtle randomization bugs, small ones can't.
- **A fading novelty effect was constructed and shown hidden by the overall average, with a genuinely honest complication** — §6 built a true effect decaying from 0.060 to 0.0034 over 21 days; the overall average (0.0213) looked like a stable, modest effect, and comparing the first-5-day average (0.0488) to the last-5-day average (0.0064) correctly revealed the fade — while a naive single-day comparison (day 0 vs. day 20) happened to coincidentally land on the same noisy value (0.038) both times, a real reminder from this notebook's own results that daily A/B metrics are noisy enough that trends need multi-day averaging, not single-point comparisons.
- **The multiple-metrics problem was measured with a striking number** — §7 found that tracking 10 metrics with **zero true effect on any of them** still produced at least one false "significant" result **39.25%** of the time — a coin-flip's chance of a false winner purely from tracking too many things at once.
- **A complete design-to-decision walkthrough tied every piece together** — §8 designed a test (n=6,687 per variant for an expected 10%→11.5% lift), ran it, passed the SRM check (p=0.979), and reached a genuinely significant, trustworthy result (p=0.0066, 14.7% relative lift) — the full real-world workflow in one place.

---

## 🧮 Mathematical Explanation

### 1. Two-proportion z-test

$$z = \frac{\hat{p}_1-\hat{p}_2}{\sqrt{\hat{p}(1-\hat{p})(1/n_1+1/n_2)}}, \qquad \hat{p}=\frac{x_1+x_2}{n_1+n_2}$$

The standard test for comparing two conversion rates, pooling both groups' data under $H_0$ to estimate the shared proportion used in the standard error.

### 2. Why peeking inflates Type I error

Each look at the data is itself a hypothesis test at level $\alpha$. Across $k$ independent looks, the probability that *at least one* falsely crosses the threshold is approximately $1-(1-\alpha)^k$ (exactly the family-wise error rate logic from Topic 06, applied across time rather than across groups) — §3 measured this empirically rather than relying on the approximation alone, since sequential looks on accumulating data aren't fully independent.

### 3. Sample Ratio Mismatch (SRM) check

$$\chi^2 = \sum_i \frac{(O_i-E_i)^2}{E_i}, \qquad df=1$$

A direct application of Topic 05's goodness-of-fit test to the *assignment counts themselves*, checked before ever looking at the outcome metric.

### 4. Multiple-metrics false-positive rate

$$P(\text{at least one false winner}) = 1-(1-\alpha)^m$$

where $m$ is the number of independently-tracked metrics — §7's 39.25% at $m=10$ is close to the theoretical $1-(0.95)^{10} \approx 40.1\%$, confirming the approximation holds well when metrics are genuinely independent.

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| Two-proportion z-test | $(\hat p_1-\hat p_2)/\sqrt{\hat p(1-\hat p)(1/n_1+1/n_2)}$ | §2 |
| Bonferroni-corrected alpha per check | $\alpha / k$ | §4 |
| SRM chi-square | $\sum(O_i-E_i)^2/E_i$ | §5 |
| Multi-metric false-positive rate | $1-(1-\alpha)^m$ | §7 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Checking significance daily and stopping as soon as p<0.05** — §3's measured 13.8% false-positive rate (5 checks) vs. the nominal 5% is a direct, quantified version of this extremely common real-world mistake.
2. **Trusting any result without an SRM check first** — §5's broken-randomizer example (p=1.97e-09) shows a subtle 3% relative skew is easily detectable and must be resolved before the outcome metric means anything at all.
3. **Judging a feature's long-run value from an experiment's overall average alone** — §6's novelty effect (early average 0.0488, late average 0.0064) would be completely invisible in the single overall-average number (0.0213).
4. **Drawing conclusions from a single day's data point in a time-series breakdown** — this notebook's own §6 result (day 0 and day 20 coincidentally landing on the identical noisy value 0.038) is a live demonstration of why multi-day averaging, not single-point comparison, is the right way to read daily A/B metrics.
5. **Tracking many metrics and reporting "significance" on whichever one crossed 0.05** — §7's 39.25% false-winner rate at 10 metrics (all with zero true effect) is a direct warning against this common "metric shopping" pattern; Topic 11 covers the formal fix.
6. **Assuming the standard SRM threshold of alpha=0.05 is appropriate** — §5 explicitly used a much stricter 0.001 threshold, since SRM checks run on every single experiment and a 5% false-alarm rate would flag healthy experiments constantly at scale.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `statsmodels.stats.proportion.proportions_ztest` | Two-proportion z-test reference |
| `statsmodels.stats.proportion.proportion_effectsize` | Cohen's h for the design step |
| `statsmodels.stats.power.NormalIndPower` | Sample size calculation, reused from Topic 09 |
| `scipy.stats.chisquare` | SRM check, reused from Topic 05 |
| `np.random.default_rng().binomial` | Simulated conversion counts throughout |

---

## 📝 Self-Test Exercises

1. Section 3 found peeking 5 times inflated the false-positive rate to 13.77%, and 20 times to 23.95%. Does the false-positive rate seem to be growing linearly, sub-linearly, or super-linearly with the number of checks? Use the $1-(1-\alpha)^k$ approximation from §2 of the Mathematical Explanation to check your intuition.
2. Section 4's Bonferroni correction required each individual check to clear a much stricter threshold (α/5 = 0.01 instead of 0.05). Explain why this makes it *harder* to detect a real, true effect early in the experiment, even though it correctly controls the false-positive rate.
3. Section 5 found a 0.1 percentage point skew (49,900/50,100 at n=100,000) was NOT reliably detected (p=0.527) while a 1.0 percentage point skew was. Using the chi-square formula from §3 of the Mathematical Explanation, explain why a fixed total sample size limits how small a skew can be reliably detected — what would need to change to catch the smaller skew?
4. Section 6's own example showed day 0 and day 20 coincidentally reporting the identical observed lift (0.038) despite the true effect decaying substantially between those days. Using what you know about binomial sampling variance at n=500/day, estimate roughly how much random day-to-day noise you'd expect around a true lift of ~0.03-0.06, and explain why this makes single-day comparisons unreliable.
5. Section 7 found tracking 10 independent metrics gave a 39.25% chance of at least one false "significant" result. If a team wanted to keep this overall false-winner risk at or below 5% while still tracking 10 metrics, what would the per-metric significance threshold need to be? (This is exactly the calculation Topic 11 formalizes.)

---

## 📓 Notebook

30 executed code cells: a randomization balance check across 3 covariates on 4,000 synthetic users, a two-proportion z-test built from scratch and matched exactly to `statsmodels`, a peeking simulation showing false-positive rate climbing from 4.3% (1 check) to 13.77% (5 checks) to 23.95% (20 checks), a Bonferroni-style correction bringing the 5-check rate back down to 3.65%, a Sample Ratio Mismatch check correctly passing a healthy randomizer and catching a broken one (p=1.97e-09), an SRM sensitivity sweep across four skew magnitudes at n=100,000, a 21-day fading novelty effect simulation with an honestly-reported single-day-noise complication resolved via early/late multi-day averaging, a multiple-metrics false-positive simulation (39.25% at 10 metrics), and a complete 5-step A/B test walkthrough from sample-size design through a final launch decision:

➡️ **[10_ab_testing_experimental_design.ipynb](10_ab_testing_experimental_design.ipynb)**

---

## 📚 Further Reading

- [statsmodels.stats.proportion documentation](https://www.statsmodels.org/stable/stats.html#proportion)
- [Kohavi, Tang & Xu, *Trustworthy Online Controlled Experiments*](https://www.cambridge.org/core/books/trustworthy-online-controlled-experiments/D97B26382EB0EB2DC2019A7A7B518F59)
- [Fabijan et al. (2019): Diagnosing Sample Ratio Mismatch in Online Controlled Experiments](https://exp-platform.com/Documents/2019_KDDFabijanSRM.pdf)
- [Johari et al. (2017): Peeking at A/B Tests — Why It Matters and What to Do About It](https://dl.acm.org/doi/10.1145/3097983.3097992)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.10_AB_Testing_Experimental_Design&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
