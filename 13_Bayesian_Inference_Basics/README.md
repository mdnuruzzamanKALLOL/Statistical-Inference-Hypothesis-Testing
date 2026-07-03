# 🎲 Bayesian Inference Basics

> Status: ✅ Complete — [Open the notebook →](13_bayesian_inference_basics.ipynb)

Topic 13 of the Statistical Inference & Hypothesis Testing repo. Every topic so far worked within the **frequentist** framework: parameters are fixed unknowns, and probability describes long-run frequency across repeated sampling (Topic 03's confidence interval coverage simulation was the clearest example). This notebook introduces the **Bayesian** alternative: parameters themselves have probability distributions representing belief, updated by data via Bayes' theorem. It builds the Beta-Binomial and Normal-Normal conjugate models from scratch, directly contrasts a credible interval's interpretation against Topic 03's confidence interval, and closes with a practical Bayesian A/B test revisiting Topic 10.

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
FREQUENTIST (Topics 01-12): theta is a fixed unknown constant
   probability describes long-run frequency across repeated sampling
   "if we repeated this study many times, 95% of intervals would contain theta"

BAYESIAN (this topic): theta is a random variable with a belief distribution
   PRIOR belief + DATA (likelihood) -> POSTERIOR belief, via Bayes' theorem
   "given this prior and this data, there is a 95% probability theta is in here"

Same data can be analyzed BOTH ways -- they are different, both valid,
frameworks for reasoning under uncertainty, not competing "right answers."
```

---

## 🎯 Why This Topic Matters

- **The Beta-Binomial conjugate model was built from scratch and closely tracked the frequentist MLE under a uniform prior** — §1 found a posterior mean (0.3636) differing from the raw MLE (0.3500) by just **0.0136** — with no strong prior belief, the two frameworks nearly agree numerically.
- **Prior strength was shown genuinely changing the posterior, and precisely how much** — §2 found a weak prior's posterior (0.3636) staying close to the data, while a strong, confidently-wrong prior (Beta(40,10), believing theta≈0.8) pulled its posterior all the way to **0.6714** — a real, measured difference, not a hand-wave.
- **That prior influence was shown fading as data accumulates, with an exact before/after number** — §2 found the gap between weak-prior and strong-prior posteriors shrinking from **0.3078** (n=20) to just **0.0403** (n=500) — the same underlying strong prior, overwhelmed by enough evidence.
- **Sequential updating was proven mathematically identical to one-shot updating** — §3 confirmed that updating on 3 separate data batches produced the *exact same* posterior, Beta(8,14), as updating on all 20 observations at once.
- **The credible interval's own coverage was checked under repeated sampling — a genuine cross-framework validation** — §4 found a weak-prior 95% credible interval's empirical coverage across 5,000 simulated datasets was **96.92%**, close to nominal — confirming that under an uninformative prior, the credible interval also behaves sensibly by frequentist standards, even though matching frequentist coverage isn't its defining guarantee.
- **The credible interval and Topic 03's Wilson confidence interval were shown numerically nearly identical, while remaining conceptually distinct** — §4 found [0.1811, 0.5697] (credible) versus [0.1812, 0.5671] (Wilson) — the two frameworks converging numerically under an uninformative prior, even though only the Bayesian interval licenses the direct "95% probability theta is in here" statement.
- **The Normal-Normal model was shown producing a genuine precision-weighted compromise** — §5 found a posterior mean of **107.75** sitting between a prior mean of 100 and a sample mean of 114.72, weighted by each source's relative precision.
- **Both the posterior mean and posterior uncertainty were shown converging with more data, starting from a deliberately biased prior** — §6 found the mean-gap shrinking from **7.264** (n=5) to **0.064** (n=1000), and the posterior standard deviation shrinking from **4.01** to **0.47** over the same range — more data narrows belief exactly the way it narrows a frequentist interval.
- **A practical Bayesian A/B test on Topic 10's own data produced a direct, intuitive answer** — §7 found **P(treatment > control | data) = 0.9961** and an expected uplift of **1.90 percentage points** — a different, arguably more directly useful way to report the same underlying result Topic 10 found significant at p=0.0078.
- **A mismatched strong prior was shown still distorting a real-data estimate even with hundreds of observations** — §8 found a deliberately wrong prior (believing the true rate was near 30%) pulling the real `tips` high-tipper estimate from 16.26% (weak prior, close to the 15.98% MLE) up to **20.06%** — even n=244 real observations couldn't fully overcome a confidently wrong prior belief.

---

## 🧮 Mathematical Explanation

### 1. Bayes' theorem

$$P(\theta \mid \text{data}) = \frac{P(\text{data} \mid \theta)\,P(\theta)}{P(\text{data})} \propto P(\text{data}\mid\theta)\,P(\theta)$$

Posterior belief is proportional to likelihood times prior — the entire Bayesian machinery in one line.

### 2. Beta-Binomial conjugacy

$$\theta \sim \text{Beta}(\alpha,\beta) \ \Rightarrow \ \theta \mid x, n \sim \text{Beta}(\alpha+x, \ \beta+n-x)$$

$\alpha+\beta$ acts as the prior's "effective sample size" — §2's strong prior (Beta(40,10), effective n=50) needed comparably large real data to be overwhelmed, exactly as §2's fading-influence result showed.

### 3. Normal-Normal conjugacy

$$\mu_{\text{post}} = \frac{\mu_0/\sigma_0^2 + n\bar{x}/\sigma^2}{1/\sigma_0^2 + n/\sigma^2}, \qquad
\sigma^2_{\text{post}} = \left(\frac{1}{\sigma_0^2}+\frac{n}{\sigma^2}\right)^{-1}$$

A precision-weighted average — each source (prior, data) contributes to the posterior mean in proportion to how precise (low-variance) it is.

### 4. Credible interval vs. confidence interval

A credible interval is a direct quantile of the posterior distribution: $P(\theta \in [\theta_{lo}, \theta_{hi}] \mid \text{data}) = 0.95$ — a genuine probability statement about $\theta$, licensed by treating $\theta$ as a random variable. A frequentist confidence interval's 95% instead describes the *construction procedure* across repeated sampling (Topic 03) — numerically the two can nearly coincide under an uninformative prior, but they are answering different questions.

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| Beta-Binomial posterior | $\text{Beta}(\alpha+x,\ \beta+n-x)$ | §1-§4 |
| Normal-Normal posterior mean | precision-weighted average | §5-§6 |
| Credible interval | posterior quantiles | §4 |
| Bayesian A/B P(better) | $P(\theta_T > \theta_C \mid \text{data})$, via Monte Carlo | §7 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Assuming a "prior" is a free, harmless modeling choice** — §8's real-data example shows a mismatched strong prior still distorting an estimate meaningfully (16.26% → 20.06%) even with 244 real observations; prior choice is a genuine, consequential decision.
2. **Confusing the credible interval's interpretation with the confidence interval's** — even when the two numerically nearly coincide (§4), only the credible interval licenses a direct probability statement about the parameter; conflating them misstates what each framework actually guarantees.
3. **Assuming Bayesian and frequentist point estimates always differ substantially** — §1 showed them agreeing to within 0.0136 under a uniform prior; the frameworks diverge most when the prior carries real information (or misinformation), not by default.
4. **Forgetting that sequential Bayesian updating is order-invariant** — §3 confirms this exactly, a useful property for streaming or incrementally-arriving data that a purely frequentist re-analysis-from-scratch approach doesn't offer as naturally.
5. **Believing a "non-informative" or "uniform" prior is truly assumption-free** — a Beta(1,1) still encodes an assumption (every value of theta equally likely a priori); it just happens to be a weak one, as §1-§2 demonstrate by contrast with the strong prior.
6. **Reporting a Bayesian A/B test result without an expected-uplift or loss estimate** — §7's P(treatment > control) = 0.9961 is more actionable paired with the expected uplift (1.90pp) than alone; a high probability of "better" doesn't say by how much.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `scipy.stats.beta.pdf` / `.ppf` | Beta distribution density and quantiles (credible intervals) |
| `scipy.stats.norm.ppf` | Normal quantiles used in the Wilson CI comparison |
| `np.random.default_rng().beta` | Monte Carlo posterior sampling for the Bayesian A/B test |
| `np.random.default_rng().binomial` | Simulated data for the credible-interval coverage check |

---

## 📝 Self-Test Exercises

1. Section 2 found the strong prior's posterior (0.6714) much closer to its own belief (~0.8) than to the data (0.35) at n=20, but the gap between weak- and strong-prior posteriors shrank to just 0.0403 at n=500. Using the Beta-Binomial update formula from §2 of the Mathematical Explanation, explain why a prior's "effective sample size" ($\alpha+\beta$) determines how much real data is needed to overwhelm it.
2. Section 4 found the credible interval's own frequentist coverage (96.92%) close to nominal under a weak prior. Would you expect this same close-to-95% coverage property to hold if the *strong, mismatched* prior from §2 were used instead? Explain your reasoning.
3. Section 6 showed both the posterior mean gap and the posterior standard deviation shrinking as n grew from 5 to 1000. Are these two convergences (mean-gap shrinking, and uncertainty/std shrinking) the same phenomenon, or two logically separate consequences of accumulating data? Justify your answer.
4. Section 7 reported both P(treatment > control | data) = 0.9961 and an expected uplift of 1.90 percentage points. Construct a hypothetical scenario where a very high P(better) (say, 0.99) is reported alongside a very small expected uplift (say, 0.05 percentage points) — what would that combination tell a decision-maker that P(better) alone would not?
5. Section 8 found the real-data weak-prior posterior mean (0.1626) close to the MLE (0.1598), while the mismatched-prior posterior (0.2006) stayed noticeably distorted. Using the Beta-Binomial formula, estimate roughly how many additional real observations (beyond the 244 already used) would be needed to pull the mismatched-prior posterior mean within 0.005 of the true MLE.

---

## 📓 Notebook

30 executed code cells: the Beta-Binomial conjugate model built from scratch and validated against the frequentist MLE, a direct weak-vs-strong-prior comparison showing the prior's fading influence as data accumulates, a proof that sequential and one-shot Bayesian updating produce identical posteriors, a credible interval built and directly contrasted against Topic 03's Wilson confidence interval (both numerically and in a 5,000-simulation coverage check), the Normal-Normal conjugate model showing a precision-weighted posterior compromise, a convergence demonstration for both the posterior mean and posterior uncertainty as sample size grows from 5 to 1000, a full Bayesian A/B test revisiting Topic 10's exact data with a Monte Carlo-estimated P(treatment > control), and a real-data application to the `tips` dataset showing a mismatched strong prior still distorting an estimate at n=244:

➡️ **[13_bayesian_inference_basics.ipynb](13_bayesian_inference_basics.ipynb)**

---

## 📚 Further Reading

- [scipy.stats.beta documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.beta.html)
- [Gelman et al., *Bayesian Data Analysis*](http://www.stat.columbia.edu/~gelman/book/)
- [Kruschke, *Doing Bayesian Data Analysis*](https://sites.google.com/site/doingbayesiandataanalysis/)
- [VanderPlas (2014): Frequentism and Bayesianism — A Practical Introduction](https://jakevdp.github.io/blog/2014/03/11/frequentism-and-bayesianism-a-practical-intro/)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.13_Bayesian_Inference_Basics&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
