# 🔔 Central Limit Theorem & Point Estimation

> Status: ✅ Complete — [Open the notebook →](02_central_limit_theorem_point_estimation.ipynb)

Topic 02 of the Statistical Inference & Hypothesis Testing repo. Topic 01 previewed the Central Limit Theorem informally on one population shape. This notebook formalizes it — testing convergence across four genuinely different population shapes, quantifying *how* Normal each sampling distribution becomes with a real Kolmogorov-Smirnov goodness-of-fit test, and showing convergence speed depends on skewness. The second half turns to point estimation: bias, consistency, and efficiency, each measured by Monte Carlo simulation rather than assumed from the formula alone.

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
Population, ANY shape (flat, skewed, two-point, bimodal)
                v
Draw a sample of size n, compute its mean, standardize:  Z = (X_bar - mu) / (sigma/sqrt(n))
                v
Repeat thousands of times  ->  distribution of Z approaches N(0,1) as n grows
                v
HOW FAST this happens depends on the population's skewness -- not the same for every shape
                v
Meanwhile: an ESTIMATOR (mean, variance, median...) is judged by three separate properties:
   BIAS         = E[estimator] - true parameter        (systematic error)
   CONSISTENCY  = does the estimator -> truth as n grows?
   EFFICIENCY   = how small is the estimator's variance / MSE?
```

---

## 🎯 Why This Topic Matters

- **CLT was tested on four population shapes that share nothing in common** — §2 built Uniform(0,1) (flat, skewness ≈ 0.001), Exponential(λ=1) (skewness 2.01), Bernoulli(p=0.1) (two-point discrete, skewness 2.67), and a bimodal Normal mixture (two separated humps, skewness ≈ 0) — none of them look remotely Normal, yet every one of their standardized sampling distributions moved toward N(0,1).
- **Convergence speed was directly measured, not assumed uniform** — §3 swept sample sizes 5→300 across all four shapes and found a **0.833 correlation** between a population's absolute skewness and how far its sampling distribution's KS statistic was from Normal at n=30 — confirming skewed populations need more data, not just "more data helps in general."
- **A concrete convergence threshold was found per shape** — §3 tracked the smallest tested n where each population's KS statistic dropped below 0.02: Uniform got there at **n=5**, the bimodal mixture at **n=30**, while Exponential and Bernoulli(p=0.1) *never* reached that threshold even at n=300 in this run — a genuinely useful, non-obvious result about which population shapes are "CLT-friendly" at practical sample sizes.
- **The ddof=0 vs ddof=1 variance bias was measured, and it matched theory exactly** — §5 ran 20,000 repeated n=5 samples from a Normal(50, 10) population and found the ddof=0 estimator averaged **79.11** against a true variance of 100 — matching the theoretical $(n-1)/n \times \sigma^2 = 80.0$ prediction (difference: 0.89, within Monte Carlo noise) — while ddof=1 averaged **98.89**, essentially unbiased.
- **The bias-shrinks-with-n table quantified exactly how much this matters at different sample sizes** — §5 showed the ddof=0 bias shrinking from **-32.4** at n=3 to **-1.1** at n=100, closely tracking the theoretical $-\sigma^2/n$ prediction at every step.
- **Consistency was watched happening in real time on one growing sample** — §6 tracked a single sample growing from n=10 to n=5000: the mean's error shrank from **1.48 to 0.05**, and the variance estimate's error shrank from **45.3 to 1.18** — direct visual and numeric proof of the Law of Large Numbers, not a citation of it.
- **Mean vs. Median efficiency genuinely flipped depending on the population, exactly as theory predicts** — §7 found the Mean had lower MSE (0.032 vs 0.051) under a clean Normal population, but under a 5%-outlier-contaminated population (kurtosis 23,224!) the Mean's MSE exploded to **96.49** while the Median barely moved to **0.055** — a real, measured efficiency-vs-robustness trade-off, not a textbook assertion.
- **Method of Moments was validated against MLE exactly** — §8 estimated an Exponential rate parameter two independent ways (MoM: $1/\bar{X}$, and `scipy.stats.expon.fit`'s MLE) and both converged to **λ̂ = 1.9296**, a difference of `0.0000` — confirming the known special property that MoM and MLE coincide for this one-parameter exponential family.
- **Everything was re-confirmed on real, messy data** — §9 repeated the CLT check on 53,940 real diamond prices (skewness 1.62, nothing like Normal) and got KS=0.0158 (p=0.270, consistent with Normality) at n=40, plus an 11.1% ddof=0-vs-ddof=1 variance gap on a real n=10 sample — matching the synthetic-data story from §5 on genuine data.

---

## 🧮 Mathematical Explanation

### 1. The Central Limit Theorem

For i.i.d. $X_1, \dots, X_n$ with finite mean $\mu$ and variance $\sigma^2$:

$$Z_n = \frac{\bar{X}_n - \mu}{\sigma/\sqrt{n}} \xrightarrow{d} N(0, 1) \quad \text{as } n \to \infty$$

The theorem places **no restriction on the population's shape** — only that it has finite mean and variance. §1-§3 test this directly rather than taking it on faith, across shapes ranging from perfectly flat to a two-point discrete distribution.

### 2. Convergence rate and skewness

The Berry-Esseen theorem gives a formal bound on how far $Z_n$'s CDF can be from $\Phi(z)$ (the standard Normal CDF) as a function of the population's third absolute moment — informally, **more skewed populations converge more slowly**. §3's correlation of 0.833 between skewness and KS-distance-from-Normal at fixed n is an empirical echo of this formal bound, not a derivation of it.

### 3. Bias of an estimator

$$\text{Bias}(\hat{\theta}) = E[\hat{\theta}] - \theta$$

For the two common variance estimators:

$$E\left[\frac{1}{n}\sum(X_i - \bar{X})^2\right] = \frac{n-1}{n}\sigma^2 \quad \text{(biased, ddof=0)} \qquad
E\left[\frac{1}{n-1}\sum(X_i - \bar{X})^2\right] = \sigma^2 \quad \text{(unbiased, ddof=1)}$$

The ddof=0 formula divides by $n$ but measures distance to $\bar{X}$ (which itself was fit *to* that same sample), losing exactly one degree of freedom — Bessel's correction ($n-1$) restores unbiasedness. §5 confirmed the $(n-1)/n$ factor numerically at five different sample sizes.

### 4. Consistency

An estimator $\hat{\theta}_n$ is consistent if $\hat{\theta}_n \xrightarrow{p} \theta$ as $n \to \infty$ — it may be biased at small $n$ (like ddof=0 above) yet still be consistent, since the bias itself shrinks to 0. §6 demonstrated this is not just a limit-point claim but an observable trend on a single growing sample.

### 5. Mean Squared Error and efficiency

$$MSE(\hat{\theta}) = \text{Bias}(\hat{\theta})^2 + \text{Var}(\hat{\theta})$$

An estimator is more **efficient** than another if it has lower MSE at the same $n$. §7 showed efficiency is not an absolute property of an estimator — it depends on the population: the Mean is more efficient under Normal data (lower variance, same near-zero bias) but far less efficient under heavy-tailed contaminated data, where its variance explodes while the Median's stays stable.

### 6. Method of Moments

Set $k$ theoretical moments equal to $k$ sample moments and solve for the parameters. For Exponential($\lambda$), $E[X] = 1/\lambda$, giving $\hat{\lambda}_{\text{MoM}} = 1/\bar{X}$ — §8 confirmed this equals the MLE exactly for this distribution family, though that equivalence does **not** hold in general (explored fully in Topic 16).

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| Standardized sample mean | $(\bar{X}-\mu)/(\sigma/\sqrt{n})$ | §1-§3 |
| Bias | $E[\hat{\theta}] - \theta$ | §5 |
| ddof=0 variance expectation | $\frac{n-1}{n}\sigma^2$ | §5 |
| MSE | $\text{Bias}^2 + \text{Variance}$ | §7 |
| Method of Moments (Exponential) | $\hat{\lambda} = 1/\bar{X}$ | §8 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Assuming CLT convergence speed is the same for every population** — §3's 0.833 skewness-vs-KS correlation and the never-reached n=300 threshold for Exponential/Bernoulli(p=0.1) both directly contradict the naive "n=30 is always enough" rule of thumb.
2. **Treating bias and consistency as the same property** — the ddof=0 variance estimator (§5) is biased at every finite $n$ yet still consistent (§6's error-shrinking trend); an estimator can be one without being the other.
3. **Assuming the mean is always the "best" estimator of center** — §7 showed the Mean's MSE exploding to 96.49 under a contaminated population where the Median stayed at 0.055; "best" depends entirely on the population's tail behavior.
4. **Assuming Method of Moments and MLE always agree** — §8's exact match is a special property of the Exponential family, not a general rule; the two methods diverge for most other distributions.
5. **Ignoring the ddof=0/ddof=1 gap "because n is large"** — §9 showed an 11.1% relative gap on a real n=10 sample of diamond prices; the gap only becomes negligible once $n$ is in the hundreds, not for typical small-sample real-world work.
6. **Non-ASCII characters inside executed `print()` statements** — same Windows kernel-encoding bug documented in Topic 01; this notebook's code cells use plain ASCII (`--` instead of `—`) throughout for exactly that reason.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `scipy.stats.kstest` | Kolmogorov-Smirnov goodness-of-fit test used to quantify closeness to N(0,1) |
| `scipy.stats.skew` / `.kurtosis` | Third/fourth moment diagnostics for each population |
| `scipy.stats.expon.fit` | MLE fit used to cross-check the Method of Moments estimate |
| `np.random.default_rng` | Seeded RNG used for every Monte Carlo simulation |
| `pandas.DataFrame.pivot` | Reshaped the n-by-population KS-statistic sweep into a comparison table |
| `seaborn.load_dataset("diamonds")` | Real, heavily right-skewed dataset used for the Section 9 validation |

---

## 📝 Self-Test Exercises

1. Section 3 found the Bernoulli(p=0.1) population needed the most samples to look Normal (KS statistic still 0.055 at n=300). Using the skewness values from §2 (Bernoulli's skewness was 2.67, the highest of the four), explain why a highly imbalanced discrete distribution is a particularly hard case for the CLT at small-to-moderate $n$.
2. Using §3's theoretical bias formula $-\sigma^2/n$, predict the ddof=0 estimator's bias at $n=20$ for a population with $\sigma^2=100$, then compare your prediction to the pattern in §5's bias-vs-n table.
3. Section 7 found the Mean's MSE was 96.49 under the contaminated population versus 0.032 under the clean one — a roughly 3000x increase — while the Median only moved from 0.051 to 0.055. What property of the sample mean's formula (as opposed to the median's) makes it so sensitive to a small fraction of extreme values?
4. Section 8 showed MoM and MLE coinciding exactly for the Exponential distribution. Without running any code, guess whether you'd expect the same exact agreement for a Normal distribution's mean and variance parameters, and justify your guess using the fact that the Normal's MLE for $\mu$ is also $\bar{X}$.
5. Section 9's real diamond-price sample showed an 11.1% relative gap between ddof=0 and ddof=1 variance at $n=10$. Using the $(n-1)/n$ formula from §5's Mathematical Explanation, calculate what sample size would be needed to shrink that gap to under 1%.

---

## 📓 Notebook

30 executed code cells: a first CLT demonstration on a Uniform(0,1) population with a Kolmogorov-Smirnov test confirming near-Normality by n=100, four population shapes (Uniform, Exponential, Bernoulli, Bimodal) tested for CLT convergence with a full n=5→300 KS-statistic sweep, a direct correlation measurement (0.833) between population skewness and convergence speed, a per-shape "smallest n to reach KS<0.02" comparison, a 20,000-repeat Monte Carlo demonstration of the ddof=0 vs ddof=1 variance bias matching theory to within 0.89, a bias-vs-sample-size sweep from n=3 to n=100, a single growing sample (n=10→5000) tracked for mean/variance consistency, a Mean-vs-Median efficiency comparison under clean and 5%-outlier-contaminated populations, a Method of Moments Exponential rate estimator validated exactly against scipy's MLE fit, and a full CLT + bias re-validation on the real, 53,940-row `diamonds` price dataset:

➡️ **[02_central_limit_theorem_point_estimation.ipynb](02_central_limit_theorem_point_estimation.ipynb)**

---

## 📚 Further Reading

- [scipy.stats.kstest documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kstest.html)
- [Wasserman, *All of Statistics*, Ch. 5: The Central Limit Theorem](https://link.springer.com/book/10.1007/978-0-387-21736-9)
- [Casella & Berger, *Statistical Inference*, Ch. 7: Point Estimation](https://www.taylorfrancis.com/books/mono/10.1201/9781003456978)
- [Berry-Esseen theorem overview](https://en.wikipedia.org/wiki/Berry%E2%80%93Esseen_theorem)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.02_Central_Limit_Theorem_Point_Estimation&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
