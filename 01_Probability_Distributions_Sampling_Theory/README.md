# 🎲 Probability Distributions & Sampling Theory

> Status: ✅ Complete — [Open the notebook →](01_probability_distributions_sampling_theory.ipynb)

Topic 01 of the Statistical Inference & Hypothesis Testing repo. Every hypothesis test built in later topics reduces to a question about a probability distribution — this notebook builds the three core distributions (Binomial, Normal, Poisson) from scratch, validates each one exactly against `scipy.stats`, and then makes the conceptual jump every test in this repo depends on: from *a* distribution to a *sampling distribution*, the object every p-value is actually computed from.

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
A single random variable's distribution        ->  Binomial, Normal, Poisson (Sections 1-4)
                v
One sample drawn from a population              ->  a single noisy point estimate (Section 6)
                v
MANY samples, each summarized to one number      ->  the SAMPLING DISTRIBUTION (Section 7)
                v
Standard Error = spread of that sampling         ->  SE = sigma / sqrt(n)  (Section 8)
distribution
                v
Every t-test, chi-square test, ANOVA in later
topics is a statement about WHERE a sample
statistic falls in its own sampling distribution
```

---

## 🎯 Why This Topic Matters

- **All three distributions were validated exactly, not approximately** — §1, §2, §3 show a from-scratch Binomial PMF, Normal PDF, and Poisson PMF each matching `scipy.stats` to within machine precision (max absolute differences of `3.33e-16`, `3.47e-18`, and `8.33e-17` respectively).
- **The Binomial → Poisson convergence was demonstrated directly, not just cited** — §3 holds `n·p = λ = 4.0` fixed while shrinking `p` and growing `n` from 20 to 100,000, and the max deviation from the true Poisson(4) PMF shrinks monotonically from `2.28e-02` down to `3.91e-06`.
- **A numerically-integrated CDF was cross-checked against the analytical one** — §4 builds the Normal CDF via trapezoidal integration (no closed form used) and matches `scipy.stats.norm.cdf` to within `8.5e-09` at five test points, including exactly `0.500000` at the mean.
- **A single sample's instability was shown concretely before the fix (more samples) was introduced** — §6 draws one sample of size 30 from a population with true mean `2.0006` and gets an estimate of `2.5182` — a 25.9% relative error — then draws 5 more single samples that swing from `1.3698` to `2.3959`, motivating why Section 7 exists at all.
- **The sampling distribution's shrinking skewness was measured, not assumed** — §7 simulates the sampling distribution of the mean at n=5, 30, 100, 500 from a deliberately skewed Exponential population (skewness 1.99); the sampling distribution's own skewness shrinks monotonically from `0.899` (n=5) down to `0.054` (n=500) — a direct, numeric preview of the Central Limit Theorem, formalized in Topic 02.
- **SE = σ/√n was checked against a real Monte Carlo simulation, not just stated** — §8 compares the closed-form standard error to the empirical standard deviation of 5,000 simulated sample means at each of 4 sample sizes; the largest relative error across all four was only `2.27%`.
- **Stratified sampling's variance reduction was measured on a constructed two-strata population** — §9 built two Normal strata (means 50 and 90) and found stratified sampling cut estimator variance by **79.1%** compared to simple random sampling at the same sample size (n=40), while both remained essentially unbiased (bias `0.084` vs `0.036`) — proof the benefit is in variance, not bias.
- **The whole pipeline was re-run on a real, messy dataset** — §10 repeats the sampling-distribution logic on seaborn's real `tips` dataset (skewness 1.46, genuinely non-Normal), and the theoretical SE (`0.3088`) matched the empirical SE from 4,000 resamples (`0.3058`) to within `0.97%`.

---

## 🧮 Mathematical Explanation

### 1. Binomial distribution

$$P(X=k) = \binom{n}{k} p^k (1-p)^{n-k}, \quad E[X] = np, \quad \text{Var}(X) = np(1-p)$$

Counts successes in $n$ independent Bernoulli($p$) trials. §1's manual implementation directly evaluates this formula per $k$ using `math.comb`.

### 2. Normal (Gaussian) distribution

$$f(x) = \frac{1}{\sigma\sqrt{2\pi}} \exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)$$

Unlike the Binomial's PMF, $f(x)$ is a *density*, not a probability — it can exceed 1, and probabilities only come from area under the curve (§4). §2's numerical-integration check confirmed $\int_{-\infty}^{\infty} f(x)\,dx = 1.0000000000$ over $\pm 10\sigma$.

### 3. Poisson distribution and the Binomial limit

$$P(X=k) = \frac{\lambda^k e^{-\lambda}}{k!}$$

Classical limiting result: $\text{Binomial}(n, p) \to \text{Poisson}(\lambda)$ as $n \to \infty$, $p \to 0$, with $np = \lambda$ held fixed. §3 confirms this is not just a symbolic limit but a genuine numerical convergence — the max PMF deviation from true Poisson(4) shrank by four orders of magnitude across the tested range of $n$.

### 4. The CDF as an integral

$$F(x) = P(X \le x) = \int_{-\infty}^{x} f(t)\,dt$$

§4 computes this via the trapezoid rule (`np.trapezoid`) rather than any closed form, confirming a purely numerical approach converges to the same answer as `scipy`'s analytical implementation.

### 5. The sampling distribution and the standard error

If $X_1, \dots, X_n$ are drawn i.i.d. from a population with mean $\mu$ and standard deviation $\sigma$, the sample mean $\bar{X} = \frac{1}{n}\sum X_i$ is itself a random variable with

$$E[\bar{X}] = \mu, \qquad SE(\bar{X}) = \frac{\sigma}{\sqrt{n}}$$

This is the single most important formula in the whole notebook: it says the *spread* of the sampling distribution shrinks at rate $1/\sqrt{n}$, not $1/n$ — doubling precision requires *quadrupling* the sample size, a fact §8's log-log-shaped curve makes visually obvious.

### 6. Stratified sampling variance

For a population split into $L$ strata with within-stratum variances $\sigma_l^2$, proportional-allocation stratified sampling gives

$$\text{Var}(\bar{X}_{\text{strat}}) = \sum_{l=1}^{L} w_l^2 \frac{\sigma_l^2}{n_l} \le \text{Var}(\bar{X}_{\text{SRS}})$$

whenever strata differ in mean (the between-strata variance that simple random sampling has to "absorb" as noise gets removed entirely by stratifying) — §9 measured this reduction directly at **79.1%** for two well-separated Normal strata.

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| Binomial PMF | $\binom{n}{k}p^k(1-p)^{n-k}$ | §1 |
| Normal PDF | $\frac{1}{\sigma\sqrt{2\pi}}e^{-(x-\mu)^2/2\sigma^2}$ | §2 |
| Poisson PMF | $\frac{\lambda^k e^{-\lambda}}{k!}$ | §3 |
| CDF | $\int_{-\infty}^{x} f(t)\,dt$ | §4 |
| Standard Error | $\sigma / \sqrt{n}$ | §8 |
| Stratified variance | $\sum_l w_l^2 \sigma_l^2 / n_l$ | §9 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Treating a single sample's statistic as the population truth** — §6 showed a single n=30 sample landing 25.9% off the true mean; always think in terms of the sampling *distribution*, not one draw.
2. **Assuming a Normal-shaped sampling distribution requires a Normal population** — §7 and §10 both used deliberately skewed populations (Exponential, and the real right-skewed `tips` data) and still watched the sampling distribution's skewness shrink toward 0 as $n$ grew — this is a preview, not a full proof, of the CLT (Topic 02 formalizes exactly when and why this holds).
3. **Confusing "unbiased" with "low variance"** — §9's bias check confirmed both simple random and stratified sampling were nearly unbiased (`0.084` vs `0.036`); the real difference stratification bought was a 79.1% variance reduction, a completely separate property from bias.
4. **Using `np.trapz`** — this function was removed in NumPy 2.x; use `np.trapezoid` instead (a real `AttributeError` hit during this notebook's first build, now fixed at §4).
5. **Non-ASCII characters (em-dashes, etc.) inside `print()` statements on Windows** — a real encoding bug was hit while building this notebook: Windows' kernel-subprocess stdout capture mangled `—` into `�` inside stored cell outputs (confirmed via `repr()` on the raw notebook JSON, not just a terminal display issue). Fixed by keeping all code-cell `print()` output plain ASCII; non-ASCII characters remain safe in markdown cells since those aren't executed through a subprocess.
6. **Forgetting the population std vs. sample std distinction** — the SE formula uses the *population* $\sigma$; §5 and §10 compute the "true" std with `ddof=0` since the full synthetic population / full `tips` dataset stood in for a known population, whereas a real sample estimate of $\sigma$ would need `ddof=1`.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `scipy.stats.binom.pmf` / `.norm.pdf` / `.poisson.pmf` | Reference implementations validated against |
| `scipy.stats.norm.cdf` | Analytical Normal CDF, validated against manual trapezoidal integration |
| `scipy.integrate.quad` | Numerical integration used to confirm the Normal PDF integrates to 1 |
| `scipy.stats.skew` | Skewness of a sample or simulated distribution |
| `np.trapezoid` | Trapezoidal-rule numerical integration (NumPy 2.x name; replaces removed `np.trapz`) |
| `np.random.default_rng` | Seeded random number generator used for all Monte Carlo simulations |
| `seaborn.load_dataset("tips")` | Real dataset used for the Section 10 real-data validation |

---

## 📝 Self-Test Exercises

1. Using §1's formula, compute by hand $P(X=3)$ for a Binomial(n=5, p=0.4) and compare it to what `scipy.stats.binom.pmf(3, 5, 0.4)` would return.
2. Section 3 held $n \cdot p = 4$ fixed while increasing $n$ from 20 to 100,000. Explain in your own words why the Binomial's max deviation from the Poisson shrinks as $n$ grows under this constraint, using the shape of the Binomial variance formula $np(1-p)$ as a hint.
3. Section 7 found the sampling distribution's skewness at n=500 (`0.054`) was much closer to 0 than at n=5 (`0.899`), even though the underlying population's skewness (`1.99`) never changed. What does this tell you about the difference between "the distribution of one observation" and "the distribution of an average of many observations"?
4. Using §5's SE formula, if a sample size of $n=100$ gives $SE=2.0$, what sample size would be needed to cut the SE in half? (This is a deliberately non-linear question — don't guess $n=50$.)
5. Section 9 found stratified sampling reduced variance by 79.1% when the two strata had means 50 and 90 (a 40-point gap). Predict, without running the code, whether that variance reduction would be larger or smaller if the two strata instead had means 65 and 75 (a 10-point gap) — then justify your answer using the variance formula in §6 of the Mathematical Explanation.

---

## 📓 Notebook

31 executed code cells: a from-scratch Binomial PMF validated against `scipy.stats.binom` to within `3.33e-16`, a from-scratch Normal PDF validated against `scipy.stats.norm` to within `3.47e-18` plus a numerical-integration check that it integrates to exactly 1.0, a from-scratch Poisson PMF validated against `scipy.stats.poisson` to within `8.33e-17`, a direct numerical demonstration of Binomial→Poisson convergence across four orders of magnitude in $n$, a trapezoidal-integration Normal CDF matching the analytical one to within `8.5e-09`, a 200,000-point synthetic Exponential population with known ground-truth mean/variance, a single-sample instability demonstration (25.9% error on one draw), a full Monte Carlo simulation of the sampling distribution of the mean at four sample sizes with skewness measured at each, a standard-error formula check against Monte Carlo simulation (max relative error 2.27%), a two-strata simple-random-vs-stratified-sampling variance comparison (79.1% variance reduction), and a full repeat of the sampling-distribution/standard-error analysis on the real `tips` dataset:

➡️ **[01_probability_distributions_sampling_theory.ipynb](01_probability_distributions_sampling_theory.ipynb)**

---

## 📚 Further Reading

- [scipy.stats: Statistical functions](https://docs.scipy.org/doc/scipy/reference/stats.html)
- [Wasserman, *All of Statistics*, Ch. 2-3: Random Variables & Expectation](https://link.springer.com/book/10.1007/978-0-387-21736-9)
- [Casella & Berger, *Statistical Inference*, Ch. 5: Properties of a Random Sample](https://www.taylorfrancis.com/books/mono/10.1201/9781003456978)
- [seaborn `tips` dataset documentation](https://github.com/mwaskom/seaborn-data)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.01_Probability_Distributions_Sampling_Theory&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
