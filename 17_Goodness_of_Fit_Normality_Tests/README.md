# 📐 Goodness-of-Fit & Normality Tests

> Status: ✅ Complete — [Open the notebook →](17_goodness_of_fit_normality_tests.ipynb)

Topic 17 of the Statistical Inference & Hypothesis Testing repo. The t-test (Topic 04), ANOVA (Topic 06), Pearson correlation (Topic 08), and the Normal MLE (Topic 16) all lean on a Normality assumption that has never, until now, been formally tested. This notebook builds a Q-Q plot from scratch, then Shapiro-Wilk, Kolmogorov-Smirnov, and Anderson-Darling normality tests, validates each against `scipy.stats`, and directly compares which test catches which *kind* of non-normality best — closing with the honest practical question these tests all eventually raise: what to actually do when normality fails.

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
Q-Q plot         -> VISUAL check: sample quantiles vs theoretical Normal quantiles
Shapiro-Wilk      -> formal test, generally the most powerful for small-to-moderate n
Kolmogorov-Smirnov -> compares empirical CDF to theoretical CDF, sup-norm distance
Anderson-Darling  -> like KS but weights the TAILS more heavily

    All of these get MORE powerful as n grows -- which cuts both ways:
    genuinely useful for detecting real problems, but also capable of flagging
    trivial, practically-irrelevant deviations as "significant" at large n.

    When normality genuinely fails: non-parametric tests (Topic 07) or
    permutation tests (Topic 12) are the principled fallback.
```

---

## 🎯 Why This Topic Matters

- **Shapiro-Wilk correctly distinguished genuinely Normal from genuinely skewed data with a dramatic contrast** — §2 found p=**0.8290** for real Normal data versus p=**2.17e-12** for Exponential-derived data — a difference of many orders of magnitude, not a borderline call.
- **The Kolmogorov-Smirnov statistic was built from scratch and matched `scipy.stats.kstest` to machine precision** — §3 found a manual D-statistic of 0.038596 agreeing with `scipy` to the sixth decimal place.
- **The Anderson-Darling statistic was also built from scratch and matched `scipy.stats.anderson` exactly** — §4 found a manual A² of 0.274703 matching `scipy`'s computation identically, with `scipy`'s own critical value table ([0.559 at 15%, 0.749 at 5%, 1.031 at 1%]) giving a direct way to read off significance without a p-value.
- **Four genuinely distinct types of non-normality were tested side by side, revealing each test's differing sensitivity profile** — §5 found the heavy-tailed distribution (t, df=3) producing the largest KS statistic (0.2068) while the skewed Exponential produced the largest Anderson-Darling statistic (7.4713) — different tests emphasizing different parts of the distribution's shape, not interchangeable measurements of "how non-normal."
- **A genuine bug was caught, diagnosed, and fixed before publication: an initial "trivial deviation" simulation produced the exact opposite of its intended finding** — the first version showed Shapiro-Wilk's rejection rate *decreasing* from 6% to 2% as sample size grew from 20 to 5,000, when a real skew injection should make rejection rate *increase* with n. Direct investigation found the injected skew's population value was effectively zero (~-0.01) — the exponential component's variance was drowned out by the Normal noise. Recalibrating the injection to a population skewness of ~0.18 produced the correct pattern.
- **After that fix, the corrected simulation showed exactly the intended, now-verified pattern** — §6 found a fixed, small population skewness (~0.18) rejected only **6.00%** of the time at n=20, but **99.33%** of the time at n=5,000 — the same real-but-modest deviation becoming almost universally "significant" purely from sample size, the normality-testing analogue of Topic 04's large-n-inflates-significance finding.
- **When normality genuinely failed on real skewed data, three different response strategies were compared side by side** — §7 found the t-test (p=0.2137), Mann-Whitney U from Topic 07 (p=0.3193), and a permutation test from Topic 12 (p=0.2112) all reaching the same qualitative conclusion on the same skewed data, reassuring agreement across three philosophically different approaches.
- **Multiple tests converged on the same real-data conclusion, strengthening confidence in the finding** — §8 found diamond price (skewness=1.610) showing dramatically larger KS (0.1722 vs 0.0784) and Anderson-Darling (23.47 vs 4.84) statistics than diamond depth (skewness=0.228) — Shapiro-Wilk, KS, and Anderson-Darling all telling the same consistent story about which variable's departure from normality is practically, not just statistically, real.

---

## 🧮 Mathematical Explanation

### 1. The Q-Q plot

Plots the $i$-th order statistic of the sample against $\Phi^{-1}\left(\frac{i-0.5}{n}\right)$, the theoretical Normal quantile at the same plotting position. A straight $y=x$ line indicates normality; systematic curvature reveals the specific kind of departure (concave/convex for skew, S-shaped for heavy/light tails).

### 2. Kolmogorov-Smirnov

$$D = \sup_x |F_n(x) - F(x)|$$

The largest vertical gap between the empirical CDF $F_n$ and the theoretical CDF $F$ — a single worst-case measure, giving equal weight to every part of the distribution.

### 3. Anderson-Darling

$$A^2 = -n - \frac{1}{n}\sum_{i=1}^n (2i-1)\left[\ln F(x_i) + \ln(1-F(x_{n+1-i}))\right]$$

A weighted integral of the squared CDF difference that up-weights the tails — often more powerful than KS specifically for detecting heavy- or light-tailed departures.

### 4. Why power grows (sometimes too much) with n

Every consistent hypothesis test's power increases with sample size for a fixed true effect (a theme running throughout this repo since Topic 04). For normality tests specifically, this means even a population skewness as mild as 0.18 — often practically negligible — becomes detected essentially every time once $n$ is large enough, exactly what §6 measured directly.

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| Q-Q plotting position | $\Phi^{-1}((i-0.5)/n)$ | §1 |
| KS statistic | $\sup_x\vert F_n(x)-F(x)\vert$ | §3 |
| Anderson-Darling | $-n-\frac1n\sum(2i-1)[\ln F(x_i)+\ln(1-F(x_{n+1-i}))]$ | §4 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Trusting a normality test's p-value without also looking at a Q-Q plot (or vice versa)** — §8's diamond price and depth comparison shows both agreeing, but the Q-Q plot makes the *magnitude* of the departure immediately visually obvious in a way a p-value alone doesn't.
2. **Assuming all three tests (Shapiro, KS, Anderson-Darling) are interchangeable** — §5 showed each test's statistic responding differently across four distinct non-normality types; a test that misses one kind of deviation may catch another.
3. **Treating a "significant" normality test result at large n as automatically practically important** — §6's corrected finding (6% → 99.33% rejection rate for a fixed, modest 0.18 skewness as n grows from 20 to 5,000) is the clearest possible warning against this.
4. **Assuming a small p-value from a normality test means the parametric method must be abandoned** — §7 found the t-test, Mann-Whitney U, and a permutation test all reaching the *same* conclusion despite Shapiro-Wilk formally rejecting normality (p=0.0001); moderate non-normality doesn't always change the substantive conclusion.
5. **Constructing a "small deviation" simulation without directly checking its actual population-level magnitude** — this notebook's own initial attempt at §6 had an injected skew that was, by direct calculation, essentially zero (~-0.01) despite looking like a real skew injection in the code; always verify a simulated effect's true size numerically before trusting the simulation's results.
6. **Comparing raw Anderson-Darling or KS statistics across datasets without accounting for how each test's scale differs** — §8's AD statistics (23.47 vs 4.84) and KS statistics (0.17 vs 0.08) live on completely different numeric scales; compare each test's own statistic across datasets, not one test's statistic to another's.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `scipy.stats.probplot` | Reference Q-Q plot construction |
| `scipy.stats.shapiro` | Shapiro-Wilk normality test |
| `scipy.stats.kstest` | Kolmogorov-Smirnov goodness-of-fit test |
| `scipy.stats.anderson` | Anderson-Darling test with critical value table |
| `scipy.stats.mannwhitneyu` | Non-parametric alternative (Topic 07), reused for comparison |
| `np.random.default_rng().permutation` | Permutation test (Topic 12), reused for comparison |

---

## 📝 Self-Test Exercises

1. Section 5 found the heavy-tailed t-distribution producing the largest KS statistic while the skewed Exponential produced the largest Anderson-Darling statistic. Using the formulas in §2-§3 of the Mathematical Explanation, explain why a distribution with heavy tails but a symmetric center might affect KS's single worst-case gap differently than it affects Anderson-Darling's tail-weighted sum.
2. Section 6 initially showed a *decreasing* rejection rate (6% → 2%) due to an accidentally near-zero injected skew, then a corrected version showed the intended increasing pattern (6% → 99.33%) after recalibrating the injection to a population skewness of ~0.18. If you were building a similar simulation yourself, what specific check would you add to your code, before running the full simulation, to avoid this exact class of bug?
3. Section 7 found the t-test, Mann-Whitney U, and permutation test all agreeing (p≈0.21-0.32) on data where Shapiro-Wilk had rejected normality with p=0.0001. Does this agreement mean the normality rejection was a "false alarm," or are both facts (normality rejected, AND all three test methods agreeing) simultaneously true and non-contradictory? Explain your reasoning.
4. Section 8 found diamond price's Shapiro-Wilk p-value (6.01e-23) far smaller than depth's (1.25e-09), even though both were technically "rejected" at n=500. Using the skewness values reported (1.610 vs 0.228), argue for which of these two rejections represents a more practically meaningful departure from normality, and why a practitioner might treat them differently despite both being "statistically significant."
5. Using the Q-Q plot plotting-position formula from §1 of the Mathematical Explanation, explain what would happen to the theoretical quantile assigned to the single largest observation in a sample as $n$ grows very large — does it approach a finite limit or diverge, and why does that make sense given the Normal distribution's unbounded support?

---

## 📓 Notebook

30 executed code cells: a Q-Q plot built from scratch and matched to `scipy.stats.probplot`'s plotting-position convention, Shapiro-Wilk correctly distinguishing Normal from skewed data by many orders of magnitude, the Kolmogorov-Smirnov statistic built from scratch and matched to `scipy` exactly, the Anderson-Darling statistic built from scratch and matched to `scipy` exactly with its critical-value table, four genuinely distinct non-normality types (skewed, heavy-tailed, light-tailed, bimodal) compared across all three tests, an honestly-diagnosed-and-fixed simulation bug (an initially near-zero injected skew corrected to a real ~0.18 population skewness) demonstrating rejection rate climbing from 6.00% to 99.33% as sample size grows from 20 to 5,000, a three-way comparison of the t-test, Mann-Whitney U, and a permutation test on the same skewed data, and a full real-data normality investigation contrasting diamond price against diamond depth across all three formal tests plus Q-Q plots:

➡️ **[17_goodness_of_fit_normality_tests.ipynb](17_goodness_of_fit_normality_tests.ipynb)**

---

## 📚 Further Reading

- [scipy.stats.shapiro / kstest / anderson documentation](https://docs.scipy.org/doc/scipy/reference/stats.html)
- [Shapiro & Wilk (1965): An Analysis of Variance Test for Normality](https://www.jstor.org/stable/2333709)
- [Anderson & Darling (1954): A Test of Goodness of Fit](https://www.jstor.org/stable/2281537)
- [D'Agostino & Stephens, *Goodness-of-Fit Techniques*](https://www.routledge.com/Goodness-of-Fit-Techniques/DAgostino-Stephens/p/book/9780824774875)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.17_Goodness_of_Fit_Normality_Tests&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
