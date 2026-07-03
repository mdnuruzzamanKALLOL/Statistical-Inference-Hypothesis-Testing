# 🧪 One-Sample & Two-Sample t-Tests

> Status: ✅ Complete — [Open the notebook →](04_one_sample_two_sample_t_tests.ipynb)

Topic 04 of the Statistical Inference & Hypothesis Testing repo. Topic 03 built a confidence interval for the difference of two means and noted informally that "if it excludes 0, that's evidence of a real difference." This notebook formalizes that into a proper hypothesis test — built from scratch, validated against `scipy.stats`, and checked the way Topic 03 checked CI coverage: by directly simulating the test's actual false-positive rate and actual power, rather than assuming the theory holds.

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
H0: no real effect  vs  H1: a real effect exists
                    |
     t = (observed difference) / (standard error of that difference)
                    |
     p-value = P(seeing a t this extreme, or more, IF H0 were actually true)
                    |
     p < alpha  ->  reject H0        p >= alpha  ->  fail to reject H0
                    |
   TWO WAYS TO BE WRONG:
   Type I error (alpha)  = reject H0 when it's actually true   (false positive)
   Type II error (beta)  = fail to reject H0 when it's false   (false negative, power = 1-beta)
```

---

## 🎯 Why This Topic Matters

- **The one-sample, two-sample, and paired t-tests were all built from scratch and matched `scipy.stats` to machine precision** — §1, §4, §6 each show manual t-statistics and p-values agreeing with `scipy.stats.ttest_1samp` / `.ttest_ind` / `.ttest_rel` to within `0.0e+00` difference.
- **The Type I error rate was directly measured across three significance levels, not just cited** — §2 simulated 10,000 samples where $H_0$ was true by construction and found an empirical false-positive rate of **4.97%** against a 5% nominal target; a follow-up sweep confirmed calibration at α=0.01 (empirical 0.96%), α=0.05 (5.26%), and α=0.10 (9.88%).
- **A full power curve was simulated, not just described** — §3 found power climbing from **4.60%** at effect=0 (correctly equal to α, since $H_0$ is then true) up to **99.12%** at effect=12, and separately showed a fixed, modest effect (size 3) becoming more detectable as n grew: power **9.10% → 19.10% → 49.83%** at n=10, 30, 100.
- **Student's pooled-variance test was caught genuinely miscalibrated under unequal variances and unequal sample sizes** — §5 built a case with $\sigma_a=5, n_a=10$ vs $\sigma_b=20, n_b=40$ (equal true means) and found Student's empirical Type I error rate was only **0.11%** (far too conservative) against Welch's **5.51%** (correctly calibrated) — a real, measured failure of the "just use pooled variance" default, and notably in the *conservative* direction here rather than the more commonly cited liberal direction, underscoring that the failure mode depends on which group has the smaller n **and** the smaller variance.
- **The cost of ignoring paired structure was demonstrated on identical data** — §6 built 25 subjects with a true within-subject treatment effect of +4, found the correct paired test's p-value effectively 0, then re-ran the *same* before/after values through an (incorrect) independent-samples test and got p=**0.274** — completely failing to detect a real, substantial effect purely because the within-subject correlation was thrown away.
- **The one-tailed/two-tailed p-value relationship was verified exactly** — §7 found the two-tailed p-value (0.1994) was almost exactly double the one-tailed p-value (0.0997), a ratio of **2.0000**.
- **A large-n "significant but meaningless" result was constructed and correctly flagged** — §8 used n=50,000 per group with a true difference of just 0.3 (against a population std of 15) and got p=**7.48e-06** — "highly significant" — while Cohen's d of **-0.0283** correctly identified the effect as negligible, directly illustrating that statistical and practical significance answer different questions.
- **Every method was re-applied to real data** — §9 tested whether lunch and dinner tip percentages genuinely differ (Welch's t=0.654, p=0.5138, d=0.075 — not significant, negligible effect) and whether the average tip percentage differs from the commonly cited 15% benchmark (t=2.763, p=0.0062 — significant, real dataset of 244 parties).

---

## 🧮 Mathematical Explanation

### 1. One-sample t-test

$$t = \frac{\bar{X} - \mu_0}{s/\sqrt{n}}, \qquad df = n-1$$

Tests whether a sample's mean plausibly came from a population with hypothesized mean $\mu_0$.

### 2. Type I and Type II error

|  | $H_0$ true | $H_0$ false |
|---|---|---|
| **Reject $H_0$** | Type I error (rate = $\alpha$) | Correct (power = $1-\beta$) |
| **Fail to reject** | Correct ($1-\alpha$) | Type II error (rate = $\beta$) |

§2 measured the left column empirically (does the real false-positive rate match $\alpha$?); §3 measured the right column (does power actually climb as theory predicts?).

### 3. Two-sample t-tests: Student's vs. Welch's

$$t_{\text{Student}} = \frac{\bar{X}_1-\bar{X}_2}{s_p\sqrt{1/n_1+1/n_2}}, \ s_p^2=\frac{(n_1-1)s_1^2+(n_2-1)s_2^2}{n_1+n_2-2}
\qquad
t_{\text{Welch}} = \frac{\bar{X}_1-\bar{X}_2}{\sqrt{s_1^2/n_1+s_2^2/n_2}}$$

Student's pools the two variances into one, assuming they're equal; Welch's does not, using the Welch-Satterthwaite df from Topic 03. §5 showed why "just assume equal variances" is a risky default.

### 4. Paired t-test

$$t = \frac{\bar{D}}{s_D/\sqrt{n}}, \quad D_i = X_{i,\text{after}} - X_{i,\text{before}}$$

Simply a one-sample t-test run on the within-subject differences — §6's function reuse (`one_sample_ttest_scratch(after - before, mu0=0)`) makes this structural relationship explicit in code, not just in the formula.

### 5. One-tailed vs. two-tailed p-values

$$p_{\text{two-tailed}} = 2 \times P(T \ge |t|) \qquad p_{\text{one-tailed}} = P(T \ge t) \text{ (or } P(T \le t) \text{)}$$

The direction must be committed to *before* seeing the data — choosing the tail after seeing which way the sample leans silently doubles the true false-positive rate.

### 6. Cohen's d

$$d = \frac{\bar{X}_1-\bar{X}_2}{s_{\text{pooled}}}$$

A standardized, sample-size-independent measure of effect magnitude — §8 is the clearest possible demonstration of why $p$-value and $d$ answer different questions: a huge $n$ can make a $d \approx 0$ effect arbitrarily significant.

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| One-sample t | $(\bar{X}-\mu_0)/(s/\sqrt{n})$ | §1 |
| Welch's t | $(\bar{X}_1-\bar{X}_2)/\sqrt{s_1^2/n_1+s_2^2/n_2}$ | §4-5 |
| Student's t (pooled) | $(\bar{X}_1-\bar{X}_2)/(s_p\sqrt{1/n_1+1/n_2})$ | §4-5 |
| Paired t | $\bar{D}/(s_D/\sqrt{n})$ | §6 |
| Cohen's d | $(\bar{X}_1-\bar{X}_2)/s_{\text{pooled}}$ | §8 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Defaulting to Student's pooled test without checking variance equality** — §5's measured 0.11% vs 5.51% Type I error gap is a real, not hypothetical, consequence.
2. **Running an independent-samples test on paired/repeated-measures data** — §6 showed a genuine effect (p≈0) becoming statistically invisible (p=0.274) purely from ignoring the pairing structure.
3. **Choosing a one-tailed test's direction after peeking at the data** — this silently doubles the effective false-positive rate; §7's exact 2.0000 ratio is the mechanism being exploited if the direction is chosen post hoc.
4. **Treating a tiny p-value as proof of a meaningful effect** — §8's p=7.48e-06 alongside d=-0.0283 is the clearest possible counterexample; always report an effect size alongside a p-value.
5. **Assuming "not statistically significant" means "no effect"** — §9's lunch-vs-dinner result (p=0.5138) only means this dataset didn't provide strong evidence of a difference, not proof there is none; absence of evidence isn't evidence of absence.
6. **Assuming low power only happens with tiny samples** — §3 showed even n=100 only reached 49.83% power for a modest effect size of 3, well under the conventional 80% power target; power depends on effect size and $\sigma$ too, not n alone.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `scipy.stats.ttest_1samp` | One-sample and paired (via difference) t-test reference |
| `scipy.stats.ttest_ind(equal_var=...)` | Student's (True) and Welch's (False) two-sample t-test reference |
| `scipy.stats.ttest_rel` | Paired t-test reference |
| `scipy.stats.t.cdf` / `.ppf` | Manual p-value and critical-value computation |
| `np.random.default_rng().normal` | Simulated samples for Type I error and power simulations |

---

## 📝 Self-Test Exercises

1. Section 2 found the empirical Type I error rate at α=0.05 was 5.26% (out of 5,000 simulations) — slightly above nominal. Is this difference large enough to suggest a real miscalibration, or is it consistent with expected Monte Carlo noise? (Hint: compare it to §2's first check, which used 10,000 simulations and got 4.97%.)
2. Using §3's power table, explain why power at effect=0 (4.60%) should theoretically equal the significance level α=0.05, and why the small discrepancy is expected.
3. Section 5 found Student's test was *too conservative* (0.11% Type I error, well below the 5% nominal) rather than too liberal. Using the pooled-variance formula from §3 of the Mathematical Explanation, explain how pooling a small, low-variance group ($\sigma_a=5, n_a=10$) with a larger, high-variance group ($\sigma_b=20, n_b=40$) could inflate the standard error estimate and make the test overly cautious rather than overly eager to reject.
4. Section 6 showed a paired test's p-value near 0 versus an incorrect independent-samples p-value of 0.274 on the *same* underlying data. Using the paired-test formula (which operates on differences $D_i$), explain in your own words why removing each subject's own baseline variability (via subtraction) can reveal an effect that gets buried when treated as two separate, unrelated groups.
5. Section 8 used n=50,000 per group to make a d=-0.03 effect "significant" at p=7.48e-06. Predict, without running any code, roughly what would happen to the p-value if n were reduced to 500 per group instead, keeping the same true 0.3-unit difference — would you still expect statistical significance?

---

## 📓 Notebook

30 executed code cells: a one-sample t-test built from scratch and matched exactly to `scipy.stats.ttest_1samp`, a 10,000-repeat Type I error rate simulation confirming 4.97% empirical vs 5% nominal (plus a 3-level alpha sweep), a power-curve simulation across 5 effect sizes (4.60% to 99.12%) and 3 sample sizes at a fixed modest effect, Student's and Welch's two-sample t-tests each matched exactly to `scipy.stats.ttest_ind`, a Type I error calibration test under deliberately unequal variances and unequal sample sizes exposing Student's test's real miscalibration (0.11% vs Welch's 5.51%), a paired t-test matched exactly to `scipy.stats.ttest_rel` plus a direct demonstration of statistical power lost when paired data is incorrectly treated as independent, a one-tailed vs two-tailed p-value comparison confirming the exact 2x relationship, a Cohen's d demonstration of a large-n "significant but meaningless" effect, and a full application of one-sample, two-sample, and effect-size methods to the real 244-row `tips` dataset:

➡️ **[04_one_sample_two_sample_t_tests.ipynb](04_one_sample_two_sample_t_tests.ipynb)**

---

## 📚 Further Reading

- [scipy.stats.ttest_1samp / ttest_ind / ttest_rel documentation](https://docs.scipy.org/doc/scipy/reference/stats.html)
- [Wasserman, *All of Statistics*, Ch. 10: Hypothesis Testing and p-values](https://link.springer.com/book/10.1007/978-0-387-21736-9)
- [Welch (1947): The Generalization of "Student's" Problem](https://www.jstor.org/stable/2332510)
- [Cohen (1988): Statistical Power Analysis for the Behavioral Sciences](https://www.routledge.com/Statistical-Power-Analysis-for-the-Behavioral-Sciences/Cohen/p/book/9780805802832)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.04_One_Sample_Two_Sample_t_Tests&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
