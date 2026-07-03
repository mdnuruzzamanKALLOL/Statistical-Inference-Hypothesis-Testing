# 🎯 Maximum Likelihood Estimation

> Status: ✅ Complete — [Open the notebook →](16_maximum_likelihood_estimation.ipynb)

Topic 16 of the Statistical Inference & Hypothesis Testing repo. Nearly every estimator built across this repo — the sample mean, the sample proportion, the Cox model's coefficients in Topic 15 — is a maximum likelihood estimate, without that principle ever being derived explicitly. This notebook builds it from first principles: the likelihood function, closed-form MLE derivations, numerical optimization when no closed form exists, Fisher Information and the Cramér-Rao bound, and the likelihood ratio test that generalizes many of this repo's earlier hypothesis tests.

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
L(theta) = P(observed data | theta)  -- how probable is what we saw, as a function of theta?
MLE: theta_hat = argmax L(theta)     -- the parameter value making the data MOST probable

Closed form exists  -> algebra gives theta_hat directly (mean, proportion, ...)
No closed form       -> numerical optimization of -log L(theta) (Gamma, Weibull, ...)

Fisher Information -> how sharply peaked the likelihood is -> how precise the MLE can be
Likelihood Ratio Test -> compares restricted vs full model likelihoods -> a general test framework
AIC/BIC -> penalize likelihood by model complexity -> choose between competing models
```

---

## 🎯 Why This Topic Matters

- **The familiar "sample proportion" formula was shown to genuinely BE the maximum likelihood estimate, not just a convenient plug-in** — §1 found a brute-force grid search over the Bernoulli likelihood surface (0.6993) landing within **0.0007** of the closed-form MLE formula (0.7000).
- **Topic 02's ddof=0 variance bias was explained as a direct mathematical consequence of MLE, not an arbitrary convention** — §2 found the MLE variance (78.30) genuinely smaller than the unbiased ddof=1 estimator (81.00) on the same data — MLE maximizes likelihood, which is a fundamentally different criterion from unbiasedness, and this is the concrete numeric proof.
- **Numerical optimization was verified to recover the exact same answer as closed-form calculus** — §2 found `scipy.optimize.minimize` landing within **0.000038** of the analytically-derived Normal MLE — confirming the general-purpose numerical approach is trustworthy even when a closed form is available to check against.
- **A genuinely closed-form-free MLE (the Gamma distribution) was solved numerically and matched `scipy`'s built-in fitter almost exactly** — §3 found a manual Nelder-Mead optimization (shape=3.5808) agreeing with `scipy.stats.gamma.fit` (shape=3.5807) to four decimal places — both landing meaningfully away from the true shape=3.0 due to genuine sampling noise at n=200, not a fitting error.
- **The Cramér-Rao lower bound was verified by direct simulation, not just cited as theory** — §4 found the empirical MLE variance staying within a ratio of **0.99 to 1.03** of the theoretical bound across five sample sizes (n=10 to n=2000) — the Bernoulli MLE achieves the bound essentially exactly, not just asymptotically.
- **Standard errors were obtained for a parameter with no analytical variance formula, using only the log-likelihood's curvature** — §5 computed a numerical Hessian at the Gamma MLE and derived 95% confidence intervals ([2.91, 4.25] for shape, [1.33, 2.00] for scale) that correctly contained the true parameters (3.0 and 2.0) — a genuinely general-purpose technique, usable for any MLE problem.
- **The likelihood ratio test was cross-validated against the exact binomial test, and their genuine difference was measured, not hidden** — §6 found the LRT's chi-square-approximation p-value (0.0696) differing from the exact binomial test's p-value (0.1153) — both correctly failing to reject $H_0$, but the gap between an asymptotic approximation and an exact finite-sample calculation was left visible rather than glossed over.
- **AIC correctly identified the true data-generating distribution among three candidates on synthetic data, using only fitted likelihoods** — §7 found the true Gamma model's AIC (990.40) decisively lower than the wrong Normal (1025.39) and Exponential (1115.85) models — with no information about the true generating process given to the selection procedure beyond the fitted log-likelihoods and parameter counts.
- **The same MLE + AIC machinery answered a genuine open question on real data** — §8 found a Lognormal distribution (AIC=2001.0) fitting real diamond carat weights meaningfully better than Gamma (2046.9) or Normal (2579.9) — a real, data-driven, non-obvious modeling decision.

---

## 🧮 Mathematical Explanation

### 1. Likelihood and log-likelihood

$$L(\theta) = \prod_{i=1}^n f(x_i\mid\theta), \qquad \ell(\theta) = \sum_{i=1}^n \log f(x_i\mid\theta)$$

Taking logs turns a product (numerically unstable for large $n$) into a sum, without changing where the maximum occurs.

### 2. Closed-form MLEs

$$\hat p_{MLE} = \bar{X} \text{ (Bernoulli)}, \qquad \hat\mu_{MLE}=\bar X, \ \hat\sigma^2_{MLE}=\frac{1}{n}\sum(x_i-\bar X)^2 \text{ (Normal)}$$

Derived by setting $\partial\ell/\partial\theta = 0$ and solving — §2 showed the Normal variance MLE is exactly Topic 02's biased ddof=0 estimator.

### 3. Fisher Information and the Cramér-Rao bound

$$I(\theta) = -E\left[\frac{\partial^2\ell}{\partial\theta^2}\right], \qquad \text{Var}(\hat\theta) \ge \frac{1}{nI(\theta)}$$

A theoretical floor on any unbiased estimator's variance — §4 confirmed the MLE achieves this floor almost exactly by direct simulation.

### 4. Observed Fisher Information (for numerical MLEs)

$$\widehat{\text{Cov}}(\hat\theta) \approx \left[-\nabla^2\ell(\hat\theta)\right]^{-1}$$

The inverse of the negative Hessian of the log-likelihood *at the MLE itself* — a general recipe for standard errors that works even when no analytical Fisher Information formula exists, as §5 demonstrated for the Gamma distribution.

### 5. Likelihood ratio test

$$\Lambda = -2\log\frac{L(\hat\theta_0)}{L(\hat\theta)} \xrightarrow{d} \chi^2_{df}$$

Compares a null-restricted model's maximized likelihood to the full model's — asymptotically chi-square distributed, an approximation §6 checked directly against the exact binomial test.

### 6. AIC and BIC

$$AIC = 2k - 2\ell(\hat\theta), \qquad BIC = k\ln(n) - 2\ell(\hat\theta)$$

Reward higher likelihood, penalize more parameters ($k$) — lower is better for both criteria.

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| Log-likelihood | $\sum\log f(x_i\mid\theta)$ | §1 |
| Normal MLE | $\bar X$, $\frac{1}{n}\sum(x_i-\bar X)^2$ | §2 |
| Cramér-Rao bound | $1/(nI(\theta))$ | §4 |
| Observed Fisher Info | $[-\nabla^2\ell(\hat\theta)]^{-1}$ | §5 |
| LRT statistic | $-2\log[L(\hat\theta_0)/L(\hat\theta)]$ | §6 |
| AIC / BIC | $2k-2\ell$ / $k\ln n - 2\ell$ | §7-§8 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Assuming MLE always produces unbiased estimates** — §2's ddof=0 variance MLE is a direct counterexample; MLE maximizes likelihood, a genuinely different criterion from unbiasedness.
2. **Trusting a numerical optimizer without a sanity check against a known closed form when one exists** — §2's cross-validation (numerical vs. analytical Normal MLE, agreeing to 0.00004) is exactly the kind of check worth doing before trusting a numerical MLE where no closed form is available.
3. **Treating the likelihood ratio test's p-value as exact** — §6's 0.0696 (LRT) vs. 0.1153 (exact binomial) gap is real; the chi-square reference distribution is only an asymptotic approximation, most accurate at larger sample sizes than this notebook's n=20 example.
4. **Comparing raw log-likelihoods across models with different numbers of parameters** — a model with more parameters can always achieve a higher (less negative) log-likelihood by construction; §7's AIC/BIC comparison is specifically designed to penalize that, which is why it — not raw log-likelihood — correctly selected the true model.
5. **Forgetting that Fisher Information (and thus the Cramér-Rao bound) depends on the true parameter value** — §4's bound was computed at the specific true $p=0.3$; the same formula gives a different bound at a different true $p$.
6. **Assuming the observed Fisher Information technique only works for "nice" distributions** — §5 applied it successfully to the Gamma distribution specifically because it has *no* closed-form variance formula, which is the whole point of the technique's generality.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `scipy.optimize.minimize` | Numerical optimization of the negative log-likelihood |
| `scipy.stats.norm.logpdf` / `.gamma.logpdf` | Log-density functions used to build custom likelihoods |
| `scipy.stats.gamma.fit` / `.lognorm.fit` / `.norm.fit` / `.expon.fit` | Reference MLE fitters for cross-validation and model comparison |
| `scipy.stats.binomtest` | Exact finite-sample test, cross-validated against the LRT |
| `np.linalg.inv` | Inverting the observed Fisher Information matrix for standard errors |

---

## 📝 Self-Test Exercises

1. Section 2 found the MLE variance (78.30) smaller than the unbiased estimator (81.00) on the same n=30 sample. Using the formulas in §2 of the Mathematical Explanation, compute the exact theoretical ratio $\hat\sigma^2_{MLE}/S^2_{\text{unbiased}}$ for n=30, and check it against the ratio of the two printed values.
2. Section 4 found the empirical MLE variance staying within a ratio of 0.99-1.03 of the Cramér-Rao bound across all tested sample sizes, not just asymptotically. Using the Bernoulli Fisher Information formula, explain why this particular MLE (the sample proportion) is special in achieving the bound exactly, rather than only in the large-n limit.
3. Section 6 found the LRT p-value (0.0696) noticeably smaller than the exact binomial test's p-value (0.1153), though both stay above 0.05. Using the LRT formula from §5 of the Mathematical Explanation, would you expect this gap between the chi-square approximation and the exact test to shrink or grow if the sample size were increased from n=20 to n=200, holding the same observed proportion?
4. Section 7 found AIC correctly preferring the true Gamma model, with all three models' AIC differing by tens of points. Using the AIC formula, explain why a difference of even 2-3 AIC points is often considered practically meaningful, while a difference of 0.1-0.2 points usually isn't.
5. Section 8 found Lognormal fitting real diamond carat data better than Gamma or Normal by AIC. Given that carat weight is a strictly positive, right-skewed measurement, propose a physical or process-based reason why a lognormal (rather than normal) distribution might be a genuinely reasonable model for this kind of data, beyond just "it fit better."

---

## 📓 Notebook

30 executed code cells: the Bernoulli likelihood surface visualized and its MLE found by both grid search and closed form, the Normal distribution's MLE derived in closed form and matched to numerical optimization (explaining Topic 02's ddof=0 bias as a genuine MLE consequence), a Gamma distribution MLE solved numerically with no closed form available and matched to `scipy.stats.gamma.fit`, the Cramér-Rao bound verified by direct simulation across five sample sizes, standard errors derived from a numerically-computed observed Fisher Information matrix for the closed-form-free Gamma MLE, a likelihood ratio test built from scratch and cross-validated against the exact binomial test, an AIC/BIC model-selection demonstration correctly identifying the true generating distribution among three candidates, and a full application of MLE-based model selection to real diamond carat weight data:

➡️ **[16_maximum_likelihood_estimation.ipynb](16_maximum_likelihood_estimation.ipynb)**

---

## 📚 Further Reading

- [scipy.optimize.minimize documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)
- [Casella & Berger, *Statistical Inference*, Ch. 7: Point Estimation](https://www.taylorfrancis.com/books/mono/10.1201/9781003456978)
- [Fisher (1922): On the Mathematical Foundations of Theoretical Statistics](https://royalsocietypublishing.org/doi/10.1098/rsta.1922.0009)
- [Akaike (1974): A New Look at the Statistical Model Identification](https://ieeexplore.ieee.org/document/1100705)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.16_Maximum_Likelihood_Estimation&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
