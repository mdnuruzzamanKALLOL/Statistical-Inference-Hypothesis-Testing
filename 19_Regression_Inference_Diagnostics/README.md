# 📉 Regression Inference & Diagnostics

> Status: ✅ Complete — [Open the notebook →](19_regression_inference_diagnostics.ipynb)

Topic 19 of the Statistical Inference & Hypothesis Testing repo. The Classical ML repo's Regression category built regression models to predict; this notebook builds the formal **inference machinery** underneath them — the F-test, coefficient t-tests, and confidence/prediction intervals — from scratch, then covers the diagnostic checks (heteroscedasticity, autocorrelation, multicollinearity, influential points) that determine whether those inferences can actually be trusted.

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
OLS: beta_hat = (X'X)^-1 X'y           -- the fit itself
F-test: does the model explain ANYTHING at all?
t-tests: is EACH coefficient individually non-zero?
CI (mean response)  vs  PI (one new point) -- PI always wider

Is the fit trustworthy?  ->  DIAGNOSTICS:
  Breusch-Pagan   -> constant error variance?
  Durbin-Watson   -> independent errors (no autocorrelation)?
  Shapiro-Wilk    -> Normal residuals?         (Topic 17's tool, reused)
  VIF             -> predictors not redundant with each other?
  Cook's distance -> no single point dominating the fit?
```

---

## 🎯 Why This Topic Matters

- **OLS, the F-test, and coefficient t-tests were all built from scratch and matched `statsmodels` exactly** — §1-§3 found manual coefficients, an F-statistic (305.85), and standard errors all agreeing with `statsmodels` to the sixth decimal place — including a real numerical-precision fix along the way: the naive `1 - CDF` p-value formula underflowed to exactly `0.0` in the far tail, while switching to the survival function `stats.f.sf` recovered the true value (1.29e-42), matching `statsmodels` exactly.
- **Confidence and prediction intervals were shown to be genuinely different widths, not just conceptually different** — §4 found a 95% confidence interval for the mean response at width **0.4688**, versus a 95% prediction interval for one new observation at width **3.7092** — nearly 8x wider, because the prediction interval must also account for the irreducible noise term that a confidence interval doesn't.
- **The Breusch-Pagan test correctly distinguished genuinely heteroscedastic from homoscedastic data** — §6 found p=**9.04e-09** on data with variance deliberately growing with the predictor, versus p=**0.2693** on data with genuinely constant variance — a clean, decisive contrast.
- **The Durbin-Watson test correctly distinguished independent from autocorrelated errors** — §7 found a statistic of **2.1424** (near the "no autocorrelation" value of 2) for independent errors, versus **0.3744** (near 0, indicating strong positive autocorrelation) for AR(1)-style dependent errors.
- **A genuinely surprising, honestly-investigated finding: heteroscedasticity and non-normality are separate violations, and this data proved it** — §8 initially assumed (in a first-draft narrative) that the heteroscedastic model's residuals would also fail a normality test, but the actual Shapiro-Wilk result was p=**0.1443** — *not* rejecting normality, despite Breusch-Pagan clearly catching the variance violation. Rather than force the expected narrative, the notebook was corrected to report this honestly: the two assumption violations are genuinely distinct, and this specific dataset happened to violate one without the other.
- **VIF correctly flagged a deliberately engineered collinear predictor pair with a dramatic number** — §9 found two predictors correlated at **0.9984** producing VIF values of **321.2** each — vastly above the conventional concern threshold of 10, exactly as designed.
- **Cook's distance correctly identified a single injected influential point by a factor of nearly 2,000** — §10 found one deliberately extreme point's Cook's distance at **18.76**, roughly **1,983x** larger than the average of the other 50 well-behaved points.
- **Every diagnostic, applied to real automotive data, found genuine, uncontrived violations** — §11 found the real `mpg` regression violating both homoscedasticity (Breusch-Pagan p=3.31e-06) and residual normality (Shapiro-Wilk p=1.31e-06), with a Durbin-Watson statistic of just **0.86** — a reminder that textbook-clean assumption satisfaction is the exception, not the rule, in real data.

---

## 🧮 Mathematical Explanation

### 1. OLS

$$\hat\beta = (X^TX)^{-1}X^Ty$$

The closed-form minimizer of squared residuals — a special case of the MLE machinery from Topic 16 under Normal errors.

### 2. F-test for overall significance

$$F = \frac{(SS_{\text{total}}-SS_{\text{residual}})/k}{SS_{\text{residual}}/(n-k-1)}$$

Tests $H_0: \beta_1=\dots=\beta_k=0$ jointly — the regression analogue of ANOVA's omnibus test (Topic 06).

### 3. Coefficient t-tests

$$t_j = \frac{\hat\beta_j}{SE(\hat\beta_j)}, \qquad SE(\hat\beta) = \sqrt{\hat\sigma^2\,[(X^TX)^{-1}]_{jj}}$$

### 4. Confidence vs. prediction interval

$$SE_{\text{mean}} = \hat\sigma\sqrt{x_0^T(X^TX)^{-1}x_0}, \qquad SE_{\text{individual}} = \hat\sigma\sqrt{1+x_0^T(X^TX)^{-1}x_0}$$

The extra "+1" inside the prediction interval's square root accounts for a new observation's own irreducible noise $\epsilon$, on top of the estimation uncertainty a confidence interval already captures — the source of §4's ~8x width gap.

### 5. Breusch-Pagan, Durbin-Watson, and VIF

$$DW \approx 2(1-\hat\rho), \qquad VIF_j = \frac{1}{1-R_j^2}$$

Breusch-Pagan regresses squared residuals on the predictors to test for variance patterns; Durbin-Watson measures first-order residual autocorrelation; $VIF_j$ comes from regressing predictor $j$ on all other predictors, exploding as $j$ becomes redundant with the others.

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| OLS | $(X^TX)^{-1}X^Ty$ | §1 |
| F-test | $\frac{(SS_t-SS_r)/k}{SS_r/(n-k-1)}$ | §2 |
| Coefficient t | $\hat\beta_j/SE(\hat\beta_j)$ | §3 |
| Prediction interval SE | $\hat\sigma\sqrt{1+x_0^T(X^TX)^{-1}x_0}$ | §4 |
| VIF | $1/(1-R_j^2)$ | §9 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Computing an extreme p-value via `1 - CDF` instead of the survival function** — this notebook's own F-test initially underflowed to exactly `0.0` for a genuinely tiny (1.29e-42) p-value; `scipy.stats.<dist>.sf` is the numerically stable alternative, always preferred for far-tail p-values.
2. **Reporting a confidence interval when a prediction interval is what's actually needed (or vice versa)** — §4's ~8x width difference means using the wrong one either wildly understates or overstates the real uncertainty for the question being asked.
3. **Assuming heteroscedasticity implies non-normal residuals, or vice versa** — §8's honestly-corrected finding is a direct counterexample; always check each assumption with its own dedicated test rather than assuming one violation implies another.
4. **Ignoring high VIF values** — §9's VIF of 321 on a deliberately collinear pair shows how dramatically multicollinearity can inflate coefficient standard errors, making individual coefficients' significance unreliable even when the model's overall fit looks fine.
5. **Not checking for influential points before trusting a fitted line** — §10's single injected point had a Cook's distance nearly 2,000x the average; one such point can meaningfully distort a small-to-moderate sample's fitted coefficients.
6. **Assuming a real dataset will cleanly satisfy OLS's assumptions** — §11's `mpg` model violated both homoscedasticity and normality; running the diagnostic suite is a standard, expected step for any real regression, not an optional afterthought.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `statsmodels.api.OLS` | Reference OLS fitting |
| `statsmodels.stats.diagnostic.het_breuschpagan` | Heteroscedasticity test |
| `statsmodels.stats.stattools.durbin_watson` | Autocorrelation test |
| `statsmodels.stats.outliers_influence.variance_inflation_factor` | Multicollinearity diagnostic |
| `statsmodels.stats.outliers_influence.OLSInfluence` | Cook's distance and other influence measures |
| `scipy.stats.f.sf` | Numerically stable F-distribution p-value |

---

## 📝 Self-Test Exercises

1. Section 1 found the manual F-test's p-value underflowing to exactly 0.0 using `1 - stats.f.cdf(...)`, fixed by switching to `stats.f.sf(...)`. Using floating-point precision concepts, explain why `1 - x` loses precision when `x` is extremely close to 1, and why a dedicated survival function avoids this.
2. Section 4 found the prediction interval about 8x wider than the confidence interval. Using the SE formulas in §4 of the Mathematical Explanation, would you expect this ratio to grow, shrink, or stay the same as the sample size n grows very large, holding the noise variance fixed?
3. Section 8's corrected finding showed heteroscedasticity present (Breusch-Pagan rejected) but residual normality NOT rejected (Shapiro-Wilk p=0.1443). Propose a plausible reason why a variance that grows with x might not necessarily produce a residual distribution shape extreme enough for Shapiro-Wilk to flag, even though the variance pattern itself is real.
4. Section 9 found VIF=321 for two predictors correlated at 0.9984. Using the VIF formula, calculate what $R_j^2$ value (from regressing one predictor on the other) would produce a VIF of exactly 321, and compare it to the squared correlation coefficient (0.9984²) reported in the notebook.
5. Section 11 found the real `mpg` model's Durbin-Watson statistic at 0.86, indicating meaningful positive autocorrelation. Given the `mpg` dataset isn't inherently time-ordered data (unlike a true time series), propose one plausible reason real cross-sectional automotive data might still show apparent residual autocorrelation if the rows happen to be sorted by a variable correlated with the outcome (e.g., model year).

---

## 📓 Notebook

30 executed code cells: OLS, the F-test, and coefficient t-tests all built from scratch and matched exactly to `statsmodels` (including a real numerical-precision fix for the F-test's p-value), confidence and prediction intervals shown at genuinely different widths (0.47 vs 3.71), Breusch-Pagan correctly distinguishing heteroscedastic from homoscedastic data, Durbin-Watson correctly distinguishing independent from autocorrelated errors, an honestly-corrected finding that heteroscedasticity didn't imply non-normal residuals on this particular dataset, VIF correctly flagging a deliberately collinear predictor pair (VIF≈321), Cook's distance correctly identifying a single injected influential point (nearly 2,000x the average), and a full real-data diagnostic workup on the `mpg` dataset revealing genuine violations of homoscedasticity and normality:

➡️ **[19_regression_inference_diagnostics.ipynb](19_regression_inference_diagnostics.ipynb)**

---

## 📚 Further Reading

- [statsmodels: Regression Diagnostics](https://www.statsmodels.org/stable/diagnostic.html)
- [Breusch & Pagan (1979): A Simple Test for Heteroscedasticity and Random Coefficient Variation](https://www.jstor.org/stable/1911963)
- [Durbin & Watson (1951): Testing for Serial Correlation in Least Squares Regression](https://www.jstor.org/stable/2332391)
- [Cook (1977): Detection of Influential Observation in Linear Regression](https://www.jstor.org/stable/1268249)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.19_Regression_Inference_Diagnostics&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
