# 🔗 Causal Inference Basics

> Status: ✅ Complete — [Open the notebook →](20_causal_inference_basics.ipynb)

Topic 20 of the Statistical Inference & Hypothesis Testing repo — **the final topic**. Every method built across this series measured *association*. This notebook confronts the harder question directly: when can an association be interpreted *causally*? It builds confounding bias, collider bias, propensity score matching, inverse propensity weighting, and regression adjustment from scratch, then closes with a bootstrap confidence interval around a real-data estimate — the standard toolkit for reasoning about causal effects from observational (non-randomized) data.

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
CONFOUNDER (biases naive comparisons):
   confounder ----> treatment
        |               |
        v               v
     outcome <----------+
   -> adjust for it (matching / IPW / regression)

COLLIDER (the opposite trap):
   treatment ----> collider <---- outcome
   -> do NOT adjust for it -- conditioning on it CREATES spurious association

Three ways to adjust for a measured confounder, cross-checked against each other:
   1. Matching       -- pair treated/control units with similar propensity scores
   2. IPW             -- reweight every unit by 1/P(received treatment they got)
   3. Regression       -- put the confounder directly in the outcome model (Topic 18's ANCOVA logic)

None of them can fix bias from an UNMEASURED confounder.
```

---

## 🎯 Why This Topic Matters

- **A known, engineered causal effect showed the naive comparison badly overestimating it** — §1 built data with a true causal effect of exactly **2.0**, but the naive treated-vs-control difference in means came out at **4.1457** — more than double, because the confounder's group means differed sharply (**-0.392** in controls vs **0.350** in treated), and the confounder independently raises the outcome.
- **Collider bias was demonstrated as the mirror-image trap to confounding** — §1b built treatment and outcome as genuinely independent (unconditional correlation **0.0059**), then showed that conditioning on a collider variable caused by both manufactured a spurious correlation of **-0.4252** out of nothing — a concrete warning against reflexively "controlling for everything available."
- **The positivity assumption was checked directly, not assumed** — §2 found the estimated propensity scores ranging from **0.0497 to 0.9623**, with only **0.1%** of units at extreme values — confirming essentially every unit had a realistic chance of either treatment group, a precondition for matching/weighting to be trustworthy.
- **Matching dramatically reduced bias and was verified via covariate balance, not just the effect estimate** — §3-§4 found the matched-sample estimate at **1.9730** (bias shrunk from **2.1457** to **0.0270**), while the confounder's standardized mean difference dropped from **0.7925** before matching to **0.0056** after — well inside the conventional |SMD| < 0.1 threshold for good balance.
- **Three independently-derived adjustment methods converged on the true effect** — §5-§5b found IPW at **1.8706**, regression adjustment at **1.9239** (SE 0.0481), and matching at **1.9730** — all far closer to the true **2.0** than the naive **4.1457**, with mutual agreement across three different methods serving as reassuring cross-validation.
- **A dedicated section proved the central limitation honestly, not glossed over**: adjusting for a *measured* confounder while an *unmeasured* one remained left the IPW estimate still badly biased at **3.6704** (error of **1.6704** vs the true 2.0, barely better than the naive **4.3232**) — demonstrating that no adjustment method, however sophisticated, can correct for a confounder that was never observed.
- **Applied to real Titanic data, two independent adjustment methods agreed with each other and shifted the naive estimate** — §7 found the naive first-class survival advantage at **0.3377**, shrinking to **0.3323** (IPW) and **0.3276** (regression) once age and sex were adjusted for — the two adjusted estimates only **0.0047** apart from each other, consistent with age and sex (tied to the "women and children first" protocol from Topic 14) acting as genuine, measured confounders.
- **A bootstrap confidence interval (reusing Topic 12's resampling machinery) attached honest uncertainty to the final real-data estimate** — §8 found a 95% CI of **[0.2571, 0.4053]** around the regression-adjusted estimate of 0.3276, excluding zero — supporting a genuine adjusted association, while the README still stops short of calling it definitively causal, since Titanic data cannot rule out unmeasured confounders like exact deck location.

---

## 🧮 Mathematical Explanation

### 1. The confounding bias formula

$$E[Y\mid T=1] - E[Y\mid T=0] = \underbrace{\tau}_{\text{true causal effect}} + \underbrace{\left(E[X\mid T=1]-E[X\mid T=0]\right)\beta_X}_{\text{confounding bias}}$$

The naive difference in means equals the true effect *plus* a bias term that is zero only when the confounder $X$ is balanced across treatment groups — exactly what randomization guarantees and what §1's naive 4.1457-vs-true-2.0 gap illustrates when it isn't.

### 2. Propensity score

$$e(X) = P(T=1 \mid X)$$

Rosenbaum & Rubin (1983): conditioning on this single scalar is asymptotically as good as conditioning on the full covariate vector $X$, provided the propensity model is correctly specified.

### 3. Propensity score matching estimator

$$\hat\tau_{\text{match}} = \frac{1}{n_1}\sum_{i \in \text{treated}} \left(Y_i - Y_{j(i)}\right), \qquad j(i) = \arg\min_{j \in \text{control}} |e(X_i) - e(X_j)|$$

### 4. Inverse propensity weighting estimator

$$\hat\tau_{IPW} = \frac{1}{n}\sum_i \frac{T_iY_i}{e(X_i)} - \frac{1}{n}\sum_i\frac{(1-T_i)Y_i}{1-e(X_i)}$$

### 5. Standardized mean difference (balance check)

$$SMD = \frac{\bar X_{\text{treated}} - \bar X_{\text{control}}}{\sqrt{(s^2_{\text{treated}}+s^2_{\text{control}})/2}}$$

Unlike a raw mean difference, SMD is scale-free, making the conventional < 0.1 "good balance" threshold comparable across any covariate.

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| Confounding bias decomposition | $\tau + (\bar X_1-\bar X_0)\beta_X$ | §1 |
| Propensity score | $P(T=1\mid X)$ | §2 |
| Matching estimator | nearest-neighbor on $e(X)$ | §3 |
| IPW estimator | $\frac1n\sum \frac{TY}{e(X)} - \frac{(1-T)Y}{1-e(X)}$ | §5 |
| Standardized mean difference | $(\bar X_1-\bar X_0)/\sqrt{(s_1^2+s_0^2)/2}$ | §4 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Adjusting for a collider instead of a confounder** — §1b's manufactured **-0.4252** correlation out of two genuinely independent variables is a direct, quantified illustration; a causal DAG (or at least careful reasoning about what causes what) is required before deciding which variables to "control for" — more adjustment is not automatically better.
2. **Skipping the positivity check** — §2's propensity range [0.0497, 0.9623] with only 0.1% extreme units was verified directly; a dataset with many units near 0 or 1 would make matching/weighting unreliable, since those units have essentially no comparable counterparts.
3. **Trusting a matched or weighted estimate without checking covariate balance** — §4's explicit before/after SMD comparison (0.7925 → 0.0056) is what actually justifies believing the matched estimate, not the estimate's proximity to expectations alone.
4. **Treating agreement between adjustment methods as proof of correctness** — §5b's three methods converging (1.87-1.97) is reassuring but not conclusive, since all three share the exact same blind spot demonstrated in §6.
5. **Assuming any adjustment removes all confounding** — §6 is the notebook's central lesson: adjusting only for the *measured* confounder left the IPW estimate still off by 1.6704 when an unmeasured confounder was also at work; every observational causal claim carries this same untestable assumption.
6. **Reporting a causal-sounding point estimate without an uncertainty interval** — §8's bootstrap CI [0.2571, 0.4053] turns the bare Titanic estimate (0.3276) into a properly hedged statistical claim, consistent with this series' standard of never reporting a number without its uncertainty.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `sklearn.linear_model.LogisticRegression` | Estimates the propensity score $e(X)$ |
| `sklearn.neighbors.NearestNeighbors` | Finds nearest-propensity matches for matching |
| `statsmodels.api.OLS` | Regression adjustment (treatment + confounder model) |
| Custom `ipw_estimate()` | Inverse propensity weighted effect estimator |
| Custom `standardized_mean_difference()` | Covariate balance diagnostic |
| Custom `bootstrap_regression_effect()` | Resampling-based confidence interval (reuses Topic 12's logic) |

---

## 📝 Self-Test Exercises

1. Section 1 found the naive estimate (4.1457) roughly double the true effect (2.0). Using the confounding bias decomposition formula in §1 of the Mathematical Explanation, compute the implied bias term from the reported confounder group means (-0.3918 for controls, 0.3502 for treated) and the outcome's confounder coefficient (3), and check it against the reported bias of 2.1457.
2. Section 1b showed conditioning on a collider manufacturing a -0.4252 correlation between two independent variables. Without running any code, explain in words why restricting a sample to "collider above its median" — such as only studying hospitalized patients, where hospitalization depends on both an exposure and an outcome — could induce a similar spurious pattern in real data.
3. Section 4 found the confounder's standardized mean difference dropping from 0.7925 to 0.0056 after matching. If a second, unmeasured confounder existed in this same dataset, would checking the *measured* confounder's balance tell you anything about whether the unmeasured one was also balanced? Justify your answer using §6's finding.
4. Section 5b found three adjustment methods (matching 1.9730, IPW 1.8706, regression 1.9239) all landing within about 0.13 of the true effect (2.0), while agreeing with each other to within roughly 0.1 of one another. Propose one plausible reason IPW's estimate might differ slightly more from the true effect than matching's, given IPW uses every observation (including those with extreme propensity scores) while matching only uses matched pairs.
5. Section 8 found a 95% bootstrap CI of [0.2571, 0.4053] around the Titanic regression-adjusted estimate (0.3276). Explain why this interval excluding zero supports an adjusted association between class and survival, but still falls short of proving a causal claim — referencing §6's unmeasured-confounder result specifically.

---

## 📓 Notebook

32 executed code cells: a known causal effect (true=2.0) recovered from a badly biased naive estimate (4.1457) via matching (1.9730), IPW (1.8706), and regression adjustment (1.9239) — all three converging on the truth and cross-validating each other; a collider-bias demonstration manufacturing a spurious -0.4252 correlation between genuinely independent variables; a positivity diagnostic confirming propensity scores stayed in a matchable range; covariate balance directly verified before and after matching (SMD 0.7925 → 0.0056); an honest demonstration that adjustment fails silently when a confounder goes unmeasured (IPW still off by 1.6704 despite "adjusting"); and a full real-data application to the Titanic dataset comparing naive, IPW, and regression-adjusted first-class survival estimates, closing with a bootstrap 95% confidence interval [0.2571, 0.4053] around the final adjusted estimate. This is the closing notebook of the 20-topic Statistical Inference & Hypothesis Testing series:

➡️ **[20_causal_inference_basics.ipynb](20_causal_inference_basics.ipynb)**

---

## 📚 Further Reading

- [Rosenbaum & Rubin (1983): The Central Role of the Propensity Score in Observational Studies for Causal Effects](https://www.jstor.org/stable/2335942)
- [Pearl, *Causality: Models, Reasoning, and Inference*](https://bayes.cs.ucla.edu/BOOK-2K/)
- [Hernán & Robins, *Causal Inference: What If*](https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/)
- [Cinelli, Forney & Pearl (2022): A Crash Course in Good and Bad Controls](https://journals.sagepub.com/doi/10.1177/00491241221099552)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.20_Causal_Inference_Basics&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
