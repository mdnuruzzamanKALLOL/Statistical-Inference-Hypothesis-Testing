# ⏳ Survival Analysis Basics

> Status: ✅ Complete — [Open the notebook →](15_survival_analysis_basics.ipynb)

Topic 15 of the Statistical Inference & Hypothesis Testing repo. Every prior topic assumed a fully observed outcome for every subject. Survival analysis handles a genuinely different data structure: **time-to-event** data where some subjects' true event time is never observed — only that it exceeds some cutoff (**censoring**). This notebook builds the Kaplan-Meier estimator, the log-rank test, the Nelson-Aalen cumulative hazard estimator, and the Cox Proportional Hazards model from scratch, validating each against `lifelines`, and demonstrates directly why ignoring censoring biases naive estimates.

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
Time-to-event data: for each subject, either
  (a) the event happened, at a known time  -- e.g. "rearrested at week 20"
  (b) the event had NOT happened when observation ended -- CENSORED, e.g. "still free at week 52"

Naive approaches (mean of observed times, or excluding censored subjects) are BOTH biased.

Kaplan-Meier    -> the survival curve S(t), correctly accounting for censoring
Log-rank test   -> compares two survival curves
Nelson-Aalen    -> the cumulative hazard H(t), showing WHEN risk concentrates
Cox PH model    -> extends survival analysis to include covariates (like regression)
```

---

## 🎯 Why This Topic Matters

- **Naive approaches to censored data were shown genuinely biased, with real numbers** — §1 found the true mean event time was **20.99**, while the naive mean of all observed/censored times was **13.11** and the naive mean using only subjects with an observed event was **8.70** — both substantially underestimating the truth, and in opposite ways (truncation vs. survivor exclusion).
- **The Kaplan-Meier estimator was built from scratch and matched `lifelines` to floating-point precision** — §2 found a maximum difference of **7.8e-16** across 123 event time points — as exact a match as double-precision arithmetic allows.
- **A genuine, honestly-investigated discrepancy in median survival time revealed a real convention difference, not a bug** — §3 found the manual median (15.0425, using the standard textbook definition $\inf\{t : S(t) \le 0.5\}$) differing from `lifelines`' reported median (15.0710) by 0.0285 — traced directly to `lifelines` using a strict-inequality convention that only disagrees with the textbook definition at the exact point a KM curve touches 0.5 precisely, which is what happened in this specific run.
- **The log-rank test correctly detected a real, constructed difference between two survival distributions** — §4 found p=0.000001 comparing two groups with genuinely different true survival distributions (mean 25 vs. mean 15), a decisive and correct rejection.
- **The Nelson-Aalen cumulative hazard estimator also matched `lifelines` to floating-point precision** — §5 found a maximum difference of **2.8e-15** across the same 123 time points.
- **A real Cox Proportional Hazards model was fit on the classic Rossi recidivism dataset (432 released prisoners, 114 rearrests)** — §6 found financial aid (`fin`) reducing rearrest hazard by **31.6%** (HR=0.6843, p=0.0474), age reducing it by **5.6%** per year (p=0.0090), and prior convictions (`prio`) increasing it by **9.6%** per prior (p=0.0014) — all real, checkable coefficients from a genuine dataset, not synthetic data.
- **The proportional hazards assumption was actually tested per-covariate, not just assumed** — §6 found `age` (p=0.0007) and `wexp` (p=0.0068) failing the proportional-hazards test — a real, honest finding that the constant-hazard-ratio assumption doesn't hold uniformly for every covariate in this model, exactly the kind of check a from-scratch or off-the-shelf Cox analysis should never skip.
- **A genuine borderline disagreement between the unadjusted and covariate-adjusted analyses was found and correctly explained** — §7 found the raw log-rank test on `fin` alone landing at p=**0.0501** (just barely *not* significant), while the Cox model — adjusting for age, prior record, and other covariates simultaneously — brought the same effect to p=**0.0474** (significant) — a real, concrete illustration of why covariate adjustment matters, not a contrived example.

---

## 🧮 Mathematical Explanation

### 1. The Kaplan-Meier estimator

$$\hat{S}(t) = \prod_{t_i \le t}\left(1 - \frac{d_i}{n_i}\right)$$

At each observed event time $t_i$, multiply the running survival estimate by the fraction of the at-risk set that survived that specific time point — censored subjects contribute to $n_i$ (the risk set) up until their censoring time, then drop out without ever counting as a $d_i$ (event).

### 2. Median survival time

$$\text{median} = \inf\{t : \hat{S}(t) \le 0.5\}$$

The textbook definition used in §3's manual implementation. `lifelines` uses a strict-inequality variant that agrees everywhere except at the exact point the curve touches 0.5 precisely — both are legitimate conventions for the same underlying concept.

### 3. The log-rank test

Compares, at every observed event time across both groups, the actual number of events in each group to the number *expected* under the null that both groups share the same hazard — summed into a chi-square statistic, the survival-analysis analogue of the two-sample tests built throughout this repo.

### 4. Nelson-Aalen cumulative hazard

$$\hat{H}(t) = \sum_{t_i \le t} \frac{d_i}{n_i}$$

The running sum (rather than product) of the same per-event-time hazard contributions — directly interpretable as the expected number of events by time $t$ under a Poisson-process approximation, and useful for visualizing exactly *when* risk concentrates.

### 5. Cox Proportional Hazards model

$$h(t \mid X) = h_0(t)\exp(\beta_1 X_1 + \beta_2 X_2 + \dots)$$

The baseline hazard $h_0(t)$ is left completely unspecified (semi-parametric); only the *relative* hazard ratios $\exp(\beta_j)$ are estimated via partial likelihood. The core assumption — that each $\exp(\beta_j)$ stays constant over time — is directly testable, as §6 demonstrated.

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| Kaplan-Meier | $\prod_{t_i\le t}(1-d_i/n_i)$ | §2 |
| Median survival | $\inf\{t:\hat S(t)\le 0.5\}$ | §3 |
| Nelson-Aalen | $\sum_{t_i\le t} d_i/n_i$ | §5 |
| Cox hazard | $h_0(t)\exp(\beta^TX)$ | §6-§7 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Computing the mean or median of observed times while ignoring censoring** — §1's 20.99-vs-13.11-vs-8.70 gap is a real, substantial bias, not a rounding artifact.
2. **Excluding censored subjects entirely instead of using proper survival methods** — §1's events-only naive mean (8.70) was the *most* biased of the two naive approaches, since it systematically drops exactly the longest-surviving subjects.
3. **Assuming two "median survival time" implementations must agree to the decimal** — §3's 0.0285 gap is a genuine, well-understood convention difference (strict vs. non-strict inequality at the 0.5 crossing), not evidence either implementation is wrong.
4. **Fitting a Cox model without checking the proportional hazards assumption** — §6 found 2 of 7 covariates (`age`, `wexp`) failing this test; taking every hazard ratio at face value without checking would silently misstate those two covariates' true, likely time-varying effects.
5. **Trusting only an unadjusted (log-rank) comparison when covariates plausibly confound the relationship** — §7's borderline log-rank p=0.0501 vs. Cox-adjusted p=0.0474 is a concrete illustration of exactly why: the raw two-group comparison alone would have led to the (wrong, in this case) "not significant" conclusion.
6. **Interpreting a Cox hazard ratio as an absolute risk statement** — hazard ratios describe *relative* risk between covariate levels at any given instant, not the probability of the event ever occurring; §6's C-index of 0.6403 is a better summary of how well the model's relative risk *ranking* performs overall.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `lifelines.KaplanMeierFitter` | Kaplan-Meier survival curve reference |
| `lifelines.NelsonAalenFitter` | Nelson-Aalen cumulative hazard reference |
| `lifelines.CoxPHFitter` | Cox Proportional Hazards model fitting |
| `lifelines.statistics.logrank_test` | Two-group survival comparison |
| `lifelines.statistics.proportional_hazard_test` | Per-covariate proportional-hazards assumption check |
| `lifelines.datasets.load_rossi` | Real recidivism dataset (432 released prisoners) |

---

## 📝 Self-Test Exercises

1. Section 1 found the events-only naive mean (8.70) was *more* biased than the all-observed naive mean (13.11), even though both underestimate the true mean (20.99). Using the definition of right censoring, explain why systematically excluding censored subjects tends to produce a *worse* bias than merely truncating their observed times.
2. Section 3 found a 0.0285 gap between the manual and `lifelines` median survival times. Using the definitions in §2 of the Mathematical Explanation, would you expect this same kind of gap to appear in a dataset where the KM curve never lands on *exactly* 0.5 at any observed event time? Explain your reasoning.
3. Section 6 found `age` and `wexp` failing the proportional-hazards test while `fin`, `mar`, `paro`, and `prio` did not. Propose one plausible real-world reason why a covariate like `age` might have a hazard effect that changes over the follow-up period, rather than staying constant.
4. Section 7 found the log-rank test (p=0.0501) and the Cox-adjusted test (p=0.0474) landing on opposite sides of the 0.05 threshold for the exact same `fin` effect. Using what you know about confounding from Topic 14 (Simpson's Paradox), explain how an unmeasured or unadjusted covariate correlated with both `fin` and the outcome could produce this kind of shift.
5. Section 6 reported a concordance index of 0.6403. Using the definition given in the Common Pitfalls section, explain why a C-index close to 0.5 would suggest the model's covariates carry almost no useful information for *ranking* subjects by risk, even if some individual hazard ratios in the summary table are statistically significant.

---

## 📓 Notebook

30 executed code cells: a direct demonstration of naive-estimate bias under censoring, the Kaplan-Meier estimator built from scratch and matched to `lifelines` to floating-point precision, an honestly-investigated median-survival-time convention difference, a log-rank test correctly detecting a genuine two-group survival difference, the Nelson-Aalen cumulative hazard estimator also matched to floating-point precision, a full Cox Proportional Hazards model fit on the real 432-subject Rossi recidivism dataset with hazard ratio interpretation and concordance index, a per-covariate proportional-hazards assumption test (finding 2 of 7 covariates genuinely violating it), and a closing cross-validation of the financial-aid effect via both an unadjusted log-rank test and a covariate-adjusted Cox model, landing on a genuine borderline disagreement that illustrates exactly why adjustment matters:

➡️ **[15_survival_analysis_basics.ipynb](15_survival_analysis_basics.ipynb)**

---

## 📚 Further Reading

- [lifelines documentation](https://lifelines.readthedocs.io/)
- [Kaplan & Meier (1958): Nonparametric Estimation from Incomplete Observations](https://www.jstor.org/stable/2281868)
- [Cox (1972): Regression Models and Life-Tables](https://www.jstor.org/stable/2985181)
- [Klein & Moeschberger, *Survival Analysis: Techniques for Censored and Truncated Data*](https://link.springer.com/book/10.1007/b97377)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.15_Survival_Analysis_Basics&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
