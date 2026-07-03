# 🧮 Multiple Testing Correction

> Status: ✅ Complete — [Open the notebook →](11_multiple_testing_correction.ipynb)

Topic 11 of the Statistical Inference & Hypothesis Testing repo. Topics 06, 07, and 10 each previewed this problem informally: pairwise t-tests after ANOVA, pairwise Mann-Whitney after Kruskal-Wallis, and tracking many A/B test metrics all inflate the false-positive rate. This notebook formalizes the fix — building Bonferroni, Holm, Sidak, and Benjamini-Hochberg from scratch, validating each against `statsmodels`, and directly measuring the actual tradeoff between the two fundamentally different error-control philosophies: family-wise error rate (FWER) and false discovery rate (FDR).

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
Run m hypothesis tests  ->  each has its own p-value
                |
   FWER methods (Bonferroni, Holm, Sidak):
   "what's the chance of EVEN ONE false positive among all rejections?"
   -> strict, conservative, good when any single mistake is costly
                |
   FDR methods (Benjamini-Hochberg):
   "what FRACTION of my rejections do I expect to be wrong?"
   -> more permissive, more powerful, good when many true effects
      are expected and some mistakes among many discoveries is acceptable
                |
   This is the SAME family-wise error rate problem previewed informally
   in Topics 06 (pairwise t-tests), 07 (Kruskal-Wallis post-hoc), and
   10 (A/B test multiple metrics) -- now formalized and fixed properly.
```

---

## 🎯 Why This Topic Matters

- **All four correction methods were built from scratch and matched `statsmodels` exactly** — §2's Bonferroni, §3's Holm, and §5's Benjamini-Hochberg each produced identical rejection decisions to `statsmodels.stats.multitest.multipletests` on the same 10 demo p-values.
- **The three methods gave three genuinely different, principled answers on identical data** — §6 found Bonferroni and Holm each rejecting just **1 of 10** tests, while Benjamini-Hochberg rejected **3 of 10** — the same p-values, different error-control philosophy, different conclusion.
- **Sidak's relative advantage over Bonferroni was measured precisely across scale, and it doesn't behave the way intuition might suggest** — §4 found Sidak's corrected alpha was only **1.28%** larger (relatively) than Bonferroni's at m=2 tests, growing slightly to **2.56%** larger at m=100 — a small but real, measured increase, not the "advantage shrinks with scale" pattern a first guess might predict. In absolute terms both correction methods still converge toward the same tiny alpha as m grows, which is the practically relevant takeaway: Sidak's edge over Bonferroni stays negligible either way once m is large.
- **A single mixed scenario (80 true nulls, 20 true effects) showed the correction working exactly as designed** — §6 found the uncorrected approach falsely rejecting **1 of 80** true nulls while only detecting **8 of 20** true effects (40% power), while all three corrected methods eliminated the false positive but also lost power, dropping to just 1 true positive detected (5% power) in that single run.
- **A 200-repeat averaged simulation revealed the real, reliable power difference between FWER and FDR methods** — §7 found Benjamini-Hochberg's average power (**9.8%**) nearly **double** Bonferroni's and Holm's (**5.2%** each) — while BH's realized FDR (**5.1%**) stayed appropriately controlled near the nominal 5% target, versus the uncorrected approach's realized FDR of **28.2%**.
- **Proper correction was applied to Topic 10's own multiple-metrics scenario and independently reproduced its finding almost exactly** — §8 found the uncorrected approach's false-"winner" rate was **39.75%** across 2,000 repeated null experiments, matching Topic 10's independently-simulated **39.25%** to within 0.5 percentage points — strong cross-validation between two separately-written notebooks.
- **The fix was shown actually working on that same scenario** — §8 found Bonferroni and Holm both bringing the false-winner rate down to **4.80%**, and Benjamini-Hochberg to **5.00%** — all landing right at the nominal 5% target, resolving the exact problem Topic 10 flagged but left unresolved.

---

## 🧮 Mathematical Explanation

### 1. Bonferroni

$$\text{reject } H_{0,i} \text{ if } p_i \le \alpha/m$$

Guaranteed FWER control via the union bound, regardless of dependence structure between tests — the most conservative of the FWER methods covered here.

### 2. Holm-Bonferroni (step-down)

Sort $p_{(1)} \le p_{(2)} \le \dots \le p_{(m)}$. Reject $H_{0,(i)}$ for all $i$ up to the first $i$ where $p_{(i)} > \frac{\alpha}{m-i+1}$. Still guarantees FWER control, but uniformly at least as powerful as Bonferroni — §3 confirmed Holm never rejects fewer tests than Bonferroni on the same data.

### 3. Sidak

$$\alpha_{\text{Sidak}} = 1-(1-\alpha)^{1/m}$$

Exact FWER control under independence — slightly less conservative than Bonferroni's more general (dependence-agnostic) bound, with the gap shrinking as $m$ grows (§4).

### 4. Benjamini-Hochberg (step-up)

Sort ascending, find the largest $i$ such that $p_{(i)} \le \frac{i}{m}\alpha$, reject all tests up to and including that $i$. Controls

$$FDR = E\left[\frac{V}{\max(R,1)}\right] \le \alpha$$

where $V$ is the number of false rejections and $R$ the total number of rejections — a fundamentally different (and, under many true effects, more powerful) guarantee than FWER.

---

## 📋 Formula Reference

| Method | Correction | Guarantee | Notebook Section |
|---|---|---|---|
| Bonferroni | $\alpha/m$ | FWER $\le \alpha$ | §2 |
| Holm | step-down, $\alpha/(m-i+1)$ | FWER $\le \alpha$ | §3 |
| Sidak | $1-(1-\alpha)^{1/m}$ | FWER $\le \alpha$ (independent tests) | §4 |
| Benjamini-Hochberg | step-up, $(i/m)\alpha$ | FDR $\le \alpha$ | §5 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Assuming one "correct" multiple-testing method exists** — §6's same 10 p-values produced 1, 1, and 3 rejections under Bonferroni, Holm, and BH respectively; the right choice depends on which error rate genuinely matters for the decision being made.
2. **Using Bonferroni by default without considering Holm** — §3 confirmed Holm is never less powerful than Bonferroni on identical data (it can only reject the same or more), making it a strictly better choice whenever plain FWER control is the goal.
3. **Assuming FDR control means "no error control"** — §7's realized FDR for Benjamini-Hochberg (5.1%) stayed right at the nominal target across 200 repeats; FDR control is a real, quantified guarantee, just a different one than FWER.
4. **Reaching for Sidak expecting a meaningful power boost over Bonferroni** — §4 showed the relative gap between the two staying under 3% across the entire tested range (m=2 to m=100), with both corrected-alpha values converging toward the same tiny number regardless.
5. **Trusting a single simulation run's power numbers** — §6's single mixed-scenario run showed all three corrected methods finding just 1 true positive (5% power), which looks identical across methods; only the 200-repeat average in §7 revealed BH's real, nearly 2x power advantage — the same single-run-vs-repeated-trials lesson from Topic 07.
6. **Ignoring the multiple-testing problem in A/B test metric tracking** — §8's exact replication of Topic 10's ~39% false-winner finding, followed by the correction bringing it back to ~5%, is the concrete, actionable fix for that topic's flagged-but-unresolved problem.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `statsmodels.stats.multitest.multipletests` | Reference implementation for Bonferroni, Holm, and Benjamini-Hochberg |
| `scipy.stats.ttest_ind` | Generates the raw p-values used throughout the simulations |
| `np.argsort` | Core sorting step underlying the Holm and BH step-down/step-up procedures |

---

## 📝 Self-Test Exercises

1. Section 6 found Bonferroni and Holm both rejecting exactly 1 test on the 10 demo p-values, while BH rejected 3. Using the comparison table in §6, identify which 2 additional p-values BH rejected that Bonferroni/Holm did not, and explain using the BH threshold formula from §4 of the Mathematical Explanation why those specific p-values cleared BH's bar but not the stricter FWER bar.
2. Section 4 found Sidak's relative advantage over Bonferroni growing slightly, from 1.28% (m=2) to 2.56% (m=100), rather than shrinking. Using the Sidak formula from §3 of the Mathematical Explanation, is this small growth in the *relative* gap consistent with the formula's behavior as m increases, and why does the *absolute* gap between the two corrected alphas still become negligible at large m even as the relative gap grows?
3. Section 7 found BH's average power (9.8%) nearly double Bonferroni's (5.2%) over 200 repeats, but §6's single run showed both methods finding just 1 true positive. Explain why a single simulation run can fail to reveal a real, reliable difference that only emerges after averaging many runs.
4. Section 8 independently reproduced Topic 10's ~39% false-winner finding to within 0.5 percentage points using a completely separate simulation written for this notebook. What does this close agreement suggest about the reliability of both notebooks' underlying simulation methodology?
5. Using the definitions in §1 of the Mathematical Explanation, construct a hypothetical scenario (in words) where you would clearly prefer FWER control over FDR control, and a second scenario where you'd clearly prefer FDR control over FWER — justify each choice by what a false positive would actually cost in that scenario.

---

## 📓 Notebook

30 executed code cells: Bonferroni, Holm-Bonferroni, and Benjamini-Hochberg each built from scratch and matched exactly to `statsmodels.stats.multitest.multipletests`, a Sidak-vs-Bonferroni corrected-alpha comparison across 5 values of m, a three-way rejection comparison on identical demo p-values (1, 1, and 3 rejections respectively), a mixed-scenario simulation (80 true nulls, 20 true effects) showing the raw power/FDR tradeoff, a 200-repeat averaged simulation revealing Benjamini-Hochberg's real ~2x power advantage over Bonferroni/Holm while keeping realized FDR controlled near nominal, and a full re-application to Topic 10's multiple-metrics A/B testing scenario, independently reproducing its ~39% false-winner finding and then resolving it back down to ~5% with proper correction:

➡️ **[11_multiple_testing_correction.ipynb](11_multiple_testing_correction.ipynb)**

---

## 📚 Further Reading

- [statsmodels.stats.multitest documentation](https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html)
- [Holm (1979): A Simple Sequentially Rejective Multiple Test Procedure](https://www.jstor.org/stable/4615733)
- [Benjamini & Hochberg (1995): Controlling the False Discovery Rate](https://www.jstor.org/stable/2346101)
- [Šidák (1967): Rectangular Confidence Regions for the Means of Multivariate Normal Distributions](https://www.jstor.org/stable/2283989)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.11_Multiple_Testing_Correction&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
