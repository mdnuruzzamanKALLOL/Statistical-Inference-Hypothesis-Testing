# 📊 ANOVA (One-way & Two-way)

> Status: ✅ Complete — [Open the notebook →](06_anova_one_way_two_way.ipynb)

Topic 06 of the Statistical Inference & Hypothesis Testing repo. Topic 04 compared exactly two groups. The obvious-looking shortcut for three or more groups — just run every pairwise t-test — has a hidden cost this notebook measures directly: the false-positive rate inflates well past the nominal alpha. ANOVA fixes this with a single omnibus F-test, built from scratch here and validated against `scipy`/`statsmodels`, followed by post-hoc tests, effect size, assumption checking, and a two-way extension with interaction effects.

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
k groups, comparing means
                |
   BAD: run every pairwise t-test  ->  false-positive rate INFLATES badly as k grows
                |
   GOOD: one ANOVA F-test  =  MS_between / MS_within  ->  stays correctly calibrated
                |
   Significant F  ->  "at least one group differs" (not WHICH one)
                |
   Post-hoc (Tukey's HSD)  ->  properly-corrected pairwise comparisons, tells you WHICH
                |
   TWO-WAY ANOVA: two factors at once, PLUS an interaction term --
   does factor A's effect on the outcome depend on factor B's level?
```

---

## 🎯 Why This Topic Matters

- **One-way ANOVA was built from scratch and matched `scipy.stats.f_oneway` exactly** — §1 shows a manual F-statistic and p-value agreeing to a difference of `7.1e-15`.
- **The pairwise-t-test shortcut's real cost was measured directly, and it's dramatic** — §2 simulated groups that all shared the *same* true mean (so $H_0$ was true everywhere) and found the family-wise error rate climbing from **5.13%** at k=2 groups to **12.60%** at k=3, **28.20%** at k=5, and **51.20%** at k=8 — meaning with 8 equal-mean groups, running every pairwise t-test at α=0.05 has roughly a coin-flip's chance of producing at least one false "significant" result, purely from the number of comparisons.
- **A single ANOVA F-test was confirmed to stay correctly calibrated on the exact same 5-group scenario** — §2 found ANOVA's empirical Type I error rate was **5.20%**, versus the pairwise approach's **28.20%** on identical data — a direct, matched comparison, not two separate claims.
- **Power was simulated across increasing effect sizes** — §3 found power climbing from **4.97%** (zero true difference, correctly equal to α) up to **44.13%**, **96.70%**, and **100%** as the group means separated further apart.
- **Tukey's HSD was shown as the principled alternative to naive pairwise testing** — §4 found it correctly flagged **2 of 3** pairwise comparisons as significant (A vs C, and B vs C) while properly declining to flag the closer A vs B pair (p-adj=0.057) — controlling the family-wise error rate the way §2 showed naive pairwise t-tests fail to.
- **Eta-squared was shown behaving like Cohen's d / Cramer's V from Topics 04-05** — §5 found raw F growing roughly 7x (14.24 → 97.64) as sample size scaled 10x, while eta-squared stayed in a comparable range (0.247 → 0.166 → 0.179) — measuring effect strength largely independent of the sample-size-driven F inflation.
- **The homogeneity-of-variance assumption was checked, and a real violation's consequence was measured** — §6 found Levene's test correctly passing equal-variance groups (p=0.749) and correctly flagging deliberately unequal-variance groups (p<0.0001), where standard ANOVA (p=0.011) and the Welch-type Alexander-Govern test (p=0.0003) genuinely disagreed by over a full order of magnitude.
- **Two-way ANOVA correctly found NO interaction when the true data-generating process was purely additive** — §7 built two factors that acted independently and found the interaction term's p-value at **0.916** — correctly not-significant.
- **A genuine interaction effect was then constructed and correctly detected** — §8 built a case where one factor's effect (fertilizer B) flips sign depending on the other factor's level (helps under high sunlight, hurts under low sunlight) — a pattern neither factor's main effect alone can describe — and found the interaction term at p=**1.10e-14**, unambiguously significant, with the non-parallel-lines plot visually confirming it.
- **Every method was re-applied to real data** — §9 found total bill significantly differs across days of the week in the real `tips` dataset (F=2.77, p=0.043, small effect), while a two-way ANOVA of day × smoker status on tip percentage found no significant interaction (p=0.326).

---

## 🧮 Mathematical Explanation

### 1. One-way ANOVA

$$F = \frac{MS_{\text{between}}}{MS_{\text{within}}} = \frac{SS_{\text{between}}/(k-1)}{SS_{\text{within}}/(N-k)}$$

$$SS_{\text{between}} = \sum_{i=1}^{k} n_i(\bar{X}_i - \bar{X}_{\text{grand}})^2, \qquad SS_{\text{within}} = \sum_{i=1}^{k}\sum_{j=1}^{n_i}(X_{ij}-\bar{X}_i)^2$$

A large $F$ means between-group variation dominates within-group noise — evidence the group means genuinely differ.

### 2. Family-wise error rate

$$P(\text{at least one false positive across } m \text{ tests}) = 1-(1-\alpha)^m \quad \text{(if independent)}$$

With $m=\binom{k}{2}$ pairwise comparisons growing quadratically in $k$, this rate climbs fast — §2 measured it empirically rather than relying on the (only approximate, since pairwise t-tests aren't independent) formula above.

### 3. Tukey's HSD

Controls the family-wise error rate across all pairwise comparisons simultaneously by using the distribution of the *studentized range* rather than a per-pair t-distribution — the principled fix for §2's problem.

### 4. Eta-squared

$$\eta^2 = \frac{SS_{\text{between}}}{SS_{\text{total}}}$$

The proportion of total variance explained by group membership — bounded in $[0,1]$, and (unlike raw $F$) not mechanically inflated by sample size alone.

### 5. Two-way ANOVA with interaction

$$Y_{ijk} = \mu + \alpha_i + \beta_j + (\alpha\beta)_{ij} + \epsilon_{ijk}$$

$\alpha_i$ and $\beta_j$ are the two factors' main effects; $(\alpha\beta)_{ij}$ is the interaction term — nonzero exactly when one factor's effect on the outcome depends on the other factor's level, as §8 constructed directly.

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| One-way F | $MS_{\text{between}}/MS_{\text{within}}$ | §1 |
| Eta-squared | $SS_{\text{between}}/SS_{\text{total}}$ | §5 |
| Two-way model | $\mu+\alpha_i+\beta_j+(\alpha\beta)_{ij}+\epsilon$ | §7-8 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Running every pairwise t-test instead of one ANOVA** — §2's measured 51.20% family-wise error rate at k=8 groups (all sharing the same true mean) is the direct, quantified cost of this shortcut.
2. **Treating a significant ANOVA F-test as identifying which groups differ** — it only says *at least one* group differs; §4's Tukey HSD is required to identify which specific pairs.
3. **Comparing raw F-statistics across studies with different sample sizes** — §5 showed F growing substantially from sample size alone at a fixed true effect; always report eta-squared (or a similar normalized effect size) alongside F.
4. **Assuming equal variances without checking** — §6's Levene's test flagged a real violation, and standard vs. Welch-type ANOVA disagreed by an order of magnitude in p-value on that same data.
5. **Testing two factors with two separate one-way ANOVAs instead of one two-way ANOVA** — this structurally cannot detect an interaction effect at all; §8's interaction (fertilizer B helping under high sunlight, hurting under low sunlight) would be invisible to any one-factor-at-a-time analysis.
6. **Assuming an interaction always exists in a two-factor design** — §7 correctly found no interaction (p=0.916) when the true data-generating process was purely additive; don't over-interpret a two-way ANOVA table without checking whether the interaction term is actually significant.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `scipy.stats.f_oneway` | One-way ANOVA reference |
| `scipy.stats.f.cdf` | Manual p-value computation from the F-statistic |
| `scipy.stats.levene` | Homogeneity-of-variance assumption check |
| `scipy.stats.alexandergovern` | Welch-type ANOVA robust to unequal variances |
| `statsmodels.stats.multicomp.pairwise_tukeyhsd` | Tukey's HSD post-hoc test |
| `statsmodels.formula.api.ols` + `statsmodels.api.stats.anova_lm` | Two-way ANOVA via OLS regression decomposition |

---

## 📝 Self-Test Exercises

1. Section 2 found the family-wise error rate climbing from 5.13% (k=2) to 51.20% (k=8). Using the formula $1-(1-\alpha)^m$ from §2 of the Mathematical Explanation, compute the *approximate* predicted FWER at k=8 groups ($m=28$ pairwise tests) and compare it to the empirically measured 51.20% — are they close, and if not, what does that suggest about the independence assumption behind the formula?
2. Section 4's Tukey HSD found A vs B not significant (p-adj=0.057) while the direct one-way ANOVA (§1) rejected the overall null at p=0.000021. Explain why a significant omnibus F-test does not guarantee every pairwise comparison will also be significant.
3. Section 5 found eta-squared varying from 0.247 to 0.166 to 0.179 across three sample-size scales, rather than staying perfectly constant. Is this degree of variation more likely a sign of a real bug, or expected Monte-Carlo-style sampling noise given each scale used a freshly drawn random sample? Justify your answer using the sample sizes involved (30, 90, 300 per group).
4. Section 8 built an interaction where fertilizer B's effect flips sign depending on sunlight level. Using the two-way model formula from §5 of the Mathematical Explanation, describe in words what the interaction term $(\alpha\beta)_{ij}$ would need to look like numerically to produce this sign-flipping pattern.
5. Section 9 found total bill significantly differs across days (p=0.043) but with only a "small" eta-squared effect size (0.033). What does combining these two facts (statistically significant, but small effect) tell you about whether a restaurant should meaningfully change its staffing or pricing strategy by day of week, based on this evidence alone?

---

## 📓 Notebook

30 executed code cells: a one-way ANOVA built from scratch and matched exactly to `scipy.stats.f_oneway`, a family-wise error rate simulation showing the real cost of pairwise t-tests climbing from 5.13% to 51.20% as group count grows from 2 to 8 (with a single ANOVA F-test staying correctly calibrated at 5.20% on the identical 5-group scenario), a power curve across four effect sizes (4.97% to 100%), a Tukey's HSD post-hoc test correctly identifying 2 of 3 significant pairwise differences, an eta-squared effect-size implementation shown behaving differently from raw F under sample-size scaling, a Levene's-test-based homogeneity-of-variance check with a standard-vs-Welch-ANOVA comparison under a genuine violation, a two-way ANOVA via `statsmodels` correctly finding no interaction on purely additive synthetic data, a constructed genuine interaction effect (with a non-parallel-lines visualization) correctly detected at p=1.10e-14, and a full one-way and two-way ANOVA application to the real 244-row `tips` dataset:

➡️ **[06_anova_one_way_two_way.ipynb](06_anova_one_way_two_way.ipynb)**

---

## 📚 Further Reading

- [scipy.stats.f_oneway / levene / alexandergovern documentation](https://docs.scipy.org/doc/scipy/reference/stats.html)
- [statsmodels: ANOVA and OLS formula interface](https://www.statsmodels.org/stable/anova.html)
- [Tukey (1949): Comparing Individual Means in the Analysis of Variance](https://www.jstor.org/stable/3001913)
- [Wasserman, *All of Statistics*, Ch. 10: Hypothesis Testing (ANOVA extensions)](https://link.springer.com/book/10.1007/978-0-387-21736-9)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.06_ANOVA_One_way_Two_way&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
