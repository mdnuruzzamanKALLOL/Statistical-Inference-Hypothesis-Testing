# 🔗 Correlation Analysis

> Status: ✅ Complete — [Open the notebook →](08_correlation_analysis.ipynb)

Topic 08 of the Statistical Inference & Hypothesis Testing repo. Every prior topic tested for a *difference* between groups. Correlation measures a *relationship* between two continuous variables instead. This notebook builds Pearson, Spearman, and Kendall correlation from scratch, validates each against `scipy.stats`, then spends real effort on the ways a single correlation number can mislead: non-linear-but-monotonic relationships, wildly different scatter patterns sharing a similar coefficient, and a constructed spurious correlation driven entirely by a hidden confounder.

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
Pearson r      -> LINEAR relationship strength (sensitive to raw magnitude, outliers)
Spearman rho   -> MONOTONIC relationship strength (rank-based, robust to outliers)
Kendall tau    -> concordant/discordant PAIR proportion (rank-based, even more robust)

                A single r number can NEVER fully describe a relationship:
                - it can be strong even when the true relationship is curved (if monotonic)
                - it can be misleadingly similar across wildly different scatter shapes
                - it can be strong and "significant" purely from a hidden confounder,
                  with NO direct causal link at all
```

---

## 🎯 Why This Topic Matters

- **All three correlation measures were built from scratch and matched `scipy.stats` exactly** — §1's Pearson r (0.9144), §3's Spearman rho (0.9093), and §4's Kendall's tau (0.7469) all agree with their `scipy` counterparts to machine precision.
- **The significance test was validated exactly, reusing Topic 04's t-distribution machinery** — §2 found the manual t-statistic (15.65, df=48) produced a p-value matching `scipy.stats.pearsonr`'s built-in p-value to `1.7e-20`.
- **A direct outlier-sensitivity comparison across all three measures showed rank-based methods are structurally more robust** — §4b injected one severe outlier and found Pearson r shifted by **0.5477** (0.914 → 0.367) while Spearman rho barely moved at **0.0025** and Kendall's tau at **0.0065** — a large, concrete gap, not a theoretical footnote.
- **Pearson's blindness to non-linear-but-monotonic structure was demonstrated directly** — §5 built a perfectly monotonic $y=x^3$ relationship and found Spearman's rho (0.9192) exceeding Pearson's r (0.8959), correctly detecting the near-perfect ranking relationship that Pearson partially misses.
- **An Anscombe-style set of datasets showed how differently "similar-looking" correlations can behave** — §6 built four datasets landing in a **0.623 to 0.985** range of Pearson r (a genuinely linear relationship, a strong curve, a linear relationship with one severe outlier inflating it, and a nearly-flat cloud rescued into "correlated" purely by one high-leverage point) — even without forcing them to an identical r, the visual comparison makes unmistakably clear that r alone says nothing about a relationship's actual shape.
- **A genuine spurious correlation was constructed with a known, verifiable ground truth** — §7 built "ice cream sales" and "drowning incidents," both driven by an unobserved confounder with **zero direct causal link** between them, and found a strong, highly "significant" correlation (r=0.9259, p=3.14e-43) anyway.
- **Partial correlation was shown actually removing the spurious relationship, not just claimed to** — §7 controlled for the hidden confounder directly and found the correlation collapsing from **0.9259** to **0.0402** — concrete proof the original relationship was confounder-driven, not direct.
- **A proper confidence interval for r was built via Fisher's z-transformation and shown narrowing with n** — §8 found a fixed r=0.5's 95% CI shrinking from a width of **1.048** at n=10 down to **0.132** at n=500 — the same $1/\sqrt{n}$ precision gain seen for the mean's CI in Topic 03.
- **Everything was re-applied to real data** — §9 found total bill and tip amount in the real `tips` dataset correlated at r=0.676 (Pearson) and rho=0.679 (Spearman) — close agreement, consistent with a relationship that's reasonably close to both linear and monotonic.

---

## 🧮 Mathematical Explanation

### 1. Pearson correlation

$$r = \frac{\sum(x_i-\bar{x})(y_i-\bar{y})}{\sqrt{\sum(x_i-\bar{x})^2}\sqrt{\sum(y_i-\bar{y})^2}}$$

Measures linear association only; sensitive to outliers since it uses raw deviations from the mean.

### 2. Significance test for r

$$t = r\sqrt{\frac{n-2}{1-r^2}}, \qquad df = n-2$$

Reuses the t-distribution from Topic 04, applied to a transformed correlation coefficient rather than a mean difference.

### 3. Spearman's rho

$$\rho_s = \text{Pearson}(\text{rank}(x), \text{rank}(y))$$

Pearson correlation computed on ranks — captures any monotonic relationship, linear or not.

### 4. Kendall's tau

$$\tau = \frac{(\text{concordant pairs}) - (\text{discordant pairs})}{\binom{n}{2}}$$

Based on pairwise ordering agreement rather than rank magnitudes — generally the most outlier-robust of the three, at higher computational cost ($O(n^2)$ pairs).

### 5. Fisher's z-transformation

$$z = \text{arctanh}(r) = \frac{1}{2}\ln\frac{1+r}{1-r}, \qquad SE_z = \frac{1}{\sqrt{n-3}}$$

$r$'s sampling distribution is skewed except near 0, so a direct Normal-based CI on $r$ is invalid; Fisher's z maps $r$ to an approximately Normal scale first, then the CI is transformed back via $\tanh$.

### 6. Partial correlation

$$r_{xy \cdot z} = \frac{r_{xy} - r_{xz}r_{yz}}{\sqrt{(1-r_{xz}^2)(1-r_{yz}^2)}}$$

The correlation between $x$ and $y$ with the linear influence of a third variable $z$ removed from both — §7's tool for unmasking a confounder-driven spurious correlation.

---

## 📋 Formula Reference

| Concept | Formula | Notebook Section |
|---|---|---|
| Pearson r | $\sum(x-\bar x)(y-\bar y)/(\ldots)$ | §1 |
| Significance t-test | $r\sqrt{(n-2)/(1-r^2)}$ | §2 |
| Spearman rho | Pearson on ranks | §3 |
| Kendall's tau | $(\text{conc}-\text{disc})/\binom{n}{2}$ | §4 |
| Fisher's z CI | $\text{arctanh}(r) \pm z_{\alpha/2}/\sqrt{n-3}$ | §8 |
| Partial correlation | $(r_{xy}-r_{xz}r_{yz})/\sqrt{(1-r_{xz}^2)(1-r_{yz}^2)}$ | §7 |

---

## ⚠️ Common Pitfalls & Gotchas

1. **Using Pearson r on a clearly non-linear (but monotonic) relationship** — §5's $y=x^3$ case showed Pearson (0.896) underselling the near-perfect monotonic relationship Spearman (0.919) correctly captured.
2. **Trusting a correlation coefficient without ever plotting the data** — §6's Anscombe-style datasets spanned r=0.623 to r=0.985 while looking nothing alike (a genuine line, a strong curve, a line dominated by one outlier, and a near-flat cloud rescued by one leverage point) — the number alone never tells the whole story.
3. **Concluding causation from a strong, significant correlation** — §7's ice cream/drowning example had r=0.9259 and p=3.14e-43 with **zero** direct causal link by construction; always ask what unobserved variable might drive both.
4. **Using a Normal-based CI directly on r** — since $r$'s sampling distribution is skewed except near 0, §8's Fisher's z-transformation is required for a valid interval, especially at small n or when $|r|$ is large.
5. **Assuming outliers affect all correlation measures equally** — §4b's single-outlier test showed Pearson shifting by 0.5477 versus Spearman's 0.0025 and Kendall's 0.0065 — over 200x more sensitive, not a marginal difference.
6. **Forgetting that partial correlation only controls for *measured* confounders** — §7's fix worked because the confounder was explicitly included; a completely unmeasured confounder can never be partialed out this way, only ruled out by careful study design.

---

## 🔑 Function Reference (used in the notebook)

| Function | Purpose |
|---|---|
| `scipy.stats.pearsonr` / `.spearmanr` / `.kendalltau` | Reference implementations for all three correlation measures |
| `scipy.stats.rankdata` | Rank computation underlying the manual Spearman implementation |
| `scipy.stats.t.cdf` / `.norm.ppf` | Manual significance test and Fisher's z critical value computation |
| `np.arctanh` / `np.tanh` | Fisher's z-transformation and its inverse |
| `seaborn.regplot` | Real-data scatter plot with fitted regression line |

---

## 📝 Self-Test Exercises

1. Section 4b found Pearson r shifted by 0.5477 from a single outlier while Spearman shifted by only 0.0025. Using the Pearson formula from §1 of the Mathematical Explanation, explain why one extreme value (with a huge deviation from the mean) can dominate the numerator's sum, while the same point in a rank-based measure can only ever become "the single highest rank" — a bounded contribution.
2. Section 5 found Spearman's rho (0.9192) exceeding Pearson's r (0.8959) on the $y=x^3$ relationship. Would you expect this gap to grow or shrink if the relationship were changed to $y=x^5$ (even more strongly curved, still monotonic)? Justify your answer.
3. Section 7 found the raw correlation (0.9259) collapsing to 0.0402 after partial correlation controlled for the confounder. Using the partial correlation formula from §6 of the Mathematical Explanation, explain what would happen to $r_{xy \cdot z}$ if $z$ had NO relationship at all with either $x$ or $y$ (i.e., $r_{xz}=r_{yz}=0$).
4. Section 8 found the 95% CI for a fixed r=0.5 narrowing from width 1.048 (n=10) to 0.132 (n=500) — roughly an 8x reduction from a 50x increase in n. Using the $1/\sqrt{n}$ relationship, is this narrowing rate consistent with what you'd expect, or slower/faster than a pure $1/\sqrt{n}$ scaling would predict?
5. Section 9 found Pearson r (0.676) and Spearman rho (0.679) very close on the real `tips` data. Propose a plausible reason why a real-world total-bill-vs-tip relationship might be expected to be closer to linear than, say, a biological growth-rate relationship might be.

---

## 📓 Notebook

30 executed code cells: Pearson, Spearman, and Kendall correlation each built from scratch and matched exactly to `scipy.stats`, a t-distribution-based significance test validated exactly against `scipy`'s built-in p-value, a direct three-way outlier-sensitivity comparison (Pearson shifted 200x+ more than the rank-based measures), a monotonic-non-linear ($y=x^3$) demonstration of Spearman correctly outperforming Pearson, four Anscombe-style datasets spanning r=0.623 to r=0.985 while looking nothing alike, a constructed spurious correlation with a known confounder (r=0.9259, collapsing to 0.0402 under partial correlation), a Fisher's z-transformation confidence interval shown narrowing from width 1.048 (n=10) to 0.132 (n=500), and a full Pearson/Spearman/Fisher's-z application to the real 244-row `tips` dataset:

➡️ **[08_correlation_analysis.ipynb](08_correlation_analysis.ipynb)**

---

## 📚 Further Reading

- [scipy.stats.pearsonr / spearmanr / kendalltau documentation](https://docs.scipy.org/doc/scipy/reference/stats.html)
- [Anscombe (1973): Graphs in Statistical Analysis](https://www.jstor.org/stable/2682899)
- [Wasserman, *All of Statistics*, Ch. 14: Association and Causation](https://link.springer.com/book/10.1007/978-0-387-21736-9)
- [Fisher (1915): Frequency Distribution of the Values of the Correlation Coefficient](https://www.jstor.org/stable/2331838)

---
[← Back to Statistical Inference & Hypothesis Testing](../README.md)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.08_Correlation_Analysis&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
