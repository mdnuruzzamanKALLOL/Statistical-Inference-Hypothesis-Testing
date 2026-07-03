# 📐 Statistical Machine Learning — Statistical Inference & Hypothesis Testing

Part of the [📘 Statistical Machine Learning for Noob](https://github.com/mdnuruzzamanKALLOL) series — a from-scratch, math-first, fully-executed ML study series designed to be read on GitHub, run locally, and actually learned from (not just skimmed).

This repo is the **statistics underneath the machine learning** — every p-value, confidence interval, and "statistically significant" claim used casually throughout the [Classical ML](https://github.com/mdnuruzzamanKALLOL/Statistical-Machine-Learning-Classical-ML) repo (VIF tests, the Evaluation Metrics topic's significance framing, A/B-testing-style model comparisons) gets its actual derivation here, from probability theory up through causal inference.

## 📑 Table of Contents

1. [Why This Repo Exists](#-why-this-repo-exists)
2. [Learning Path](#-learning-path)
3. [Topics](#-topics)
4. [Repository Structure](#-repository-structure)
5. [Prerequisites](#-prerequisites)
6. [How to Use](#-how-to-use)
7. [Datasets](#-datasets)
8. [The Full Series](#-the-full-series)
9. [Self-Test Philosophy](#-self-test-philosophy)
10. [Feedback & Corrections](#-feedback--corrections)
11. [Related](#-related)

---

## 🎯 Why This Repo Exists

Machine learning results are often reported with statistical language ("significant improvement", "95% confidence") without the underlying machinery ever being derived. This repo builds that machinery from scratch: every test statistic's formula, every distributional assumption, and every common pitfall (p-hacking, multiple comparisons, mistaking correlation for causation) — the same honesty-first, fully-executed-notebook standard as the rest of this series.

## 🗺️ Learning Path

```
01 Probability & Sampling Theory   ->  the foundation every test statistic is built on
        v
02-03 Point/Interval Estimation      ->  turning a sample into a population-level claim
        v
04-08 Classical Hypothesis Tests     ->  t-tests, chi-square, ANOVA, non-parametric, correlation
        v
09-12 Experimental Design & Rigor    ->  power analysis, A/B testing, multiple-comparison correction, resampling
        v
13-20 Advanced & Applied Inference   ->  Bayesian methods, categorical/survival analysis, MLE, ANCOVA/MANOVA,
                                          regression diagnostics, causal inference
```

Each topic explicitly says what it's setting up for — read the "Why This Topic Matters" section at the top of every README once topics are built out.

## 📚 Topics

| # | Topic | Key Concepts | Status |
|---|-------|-------------|:---:|
| 01 | [Probability Distributions & Sampling Theory](01_Probability_Distributions_Sampling_Theory/) | Normal, Binomial, Poisson distributions and the sampling distributions built from them | ✅ |
| 02 | [Central Limit Theorem & Point Estimation](02_Central_Limit_Theorem_Point_Estimation/) | Why sample means become normal, and how to estimate population parameters | ✅ |
| 03 | [Confidence Intervals](03_Confidence_Intervals/) | Quantifying uncertainty around a point estimate | ✅ |
| 04 | [One-Sample & Two-Sample t-Tests](04_One_Sample_Two_Sample_t_Tests/) | Testing whether a sample mean differs from a hypothesis, or two groups differ | ✅ |
| 05 | [Chi-Square Tests](05_Chi_Square_Tests/) | Testing goodness of fit and independence for categorical data | ✅ |
| 06 | [ANOVA (One-way & Two-way)](06_ANOVA_One_way_Two_way/) | Comparing means across three or more groups | ✅ |
| 07 | [Non-Parametric Tests](07_Non_Parametric_Tests/) | Distribution-free alternatives to t-tests and ANOVA | ✅ |
| 08 | [Correlation Analysis](08_Correlation_Analysis/) | Measuring linear and monotonic relationships between variables | ✅ |
| 09 | [Power Analysis & Sample Size Determination](09_Power_Analysis_Sample_Size_Determination/) | How many samples are enough to detect a real effect | ✅ |
| 10 | [A/B Testing & Experimental Design](10_AB_Testing_Experimental_Design/) | Designing and analyzing controlled experiments | ✅ |
| 11 | [Multiple Testing Correction](11_Multiple_Testing_Correction/) | Controlling false positives when running many hypothesis tests | ✅ |
| 12 | [Bootstrap & Permutation Tests](12_Bootstrap_Permutation_Tests/) | Resampling-based inference without distributional assumptions | ✅ |
| 13 | [Bayesian Inference Basics](13_Bayesian_Inference_Basics/) | Priors, posteriors, and credible intervals | ✅ |
| 14 | [Categorical Data Analysis](14_Categorical_Data_Analysis/) | Contingency tables, odds ratios, and Fisher's Exact Test | ✅ |
| 15 | [Survival Analysis Basics](15_Survival_Analysis_Basics/) | Time-to-event data, Kaplan-Meier curves, and the Cox model | 🚧 |
| 16 | [Maximum Likelihood Estimation](16_Maximum_Likelihood_Estimation/) | The optimization principle behind most parametric estimators | 🚧 |
| 17 | [Goodness-of-Fit & Normality Tests](17_Goodness_of_Fit_Normality_Tests/) | Shapiro-Wilk, Kolmogorov-Smirnov, and Q-Q plots | 🚧 |
| 18 | [ANCOVA & MANOVA](18_ANCOVA_MANOVA/) | Extending ANOVA with covariates and multiple outcomes | 🚧 |
| 19 | [Regression Inference & Diagnostics](19_Regression_Inference_Diagnostics/) | F-tests, coefficient significance, and residual diagnostics | 🚧 |
| 20 | [Causal Inference Basics](20_Causal_Inference_Basics/) | Confounding, causal graphs, and propensity score matching | 🚧 |

**20 topics planned** — built one at a time, each with a deep-dive `README.md` (full math derivation in LaTeX, pitfalls, self-test exercises) and a fully-executed Jupyter notebook (30+ code cells), following the exact standard established in [Foundation](https://github.com/mdnuruzzamanKALLOL/Statistical-Machine-Learning-Foundation) and [Classical ML](https://github.com/mdnuruzzamanKALLOL/Statistical-Machine-Learning-Classical-ML).

## 📁 Repository Structure

```
Statistical-Inference-Hypothesis-Testing/
├── README.md                          ← you are here
├── 01_Probability_Distributions_Sampling_Theory/
│   └── README.md
├── 02_Central_Limit_Theorem_Point_Estimation/
│   └── README.md
├── 03_Confidence_Intervals/
│   └── README.md
├── 04_One_Sample_Two_Sample_t_Tests/
│   └── README.md
├── 05_Chi_Square_Tests/
│   └── README.md
├── 06_ANOVA_One_way_Two_way/
│   └── README.md
├── 07_Non_Parametric_Tests/
│   └── README.md
├── 08_Correlation_Analysis/
│   └── README.md
├── 09_Power_Analysis_Sample_Size_Determination/
│   └── README.md
├── 10_AB_Testing_Experimental_Design/
│   └── README.md
├── 11_Multiple_Testing_Correction/
│   └── README.md
├── 12_Bootstrap_Permutation_Tests/
│   └── README.md
├── 13_Bayesian_Inference_Basics/
│   └── README.md
├── 14_Categorical_Data_Analysis/
│   └── README.md
├── 15_Survival_Analysis_Basics/
│   └── README.md
├── 16_Maximum_Likelihood_Estimation/
│   └── README.md
├── 17_Goodness_of_Fit_Normality_Tests/
│   └── README.md
├── 18_ANCOVA_MANOVA/
│   └── README.md
├── 19_Regression_Inference_Diagnostics/
│   └── README.md
└── 20_Causal_Inference_Basics/
    └── README.md
```

Every topic folder will be self-contained once built: read the `README.md` for the theory, open the `.ipynb` for the hands-on implementation.

## 🧰 Prerequisites

- Python 3.9+
- No machine learning library assumptions beyond `scipy.stats` and `statsmodels` — this repo is deliberately closer to a statistics course than an ML course.

```bash
pip install numpy pandas matplotlib seaborn scipy statsmodels scikit-learn jupyter
```

## 🚀 How to Use

**Just reading?** Every notebook (once built) renders directly on GitHub with full output — click any `.ipynb` link in the table above.

**Running it yourself:**

```bash
git clone https://github.com/mdnuruzzamanKALLOL/Statistical-Inference-Hypothesis-Testing.git
cd Statistical-Inference-Hypothesis-Testing
pip install numpy pandas matplotlib seaborn scipy statsmodels scikit-learn jupyter
jupyter notebook
```

## 📦 Datasets

Real-data topics draw from the central **[Datasets](https://github.com/mdnuruzzamanKALLOL/Datasets)** repo (162-entry catalog, 418 files) shared across this entire series.

## 🔭 The Full Series

| Repo | Covers | Status |
|---|---|---|
| [Foundation](https://github.com/mdnuruzzamanKALLOL/Statistical-Machine-Learning-Foundation) | Python/NumPy/Pandas, Visualization, Preprocessing, Feature Engineering, Math | ✅ Complete |
| [Classical ML](https://github.com/mdnuruzzamanKALLOL/Statistical-Machine-Learning-Classical-ML) | Regression, Classification, Ensembles, Unsupervised, Model Evaluation & Tuning (29 algorithms) | ✅ Complete |
| [Statistical Inference & Hypothesis Testing](https://github.com/mdnuruzzamanKALLOL/Statistical-Inference-Hypothesis-Testing) | Probability through causal inference (20 topics) | 🚧 In progress |
| [Time Series Analysis](https://github.com/mdnuruzzamanKALLOL/Time-Series-Analysis) | Decomposition through neural forecasting (20 topics) | 🚧 In progress |
| [Deep Learning Foundations](https://github.com/mdnuruzzamanKALLOL/Deep-Learning-Foundations) | Perceptron through GANs and interpretability (20 topics) | 🚧 In progress |
| [Datasets](https://github.com/mdnuruzzamanKALLOL/Datasets) | 418 datasets (162 cataloged + personal collection) backing the whole series | ✅ Complete |

## 📝 Self-Test Philosophy

Every topic README (once built) ends with a **Self-Test Exercises** section — deliberately *not* answered inline. The point is to predict the answer before running the corresponding notebook cell, not to read a solution.

## 💬 Feedback & Corrections

Found a bug, a math error, or a broken dataset link? Open an issue or a pull request — this series documents its own mistakes on purpose, so a caught error is a feature, not an embarrassment.

## 🔗 Related

- [Foundation →](https://github.com/mdnuruzzamanKALLOL/Statistical-Machine-Learning-Foundation)
- [Classical ML →](https://github.com/mdnuruzzamanKALLOL/Statistical-Machine-Learning-Classical-ML)
- [Datasets →](https://github.com/mdnuruzzamanKALLOL/Datasets)
- [Author's GitHub profile](https://github.com/mdnuruzzamanKALLOL)

<!-- page-views-badge -->
<div align="center" style="margin-top: 16px;">

![Page Views](https://visitor-badge.laobi.icu/badge?page_id=mdnuruzzamanKALLOL.StatisticalInferenceHypothesisTesting.root&left_color=%23FF6F00&right_color=%230e75b6&left_text=Page%20Views)

</div>
