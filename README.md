# Survival Analysis of Heart Failure Clinical Records

## Overview

This repository contains a coursework project in survival analysis based on the `heart_failure_clinical_records_dataset.csv` data set. The empirical problem is the modelling of time from clinical observation to death among patients with heart failure, with `time` measured in days and `DEATH_EVENT` used as the event indicator.

The final analytical sample contains 299 patients, including 96 observed deaths and 203 right-censored observations. Right censoring is therefore substantial: approximately 67.9% of patients did not experience the event during the observation window. This feature is central to the analysis, because it affects the precision of survival estimates, especially in the later part of follow-up.

The project compares nonparametric, parametric, Bayesian and semiparametric survival models. The emphasis is not only on fitting models, but also on checking their assumptions, interpreting coefficient scales correctly, and comparing the conclusions obtained from different modelling frameworks.

## Data and variables

The input file is stored in:

```text
data/heart_failure_clinical_records_dataset.csv
```

The data set contains 13 source variables and 299 observations. The main variables used in the analysis are:

* `time` — follow-up time in days;
* `DEATH_EVENT` — event indicator, where `1` denotes death and `0` denotes right-censoring;
* `age` — patient age;
* `anaemia`, `diabetes`, `high_blood_pressure`, `sex`, `smoking` — binary clinical or demographic variables;
* `ejection_fraction` — left ventricular ejection fraction;
* `serum_creatinine` — serum creatinine level;
* `serum_sodium` — serum sodium level;
* `platelets` — platelet count;
* `creatinine_phosphokinase` — CPK level.

The analysis also uses the transformation

$$
\mathrm{log_cpk} = \log(\mathrm{CPK}+1),
$$

to reduce the influence of the highly right-skewed CPK variable.

Selected descriptive statistics:

| Variable          |   Mean | Median |    Min |    Max |
| ----------------- | -----: | -----: | -----: | -----: |
| age               |  60.83 |  60.00 |  40.00 |  95.00 |
| ejection fraction |  38.08 |  38.00 |  14.00 |  80.00 |
| serum creatinine  |   1.39 |   1.10 |   0.50 |   9.40 |
| serum sodium      | 136.63 | 137.00 | 113.00 | 148.00 |
| follow-up time    | 130.26 | 115.00 |   4.00 | 285.00 |

For subgroup analysis, several clinically interpretable binary splits are used:

* `ejection_fraction < 35%`;
* `serum_creatinine > 1.5`;
* `serum_sodium < 135`;
* presence of hypertension, anaemia, diabetes, smoking, and sex category.

## Analytical workflow

The project is organised into three modelling notebooks and a final Beamer presentation.

### 1. Nonparametric survival analysis

The nonparametric notebook estimates survival and cumulative hazard without imposing a parametric distribution on event times.

The main outputs are:

* Kaplan-Meier survival curve for the full sample;
* Nelson-Aalen cumulative hazard curve;
* smoothed nonparametric hazard estimate;
* actuarial life table using 30-day intervals;
* Kaplan-Meier curves for clinically defined subgroups;
* log-rank and Wilcoxon tests for equality of survival curves.

The actuarial life table shows a clear decline in survival probability across the observation window. In the first 30-day interval, 31 deaths are observed and the estimated interval survival is approximately 0.896. By the final interval, the life-table survival estimate is approximately 0.559.

The subgroup tests show the strongest univariate differences for:

| Variable                | Log-rank statistic | p-value |
| ----------------------- | -----------------: | ------: |
| ejection fraction < 35% |              30.63 | < 0.001 |
| serum creatinine > 1.5  |              39.61 | < 0.001 |
| serum sodium < 135      |              15.65 | < 0.001 |
| high blood pressure     |               4.41 |   0.036 |

Anaemia is borderline in the univariate tests, while diabetes, sex and smoking do not produce meaningful separation of Kaplan-Meier curves in this sample.

These nonparametric results are treated as exploratory evidence. They identify variables worth testing in regression models, but they do not control for multiple patient characteristics simultaneously.

### 2. Parametric AFT models

The parametric part compares accelerated failure time models fitted under several distributional assumptions:

* Weibull AFT;
* log-normal AFT;
* log-logistic AFT;
* exponential model without covariates;
* reduced Weibull AFT model.

The AFT specification is written as

$$
\log(T_i)
=========

\beta_0 + x_i^\top \beta + \sigma \varepsilon_i.
$$

In an AFT model, the quantity

$$
\exp(\beta_j)
$$

is a time ratio, not a hazard ratio. Values above 1 lengthen the predicted survival time, while values below 1 shorten it, holding the other variables fixed.

The likelihood contribution for an observation with time (t_i), covariates (x_i), and event indicator (\delta_i) is

$$
L_i
===

f(t_i \mid x_i)^{\delta_i}
S(t_i \mid x_i)^{1-\delta_i}.
$$

This explicitly accounts for right-censored observations through the survival term (S(t_i \mid x_i)).

The full parametric model comparison is:

| Model                          |  logLik |     AIC |     BIC | Parameters |
| ------------------------------ | ------: | ------: | ------: | ---------: |
| Weibull AFT full               | -631.52 | 1287.04 | 1331.45 |         12 |
| Log-normal AFT                 | -634.36 | 1292.72 | 1337.12 |         12 |
| Log-logistic AFT               | -633.04 | 1290.09 | 1334.49 |         12 |
| Weibull AFT reduced            | -632.13 | 1280.27 | 1309.87 |          8 |
| Exponential without covariates | -672.54 | 1347.08 | 1350.78 |          1 |

The full Weibull AFT model performs best among the full AFT candidates by AIC. The reduced Weibull AFT model improves the information criteria further by removing weak predictors.

In the reduced Weibull AFT model, the main time-ratio effects are:

| Variable            |  Coef. | Time ratio | p-value |
| ------------------- | -----: | ---------: | ------: |
| age                 | -0.045 |      0.956 | < 0.001 |
| anaemia             | -0.405 |      0.667 |   0.055 |
| ejection fraction   |  0.048 |      1.049 | < 0.001 |
| high blood pressure | -0.496 |      0.609 |   0.021 |
| serum creatinine    | -0.310 |      0.733 | < 0.001 |
| serum sodium        |  0.044 |      1.045 |   0.068 |

The estimated Weibull shape parameter is close to one, approximately

$$
\widehat{\rho} \approx 0.96,
$$

which suggests a hazard close to constant, or only mildly decreasing, over time.

The probability plots show good marginal fit for all three candidate distributions. The log-normal plot has the highest graphical (R^2), but the regression-level information criteria favour Weibull AFT. This is an important modelling tension: graphical marginal fit and regression AIC do not have to select the same model.

### 3. Bayesian Weibull extension

The Bayesian part is included as a methodological extension rather than the main empirical model. It uses a simplified Weibull model with three key predictors:

* `age`;
* `ejection_fraction`;
* `serum_creatinine`.

The model is specified as

$$
T_i \sim \mathrm{Weibull}(\rho, \lambda_i),
\qquad
\log(\lambda_i)
===============

\beta_0 + x_i^\top \beta.
$$

The priors are

$$
\beta_j \sim N(0, 2.5^2),
\qquad
\log(\rho) \sim N(0,1).
$$

Right-censoring is handled by using the survival contribution (S(t_i \mid x_i)) for observations without an event.

The posterior summary supports the same qualitative direction as the classical AFT models:

| Parameter         |   Mean |   2.5% |  97.5% | (P(\beta>0)) |
| ----------------- | -----: | -----: | -----: | -----------: |
| age               | -0.592 | -0.917 | -0.379 |        0.000 |
| ejection fraction |  0.701 |  0.381 |  1.114 |        1.000 |
| serum creatinine  | -0.404 | -0.562 | -0.246 |        0.000 |
| log rho           | -0.124 | -0.326 |  0.057 |        0.106 |

The Bayesian model is deliberately narrower than the full classical models. Its role is to illustrate posterior uncertainty and probabilistic interpretation of effects, not to replace the full AFT or Cox analyses.

### 4. Cox proportional hazards model and diagnostics

The semiparametric part fits a Cox proportional hazards model:

$$
h(t \mid x)
===========

h_0(t)\exp(x^\top \beta).
$$

Here (h_0(t)) is left unspecified, while (\exp(\beta_j)) is interpreted as a hazard ratio for a one-unit increase in (x_j), holding the other covariates fixed.

The main Cox model estimates are:

| Variable            |  Coef. |    HR | p-value |
| ------------------- | -----: | ----: | ------: |
| age                 |  0.038 | 1.039 | < 0.001 |
| ejection fraction   | -0.039 | 0.962 | < 0.001 |
| serum creatinine    |  0.294 | 1.342 | < 0.001 |
| serum sodium        | -0.041 | 0.960 |   0.058 |
| anaemia             |  0.337 | 1.401 |   0.086 |
| diabetes            |  0.072 | 1.075 |   0.721 |
| high blood pressure |  0.401 | 1.493 |   0.045 |
| sex                 | -0.115 | 0.891 |   0.600 |
| smoking             |  0.085 | 1.088 |   0.705 |
| log(CPK+1)          |  0.048 | 1.050 |   0.588 |

The strongest and most stable predictors in the Cox model are age, ejection fraction and serum creatinine. Hypertension remains significant at the 5% level, while anaemia and sodium are borderline.

The proportional hazards assumption is tested using Schoenfeld-residual-based diagnostics. The test identifies a violation for `ejection_fraction`:

| Variable          | Chi-square | p-value |
| ----------------- | ---------: | ------: |
| ejection fraction |      4.227 |   0.040 |

This means that the Cox hazard ratio for ejection fraction should be interpreted cautiously as an average effect over time rather than as a strictly constant hazard ratio.

A time-interaction extension is fitted using an interaction between ejection fraction and a log-time transformation. The comparison is:

| Model                    |  logLik | Partial AIC | Parameters | C-index |
| ------------------------ | ------: | ----------: | ---------: | ------: |
| Cox baseline             | -471.28 |      962.57 |         10 |   0.734 |
| Cox with EF × log(1 + t) | -490.46 |     1002.92 |         11 |   0.732 |

The time-interaction model does not improve the overall fit according to partial AIC or concordance. The baseline Cox model is therefore retained as the main semiparametric model, with an explicit caveat about the time-varying effect of ejection fraction.

The Cox notebook also includes:

* coefficient plots on the log hazard-ratio scale;
* Schoenfeld residual diagnostics;
* delta-beta influence analysis;
* predicted survival curves for low-, medium- and high-risk patient profiles;
* predicted hazard profiles under the fitted Cox model.

## Main empirical conclusions

The different modelling approaches produce a coherent risk ranking.

The variables most consistently associated with worse survival are:

* higher age;
* lower ejection fraction;
* higher serum creatinine;
* hypertension;
* anaemia, with weaker statistical support;
* lower serum sodium, with borderline support.

The nonparametric analysis shows strong univariate separation for low ejection fraction, high creatinine and low sodium. The AFT models translate these differences into survival-time ratios. The Cox model confirms the same broad structure using hazard ratios, while also requiring diagnostic qualification because the proportional hazards assumption is not fully satisfied for ejection fraction.

The project should therefore be read as a methodological survival-analysis study, not as a clinical decision system. The sample is small, the number of observed deaths is limited, censoring is substantial, and no external validation is performed.

## Repository structure

```text
data/
    heart_failure_clinical_records_dataset.csv

notebooks/
    01_nonparametric_survival_analysis.ipynb
    02_parametric_and_bayesian_survival_models.ipynb
    03_cox_semiparametric_model_diagnostics.ipynb

presentation/
    main.tex
    survival_analysis_presentation.pdf

reports/
    interpretation_notes.pdf

tables/
    corrected LaTeX tables used in the presentation

figures/
    corrected model, survival, hazard and diagnostic figures

results/
    compact run summary metadata
```

## Recommended review order

1. `presentation/survival_analysis_presentation.pdf` — final presentation and high-level narrative.
2. `notebooks/01_nonparametric_survival_analysis.ipynb` — Kaplan-Meier, Nelson-Aalen, life table and subgroup tests.
3. `notebooks/02_parametric_and_bayesian_survival_models.ipynb` — AFT models, distributional comparison and Bayesian Weibull extension.
4. `notebooks/03_cox_semiparametric_model_diagnostics.ipynb` — Cox model, diagnostics, Schoenfeld residuals, influence analysis and risk profiles.
5. `reports/interpretation_notes.pdf` — supplementary interpretation notes.

## Reproducibility

Install the Python dependencies:

```bash
pip install -r requirements.txt
```

Run the notebooks from the repository root so that relative paths resolve correctly:

```text
data/
tables/
figures/
results/
```

The presentation can be rebuilt from the `presentation/` directory with a LaTeX installation that supports Beamer:

```bash
cd presentation
pdflatex main.tex
pdflatex main.tex
```

## Outputs

The repository contains the following main outputs:

* Kaplan-Meier survival curves for the full sample and clinical subgroups;
* Nelson-Aalen cumulative hazard estimates;
* actuarial life table for 30-day intervals;
* log-rank and Wilcoxon tests for subgroup comparisons;
* AFT model comparison tables;
* Weibull AFT full and reduced model estimates;
* probability plots for Weibull, log-normal and log-logistic distributions;
* predicted survival and hazard profiles from AFT models;
* Bayesian Weibull posterior summary and traceplots;
* Cox model coefficient table and log hazard-ratio plot;
* proportional hazards diagnostics based on Schoenfeld residuals;
* Cox model with time interaction for ejection fraction;
* delta-beta influence diagnostics;
* predicted Cox survival and hazard profiles for patient risk classes.

## Limitations

The main limitations are:

* only 299 observations are available;
* only 96 deaths are observed;
* 203 observations are right-censored;
* the analysis assumes independent censoring;
* later survival estimates are less precise because the risk set becomes smaller;
* parametric AFT conclusions depend on the selected event-time distribution;
* the Cox model shows a proportional hazards violation for `ejection_fraction`;
* the Bayesian model is simplified and should be treated as demonstrative;
* no external validation sample is used.

## References

* Cox, D. R. (1972). Regression models and life-tables. *Journal of the Royal Statistical Society: Series B*.
* Kaplan, E. L., & Meier, P. (1958). Nonparametric estimation from incomplete observations. *Journal of the American Statistical Association*.
* Nelson, W. (1972). Theory and applications of hazard plotting for censored failure data. *Technometrics*.
* Davidson-Pilon, C. et al. `lifelines` Python survival analysis library documentation.
