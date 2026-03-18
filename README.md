# HTA Markov Model: Telemonitoring for Heart Failure in Singapore

A cost-effectiveness analysis comparing home telemonitoring versus structured telephone support for heart failure patients aged 65+ following hospital discharge in Singapore.

## Background

Developed for SPH5412 (Economic Methods in Health Technology Assessment) at NUS Saw Swee Hock School of Public Health. Based on the Changi General Hospital (CGH) telemonitoring programme.

## Model Structure

- **Type**: Cohort-based state-transition Markov model
- **Health states**: NYHA Class I, II, III, IV, Death
- **Cycle length**: 1 month
- **Time horizon**: Lifetime (35 years); scenario analyses at 5 and 10 years
- **Perspective**: Healthcare system (ACE Singapore)
- **Discount rate**: 3%

## Key Parameters

| Parameter | Value | Source |
|---|---|---|
| Transition probabilities | NYHA-stratified monthly | Liang et al. (2017) |
| Hospitalisation costs | SGD 8,500/admission | Senanayake et al. (2025) |
| TM effectiveness (OR) | 0.71 hospitalisation | Kotb et al. (2015) |
| Utility values | 0.432-0.848 by NYHA | Liang et al. (2017) |
| TM programme cost | SGD 380/month | Leng Chow et al. (2020) |

## Analyses

- Base-case ICER (SGD per QALY)
- One-way deterministic sensitivity analysis (tornado diagram)
- Probabilistic sensitivity analysis (1,000 Monte Carlo iterations)
- Cost-effectiveness plane and CEAC
- Scenario analyses (5-year and 10-year horizons)

## How to Run

1. Open R/guardian_hf_markov.R in RStudio
2. Install dependencies: install.packages(c("dampack", "ggplot2", "dplyr", "scales"))
3. Run the full script - results print to console, 4 figures saved as PNG

## Tools

- R 4.x
- Packages: dampack, ggplot2, dplyr, scales
