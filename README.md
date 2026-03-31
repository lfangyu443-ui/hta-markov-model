# HTA Markov Model: Telemonitoring for Heart Failure in Singapore

A cost-effectiveness analysis comparing home telemonitoring (TM) versus structured telephone support (STS) for heart failure patients aged 65+ following hospital discharge in Singapore.

## Background

Developed for SPH5412 (Economic Methods in Health Technology Assessment) at NUS Saw Swee Hock School of Public Health. Based on the Changi General Hospital (CGH) telemonitoring programme (Leng Chow et al. 2020), evaluated against STS as the active comparator reflecting current CGH clinical practice.

Heart failure affects 26 million people worldwide, with Singapore reporting a mean admission age of ~66 years and a 30-day readmission rate of ~23% (Senanayake et al. 2025). Despite clinical evidence from CGH showing TM superiority over STS, no formal Markov-based cost-effectiveness analysis had been conducted for the Singapore setting — the gap this model addresses.

## PICO

| Element | Description |
|---|---|
| **Population** | Adults ≥65 with HFrEF, recently discharged from a Singapore public hospital |
| **Intervention** | Home telemonitoring (TM) — CGH programme: weight, BP, HR via tablet + nurse tele-support for 12 months |
| **Comparator** | Structured telephone support (STS) — daily calls for patients at risk of 30-day readmission |
| **Outcomes** | Costs (SGD), QALYs, ICER |

## Model Structure

- **Type:** Cohort state-transition Markov model (hybrid: Jiang et al. 2020 + Caillon et al. 2022)
- **Health states:** NYHA Class I, II, III, IV, Dead (absorbing)
- **Cycle length:** 1 month
- **Time horizon:** 10 years (base case); 5 years (scenario)
- **Perspective:** Healthcare system (ACE Singapore)
- **Discount rate:** 3% per annum (Liang et al. 2018)
- **WTP threshold:** SGD 80,000/QALY (consistent with prior Singapore CEA)

### Treatment Effect Timeline

| Phase | Months | Description |
|---|---|---|
| Phase 1 | 1–12 | Programme cost active: TM SGD 150/month, STS SGD 60/month. HR 0.70 · OR 0.64 active |
| Phase 2 | 13–60 | No programme cost. HR 0.70 · OR 0.64 still active (Caillon et al. 2022) |
| Phase 3 | 61–120 | No programme cost. TM effect wanes to zero — both arms use STS matrix |

The TM benefit is modelled exclusively through:
- **Mortality:** HR 0.70 (95% CI 0.50–0.96) from TIM-HF2 trial (cited in Jiang et al. 2020)
- **Hospitalisation:** OR 0.64 (95% CI 0.39–0.95) from Kotb et al. (2015) NMA (cited in Jiang et al. 2020)

Both applied via logistic transformation. NYHA transition probabilities are identical across both arms.

## Key Parameters

| Parameter | Value | DSA Range | Source |
|---|---|---|---|
| Mortality NYHA I–IV (6-month) | 0.65%, 3.56%, 11.68%, 11.68% | — | Jiang et al. 2020 (trial n=2,737) |
| TM mortality effect (HR) | 0.70 | 0.50–0.96 (95% CI) | TIM-HF2 trial, cited in Jiang et al. 2020 |
| TM hospitalisation effect (OR) | 0.64 | 0.39–0.95 (95% CI) | Kotb et al. 2015 NMA, cited in Jiang et al. 2020 |
| Hospitalisation prob. NYHA I–IV | 4%, 8%, 14%, 22%/month | — | Assumption; Jiang et al. 2020 + Senanayake et al. 2025 |
| NYHA transition matrix | Monthly probabilities | — | Liang et al. 2018; King et al. 2016; CARE-HF |
| Starting cohort NYHA I–IV | 10%, 40%, 35%, 15% | — | Clinical assumption |
| Utility NYHA I/II | 0.771 | ±20% | Caillon et al. 2022 Table 2; SHIFT trial |
| Utility NYHA III/IV | 0.639 | ±20% | Caillon et al. 2022 Table 2; SHIFT trial |
| Hospitalisation disutility | −0.212/event | — | Caillon et al. 2022 Table 2 |
| TM programme cost | SGD 150/month (M1–12) | ±20% | Estimated; Leng Chow et al. 2020 |
| STS programme cost | SGD 60/month (M1–12) | ±20% | Estimated; Leng Chow et al. 2020 |
| HF hospitalisation cost | SGD 5,900/event | ±20% | Derived; Senanayake et al. 2025 Table 2 |
| Palliative care cost | SGD 3,500/death | — | Assumption |
| Outpatient + Rx NYHA I–IV | SGD 280–600/month | — | Assumption; Senanayake et al. 2025 + Jiang et al. 2020 |
| Treatment effect duration | 60 months | 36–120 months | Caillon et al. 2022 |
| Discount rate | 3% p.a. | 0–5% | Liang et al. 2018 |
| WTP threshold | SGD 80,000/QALY | — | Liang et al. 2018 |

## Results

### Base Case (10-Year Horizon)

| Strategy | Total Cost (SGD) | Total QALYs | Incr. Cost (SGD) | Incr. QALYs | ICER (SGD/QALY) |
|---|---|---|---|---|---|
| STS (Comparator) | 5,629 | 3.14 | — | — | — |
| Telemonitoring (TM) | 5,526 | 3.62 | −103 | 0.48 | **−214 (Dominant)** |

TM is **dominant** — cheaper and more effective than STS. Cost savings arise from averted hospitalisations (OR 0.64) and reduced mortality (HR 0.70) during years 1–5, more than offsetting the SGD 90/month incremental programme cost.

### Scenario: 5-Year Horizon

TM remains dominant: **ICER −SGD 1,968/QALY** (full treatment effect window captured).

### Deterministic Sensitivity Analysis (DSA)

One-way DSA varied: mortality HR (0.50–0.96), hospitalisation OR (0.39–0.95), effect duration (36–120 months), costs (±20%), discount rate (0–5%), and utilities (±20%).

- Most influential parameter: **mortality HR**, followed by hospitalisation OR and effect duration
- TM remained **dominant or cost-effective across all DSA scenarios**
- Worst case (HR=0.96): ICER ~SGD 514/QALY — still well below SGD 80,000/QALY

### Probabilistic Sensitivity Analysis (PSA, n=1,000 Monte Carlo)

| Metric | Value |
|---|---|
| Median ICER | −SGD 232/QALY |
| Mean ICER | −SGD 2,057/QALY (skewed; median preferred) |
| 95% Credible Interval | −SGD 8,437 to +SGD 2,168 |
| Probability CE at SGD 80,000/QALY | **99.8%** |

Distributions: HR/OR — log-normal (SDs from 95% CIs); utilities — beta (truncated at disutility floor 0.212); costs — gamma.

## Analyses

- Base-case ICER (SGD per QALY)
- One-way deterministic sensitivity analysis (tornado diagram — two-panel)
- Probabilistic sensitivity analysis (1,000 Monte Carlo iterations)
- Cost-effectiveness plane and CEAC
- Scenario analyses (5-year horizon)

## Key Assumptions

- TM benefit applied via HR 0.70 (mortality) and OR 0.64 (hospitalisation) only; NYHA transitions identical across arms
- Programme delivery cost applies months 1–12 only; clinical benefit retained months 1–60, then wanes to zero (Caillon et al. 2022)
- Palliative care cost (SGD 3,500) applied as a one-off event-based cost at death transition
- After month 60, both arms revert to identical matrices (Markov memoryless property)
- Non-adherence not explicitly modelled (TM no longer preferred if adherence falls below ~20.9%, per Jiang et al. 2020)
- Utility values truncated at hospitalisation disutility (0.212) in PSA to prevent negative QALYs; median ICER reported as primary PSA statistic

## Limitations

- **Treatment effect duration:** 60-month waning assumption from French cohort (Caillon et al. 2022) — Singapore-specific data unavailable
- **Adherence not modelled:** Real-world adherence among older Singaporeans (mean age 66, varied digital literacy) unknown
- **Programme costs assumed:** No published unit costs for TM/STS in Singapore; SGD 150/SGD 60 per month estimated from CGH programme descriptions, tested with ±20% DSA
- **Utility values:** Derived from SHIFT trial (European population); Singapore-specific EQ-5D data by NYHA class and ethnicity would improve precision
- **Starting NYHA distribution assumed:** Senanayake et al. (2025) does not report NYHA distribution at discharge

## Policy Implications

Given TM dominance and 99.8% PSA probability of cost-effectiveness, there is a strong economic rationale for policy intervention. MediFund and MediShield Life do not currently cover remote monitoring or nurse-led digital care. Recommended approach for MOH/ACE:

1. Structured pilot of the CGH TM model across ≥2 public hospital clusters with outcomes evaluation
2. Use Singapore-collected long-term effectiveness data to submit for ACE subsidy consideration
3. If continued evidence of effectiveness, extend coverage under the CDMP umbrella as means-tested post-discharge bundles

## How to Run

1. Open `R/markov_hf_telemedicine.R` in RStudio
2. Install dependencies: `install.packages(c("ggplot2", "gridExtra"))`
3. Run the full script — results print to console, 4 figures saved as PNG:
   - `trace_plot.png` — Markov cohort trace (TM vs STS, 10-year)
   - `tornado_plot.png` — One-way DSA tornado diagram (two-panel)
   - `ce_plane.png` — PSA cost-effectiveness plane
   - `ceac_plot.png` — Cost-effectiveness acceptability curve

## Tools

- R 4.x
- Packages: `ggplot2`, `gridExtra` (base R used for all reshaping — no `reshape2` or `MASS` required)
- Random seed: `set.seed(12345)` for PSA reproducibility

## References

- Caillon et al. (2022). BMC Cardiovasc Disord. doi:10.1186/s12872-022-02878-1
- Jiang et al. (2020). JMIR Mhealth Uhealth. doi:10.2196/17846
- Kotb et al. (2015). PLoS One. doi:10.1371/journal.pone.0118681
- Leng Chow et al. (2020). J Telemed Telecare. doi:10.1177/1357633X18825164
- Liang et al. (2018). J Med Econ. doi:10.1080/13696998.2017.1387119
- Senanayake et al. (2025). Value Health Reg Issues. doi:10.1016/j.vhri.2024.101037
- ACE Singapore. Guidance on economic evaluation for HTA. MOH; 2022.

## Transparency

An AI tool was used to assist in structuring, drafting, and fixing R code. All analytical decisions, parameter choices, and interpretations were reviewed and verified by the author. References were independently checked for accuracy.
