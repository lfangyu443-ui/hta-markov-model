################################################################################
# SPH5412 Individual Assignment
# Cost-Effectiveness Analysis of Telemonitoring vs Structured Telephone Support
# for Heart Failure Patients in Singapore
# Markov Model Implementation (Hybrid: Jiang et al. 2020 + Caillon et al. 2022)
################################################################################

# ---- 0. Load Libraries -------------------------------------------------------
# Auto-install ggplot2 if not present (only external package needed)
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)
# reshape2 and MASS are NOT required — all reshaping done in base R

set.seed(12345)

setwd("~/Desktop/HTA")
cat("Saving outputs to:", getwd(), "\n")

# ---- 1. MODEL PARAMETERS -----------------------------------------------------

## 1a. Time horizon & cycle
n_cycles   <- 12 * 10          # 10-year horizon, 1-month cycles (120 cycles)
cycle_len  <- 1/12             # Cycle length in years
discount_r <- 0.03             # Annual discount rate (ACE Singapore recommendation)

## 1b. Starting cohort distribution across NYHA states (post-discharge)
# Clinical assumption for Singapore post-discharge HF cohort
# (Senanayake 2025 does not report NYHA distribution; values based on clinical judgement)
# NYHA I=10%, II=40%, III=35%, IV=15% at discharge
start_dist <- c(NYHA_I = 0.10, NYHA_II = 0.40, NYHA_III = 0.35, NYHA_IV = 0.15, Dead = 0.00)
n_states   <- length(start_dist)
state_names <- names(start_dist)

## 1c. Monthly transition probabilities — Standard of Care (STS)
# Source: control group mortality from clinical trial cited in Jiang et al. (2020)
# NYHA transition matrix adapted from Liang et al. (2018), derived from CARE-HF trial via King et al. (2016)
# All-cause monthly mortality by NYHA class
# Source: Jiang et al. (2020) Table 1 — 6-month mortality from control group
#         of clinical trial of 2,737 HF patients
# 6-month mortality: NYHA I=0.65%, II=3.56%, III=11.68%, IV=11.68%
# Converted to monthly using: p_monthly = 1 - (1 - p_6month)^(1/6)
mort_6month_STS <- c(0.0065, 0.0356, 0.1168, 0.1168)
mort_monthly_STS <- 1 - (1 - mort_6month_STS)^(1/6)

# NYHA transition matrix for STS (rows = from, cols = to; last state = Dead)
# Adapted from Liang et al. (2018) Singapore ACE model (from CARE-HF via King et al. 2016)
# Diagonal = stay; off-diagonal = worsen/improve; last col = death
tp_STS <- matrix(0, nrow = 5, ncol = 5,
                  dimnames = list(state_names, state_names))

# Death probabilities placed in last column
tp_STS[1, 5] <- mort_monthly_STS[1]   # NYHA I -> Dead
tp_STS[2, 5] <- mort_monthly_STS[2]   # NYHA II -> Dead
tp_STS[3, 5] <- mort_monthly_STS[3]   # NYHA III -> Dead
tp_STS[4, 5] <- mort_monthly_STS[4]   # NYHA IV -> Dead
tp_STS[5, 5] <- 1.00                   # Dead -> Dead (absorbing)

# NYHA transitions (non-death portion); adapted from Liang et al. (2017)
# Probabilities for alive patients transitioning between NYHA classes
nyha_alive_STS <- matrix(
  c(0.70, 0.20, 0.08, 0.02,   # From NYHA I
    0.15, 0.60, 0.20, 0.05,   # From NYHA II
    0.05, 0.20, 0.55, 0.20,   # From NYHA III
    0.02, 0.05, 0.20, 0.73),  # From NYHA IV
  nrow = 4, byrow = TRUE
)

# Scale alive transitions by (1 - monthly death prob)
for (i in 1:4) {
  tp_STS[i, 1:4] <- nyha_alive_STS[i, ] * (1 - mort_monthly_STS[i])
}

# Verify rows sum to 1
stopifnot(all(abs(rowSums(tp_STS) - 1) < 1e-8))

## 1d. Telemonitoring (TM) effect
# Mortality: HR = 0.70 (95% CI 0.50-0.96) from TIM-HF2 trial
#            cited in Jiang et al. (2020), p.3 — applied as log-odds ratio
# Hospitalisation: OR = 0.64 (95% CI 0.39-0.95) from Kotb et al. (2015) NMA
#            cited in Jiang et al. (2020), p.3 — physiologic telemonitoring vs UC
OR_mort_TM   <- 0.70   # HR from TIM-HF2 (Jiang et al. 2020, citing Koehler et al. 2018)
OR_hosp_TM   <- 0.64   # OR from Kotb et al. (2015) NMA (cited in Jiang et al. 2020)

# TM mortality probabilities derived by applying OR to STS log-odds
logit   <- function(p) log(p / (1 - p))
inv_logit <- function(x) exp(x) / (1 + exp(x))

mort_monthly_TM <- inv_logit(logit(mort_monthly_STS) + log(OR_mort_TM))

# Rebuild TM transition matrix
tp_TM <- tp_STS  # Start from STS matrix

# Update death column and rescale NYHA transitions
for (i in 1:4) {
  tp_TM[i, 5] <- mort_monthly_TM[i]
  tp_TM[i, 1:4] <- nyha_alive_STS[i, ] * (1 - mort_monthly_TM[i])
}

stopifnot(all(abs(rowSums(tp_TM) - 1) < 1e-8))

## 1e. Monthly hospitalisation probabilities
# Base hospitalisation probabilities derived from Jiang et al. (2020):
#   NYHA I: 3-year rate 14.3% (DIG trial) → 6-month = 2.36% → monthly ≈ 0.004
#   NYHA II–IV: multiplied by DIG trial hazard ratios (1.16, 2.27, 3.71)
# Converted from 6-month Jiang cycle to monthly, then calibrated to
# Senanayake et al. (2025, citing Khera et al. 2019) 30-day readmission = 23% for Singapore
# NYHA stratification reflects DIG trial relative risks (Jiang et al. 2020, p.3)
hosp_monthly_STS <- c(0.04, 0.08, 0.14, 0.22)  # NYHA I–IV (assumption, see report)
hosp_monthly_TM  <- hosp_monthly_STS * OR_hosp_TM

## 1f. Utility values (EQ-5D) by NYHA class
# Source: Caillon et al. (2022) Table 2 — derived from SHIFT trial EQ-5D data
# (Griffiths et al. 2017), stratified by age and NYHA class.
# Weighted average across age groups using Caillon cohort split: 56% <70, 44% >=70
# NYHA I/II:   0.56*0.788 + 0.44*0.749 = 0.771
# NYHA III/IV: 0.56*0.669 + 0.44*0.603 = 0.639
# NYHA I and II use the NYHA I/II combined estimate (Caillon Table 2)
# NYHA III and IV use the NYHA III/IV combined estimate (Caillon Table 2)
# Note: Liang et al. (2018) reports similar values from the PARADIGM-HF trial
utils <- c(
  NYHA_I   = round(0.56*0.788 + 0.44*0.749, 3),  # = 0.771 (Caillon et al. 2022, Table 2)
  NYHA_II  = round(0.56*0.788 + 0.44*0.749, 3),  # = 0.771 (same NYHA I/II class)
  NYHA_III = round(0.56*0.669 + 0.44*0.603, 3),  # = 0.639 (Caillon et al. 2022, Table 2)
  NYHA_IV  = round(0.56*0.669 + 0.44*0.603, 3),  # = 0.639 (same NYHA III/IV class)
  Dead     = 0.00
)

# Hospitalisation disutility: -0.212 per hospitalisation event
# Source: Caillon et al. (2022) Table 2 — derived from SHIFT trial EQ-5D regression
# Applied as a one-off deduction in the cycle a patient is hospitalised
disutil_hosp <- 0.212  # Caillon et al. (2022), Table 2

## 1g. Costs (SGD, 2023 prices)
# Programme costs: Estimated from CGH programme descriptions (Leng Chow et al. 2020)
# The paper reports total cost saving of SGD 2,774/patient but no unit monthly breakdown.
# SGD 150/month (TM) and SGD 60/month (STS) are ESTIMATES based on programme scope.

# Monthly programme costs (intervention period = 12 months)
cost_TM_monthly  <- 150   # SGD/month — estimated (Leng Chow et al. 2020, programme description)
cost_STS_monthly <- 60    # SGD/month — estimated (Leng Chow et al. 2020, programme description)

# Hospitalisation cost per HF admission
# Senanayake et al. (2025) Table 2: HF-related per-admission cost
#   Year 1 mean = SGD 3,157; Years 2-4 mean range = SGD 5,698-6,566
#   SGD 5,900 used as approximate midpoint of years 2-4 (ongoing HF admissions)
# NOTE: SGD 5,900 ≠ a single reported figure; it is derived from Table 2
cost_hosp <- 5900         # SGD per HF admission (derived from Senanayake et al. 2025, Table 2)

# Palliative care cost per death — assumption, not directly reported in any source
# Senanayake et al. (2025) Table 4: final-year all-cause cost = SGD 29,479
# End-of-life care estimate: SGD 3,500 (assumption based on MOH palliative care tariffs)
cost_palliative <- 3500   # SGD — assumption (MOH palliative care context)

# Background outpatient/medication cost per month by NYHA class
# Senanayake et al. (2025) Table 3: overall annual outpatient cost SGD 3,032-5,034/year
#   = SGD 253-420/month overall (no NYHA stratification available in paper)
# NYHA-stratified values are ASSUMPTIONS informed by:
#   - Senanayake et al. (2025) overall range (SGD 253-420/month)
#   - Jiang et al. (2020) relative cost ratios across NYHA classes (US data)
# Values below are clinical assumptions; not directly from any single source
cost_outpatient <- c(NYHA_I = 280, NYHA_II = 350, NYHA_III = 450, NYHA_IV = 600)

# ---- 2. MARKOV MODEL SIMULATION ----------------------------------------------

run_markov <- function(tp_matrix, tp_matrix_base, hosp_probs, hosp_probs_base,
                       prog_cost_monthly,
                       n_cycles, start_dist, utils, cost_outpatient,
                       cost_hosp, cost_palliative, disutil_hosp,
                       discount_r, cycle_len,
                       programme_duration = 12,
                       effect_duration    = 60) {
  # effect_duration: months after which TM benefit wanes to zero
  # Based on Caillon et al. (2022): SCAD effectiveness assumed constant for
  # first 60 months, then no further effect beyond that period (p.4-5).
  # After effect_duration, the STS (base) transition matrix and hospitalisation
  # probabilities are used for the TM arm — i.e. both arms become identical.
  # For the STS arm, tp_matrix == tp_matrix_base so this has no effect.

  n_states <- length(start_dist)
  trace    <- matrix(0, nrow = n_cycles + 1, ncol = n_states,
                     dimnames = list(0:n_cycles, names(start_dist)))
  trace[1, ] <- start_dist

  total_cost  <- 0
  total_qaly  <- 0

  for (t in 1:n_cycles) {
    # Discount factor
    disc_t <- (1 / (1 + discount_r))^((t - 1) * cycle_len)

    # Current state distribution
    curr       <- trace[t, ]
    alive_prop <- curr[1:4]

    # Select active transition matrix and hospitalisation probabilities:
    # Within effect_duration → use TM matrix; after → revert to STS matrix
    active_tp   <- if (t <= effect_duration) tp_matrix   else tp_matrix_base
    active_hosp <- if (t <= effect_duration) hosp_probs  else hosp_probs_base

    # Hospitalisation events this cycle
    hosp_events <- alive_prop * active_hosp

    # Programme cost (applies only during programme_duration months)
    prog_cost <- if (t <= programme_duration) prog_cost_monthly else 0

    # Outpatient cost
    out_cost  <- sum(alive_prop * cost_outpatient)

    # Hospitalisation cost
    hosp_cost <- sum(hosp_events) * cost_hosp

    # Palliative cost — applied at transition to Dead this cycle
    new_dead   <- sum(alive_prop * active_tp[1:4, 5])
    death_cost <- new_dead * cost_palliative

    # Total discounted cost this cycle
    cycle_cost <- (prog_cost + out_cost + hosp_cost + death_cost) * disc_t * cycle_len

    # QALYs this cycle
    # Net utility per state = state utility - hospitalisation disutility if event occurs
    # Truncated at zero to prevent negative QALYs (per professor's requirement)
    # This can occur in PSA when sampled utility is very low and disutility is large
    net_utils <- pmax(utils[1:4] - disutil_hosp * active_hosp, 0)
    cycle_qaly <- sum(alive_prop * net_utils) * disc_t * cycle_len

    total_cost  <- total_cost  + cycle_cost
    total_qaly  <- total_qaly  + cycle_qaly

    # State transition
    trace[t + 1, ] <- curr %*% active_tp
  }

  return(list(trace = trace, cost = total_cost, qaly = total_qaly))
}

# ---- 3. BASE CASE RESULTS ----------------------------------------------------

cat("===== Running Base Case =====\n")

# STS arm: tp_matrix == tp_matrix_base (no TM effect), effect_duration irrelevant
result_STS <- run_markov(
  tp_matrix = tp_STS, tp_matrix_base = tp_STS,
  hosp_probs = hosp_monthly_STS, hosp_probs_base = hosp_monthly_STS,
  prog_cost_monthly = cost_STS_monthly,
  n_cycles = n_cycles, start_dist = start_dist, utils = utils,
  cost_outpatient = cost_outpatient, cost_hosp = cost_hosp,
  cost_palliative = cost_palliative, disutil_hosp = disutil_hosp,
  discount_r = discount_r, cycle_len = cycle_len
)

# TM arm: TM effect (HR 0.70, OR 0.64) applied for months 1-60, then wanes to STS
# Based on Caillon et al. (2022): SCAD effectiveness constant for 60 months,
# no further benefit assumed beyond that period (conservative assumption)
result_TM <- run_markov(
  tp_matrix = tp_TM, tp_matrix_base = tp_STS,
  hosp_probs = hosp_monthly_TM, hosp_probs_base = hosp_monthly_STS,
  prog_cost_monthly = cost_TM_monthly,
  n_cycles = n_cycles, start_dist = start_dist, utils = utils,
  cost_outpatient = cost_outpatient, cost_hosp = cost_hosp,
  cost_palliative = cost_palliative, disutil_hosp = disutil_hosp,
  discount_r = discount_r, cycle_len = cycle_len,
  effect_duration = 60   # TM benefit wanes after 60 months (Caillon et al. 2022)
)

incr_cost <- result_TM$cost - result_STS$cost
incr_qaly <- result_TM$qaly - result_STS$qaly
ICER      <- incr_cost / incr_qaly

cat("\n--- Base Case Summary ---\n")
cat(sprintf("  STS  - Total Cost: SGD %.0f | QALYs: %.4f\n", result_STS$cost, result_STS$qaly))
cat(sprintf("  TM   - Total Cost: SGD %.0f | QALYs: %.4f\n", result_TM$cost, result_TM$qaly))
cat(sprintf("  Incremental Cost: SGD %.0f\n", incr_cost))
cat(sprintf("  Incremental QALYs: %.4f\n", incr_qaly))
cat(sprintf("  ICER: SGD %.0f per QALY\n", ICER))


wtp <- 80000
cat(sprintf("\n  WTP Threshold: SGD %d/QALY\n", wtp))
cat(sprintf("  Cost-effective at WTP? %s\n", ifelse(ICER < wtp, "YES", "NO")))

# ---- 4. MARKOV TRACE PLOT ----------------------------------------------------

trace_STS_df <- as.data.frame(result_STS$trace)
trace_STS_df$Cycle <- 0:n_cycles
trace_STS_df$Strategy <- "STS"
# Base R reshape: wide -> long (no reshape2 needed)
trace_STS_long <- reshape(trace_STS_df,
                           varying   = state_names,
                           v.names   = "Proportion",
                           timevar   = "State",
                           times     = state_names,
                           direction = "long")

trace_TM_df <- as.data.frame(result_TM$trace)
trace_TM_df$Cycle <- 0:n_cycles
trace_TM_df$Strategy <- "TM"
trace_TM_long <- reshape(trace_TM_df,
                          varying   = state_names,
                          v.names   = "Proportion",
                          timevar   = "State",
                          times     = state_names,
                          direction = "long")

trace_all <- rbind(trace_STS_long, trace_TM_long)

trace_plot_data <- trace_all[trace_all$State != "Dead", ]

state_labels <- c(NYHA_I = "NYHA I", NYHA_II = "NYHA II",
                  NYHA_III = "NYHA III", NYHA_IV = "NYHA IV")
trace_plot_data$State_label    <- state_labels[trace_plot_data$State]
trace_plot_data$Strategy_label <- ifelse(trace_plot_data$Strategy == "TM",
                                         "Telemonitoring (TM)", "STS (comparator)")

p_trace <- ggplot(trace_plot_data,
                  aes(x     = Cycle / 12,
                      y     = Proportion,
                      colour = State_label,
                      linetype = Strategy_label,
                      linewidth = Strategy_label)) +
  geom_line() +
  scale_colour_manual(
    values = c("NYHA I" = "#2166ac", "NYHA II" = "#4dac26",
               "NYHA III" = "#d6604d", "NYHA IV" = "#f4a582"),
    name = "Health state") +
  scale_linetype_manual(
    values = c("Telemonitoring (TM)" = "solid", "STS (comparator)" = "dashed"),
    name = "Strategy") +
  scale_linewidth_manual(
    values = c("Telemonitoring (TM)" = 1.0, "STS (comparator)" = 0.7),
    name = "Strategy") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  scale_x_continuous(breaks = 0:10) +
  geom_vline(xintercept = 1,  linetype = "dotted", colour = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = 5,  linetype = "dotted", colour = "grey50", linewidth = 0.5) +
  annotate("text", x = 1.1, y = 0.38, label = "Programme\nends (yr 1)",
           size = 3, hjust = 0, colour = "grey40") +
  annotate("text", x = 5.1, y = 0.38, label = "Effect\nwanes (yr 5)",
           size = 3, hjust = 0, colour = "grey40") +
  labs(title    = "Markov cohort trace: telemonitoring vs STS",
       subtitle = "10-year horizon · 1-month cycles · Singapore HF cohort",
       x = "Time (years)", y = "Proportion of cohort") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.box        = "vertical",
        panel.grid.minor  = element_blank())

ggsave("trace_plot.png", p_trace, width = 9, height = 6, dpi = 150)
print(p_trace)
cat("\nTrace plot saved.\n")

# ---- 5. ONE-WAY DETERMINISTIC SENSITIVITY ANALYSIS (DSA) --------------------

cat("\n===== Running DSA =====\n")

# Parameters and ranges for DSA
# OR ranges from 95% CIs reported in source papers:
#   OR_mort_TM: HR 0.70 (95% CI 0.50-0.96) — TIM-HF2 (Jiang et al. 2020)
#   OR_hosp_TM: OR 0.64 (95% CI 0.39-0.95) — Kotb et al. 2015 (Jiang et al. 2020)
#   Utilities: ±20% applied (Caillon et al. 2022 DSA convention, Table 4)
#   Costs: ±20% applied (Caillon et al. 2022 DSA convention, Table 4)
#   Discount: 0-5% range (ACE Singapore guideline sensitivity)
dsa_params <- list(
  OR_mort_TM       = c(base = 0.70, low = 0.50, high = 0.96),  # 95% CI, TIM-HF2
  OR_hosp_TM       = c(base = 0.64, low = 0.39, high = 0.95),  # 95% CI, Kotb 2015
  effect_duration  = c(base = 60,   low = 36,   high = 120),   # months; Caillon: 60; range: 3-10 yrs
  cost_hosp        = c(base = 5900, low = 4720, high = 7080),   # ±20%
  cost_TM_monthly  = c(base = 150,  low = 120,  high = 180),    # ±20%
  cost_STS_monthly = c(base = 60,   low = 48,   high = 72),     # ±20%
  discount_r       = c(base = 0.03, low = 0.00, high = 0.05),   # ACE guideline range
  util_NYHA_II     = c(base = 0.771,low = 0.617,high = 0.925),  # ±20%, Caillon Table 2
  util_NYHA_III    = c(base = 0.639,low = 0.511,high = 0.767)   # ±20%, Caillon Table 2
)

run_dsa_scenario <- function(param_name, param_val) {
  # Rebuild parameters with the modified value
  or_m  <- if (param_name == "OR_mort_TM")   param_val else OR_mort_TM
  or_h  <- if (param_name == "OR_hosp_TM")   param_val else OR_hosp_TM
  ch    <- if (param_name == "cost_hosp")     param_val else cost_hosp
  ctm   <- if (param_name == "cost_TM_monthly")  param_val else cost_TM_monthly
  csts  <- if (param_name == "cost_STS_monthly") param_val else cost_STS_monthly
  dr    <- if (param_name == "discount_r")    param_val else discount_r
  ed    <- if (param_name == "effect_duration") param_val else 60

  u_loc <- utils
  if (param_name == "util_NYHA_II")  u_loc["NYHA_II"]  <- param_val
  if (param_name == "util_NYHA_III") u_loc["NYHA_III"] <- param_val

  # Rebuild TM matrix with new OR_mort
  mort_m_TM <- inv_logit(logit(mort_monthly_STS) + log(or_m))
  tp_TM_sc  <- tp_STS
  for (i in 1:4) {
    tp_TM_sc[i, 5]   <- mort_m_TM[i]
    tp_TM_sc[i, 1:4] <- nyha_alive_STS[i, ] * (1 - mort_m_TM[i])
  }
  hosp_TM_sc <- hosp_monthly_STS * or_h

  r_sts <- run_markov(tp_STS, tp_STS, hosp_monthly_STS, hosp_monthly_STS,
                      csts, n_cycles,
                      start_dist, u_loc, cost_outpatient, ch, cost_palliative,
                      disutil_hosp, dr, cycle_len)
  r_tm  <- run_markov(tp_TM_sc, tp_STS, hosp_TM_sc, hosp_monthly_STS,
                      ctm, n_cycles,
                      start_dist, u_loc, cost_outpatient, ch, cost_palliative,
                      disutil_hosp, dr, cycle_len, effect_duration = ed)

  icer_sc <- (r_tm$cost - r_sts$cost) / (r_tm$qaly - r_sts$qaly)
  return(icer_sc)
}

dsa_results <- data.frame(
  Parameter = character(), Low_ICER = numeric(), High_ICER = numeric(),
  stringsAsFactors = FALSE
)

for (p in names(dsa_params)) {
  icer_low  <- run_dsa_scenario(p, dsa_params[[p]]["low"])
  icer_high <- run_dsa_scenario(p, dsa_params[[p]]["high"])
  dsa_results <- rbind(dsa_results, data.frame(
    Parameter = p, Low_ICER = icer_low, High_ICER = icer_high
  ))
}

dsa_results$Range <- abs(dsa_results$High_ICER - dsa_results$Low_ICER)
dsa_results <- dsa_results[order(dsa_results$Range, decreasing = TRUE), ]

cat("\nDSA Results (ICER at low/high parameter values):\n")
print(dsa_results)

# Tornado plot — two panels because scale difference is too large for one axis
# Top panel: clinical effect parameters (large range)
# Bottom panel: cost/utility/discount parameters (small range)
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
library(gridExtra)

param_labels <- c(
  OR_mort_TM       = "Mortality HR (TM): 0.50 to 0.96",
  OR_hosp_TM       = "Hospitalisation OR (TM): 0.39 to 0.95",
  effect_duration  = "Effect duration: 36 to 120 months",
  cost_hosp        = "Hospitalisation cost: +/-20%",
  discount_r       = "Discount rate: 0 to 5%",
  cost_TM_monthly  = "TM programme cost: +/-20%",
  cost_STS_monthly = "STS programme cost: +/-20%",
  util_NYHA_II     = "Utility NYHA I/II: +/-20%",
  util_NYHA_III    = "Utility NYHA III/IV: +/-20%"
)
dsa_results$Label <- param_labels[dsa_results$Parameter]
dsa_results <- dsa_results[order(dsa_results$Range), ]
dsa_results$Label <- factor(dsa_results$Label, levels = dsa_results$Label)

# Split by range size
dsa_top <- dsa_results[dsa_results$Range > 300, ]
dsa_bot <- dsa_results[dsa_results$Range <= 300, ]

make_panel <- function(df, subtitle, show_base_icer = TRUE) {
  lo <- min(df$Low_ICER, df$High_ICER, na.rm = TRUE)
  hi <- max(df$Low_ICER, df$High_ICER, na.rm = TRUE)
  pad <- (hi - lo) * 0.15
  x_lo <- lo - pad
  x_hi <- hi + pad
  p <- ggplot(df, aes(y = Label,
                 xmin = pmin(Low_ICER, High_ICER),
                 xmax = pmax(Low_ICER, High_ICER))) +
    geom_errorbarh(height = 0.35, colour = "steelblue", linewidth = 5, alpha = 0.75) +
    scale_x_continuous(limits = c(x_lo, x_hi), labels = scales::comma) +
    labs(subtitle = subtitle, x = "ICER (SGD per QALY)", y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor  = element_blank(),
          axis.text.y       = element_text(size = 10),
          plot.subtitle     = element_text(size = 10, colour = "grey40"))
  # Only draw zero line if within range
  if (0 >= x_lo && 0 <= x_hi) {
    p <- p + geom_vline(xintercept = 0, linetype = "solid",
                        colour = "grey40", linewidth = 0.5)
  }
  # Only draw base case ICER line if within range
  if (show_base_icer && ICER >= x_lo && ICER <= x_hi) {
    p <- p + geom_vline(xintercept = ICER, linetype = "dashed",
                        colour = "red", linewidth = 0.8)
  }
  p
}

p_top <- make_panel(dsa_top,
  "Clinical effect parameters (large range — note different x-axis scale)")
p_bot <- make_panel(dsa_bot,
  "Cost, utility, and discount parameters (small range — zoomed x-axis)")

p_combined <- arrangeGrob(
  p_top, p_bot, nrow = 2,
  top = grid::textGrob(
    "Tornado diagram — one-way sensitivity analysis (red dashed = base case ICER: SGD -214/QALY)",
    gp = grid::gpar(fontsize = 12, fontface = "bold"))
)

ggsave("tornado_plot.png", p_combined, width = 10, height = 9, dpi = 150)
grid::grid.newpage(); grid::grid.draw(p_combined)
cat("Tornado plot saved.\n")

# ---- 6. PROBABILISTIC SENSITIVITY ANALYSIS (PSA) ----------------------------

cat("\n===== Running PSA (n=1000 simulations) =====\n")

n_psa <- 1000

# Parameter distributions for PSA
# OR: log-normal; utilities: beta; costs: gamma; discount: fixed

psa_results <- data.frame(
  sim = 1:n_psa, cost_STS = NA, qaly_STS = NA,
  cost_TM = NA, qaly_TM = NA, ICER = NA
)

for (s in 1:n_psa) {
  # Sample ORs from log-normal distributions
  # SD of log(OR) derived from 95% CI: sd = (log(upper) - log(lower)) / (2*1.96)
  # OR_mort: HR 0.70 (95% CI 0.50-0.96) — TIM-HF2 trial
  # OR_hosp: OR 0.64 (95% CI 0.39-0.95) — Kotb et al. 2015
  # Full log-normal distribution retained — occasional draws >1.0 reflect
  # genuine parameter uncertainty and keep the PSA unbiased
  or_m_s  <- exp(rnorm(1, log(0.70), (log(0.96)-log(0.50))/(2*1.96)))
  or_h_s  <- exp(rnorm(1, log(0.64), (log(0.95)-log(0.39))/(2*1.96)))
  
  # Beta-distributed utilities; SE = 20% of mean (Caillon et al. 2022 convention)
  # Caillon Table 2: NYHA I/II weighted = 0.771; NYHA III/IV weighted = 0.639
  # Truncated at disutil_hosp (0.212) as lower bound so net utility (u - disutil)
  # is always >= 0, preventing negative QALYs in PSA
  beta_params <- function(mu, se) {
    alpha <- mu * ((mu*(1-mu)/se^2) - 1)
    beta  <- (1-mu) * ((mu*(1-mu)/se^2) - 1)
    c(alpha, beta)
  }
  u_s <- c(
    NYHA_I   = max(rbeta(1, beta_params(0.771, 0.077)[1], beta_params(0.771, 0.077)[2]), disutil_hosp),
    NYHA_II  = max(rbeta(1, beta_params(0.771, 0.077)[1], beta_params(0.771, 0.077)[2]), disutil_hosp),
    NYHA_III = max(rbeta(1, beta_params(0.639, 0.064)[1], beta_params(0.639, 0.064)[2]), disutil_hosp),
    NYHA_IV  = max(rbeta(1, beta_params(0.639, 0.064)[1], beta_params(0.639, 0.064)[2]), disutil_hosp),
    Dead     = 0
  )
  
  # Gamma-distributed costs (shape = (mean/se)^2, rate = mean/se^2)
  ch_s   <- rgamma(1, shape = (5900/590)^2, rate = 5900/590^2)
  ctm_s  <- rgamma(1, shape = (150/20)^2,   rate = 150/20^2)
  csts_s <- rgamma(1, shape = (60/10)^2,    rate = 60/10^2)
  
  # Rebuild TM transition matrix
  mort_m_TM_s <- inv_logit(logit(mort_monthly_STS) + log(or_m_s))
  tp_TM_s     <- tp_STS
  for (i in 1:4) {
    tp_TM_s[i, 5]   <- mort_m_TM_s[i]
    tp_TM_s[i, 1:4] <- nyha_alive_STS[i, ] * (1 - mort_m_TM_s[i])
  }
  hosp_TM_s <- hosp_monthly_STS * or_h_s
  
  r_sts_s <- run_markov(tp_STS, tp_STS, hosp_monthly_STS, hosp_monthly_STS,
                        csts_s, n_cycles,
                        start_dist, u_s, cost_outpatient, ch_s, cost_palliative,
                        disutil_hosp, discount_r, cycle_len)
  r_tm_s  <- run_markov(tp_TM_s, tp_STS, hosp_TM_s, hosp_monthly_STS,
                        ctm_s, n_cycles,
                        start_dist, u_s, cost_outpatient, ch_s, cost_palliative,
                        disutil_hosp, discount_r, cycle_len, effect_duration = 60)
  
  psa_results$cost_STS[s] <- r_sts_s$cost
  psa_results$qaly_STS[s] <- r_sts_s$qaly
  psa_results$cost_TM[s]  <- r_tm_s$cost
  psa_results$qaly_TM[s]  <- r_tm_s$qaly
  psa_results$ICER[s]     <- (r_tm_s$cost - r_sts_s$cost) / (r_tm_s$qaly - r_sts_s$qaly)
}

psa_results$incr_cost <- psa_results$cost_TM - psa_results$cost_STS
psa_results$incr_qaly <- psa_results$qaly_TM - psa_results$qaly_STS

cat(sprintf("\nPSA Summary (n=%d):\n", n_psa))
cat(sprintf("  Median ICER: SGD %.0f/QALY\n", median(psa_results$ICER, na.rm = TRUE)))
cat(sprintf("  Mean ICER: SGD %.0f/QALY (note: skewed by extreme simulations; median preferred)\n",
            mean(psa_results$ICER, na.rm = TRUE)))
cat(sprintf("  95%% CI: [SGD %.0f, SGD %.0f]\n",
            quantile(psa_results$ICER, 0.025, na.rm = TRUE),
            quantile(psa_results$ICER, 0.975, na.rm = TRUE)))
pct_ce <- mean(psa_results$ICER < wtp, na.rm = TRUE) * 100
cat(sprintf("  Probability cost-effective at WTP SGD %d: %.1f%%\n", wtp, pct_ce))

# Cost-effectiveness scatter plot — no WTP line, colour by cost-effective status
cost_range <- range(psa_results$incr_cost, na.rm = TRUE)
qaly_range <- range(psa_results$incr_qaly, na.rm = TRUE)
x_lim <- c(min(qaly_range[1], -0.05), max(qaly_range[2], 0.05) * 1.1)
y_lim <- c(min(cost_range[1], -500) * 1.1, max(cost_range[2], 500) * 1.1)

psa_ce <- psa_results
psa_ce$ce_label <- ifelse(psa_ce$incr_cost < wtp * psa_ce$incr_qaly,
                          "Cost-effective", "Not cost-effective")

p_ce_plane <- ggplot(psa_ce, aes(x = incr_qaly, y = incr_cost, colour = ce_label)) +
  geom_point(alpha = 0.4, size = 0.9) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.6) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.6) +
  geom_point(data = data.frame(incr_qaly = incr_qaly, incr_cost = incr_cost),
             aes(x = incr_qaly, y = incr_cost),
             colour = "black", size = 4, shape = 18, inherit.aes = FALSE) +
  annotate("text", x = incr_qaly + 0.02, y = incr_cost + 110,
           label = sprintf("Base case\nICER = SGD %d/QALY", round(ICER)),
           size = 3, hjust = 0, colour = "grey20") +
  annotate("text", x = x_lim[1] * 0.5, y = y_lim[2] * 0.85,
           label = "More costly\nLess effective",
           size = 3, colour = "grey50", hjust = 0.5) +
  annotate("text", x = x_lim[2] * 0.75, y = y_lim[1] * 0.85,
           label = "Less costly\nMore effective",
           size = 3, colour = "grey50", hjust = 0.5) +
  annotate("text", x = x_lim[2] * 0.75, y = y_lim[2] * 0.85,
           label = sprintf("%.1f%% cost-effective\nat SGD 80,000/QALY", pct_ce),
           size = 3, colour = "steelblue", hjust = 0.5) +
  scale_colour_manual(values = c("Cost-effective"     = "steelblue",
                                  "Not cost-effective" = "#d62728"),
                      name = NULL) +
  scale_x_continuous(limits = x_lim,
                     labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(limits = y_lim, labels = scales::comma) +
  labs(title    = "Cost-effectiveness plane (PSA)",
       subtitle = sprintf("n=1,000 simulations · diamond = base case · WTP threshold = SGD 80,000/QALY"),
       x = "Incremental QALYs (TM vs STS)",
       y = "Incremental cost (SGD)") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave("ce_plane.png", p_ce_plane, width = 8, height = 6, dpi = 150)
print(p_ce_plane)
cat("CE plane plot saved.\n")

# CEAC — correct formula handles negative ICERs properly
# Uses net monetary benefit: TM is cost-effective at threshold w if
# incr_cost < w * incr_qaly (works for all signs of incr_cost and incr_qaly)
wtp_range <- seq(0, 200000, by = 2000)
ceac <- sapply(wtp_range, function(w) {
  mean(psa_results$incr_cost < w * psa_results$incr_qaly, na.rm = TRUE)
})
ceac_df <- data.frame(WTP = wtp_range, Prob_CE = ceac)
prob_at_wtp <- ceac_df$Prob_CE[ceac_df$WTP == 80000]

p_ceac <- ggplot(ceac_df, aes(x = WTP / 1000, y = Prob_CE)) +
  geom_line(colour = "steelblue", linewidth = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  geom_vline(xintercept = wtp / 1000, linetype = "dotted",
             colour = "red", linewidth = 0.8) +
  annotate("text", x = wtp / 1000 + 5, y = 0.55,
           label = sprintf("SGD 80,000/QALY\n%.1f%% cost-effective", prob_at_wtp * 100),
           colour = "red", size = 3, hjust = 0) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"),
                     limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.25)) +
  scale_x_continuous(breaks = seq(0, 200, by = 25)) +
  labs(title    = "(b) Cost-effectiveness acceptability curve (CEAC)",
       x = "Willingness-to-pay threshold (SGD thousands per QALY)",
       y = "Probability cost-effective") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

# Add panel label to CE plane
p_ce_plane <- p_ce_plane +
  labs(title = "(a) Cost-effectiveness plane (PSA)")

# Combine into one figure
p_psa_combined <- gridExtra::arrangeGrob(
  p_ce_plane, p_ceac,
  ncol = 2,
  top = grid::textGrob(
    "Figure 3. Probabilistic sensitivity analysis (PSA, n=1,000 simulations)",
    gp = grid::gpar(fontsize = 12, fontface = "bold"))
)

ggsave("psa_combined.png", p_psa_combined, width = 14, height = 6, dpi = 150)
grid::grid.newpage(); grid::grid.draw(p_psa_combined)
cat("Combined PSA figure saved.\n")

# Also save individually for reference
ggsave("ce_plane.png",  p_ce_plane, width = 7, height = 5, dpi = 150)
ggsave("ceac_plot.png", p_ceac,     width = 7, height = 5, dpi = 150)
cat("Individual PSA plots also saved.\n")

# ---- 7. SCENARIO ANALYSIS: 5-year horizon ------------------------------------

cat("\n===== Scenario: 5-year horizon =====\n")

n_cycles_5y <- 12 * 5

# 5-year horizon: TM effect is fully within 60-month window, so waning has no
# effect here — but pass consistent parameters for code integrity
r_sts_5y <- run_markov(tp_STS, tp_STS, hosp_monthly_STS, hosp_monthly_STS,
                       cost_STS_monthly, n_cycles_5y,
                       start_dist, utils, cost_outpatient, cost_hosp, cost_palliative,
                       disutil_hosp, discount_r, cycle_len)
r_tm_5y  <- run_markov(tp_TM,  tp_STS, hosp_monthly_TM,  hosp_monthly_STS,
                       cost_TM_monthly,  n_cycles_5y,
                       start_dist, utils, cost_outpatient, cost_hosp, cost_palliative,
                       disutil_hosp, discount_r, cycle_len, effect_duration = 60)

ICER_5y <- (r_tm_5y$cost - r_sts_5y$cost) / (r_tm_5y$qaly - r_sts_5y$qaly)
cat(sprintf("  5-year ICER: SGD %.0f/QALY\n", ICER_5y))

# ---- 8. SUMMARY TABLE --------------------------------------------------------

cat("\n===== Final Results Summary =====\n\n")

summary_table <- data.frame(
  Strategy = c("STS (Comparator)", "Telemonitoring (TM)"),
  Total_Cost_SGD = c(round(result_STS$cost), round(result_TM$cost)),
  Total_QALYs    = c(round(result_STS$qaly, 4), round(result_TM$qaly, 4)),
  Incr_Cost_SGD  = c(NA, round(incr_cost)),
  Incr_QALYs     = c(NA, round(incr_qaly, 4)),
  ICER_SGD_per_QALY = c(NA, round(ICER))
)
print(summary_table)

cat("\n===== Analysis Complete =====\n")
cat("Output files saved to your working directory:\n")
cat("  - trace_plot.png\n  - tornado_plot.png\n  - ce_plane.png\n  - ceac_plot.png\n")
