# =============================================================================
# SPH5412 Individual Assignment
# Cost-Effectiveness of Telemonitoring vs Structured Telephone Support
# for Heart Failure in Singapore
# Author: Lin Fang Yu
# Date: March 2026
# =============================================================================

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dampack, ggplot2, dplyr, scales, knitr)

# =============================================================================
# SECTION 1: MODEL PARAMETERS
# =============================================================================

# --- Cycle and time settings ---
n_cycles     <- 12 * 35   # lifetime = 35 years x 12 months
cycle_length <- 1/12      # 1 month in years
discount_r   <- 0.03      # 3% annual discount rate (ACE Singapore recommendation)

# --- Health states ---
states <- c("NYHA_I", "NYHA_II", "NYHA_III", "NYHA_IV", "Death")
n_states <- length(states)

# --- Starting cohort distribution (NYHA class at discharge) ---
# Based on Jiang et al. (2020): mixed severity post-discharge cohort
start_dist <- c(NYHA_I = 0.10, NYHA_II = 0.30, NYHA_III = 0.40, NYHA_IV = 0.20, Death = 0.00)

# =============================================================================
# SECTION 2: TRANSITION PROBABILITY MATRICES
# =============================================================================
# Source: Liang et al. (2017) - sacubitril/valsartan Singapore model
# Monthly transition probabilities between NYHA classes

# --- Standard of Care (Structured Telephone Support) ---
tp_soc <- matrix(c(
  # To:  I      II     III    IV     Death
       0.700, 0.180, 0.080, 0.020, 0.020,  # From NYHA I
       0.120, 0.580, 0.220, 0.050, 0.030,  # From NYHA II
       0.050, 0.150, 0.560, 0.180, 0.060,  # From NYHA III
       0.020, 0.050, 0.180, 0.620, 0.130,  # From NYHA IV
       0.000, 0.000, 0.000, 0.000, 1.000   # From Death (absorbing)
), nrow = n_states, byrow = TRUE,
dimnames = list(states, states))

# Verify rows sum to 1
stopifnot(all(abs(rowSums(tp_soc) - 1) < 1e-9))

# --- Telemonitoring (TM) ---
# Odds ratios applied from Kotb et al. (2015) meta-analysis:
# OR for hospitalisation: 0.71; OR for mortality: 0.87
# TM improves transitions: less worsening, more improvement, lower death

tp_tm <- matrix(c(
  # To:  I      II     III    IV     Death
       0.740, 0.170, 0.060, 0.015, 0.015,  # From NYHA I
       0.150, 0.600, 0.180, 0.040, 0.030,  # From NYHA II
       0.070, 0.180, 0.570, 0.140, 0.040,  # From NYHA III
       0.030, 0.070, 0.200, 0.610, 0.090,  # From NYHA IV
       0.000, 0.000, 0.000, 0.000, 1.000   # From Death (absorbing)
), nrow = n_states, byrow = TRUE,
dimnames = list(states, states))

stopifnot(all(abs(rowSums(tp_tm) - 1) < 1e-9))

# =============================================================================
# SECTION 3: UTILITY VALUES
# =============================================================================
# Source: Liang et al. (2017); disutility for hospitalisation from published literature

utilities <- c(
  NYHA_I   = 0.848,
  NYHA_II  = 0.762,
  NYHA_III = 0.644,
  NYHA_IV  = 0.432,
  Death    = 0.000
)

disutility_hosp <- 0.082  # per hospitalisation event

# =============================================================================
# SECTION 4: COST INPUTS (SGD)
# =============================================================================
# Source: Senanayake et al. (2025) - SG HFrEF costs; local programme data

costs <- c(
  NYHA_I   = 800,    # monthly outpatient/ambulatory cost
  NYHA_II  = 1200,
  NYHA_III = 2100,
  NYHA_IV  = 3800,
  Death    = 0
)

cost_hosp       <- 8500   # per HF hospitalisation (Senanayake et al., 2025)
cost_palliative <- 4200   # per death in NYHA IV (local hospital data)
cost_tm         <- 380    # monthly cost of telemonitoring programme (CGH)
cost_sts        <- 120    # monthly cost of structured telephone support

# Hospitalisation probability per cycle (monthly) by NYHA class
p_hosp <- c(NYHA_I = 0.02, NYHA_II = 0.05, NYHA_III = 0.12, NYHA_IV = 0.25, Death = 0)
p_hosp_tm <- p_hosp * 0.71  # applying OR from Kotb et al. (2015)

# =============================================================================
# SECTION 5: MARKOV MODEL FUNCTION
# =============================================================================

run_markov <- function(tp, util, cost_state, cost_intervention,
                       p_hospitalisation, n_cycles, cycle_length,
                       discount_r, start_dist) {

  n_states <- nrow(tp)
  state_names <- rownames(tp)

  # Initialise trace matrix
  trace <- matrix(0, nrow = n_cycles + 1, ncol = n_states,
                  dimnames = list(0:n_cycles, state_names))
  trace[1, ] <- start_dist

  # Discount vectors
  disc_costs  <- 1 / (1 + discount_r)^((0:n_cycles) * cycle_length)
  disc_qalys  <- 1 / (1 + discount_r)^((0:n_cycles) * cycle_length)

  total_cost  <- 0
  total_qaly  <- 0

  for (i in 1:n_cycles) {
    # State occupancy this cycle
    occ <- trace[i, ]

    # Cycle costs: state costs + intervention cost + hospitalisation cost
    hosp_cost <- sum(occ * p_hospitalisation * cost_hosp)
    state_cost <- sum(occ * cost_state)
    cycle_cost <- (state_cost + hosp_cost + cost_intervention) * disc_costs[i]

    # Cycle QALYs: adjust for hospitalisation disutility
    hosp_disutil <- sum(occ * p_hospitalisation * disutility_hosp)
    cycle_qaly   <- (sum(occ * util) - hosp_disutil) * cycle_length * disc_qalys[i]

    total_cost <- total_cost + cycle_cost
    total_qaly <- total_qaly + cycle_qaly

    # Transition to next cycle
    trace[i + 1, ] <- occ %*% tp
  }

  return(list(
    trace      = trace,
    total_cost = total_cost,
    total_qaly = total_qaly
  ))
}

# =============================================================================
# SECTION 6: BASE CASE ANALYSIS
# =============================================================================

cat("\n========== BASE CASE RESULTS ==========\n")

# Run both arms
res_soc <- run_markov(tp_soc, utilities, costs, cost_sts,
                      p_hosp, n_cycles, cycle_length, discount_r, start_dist)

res_tm  <- run_markov(tp_tm,  utilities, costs, cost_tm,
                      p_hosp_tm, n_cycles, cycle_length, discount_r, start_dist)

# ICER
delta_cost <- res_tm$total_cost - res_soc$total_cost
delta_qaly <- res_tm$total_qaly - res_soc$total_qaly
icer       <- delta_cost / delta_qaly

# Print results table
results_table <- data.frame(
  Strategy      = c("Structured Tel. Support (STS)", "Telemonitoring (TM)"),
  Total_Cost_SGD = round(c(res_soc$total_cost, res_tm$total_cost), 0),
  Total_QALYs   = round(c(res_soc$total_qaly, res_tm$total_qaly), 4),
  Inc_Cost      = c(NA, round(delta_cost, 0)),
  Inc_QALYs     = c(NA, round(delta_qaly, 4)),
  ICER          = c(NA, round(icer, 0))
)
print(results_table)
cat(sprintf("\nICER: SGD %.0f per QALY gained\n", icer))
cat("Singapore WTP threshold: SGD 50,000–100,000 per QALY\n")
if (icer < 100000) cat(">> Telemonitoring is COST-EFFECTIVE at SGD 100,000/QALY threshold\n")

# =============================================================================
# SECTION 7: MARKOV TRACE PLOT
# =============================================================================

trace_soc <- as.data.frame(res_soc$trace)
trace_soc$cycle <- 0:n_cycles
trace_soc$arm   <- "STS"

trace_tm <- as.data.frame(res_tm$trace)
trace_tm$cycle <- 0:n_cycles
trace_tm$arm   <- "TM"

trace_long <- bind_rows(
  trace_soc %>% tidyr::pivot_longer(cols = all_of(states), names_to = "State", values_to = "Proportion"),
  trace_tm  %>% tidyr::pivot_longer(cols = all_of(states), names_to = "State", values_to = "Proportion")
)

p_trace <- ggplot(trace_long %>% filter(State != "Death"),
       aes(x = cycle / 12, y = Proportion, colour = State, linetype = arm)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c("#2166AC", "#4DAC26", "#F4A582", "#D6604D")) +
  labs(title = "Markov Trace: State Occupancy Over Time",
       x = "Time (years)", y = "Proportion of cohort",
       colour = "Health State", linetype = "Strategy") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(p_trace)
ggsave("markov_trace.png", p_trace, width = 8, height = 5, dpi = 150)
cat("\nSaved: markov_trace.png\n")

# =============================================================================
# SECTION 8: DETERMINISTIC SENSITIVITY ANALYSIS (DSA)
# =============================================================================

cat("\n========== ONE-WAY DSA ==========\n")

base_icer <- icer
params_dsa <- list(
  "Cost of TM (SGD 280)"       = list(cost_tm = 280),
  "Cost of TM (SGD 480)"       = list(cost_tm = 480),
  "Cost of hospitalisation (SGD 6000)" = list(cost_hosp = 6000),
  "Cost of hospitalisation (SGD 11000)"= list(cost_hosp = 11000),
  "Discount rate 0%"           = list(discount_r = 0.00),
  "Discount rate 5%"           = list(discount_r = 0.05),
  "Utility NYHA II (0.70)"     = list(util_II = 0.700),
  "Utility NYHA II (0.82)"     = list(util_II = 0.820),
  "OR hospitalisation TM 0.60" = list(or_hosp = 0.60),
  "OR hospitalisation TM 0.82" = list(or_hosp = 0.82)
)

dsa_results <- data.frame(Parameter = character(), ICER = numeric())

for (param_name in names(params_dsa)) {
  p <- params_dsa[[param_name]]

  ct  <- if (!is.null(p$cost_tm))    p$cost_tm    else cost_tm
  ch  <- if (!is.null(p$cost_hosp))  p$cost_hosp  else cost_hosp
  dr  <- if (!is.null(p$discount_r)) p$discount_r else discount_r
  uII <- if (!is.null(p$util_II))    p$util_II    else utilities["NYHA_II"]
  orh <- if (!is.null(p$or_hosp))    p$or_hosp    else 0.71

  util_tmp       <- utilities
  util_tmp["NYHA_II"] <- uII
  p_hosp_tmp     <- p_hosp * orh

  r_soc_tmp <- run_markov(tp_soc, util_tmp, costs, cost_sts,
                           p_hosp, n_cycles, cycle_length, dr, start_dist)
  r_tm_tmp  <- run_markov(tp_tm,  util_tmp, costs, ct,
                           p_hosp_tmp, n_cycles, cycle_length, dr, start_dist)

  icer_tmp <- (r_tm_tmp$total_cost - r_soc_tmp$total_cost) /
              (r_tm_tmp$total_qaly - r_soc_tmp$total_qaly)

  dsa_results <- rbind(dsa_results, data.frame(Parameter = param_name, ICER = round(icer_tmp, 0)))
}

print(dsa_results)

# Tornado plot
dsa_results$Difference <- dsa_results$ICER - base_icer
dsa_results <- dsa_results[order(abs(dsa_results$Difference)), ]
dsa_results$Parameter <- factor(dsa_results$Parameter, levels = dsa_results$Parameter)

p_tornado <- ggplot(dsa_results, aes(x = ICER, y = Parameter,
                                      fill = Difference > 0)) +
  geom_col(width = 0.6) +
  geom_vline(xintercept = base_icer, linetype = "dashed", colour = "black") +
  scale_fill_manual(values = c("#4DAC26", "#D6604D"), guide = "none") +
  labs(title = "Tornado Diagram — One-Way DSA",
       x = "ICER (SGD per QALY)", y = "") +
  theme_minimal(base_size = 11)

print(p_tornado)
ggsave("tornado_dsa.png", p_tornado, width = 8, height = 5, dpi = 150)
cat("Saved: tornado_dsa.png\n")

# =============================================================================
# SECTION 9: PROBABILISTIC SENSITIVITY ANALYSIS (PSA)
# =============================================================================

cat("\n========== PSA (1000 iterations) ==========\n")

set.seed(42)
n_psa <- 1000

# Sample parameters from distributions
psa_cost_tm   <- rgamma(n_psa, shape = 25,  rate = 25/cost_tm)
psa_cost_hosp <- rgamma(n_psa, shape = 25,  rate = 25/cost_hosp)
psa_util_I    <- rbeta(n_psa,  shape1 = 85, shape2 = 15)
psa_util_II   <- rbeta(n_psa,  shape1 = 76, shape2 = 24)
psa_util_III  <- rbeta(n_psa,  shape1 = 64, shape2 = 36)
psa_util_IV   <- rbeta(n_psa,  shape1 = 43, shape2 = 57)
psa_or_hosp   <- rlnorm(n_psa, meanlog = log(0.71), sdlog = 0.15)

psa_results <- data.frame(
  iter       = 1:n_psa,
  cost_soc   = NA_real_,
  cost_tm    = NA_real_,
  qaly_soc   = NA_real_,
  qaly_tm    = NA_real_,
  icer       = NA_real_
)

for (i in 1:n_psa) {
  util_i <- c(psa_util_I[i], psa_util_II[i], psa_util_III[i], psa_util_IV[i], 0)
  names(util_i) <- states

  p_hosp_i    <- p_hosp
  p_hosp_tm_i <- p_hosp * psa_or_hosp[i]

  r_s <- run_markov(tp_soc, util_i, costs, cost_sts,
                    p_hosp_i, n_cycles, cycle_length, discount_r, start_dist)
  r_t <- run_markov(tp_tm, util_i, costs, psa_cost_tm[i],
                    p_hosp_tm_i, n_cycles, cycle_length, discount_r, start_dist)

  psa_results$cost_soc[i] <- r_s$total_cost
  psa_results$cost_tm[i]  <- r_t$total_cost
  psa_results$qaly_soc[i] <- r_s$total_qaly
  psa_results$qaly_tm[i]  <- r_t$total_qaly
  psa_results$icer[i]     <- (r_t$total_cost - r_s$total_cost) /
                              (r_t$total_qaly - r_s$total_qaly)
}

# Cost-effectiveness plane
p_ce_plane <- ggplot(psa_results,
       aes(x = qaly_tm - qaly_soc, y = cost_tm - cost_soc)) +
  geom_point(alpha = 0.25, colour = "#2166AC", size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_abline(slope = 100000, intercept = 0,
              linetype = "dotted", colour = "red") +
  annotate("text", x = max(psa_results$qaly_tm - psa_results$qaly_soc) * 0.7,
           y = max(psa_results$cost_tm - psa_results$cost_soc) * 0.9,
           label = "WTP = SGD 100,000/QALY", colour = "red", size = 3) +
  labs(title = "Cost-Effectiveness Plane (PSA, n=1000)",
       x = "Incremental QALYs (TM vs STS)",
       y = "Incremental Cost SGD (TM vs STS)") +
  theme_minimal(base_size = 12)

print(p_ce_plane)
ggsave("ce_plane_psa.png", p_ce_plane, width = 7, height = 6, dpi = 150)
cat("Saved: ce_plane_psa.png\n")

# CEAC
wtp_range <- seq(0, 200000, by = 5000)
ceac <- sapply(wtp_range, function(wtp) {
  nb_tm  <- psa_results$qaly_tm  * wtp - psa_results$cost_tm
  nb_soc <- psa_results$qaly_soc * wtp - psa_results$cost_soc
  mean(nb_tm > nb_soc)
})

ceac_df <- data.frame(WTP = wtp_range, Prob_CE = ceac)

p_ceac <- ggplot(ceac_df, aes(x = WTP / 1000, y = Prob_CE)) +
  geom_line(colour = "#2166AC", linewidth = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 100, linetype = "dotted", colour = "red") +
  annotate("text", x = 105, y = 0.15,
           label = "SGD 100k\nthreshold", colour = "red", size = 3) +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Cost-Effectiveness Acceptability Curve (CEAC)",
       x = "Willingness-to-Pay Threshold (SGD '000 per QALY)",
       y = "Probability TM is Cost-Effective") +
  theme_minimal(base_size = 12)

print(p_ceac)
ggsave("ceac.png", p_ceac, width = 7, height = 5, dpi = 150)
cat("Saved: ceac.png\n")

# PSA summary
cat(sprintf("\nPSA Summary (n = %d):\n", n_psa))
cat(sprintf("  Mean ICER: SGD %.0f per QALY\n", mean(psa_results$icer)))
cat(sprintf("  95%% CI:    SGD %.0f to %.0f\n",
    quantile(psa_results$icer, 0.025),
    quantile(psa_results$icer, 0.975)))
cat(sprintf("  P(CE at SGD 100,000/QALY): %.1f%%\n",
    ceac_df$Prob_CE[ceac_df$WTP == 100000] * 100))

# =============================================================================
# SECTION 10: SCENARIO ANALYSES
# =============================================================================

cat("\n========== SCENARIO ANALYSES ==========\n")

scenarios <- list(
  "Base case (35 years)"   = list(n_cycles = 12 * 35),
  "10-year horizon"        = list(n_cycles = 12 * 10),
  "5-year horizon"         = list(n_cycles = 12 * 5)
)

scen_results <- data.frame()
for (s in names(scenarios)) {
  nc <- scenarios[[s]]$n_cycles
  r_s <- run_markov(tp_soc, utilities, costs, cost_sts,
                    p_hosp, nc, cycle_length, discount_r, start_dist)
  r_t <- run_markov(tp_tm,  utilities, costs, cost_tm,
                    p_hosp_tm, nc, cycle_length, discount_r, start_dist)
  icer_s <- (r_t$total_cost - r_s$total_cost) / (r_t$total_qaly - r_s$total_qaly)
  scen_results <- rbind(scen_results, data.frame(
    Scenario   = s,
    Cost_STS   = round(r_s$total_cost, 0),
    Cost_TM    = round(r_t$total_cost, 0),
    QALY_STS   = round(r_s$total_qaly, 3),
    QALY_TM    = round(r_t$total_qaly, 3),
    ICER_SGD   = round(icer_s, 0)
  ))
}
print(scen_results)

cat("\n========== MODEL COMPLETE ==========\n")
cat("Outputs: markov_trace.png, tornado_dsa.png, ce_plane_psa.png, ceac.png\n")
