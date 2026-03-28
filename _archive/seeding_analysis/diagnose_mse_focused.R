#!/usr/bin/env Rscript
# ==============================================================================
# MSE-excess-only diagnostic: DIC vs DT across designs
# Key questions:
#   1. MSE excess distribution (not MAE) — where does each method hurt predictions?
#   2. How does curve flatness mediate failures?
#   3. Do DIC failure rates increase with n? DT failure rates?
# ==============================================================================

library(tidyverse)
EPS <- 0.7; REPS <- 5

# ==============================================================================
# LOAD + COMPUTE PER-COMPARISON MSE EXCESS
# ==============================================================================

per_comp <- function(dt_data, oracle_data) {
  oracle_sel <- oracle_data %>%
    filter(!is.na(mse_true)) %>%
    group_by(config, comp_no) %>%
    slice_min(mse_true, n=1, with_ties=FALSE) %>%
    select(config, comp_no, oracle_nbasis=nbasis, mse_oracle=mse_true)

  mse_lkp <- oracle_data %>% filter(!is.na(mse_true)) %>%
    select(config, comp_no, nbasis, mse_true)

  # MSE curve shape: width of band within X% of oracle
  curve_shape <- oracle_data %>%
    filter(!is.na(mse_true)) %>%
    left_join(oracle_sel %>% select(config, comp_no, mse_oracle),
              by=c("config","comp_no")) %>%
    mutate(rel_excess = (mse_true - mse_oracle)/mse_oracle*100) %>%
    group_by(config, comp_no) %>%
    summarise(
      band_2pct  = sum(rel_excess <= 2),   # models within 2% of oracle
      band_5pct  = sum(rel_excess <= 5),   # models within 5%
      slope_low  = (mse_true[nbasis == min(nbasis)] - mse_oracle) / mse_oracle * 100,
      .groups="drop"
    )

  dt_filt <- dt_data %>%
    filter(abs(epsilon - EPS) < 0.001, n_reps_used == REPS)

  dt_nll <- dt_filt %>% filter(loss_function == "plugin_NLL") %>%
    group_by(config, comp_no) %>%
    slice_min(metric, n=1, with_ties=FALSE) %>%
    select(config, comp_no, dt_nll_nbasis=nbasis)

  dic_sel <- oracle_data %>% filter(!is.na(DIC)) %>%
    group_by(config, comp_no) %>%
    slice_min(DIC, n=1, with_ties=FALSE) %>%
    select(config, comp_no, dic_nbasis=nbasis)

  oracle_sel %>%
    left_join(dic_sel,    by=c("config","comp_no")) %>%
    left_join(dt_nll,     by=c("config","comp_no")) %>%
    left_join(curve_shape,by=c("config","comp_no")) %>%
    left_join(mse_lkp %>% rename(dic_mse=mse_true),
              by=c("config","comp_no","dic_nbasis"="nbasis")) %>%
    left_join(mse_lkp %>% rename(dt_nll_mse=mse_true),
              by=c("config","comp_no","dt_nll_nbasis"="nbasis")) %>%
    mutate(
      dic_excess    = (dic_mse    - mse_oracle)/mse_oracle*100,
      dt_nll_excess = (dt_nll_mse - mse_oracle)/mse_oracle*100,
      flat_curve    = band_2pct >= 10,   # >= 10 of 20 models within 2% = flat
      steep_curve   = band_2pct <= 3
    )
}

cat("Loading data...\n")
pa_old <- readRDS("results_multi_config/paper_final/pa_method_comparison.RDS")
pa_new <- readRDS("results_multi_config/paper_final/pa_method_comparison_final.RDS")
eq     <- readRDS("results_multi_config/paper_final/equal_allocation_results.RDS")

pc_pa_old <- per_comp(pa_old$dt,    pa_old$oracle)
pc_pa_new <- per_comp(pa_new$dt,    pa_new$oracle)
pc_eq     <- per_comp(eq$dt_1fold,  eq$oracle)

all_pc <- bind_rows(
  pc_pa_old %>% mutate(dataset="PA-old"),
  pc_pa_new %>% mutate(dataset="PA-new"),
  pc_eq     %>% mutate(dataset="Equal")
)

# Add effective avg sample size per PUMA (approximate)
n_map <- c(
  equal_30=30, equal_40=40, equal_50=50, equal_75=75, equal_100=100, equal_125=125,
  prop_0.75pct=47, prop_1p25pct=78, prop_1p5pct=94, prop_1p75pct=110
)
all_pc <- all_pc %>% mutate(avg_n = n_map[config])

# ==============================================================================
# 1. MSE EXCESS DISTRIBUTION by design — mean, median, 90th pct, failure rate
# ==============================================================================

cat("\n=== MSE excess distribution by design (DIC | DT-NLL) ===\n")
cat("Format: mean / median / 90th pct / fail>10%\n\n")

all_pc %>%
  group_by(avg_n, config) %>%
  summarise(
    n = n(),
    dic_mean   = mean(dic_excess, na.rm=TRUE),
    dic_med    = median(dic_excess, na.rm=TRUE),
    dic_p90    = quantile(dic_excess, 0.9, na.rm=TRUE),
    dic_fail10 = mean(dic_excess > 10, na.rm=TRUE)*100,
    dtnll_mean = mean(dt_nll_excess, na.rm=TRUE),
    dtnll_med  = median(dt_nll_excess, na.rm=TRUE),
    dtnll_p90  = quantile(dt_nll_excess, 0.9, na.rm=TRUE),
    dtnll_fail10 = mean(dt_nll_excess > 10, na.rm=TRUE)*100,
    .groups="drop"
  ) %>%
  arrange(avg_n) %>%
  mutate(
    dic_str   = sprintf("%.1f / %.1f / %.1f / %.0f%%",
                        dic_mean, dic_med, dic_p90, dic_fail10),
    dtnll_str = sprintf("%.1f / %.1f / %.1f / %.0f%%",
                        dtnll_mean, dtnll_med, dtnll_p90, dtnll_fail10)
  ) %>%
  select(avg_n, config, n, DIC=dic_str, `DT-NLL`=dtnll_str) %>%
  print(n=50)

# ==============================================================================
# 2. CURVE FLATNESS: does flat curve = low MSE excess regardless of selection?
# ==============================================================================

cat("\n=== MSE excess by curve flatness ===\n")
all_pc %>%
  mutate(curve_type = case_when(
    flat_curve  ~ "flat (>=10 models within 2%)",
    steep_curve ~ "steep (<=3 models within 2%)",
    TRUE        ~ "medium"
  )) %>%
  group_by(curve_type) %>%
  summarise(
    n          = n(),
    pct_total  = n()/nrow(all_pc)*100,
    dic_mean   = mean(dic_excess, na.rm=TRUE),
    dtnll_mean = mean(dt_nll_excess, na.rm=TRUE),
    dic_fail10 = mean(dic_excess > 10, na.rm=TRUE)*100,
    dtnll_fail10 = mean(dt_nll_excess > 10, na.rm=TRUE)*100,
    .groups="drop"
  ) %>%
  print()

# ==============================================================================
# 3. DIC FAILURE RATE vs SAMPLE SIZE — does it increase with n?
# ==============================================================================

cat("\n=== DIC failure rate (excess >10%) by avg_n ===\n")
cat("   (clear monotone increase would support the 'DIC fails at large n' claim)\n\n")

all_pc %>%
  group_by(avg_n) %>%
  summarise(
    n_comps      = n(),
    dic_fail_pct = mean(dic_excess > 10, na.rm=TRUE)*100,
    dtnll_fail_pct = mean(dt_nll_excess > 10, na.rm=TRUE)*100,
    dic_mean     = mean(dic_excess, na.rm=TRUE),
    dtnll_mean   = mean(dt_nll_excess, na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(avg_n) %>%
  print()

# ==============================================================================
# 4. WHEN DIC FAILS (>10%): what does the curve look like?
#    Is DIC failure associated with steep curves?
# ==============================================================================

cat("\n=== When DIC fails: curve flatness distribution ===\n")
all_pc %>%
  filter(!is.na(dic_excess)) %>%
  mutate(dic_failed = dic_excess > 10) %>%
  group_by(dic_failed) %>%
  summarise(
    n            = n(),
    band_2pct    = mean(band_2pct),
    band_5pct    = mean(band_5pct),
    pct_flat     = mean(flat_curve)*100,
    pct_steep    = mean(steep_curve)*100,
    oracle_mean  = mean(oracle_nbasis),
    dic_sel      = mean(dic_nbasis, na.rm=TRUE),
    avg_n        = mean(avg_n),
    .groups="drop"
  ) %>%
  print()

cat("\n=== When DT-NLL fails: curve flatness distribution ===\n")
all_pc %>%
  filter(!is.na(dt_nll_excess)) %>%
  mutate(dt_failed = dt_nll_excess > 10) %>%
  group_by(dt_failed) %>%
  summarise(
    n            = n(),
    band_2pct    = mean(band_2pct),
    band_5pct    = mean(band_5pct),
    pct_flat     = mean(flat_curve)*100,
    pct_steep    = mean(steep_curve)*100,
    oracle_mean  = mean(oracle_nbasis),
    dt_sel       = mean(dt_nll_nbasis, na.rm=TRUE),
    avg_n        = mean(avg_n),
    .groups="drop"
  ) %>%
  print()

# ==============================================================================
# 5. CONDITIONAL ON STEEP CURVE: which method is more robust?
# ==============================================================================

cat("\n=== On STEEP curves only: DIC vs DT-NLL MSE excess ===\n")
all_pc %>%
  filter(steep_curve) %>%
  group_by(avg_n) %>%
  summarise(
    n            = n(),
    dic_mean     = mean(dic_excess, na.rm=TRUE),
    dtnll_mean   = mean(dt_nll_excess, na.rm=TRUE),
    dic_fail10   = mean(dic_excess > 10, na.rm=TRUE)*100,
    dtnll_fail10 = mean(dt_nll_excess > 10, na.rm=TRUE)*100,
    .groups="drop"
  ) %>%
  arrange(avg_n) %>%
  print()

# ==============================================================================
# 6. HIGH n only (avg_n >= 75): focus on the regime where DIC should fail
# ==============================================================================

cat("\n=== High n (avg_n >= 75): per-comparison excess scatter ===\n")
high_n <- all_pc %>% filter(avg_n >= 75) %>%
  select(config, avg_n, comp_no, oracle_nbasis, dic_nbasis, dt_nll_nbasis,
         dic_excess, dt_nll_excess, band_2pct, flat_curve) %>%
  arrange(avg_n, desc(dic_excess))

cat("DIC excess > 15% cases:\n")
high_n %>% filter(dic_excess > 15) %>%
  select(config, comp_no, oracle_nbasis, dic_nbasis, dic_excess,
         dt_nll_nbasis, dt_nll_excess, band_2pct) %>%
  print(n=30)

cat("\nDT-NLL excess > 15% cases:\n")
high_n %>% filter(dt_nll_excess > 15) %>%
  select(config, comp_no, oracle_nbasis, dt_nll_nbasis, dt_nll_excess,
         dic_nbasis, dic_excess, band_2pct) %>%
  print(n=30)
