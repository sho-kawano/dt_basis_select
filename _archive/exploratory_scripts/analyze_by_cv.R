#!/usr/bin/env Rscript
# ==============================================================================
# Analysis by Coefficient of Variation (Signal-to-Noise)
# ==============================================================================
# Examine how CV relates to method performance

library(tidyverse)

cat("\n=== ANALYSIS BY COEFFICIENT OF VARIATION ===\n\n")

# Load existing oracle and DT results
oracle_all <- readRDS("results_multi_config/oracle_all.RDS")
dt_all <- readRDS("results_multi_config/dt_all.RDS")

configs <- c("equal_40", "equal_50", "prop_0.5pct", "prop_1pct", "prop_2pct")

# Map config names to directory names
config_to_dir <- c(
  "equal_40" = "_results_equal40_comparison",
  "equal_50" = "_results_equal50_comparison",
  "prop_0.5pct" = "_results_prop0.5pct_comparison",
  "prop_1pct" = "_results_prop1pct_comparison",
  "prop_2pct" = "_results_prop2pct_comparison"
)

# ==============================================================================
# 1. Compute CV statistics per config
# ==============================================================================

cv_stats <- map_df(configs, function(cfg) {
  results_dir <- config_to_dir[cfg]
  if (!dir.exists(results_dir)) {
    return(NULL)
  }

  comp_dirs <- list.files(results_dir, pattern = "^comparison_", full.names = TRUE)

  comp_data <- map_df(comp_dirs, function(comp_dir) {
    comp_no <- as.integer(str_extract(basename(comp_dir), "[0-9]+"))
    z_list <- readRDS(file.path(comp_dir, "z.RDS"))
    d_list <- readRDS(file.path(comp_dir, "d.RDS"))

    z <- z_list$values
    d <- d_list$values

    tibble(
      config = cfg,
      comp_no = comp_no,
      puma_id = 1:length(z),
      z = z,
      d = d,
      se = sqrt(d),
      cv = sqrt(d) / abs(z)
    )
  })

  # Summary per config
  comp_data %>%
    group_by(config) %>%
    summarise(
      mean_cv = mean(cv, na.rm = TRUE),
      median_cv = median(cv, na.rm = TRUE),
      sd_cv = sd(cv, na.rm = TRUE),
      min_cv = min(cv, na.rm = TRUE),
      max_cv = max(cv, na.rm = TRUE),
      .groups = "drop"
    )
})

cat("=== CV STATISTICS BY CONFIG ===\n\n")
print(cv_stats %>% knitr::kable(digits = 3))

# ==============================================================================
# 2. Oracle characteristics by config
# ==============================================================================

oracle_summary <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  group_by(config) %>%
  summarise(
    oracle_mean_nbasis = mean(nbasis),
    oracle_sd_nbasis = sd(nbasis),
    .groups = "drop"
  )

# MSE variation
mse_var <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  summarise(
    mse_range = max(mse_true) - min(mse_true),
    mse_min = min(mse_true),
    .groups = "drop"
  ) %>%
  group_by(config) %>%
  summarise(
    mse_variation_pct = mean(mse_range / mse_min) * 100,
    .groups = "drop"
  )

# ==============================================================================
# 3. Method performance by config
# ==============================================================================

# Get oracle selections
oracle_selections <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, oracle_nbasis = nbasis)

# DIC performance
dic_perf <- oracle_all %>%
  filter(!is.na(DIC)) %>%
  group_by(config, comp_no) %>%
  slice_min(DIC, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, dic_nbasis = nbasis) %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(dic_deviation = dic_nbasis - oracle_nbasis) %>%
  group_by(config) %>%
  summarise(
    dic_mad = mean(abs(dic_deviation)),
    .groups = "drop"
  )

# WAIC performance
waic_perf <- oracle_all %>%
  filter(!is.na(WAIC)) %>%
  group_by(config, comp_no) %>%
  slice_min(WAIC, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, waic_nbasis = nbasis) %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(waic_deviation = waic_nbasis - oracle_nbasis) %>%
  group_by(config) %>%
  summarise(
    waic_mad = mean(abs(waic_deviation)),
    .groups = "drop"
  )

# Best DT performance (across all configs)
dt_perf <- dt_all %>%
  filter(method_name == "dt_1fold") %>%
  group_by(config, comp_no, epsilon, n_reps_used) %>%
  slice_min(metric, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, epsilon, n_reps_used, dt_nbasis = nbasis) %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(dt_deviation = dt_nbasis - oracle_nbasis) %>%
  group_by(config, epsilon, n_reps_used) %>%
  summarise(
    dt_mad = mean(abs(dt_deviation)),
    .groups = "drop"
  ) %>%
  group_by(config) %>%
  slice_min(dt_mad, n = 1, with_ties = FALSE) %>%
  select(config, best_dt_mad = dt_mad, best_dt_eps = epsilon, best_dt_nreps = n_reps_used)

# ==============================================================================
# 4. Combine everything
# ==============================================================================

combined <- cv_stats %>%
  left_join(oracle_summary, by = "config") %>%
  left_join(mse_var, by = "config") %>%
  left_join(dic_perf, by = "config") %>%
  left_join(waic_perf, by = "config") %>%
  left_join(dt_perf, by = "config") %>%
  arrange(mean_cv)

cat("\n\n=== COMBINED ANALYSIS (sorted by mean CV) ===\n\n")
print(combined %>%
        select(config, mean_cv, mse_variation_pct,
               oracle_mean_nbasis, oracle_sd_nbasis,
               dic_mad, waic_mad, best_dt_mad) %>%
        knitr::kable(digits = 2))

# ==============================================================================
# 5. Which method wins by CV regime?
# ==============================================================================

cat("\n\n=== WINNER BY CONFIG (sorted by CV) ===\n\n")

winners <- combined %>%
  mutate(
    best_mad = pmin(dic_mad, waic_mad, best_dt_mad, na.rm = TRUE),
    winner = case_when(
      best_mad == dic_mad ~ "DIC",
      best_mad == waic_mad ~ "WAIC",
      best_mad == best_dt_mad ~ sprintf("DT eps=%.1f n=%d", best_dt_eps, best_dt_nreps),
      TRUE ~ "Unknown"
    )
  ) %>%
  select(config, mean_cv, winner, best_mad, dic_mad, waic_mad, best_dt_mad)

print(winners %>% knitr::kable(digits = 2))

# ==============================================================================
# 6. Correlations
# ==============================================================================

cat("\n\n=== CORRELATIONS ===\n\n")

# Does CV predict which method wins?
cv_vs_perf <- combined %>%
  select(mean_cv, mse_variation_pct, oracle_sd_nbasis,
         dic_mad, waic_mad, best_dt_mad)

correlations <- tibble(
  cor_cv_msevar = cor(cv_vs_perf$mean_cv, cv_vs_perf$mse_variation_pct),
  cor_cv_oracle_sd = cor(cv_vs_perf$mean_cv, cv_vs_perf$oracle_sd_nbasis),
  cor_cv_dic_mad = cor(cv_vs_perf$mean_cv, cv_vs_perf$dic_mad),
  cor_cv_waic_mad = cor(cv_vs_perf$mean_cv, cv_vs_perf$waic_mad),
  cor_cv_dt_mad = cor(cv_vs_perf$mean_cv, cv_vs_perf$best_dt_mad)
)

print(correlations %>% knitr::kable(digits = 3))

cat("\nInterpretation:\n")
cat("- cor_cv_msevar: Does higher CV → higher MSE variation?\n")
cat("- cor_cv_oracle_sd: Does higher CV → more/less stable oracle?\n")
cat("- cor_cv_dic_mad: Does higher CV → worse DIC performance?\n")
cat("- cor_cv_waic_mad: Does higher CV → worse WAIC performance?\n")
cat("- cor_cv_dt_mad: Does higher CV → worse/better DT performance?\n")

# ==============================================================================
# 7. Visualize CV regimes
# ==============================================================================

cat("\n\n=== CV REGIMES ===\n\n")

cv_regimes <- combined %>%
  mutate(
    cv_regime = case_when(
      mean_cv < 0.15 ~ "Low CV (<0.15)",
      mean_cv < 0.20 ~ "Medium CV (0.15-0.20)",
      TRUE ~ "High CV (>0.20)"
    )
  ) %>%
  select(config, cv_regime, mean_cv, winner, best_mad)

print(cv_regimes %>% knitr::kable(digits = 3))

cat("\n")
