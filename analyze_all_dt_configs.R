#!/usr/bin/env Rscript
# ==============================================================================
# All DT Configurations Analysis
# ==============================================================================
# Compare DIC, WAIC, and ALL DT configs (eps × n_reps)
# ==============================================================================

library(tidyverse)

cat("\n=== ALL DT CONFIGURATIONS ANALYSIS ===\n\n")

# Load data (10 configs)
dt_all <- readRDS("results_multi_config/dt_all_10configs.RDS")
oracle_all <- readRDS("results_multi_config/oracle_all_10configs.RDS")

oracle_selections <- oracle_all %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, oracle_nbasis = nbasis)

# DIC
dic_mad <- oracle_all %>%
  filter(!is.na(DIC)) %>%
  group_by(config, comp_no) %>%
  slice_min(DIC, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, selected_nbasis = nbasis) %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(deviation = selected_nbasis - oracle_nbasis) %>%
  group_by(config) %>%
  summarise(MAD = mean(abs(deviation)), .groups = "drop") %>%
  mutate(method = "DIC")

# WAIC
waic_mad <- oracle_all %>%
  filter(!is.na(WAIC)) %>%
  group_by(config, comp_no) %>%
  slice_min(WAIC, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, selected_nbasis = nbasis) %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(deviation = selected_nbasis - oracle_nbasis) %>%
  group_by(config) %>%
  summarise(MAD = mean(abs(deviation)), .groups = "drop") %>%
  mutate(method = "WAIC")

# All DT 1-fold configs
dt1_all_configs <- dt_all %>%
  filter(method_name == "dt_1fold") %>%
  group_by(config, comp_no, epsilon, n_reps_used) %>%
  slice_min(metric, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, epsilon, n_reps_used, selected_nbasis = nbasis) %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(deviation = selected_nbasis - oracle_nbasis)

dt1_mad <- dt1_all_configs %>%
  group_by(config, epsilon, n_reps_used) %>%
  summarise(MAD = mean(abs(deviation)), .groups = "drop") %>%
  mutate(method = sprintf("DT eps=%.1f n=%d", epsilon, n_reps_used))

# Combine
all_methods <- bind_rows(dic_mad, waic_mad, dt1_mad)

# For each config, show all methods
cat("=== ALL METHODS BY CONFIG ===\n\n")

for (cfg in unique(all_methods$config)) {
  cat(sprintf("\n--- %s ---\n\n", cfg))

  cfg_data <- all_methods %>%
    filter(config == cfg) %>%
    arrange(MAD)

  print(cfg_data %>%
          select(method, MAD) %>%
          knitr::kable(digits = 2))

  winner <- cfg_data$method[1]
  cat(sprintf("\nWinner: %s (MAD=%.2f)\n", winner, cfg_data$MAD[1]))
}

# Summary table: winner per config
winners <- all_methods %>%
  group_by(config) %>%
  slice_min(MAD, n = 1, with_ties = FALSE) %>%
  select(config, winner = method, best_mad = MAD)

cat("\n\n=== WINNER SUMMARY ===\n\n")
print(winners %>% knitr::kable(digits = 2))

cat("\n\nWinner counts:\n")
print(table(winners$winner))

# Separate DT winners by epsilon and n_reps
dt_winners <- winners %>%
  filter(grepl("^DT", winner)) %>%
  mutate(
    epsilon = as.numeric(str_extract(winner, "(?<=eps=)[0-9.]+")),
    n_reps = as.numeric(str_extract(winner, "(?<=n=)[0-9]+"))
  )

if (nrow(dt_winners) > 0) {
  cat("\n\nDT winner breakdown:\n")
  cat(sprintf("  eps=0.3: %d configs\n", sum(dt_winners$epsilon == 0.3)))
  cat(sprintf("  eps=0.5: %d configs\n", sum(dt_winners$epsilon == 0.5)))
  cat(sprintf("  eps=0.7: %d configs\n", sum(dt_winners$epsilon == 0.7)))
  cat(sprintf("  n=1: %d configs\n", sum(dt_winners$n_reps == 1)))
  cat(sprintf("  n=3: %d configs\n", sum(dt_winners$n_reps == 3)))
  cat(sprintf("  n=5: %d configs\n", sum(dt_winners$n_reps == 5)))

  cat("\n\nDT winning configs:\n")
  print(dt_winners %>%
          select(config, winner, epsilon, n_reps, best_mad) %>%
          knitr::kable(digits = 2))
}

# Show top 3 methods per config + Best DT vs DIC gap
cat("\n\n=== COMPREHENSIVE SUMMARY: TOP 3 + DT vs DIC GAP ===\n\n")

for (cfg in unique(all_methods$config)) {
  cat(sprintf("\n--- %s ---\n", cfg))

  cfg_data <- all_methods %>%
    filter(config == cfg) %>%
    arrange(MAD)

  # Top 3
  top3 <- cfg_data %>% head(3)
  cat("Top 3 Methods:\n")
  for (i in 1:nrow(top3)) {
    cat(sprintf("  %d. %-20s MAD=%.2f\n", i, top3$method[i], top3$MAD[i]))
  }

  # Best DT vs DIC gap
  dic_mad <- cfg_data %>% filter(method == "DIC") %>% pull(MAD)
  dt_mads <- cfg_data %>% filter(grepl("^DT", method))

  if (nrow(dt_mads) > 0 && length(dic_mad) > 0) {
    best_dt_mad <- min(dt_mads$MAD)
    best_dt_method <- dt_mads %>% filter(MAD == best_dt_mad) %>% slice(1) %>% pull(method)
    gap <- best_dt_mad - dic_mad

    if (gap < 0) {
      cat(sprintf("Best DT vs DIC Gap: %.2f → DT wins by %.2f (%s)\n",
                  gap, abs(gap), best_dt_method))
    } else {
      cat(sprintf("Best DT vs DIC Gap: +%.2f → DIC wins by %.2f\n", gap, gap))
    }
  }
}

cat("\n")
