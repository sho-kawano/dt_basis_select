#!/usr/bin/env Rscript
# ==============================================================================
# Why did DIC/WAIC fail at prop_2pct?
# ==============================================================================

library(tidyverse)

cat("\n=== DIAGNOSTIC: Why DIC/WAIC Failed at prop_2pct ===\n\n")

# Load data
oracle_all <- readRDS("results_multi_config/oracle_all.RDS")

# Get optimal nbasis per config/comparison
optimal <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  ungroup()

# Oracle characteristics by config
oracle_stats <- optimal %>%
  group_by(config) %>%
  summarise(
    n_comparisons = n(),
    oracle_mean_nbasis = mean(nbasis),
    oracle_sd_nbasis = sd(nbasis),
    oracle_min = min(nbasis),
    oracle_max = max(nbasis),
    .groups = "drop"
  )

# MSE variation by config
mse_variation <- oracle_all %>%
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
    avg_mse_range = mean(mse_range),
    avg_mse_min = mean(mse_min),
    .groups = "drop"
  )

# Sample size info (from comp 1 of each config)
sample_info <- map_df(unique(oracle_all$config), function(cfg) {
  d_file <- sprintf("results_multi_config/%s/comparison_001/d.RDS", cfg)
  if (!file.exists(d_file)) {
    return(tibble(config = cfg, n_total = NA, avg_n_per_puma = NA, d_ratio = NA))
  }
  d <- readRDS(d_file)
  tibble(
    config = cfg,
    n_total = sum(1/d),
    n_pumas = length(d),
    avg_n_per_puma = sum(1/d) / length(d),
    d_ratio = max(d) / min(d)
  )
})

# DIC performance
dic_perf <- oracle_all %>%
  filter(!is.na(DIC)) %>%
  group_by(config, comp_no) %>%
  slice_min(DIC, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, dic_nbasis = nbasis) %>%
  left_join(optimal %>% select(config, comp_no, oracle_nbasis = nbasis),
            by = c("config", "comp_no")) %>%
  mutate(dic_deviation = dic_nbasis - oracle_nbasis) %>%
  group_by(config) %>%
  summarise(
    dic_mean_nbasis = mean(dic_nbasis),
    dic_sd_nbasis = sd(dic_nbasis),
    dic_mean_dev = mean(dic_deviation),
    dic_mad = mean(abs(dic_deviation)),
    .groups = "drop"
  )

# WAIC performance
waic_perf <- oracle_all %>%
  filter(!is.na(WAIC)) %>%
  group_by(config, comp_no) %>%
  slice_min(WAIC, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, waic_nbasis = nbasis) %>%
  left_join(optimal %>% select(config, comp_no, oracle_nbasis = nbasis),
            by = c("config", "comp_no")) %>%
  mutate(waic_deviation = waic_nbasis - oracle_nbasis) %>%
  group_by(config) %>%
  summarise(
    waic_mean_nbasis = mean(waic_nbasis),
    waic_sd_nbasis = sd(waic_nbasis),
    waic_mean_dev = mean(waic_deviation),
    waic_mad = mean(abs(waic_deviation)),
    .groups = "drop"
  )

# Combine all stats
summary_stats <- oracle_stats %>%
  left_join(mse_variation, by = "config") %>%
  left_join(sample_info, by = "config") %>%
  left_join(dic_perf, by = "config") %>%
  left_join(waic_perf, by = "config")

cat("=== CONFIG CHARACTERISTICS ===\n\n")
print(summary_stats %>%
        select(config, n_total, avg_n_per_puma, d_ratio,
               oracle_mean_nbasis, oracle_sd_nbasis,
               mse_variation_pct) %>%
        knitr::kable(digits = 1))

cat("\n\n=== DIC PERFORMANCE ===\n\n")
print(summary_stats %>%
        select(config, oracle_mean_nbasis, dic_mean_nbasis, dic_sd_nbasis,
               dic_mean_dev, dic_mad) %>%
        knitr::kable(digits = 1))

cat("\n\n=== WAIC PERFORMANCE ===\n\n")
print(summary_stats %>%
        select(config, oracle_mean_nbasis, waic_mean_nbasis, waic_sd_nbasis,
               waic_mean_dev, waic_mad) %>%
        knitr::kable(digits = 1))

# Deep dive into prop_2pct
cat("\n\n=== DEEP DIVE: prop_2pct ===\n\n")

oracle_2pct <- oracle_all %>% filter(config == "prop_2pct")
optimal_2pct <- optimal %>% filter(config == "prop_2pct")

# DIC selections
dic_2pct <- oracle_2pct %>%
  filter(!is.na(DIC)) %>%
  group_by(comp_no) %>%
  slice_min(DIC, n = 1, with_ties = FALSE) %>%
  select(comp_no, dic_nbasis = nbasis)

comparison_2pct <- optimal_2pct %>%
  select(comp_no, oracle_nbasis = nbasis, oracle_mse = mse_true) %>%
  left_join(dic_2pct, by = "comp_no") %>%
  mutate(dic_deviation = dic_nbasis - oracle_nbasis)

cat("DIC selections vs Oracle (prop_2pct):\n\n")
print(comparison_2pct %>%
        select(comp_no, oracle_nbasis, dic_nbasis, dic_deviation) %>%
        arrange(comp_no) %>%
        knitr::kable())

cat("\n\nDIC deviation distribution:\n")
print(table(comparison_2pct$dic_deviation))

cat("\n\nDIC selected nbasis distribution:\n")
print(table(comparison_2pct$dic_nbasis))

cat("\n\nOracle nbasis distribution:\n")
print(table(comparison_2pct$oracle_nbasis))

# Check DIC curve shape for prop_2pct
cat("\n\n=== HYPOTHESIS: Is DIC consistently selecting high nbasis? ===\n\n")

# For each comp, show top 5 models by DIC
top5_dic <- oracle_2pct %>%
  filter(!is.na(DIC)) %>%
  group_by(comp_no) %>%
  arrange(comp_no, desc(DIC)) %>%  # Higher DIC = better (less negative)
  slice(1:5) %>%
  summarise(
    top1_nbasis = nbasis[1],
    top2_nbasis = nbasis[2],
    top3_nbasis = nbasis[3],
    top4_nbasis = nbasis[4],
    top5_nbasis = nbasis[5],
    .groups = "drop"
  )

cat("Top 5 models by DIC (highest/least negative) for each comparison:\n\n")
print(top5_dic %>% knitr::kable())

# Frequency analysis
all_top5 <- c(top5_dic$top1_nbasis, top5_dic$top2_nbasis,
              top5_dic$top3_nbasis, top5_dic$top4_nbasis,
              top5_dic$top5_nbasis)
cat("\n\nFrequency of nbasis in top 5 DIC selections:\n")
print(sort(table(all_top5), decreasing = TRUE))

# Compare to other configs
cat("\n\n=== COMPARISON: DIC behavior across configs ===\n\n")

dic_behavior <- oracle_all %>%
  filter(!is.na(DIC)) %>%
  group_by(config, comp_no) %>%
  arrange(config, comp_no, desc(DIC)) %>%
  slice(1:3) %>%
  summarise(
    top1_nbasis = nbasis[1],
    top2_nbasis = nbasis[2],
    top3_nbasis = nbasis[3],
    .groups = "drop"
  ) %>%
  group_by(config) %>%
  summarise(
    mean_top1 = mean(top1_nbasis),
    mean_top2 = mean(top2_nbasis),
    mean_top3 = mean(top3_nbasis),
    .groups = "drop"
  )

print(dic_behavior %>% knitr::kable(digits = 1))

cat("\n")
