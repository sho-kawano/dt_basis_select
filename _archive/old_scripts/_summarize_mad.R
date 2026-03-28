library(tidyverse)

# Load metrics
metrics_raw <- readRDS("results_multi_config/methods_comparison/metrics_summary.RDS")

# Filter out ε=0.6
metrics <- metrics_raw %>%
  filter(is.na(epsilon) | abs(epsilon - 0.6) > 0.001) %>%
  mutate(
    method_label = case_when(
      method_type == "DT-MSE" ~ sprintf("DT-MSE ε=%.1f", epsilon),
      method_type == "DT-Likelihood" ~ sprintf("DT-Likelihood ε=%.1f", epsilon),
      TRUE ~ method_type
    )
  )

cat("=== MAD by Method and Design ===\n\n")
mad_by_design <- metrics %>%
  select(config, method_label, MAD) %>%
  pivot_wider(names_from = config, values_from = MAD)

print(mad_by_design)

cat("\n\n=== Average MAD and Consistency ===\n\n")
consistency <- metrics %>%
  group_by(method_label) %>%
  summarise(
    avg_MAD = mean(MAD),
    sd_MAD = sd(MAD),
    min_MAD = min(MAD),
    max_MAD = max(MAD),
    range_MAD = max(MAD) - min(MAD),
    .groups = "drop"
  ) %>%
  arrange(avg_MAD)

print(consistency)

cat("\n\n=== Rank by Design ===\n\n")
ranks <- metrics %>%
  group_by(config) %>%
  mutate(rank = rank(MAD, ties.method = "min")) %>%
  ungroup() %>%
  select(config, method_label, MAD, rank) %>%
  arrange(config, rank)

print(ranks)
