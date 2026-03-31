library(tidyverse)

# Load aggregated results (matching the notebook)
dt_all <- readRDS("results_multi_config/methods_comparison/dt_results.RDS")
oracle_all <- readRDS("results_multi_config/methods_comparison/oracle_results.RDS")
esim_all <- readRDS("results_multi_config/methods_comparison/esim_results.RDS")
metrics_raw <- readRDS("results_multi_config/methods_comparison/metrics_summary.RDS")
selections <- readRDS("results_multi_config/methods_comparison/method_selections.RDS")

# Filter out ε=0.6 for cleaner analysis
metrics <- metrics_raw %>%
  filter(is.na(epsilon) | abs(epsilon - 0.6) > 0.001)

# Get oracle selections (best model per comparison)
oracle_selections <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, oracle_model = model, oracle_nbasis = nbasis) %>%
  ungroup()

# Compute deviations from oracle for each comparison
deviations <- selections %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  filter(is.na(epsilon) | abs(epsilon - 0.6) > 0.001) %>%  # Filter out ε=0.6
  mutate(
    deviation = selected_nbasis - oracle_nbasis,
    method_type = case_when(
      str_detect(method, "DT-MSE") ~ "DT-MSE",
      str_detect(method, "DT-plugin_NLL") ~ "DT-Likelihood",
      TRUE ~ method
    ),
    method_label = case_when(
      str_detect(method, "DT-MSE") ~ sprintf("DT-MSE ε=%.1f", epsilon),
      str_detect(method, "DT-plugin_NLL") ~ sprintf("DT-Likelihood ε=%.1f", epsilon),
      TRUE ~ method
    ),
    config_label = case_when(
      config == "prop_0.75pct" ~ "0.75% PA",
      config == "prop_1p25pct" ~ "1.25% PA",
      config == "prop_1p75pct" ~ "1.75% PA"
    )
  )

cat("Filtered deviations summary by method:\n")
deviations %>%
  group_by(method_label) %>%
  summarise(
    n = n(),
    mean_dev = mean(deviation),
    median_dev = median(deviation),
    sd_dev = sd(deviation)
  ) %>%
  print()

cat("\nCreating plot...\n")

# Create the plot
p <- ggplot(deviations, aes(x = method_label, y = deviation)) +
  geom_boxplot(aes(fill = method_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~config_label, ncol = 3) +
  labs(
    title = "Distribution of Deviations from Oracle Basis",
    subtitle = "Red line = perfect agreement with oracle",
    x = "Method",
    y = "Deviation (selected - oracle basis)",
    fill = "Method Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  scale_fill_brewer(palette = "Set2")

ggsave("_test_notebook_chunk.png", p, width = 10, height = 4)
cat("Plot saved to _test_notebook_chunk.png\n")
