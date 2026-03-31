library(tidyverse)

# Load data
oracle_all <- readRDS("results_multi_config/methods_comparison/oracle_results.RDS")
metrics_raw <- readRDS("results_multi_config/methods_comparison/metrics_summary.RDS")

# Filter out ε=0.6
metrics <- metrics_raw %>%
  filter(is.na(epsilon) | abs(epsilon - 0.6) > 0.001)

cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("QUESTION 1: Oracle Basis Distribution\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

# Get oracle selections for each sample
oracle_selections <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, oracle_nbasis = nbasis, oracle_mse = mse_true) %>%
  ungroup()

# Oracle basis statistics by design
oracle_stats <- oracle_selections %>%
  group_by(config) %>%
  summarise(
    n_samples = n(),
    mean_oracle = mean(oracle_nbasis),
    sd_oracle = sd(oracle_nbasis),
    min_oracle = min(oracle_nbasis),
    max_oracle = max(oracle_nbasis),
    range_oracle = max(oracle_nbasis) - min(oracle_nbasis),
    cv_oracle = sd(oracle_nbasis) / mean(oracle_nbasis)  # Coefficient of variation
  )

print(oracle_stats)

cat("\n\nOracle basis distribution by design:\n")
for (cfg in unique(oracle_selections$config)) {
  vals <- oracle_selections %>%
    filter(config == cfg) %>%
    pull(oracle_nbasis) %>%
    table()
  cat(sprintf("\n%s:\n", cfg))
  print(vals)
}

cat("\n\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("QUESTION 2: MAD-Penalty Relationship\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

# Create method labels
metrics_labeled <- metrics %>%
  mutate(
    method_label = case_when(
      method_type == "DT-MSE" ~ sprintf("DT-MSE ε=%.1f", epsilon),
      method_type == "DT-Likelihood" ~ sprintf("DT-Likelihood ε=%.1f", epsilon),
      TRUE ~ method_type
    )
  )

# Overall correlation
overall_cor <- cor(metrics_labeled$MAD, metrics_labeled$oracle_penalty,
                   use = "complete.obs")

cat(sprintf("Overall correlation (MAD vs oracle_penalty): r = %.3f\n\n", overall_cor))

# Linear model
lm_fit <- lm(oracle_penalty ~ MAD, data = metrics_labeled)
cat("Linear model: oracle_penalty ~ MAD\n")
print(summary(lm_fit))

cat("\n\nInterpretation:\n")
slope <- coef(lm_fit)[2]
cat(sprintf("Slope: %.6f\n", slope))
cat(sprintf("Each additional unit of MAD (1 basis function off on average) \n"))
cat(sprintf("costs approximately %.2f%% in oracle penalty (MSE increase).\n\n", slope))

# By design
cat("Correlation by design:\n")
cor_by_design <- metrics_labeled %>%
  group_by(config) %>%
  summarise(
    correlation = cor(MAD, oracle_penalty, use = "complete.obs"),
    n = n()
  )
print(cor_by_design)

cat("\n\nSlope by design (oracle_penalty ~ MAD):\n")
slopes_by_design <- metrics_labeled %>%
  group_by(config) %>%
  summarise(
    slope = coef(lm(oracle_penalty ~ MAD))[2],
    intercept = coef(lm(oracle_penalty ~ MAD))[1],
    r_squared = summary(lm(oracle_penalty ~ MAD))$r.squared
  )
print(slopes_by_design)

# Scatterplot data
cat("\n\nScatterplot summary (MAD vs oracle_penalty):\n")
metrics_labeled %>%
  select(config, method_label, MAD, oracle_penalty) %>%
  arrange(config, MAD) %>%
  print(n = 21)

# Create visualization
p <- ggplot(metrics_labeled, aes(x = MAD, y = oracle_penalty)) +
  geom_point(aes(color = config), size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  facet_wrap(~config, ncol = 3, scales = "free") +
  labs(
    title = "MAD vs Oracle MSE Penalty",
    subtitle = sprintf("Overall correlation: r = %.3f | Slope: %.3f%% per MAD unit",
                      overall_cor, slope),
    x = "MAD (mean absolute deviation from oracle basis)",
    y = "Oracle penalty (%)",
    color = "Design"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("_mad_penalty_relationship.png", p, width = 10, height = 4)
cat("\n\nPlot saved to _mad_penalty_relationship.png\n")
