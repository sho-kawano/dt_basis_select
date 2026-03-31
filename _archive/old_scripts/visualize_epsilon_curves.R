#!/usr/bin/env Rscript
# Visualize DT epsilon curves vs DIC/WAIC across 10 configs
# Shows percentage increase in MSE from oracle model

library(tidyverse)

# Load data
dt_all <- readRDS("results_multi_config/dt_all_10configs.RDS")
oracle_all <- readRDS("results_multi_config/oracle_all_10configs.RDS")

# Oracle selections (minimal MSE per comparison)
oracle_selections <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, oracle_nbasis = nbasis, oracle_mse = mse_true)

# === DT Oracle Penalties ===
dt_selections <- dt_all %>%
  filter(method_name == "dt_1fold", loss_function == "MSE") %>%
  group_by(config, comp_no, epsilon, n_reps_used) %>%
  slice_min(metric, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, epsilon, n_reps_used, dt_nbasis = nbasis)

dt_oracle_penalty <- dt_selections %>%
  left_join(oracle_all %>% select(config, comp_no, nbasis, mse_true),
            by = c("config", "comp_no", "dt_nbasis" = "nbasis")) %>%
  rename(dt_mse = mse_true) %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(pct_increase = (dt_mse - oracle_mse) / oracle_mse * 100) %>%
  group_by(config, epsilon, n_reps_used) %>%
  summarise(
    mean_pct_increase = mean(pct_increase),
    se_pct_increase = sd(pct_increase) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(method = "DT")

# === DIC Oracle Penalties ===
dic_selections <- oracle_all %>%
  filter(!is.na(DIC)) %>%
  group_by(config, comp_no) %>%
  slice_min(DIC, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, dic_nbasis = nbasis)

dic_oracle_penalty <- dic_selections %>%
  left_join(oracle_all %>% select(config, comp_no, nbasis, mse_true),
            by = c("config", "comp_no", "dic_nbasis" = "nbasis")) %>%
  rename(dic_mse = mse_true) %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(pct_increase = (dic_mse - oracle_mse) / oracle_mse * 100) %>%
  group_by(config) %>%
  summarise(
    mean_pct_increase = mean(pct_increase),
    se_pct_increase = sd(pct_increase) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(method = "DIC", epsilon = NA, n_reps_used = NA)

# === WAIC Oracle Penalties ===
waic_selections <- oracle_all %>%
  filter(!is.na(WAIC)) %>%
  group_by(config, comp_no) %>%
  slice_min(WAIC, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, waic_nbasis = nbasis)

waic_oracle_penalty <- waic_selections %>%
  left_join(oracle_all %>% select(config, comp_no, nbasis, mse_true),
            by = c("config", "comp_no", "waic_nbasis" = "nbasis")) %>%
  rename(waic_mse = mse_true) %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(pct_increase = (waic_mse - oracle_mse) / oracle_mse * 100) %>%
  group_by(config) %>%
  summarise(
    mean_pct_increase = mean(pct_increase),
    se_pct_increase = sd(pct_increase) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(method = "WAIC", epsilon = NA, n_reps_used = NA)

# Combine all
all_penalties <- bind_rows(
  dt_oracle_penalty,
  dic_oracle_penalty,
  waic_oracle_penalty
)

# === MAD Data ===
dt_mad_by_epsilon <- dt_selections %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(deviation = abs(dt_nbasis - oracle_nbasis)) %>%
  group_by(config, epsilon, n_reps_used) %>%
  summarise(MAD = mean(deviation), .groups = "drop") %>%
  mutate(method = "DT")

dic_mad_data <- dic_selections %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(deviation = abs(dic_nbasis - oracle_nbasis)) %>%
  group_by(config) %>%
  summarise(MAD = mean(deviation), .groups = "drop") %>%
  mutate(method = "DIC", epsilon = NA, n_reps_used = NA)

waic_mad_data <- waic_selections %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(deviation = abs(waic_nbasis - oracle_nbasis)) %>%
  group_by(config) %>%
  summarise(MAD = mean(deviation), .groups = "drop") %>%
  mutate(method = "WAIC", epsilon = NA, n_reps_used = NA)

all_mad <- bind_rows(dt_mad_by_epsilon, dic_mad_data, waic_mad_data)

# Compute sample sizes for ordering
ca_pop <- readRDS("data/ca_pums_population.rds")
puma_sizes <- ca_pop %>%
  group_by(PUMA) %>%
  summarise(pop = n(), .groups = "drop")

prop_configs <- tibble(
  config = c("prop_0.5pct", "prop_1pct", "prop_1p25pct", "prop_1p5pct", "prop_1p75pct", "prop_2pct"),
  samp_frac = c(0.005, 0.01, 0.0125, 0.015, 0.0175, 0.02)
)

sample_sizes <- puma_sizes %>%
  crossing(prop_configs) %>%
  mutate(sample_n = pmax(10, round(pop * samp_frac))) %>%
  group_by(config) %>%
  summarise(mean_n = mean(sample_n), .groups = "drop")

equal_configs <- tibble(
  config = c("equal_40", "equal_50", "equal_75", "equal_100"),
  mean_n = c(40, 50, 75, 100)
)

all_sample_sizes <- bind_rows(sample_sizes, equal_configs) %>%
  arrange(mean_n)

config_levels <- all_sample_sizes %>% pull(config)

# === Plot 1: MAD - DT n=5 only ===
p1 <- all_mad %>%
  filter(method == "DT" & n_reps_used == 5 | method %in% c("DIC", "WAIC")) %>%
  mutate(config = factor(config, levels = config_levels)) %>%
  ggplot(aes(x = epsilon, y = MAD)) +
  geom_hline(data = all_mad %>% filter(method == "DIC") %>%
               mutate(config = factor(config, levels = config_levels)),
             aes(yintercept = MAD, color = "DIC"),
             linetype = "dashed", linewidth = 0.8) +
  geom_hline(data = all_mad %>% filter(method == "WAIC") %>%
               mutate(config = factor(config, levels = config_levels)),
             aes(yintercept = MAD, color = "WAIC"),
             linetype = "dashed", linewidth = 0.8) +
  geom_line(data = all_mad %>% filter(method == "DT", n_reps_used == 5) %>%
              mutate(config = factor(config, levels = config_levels)),
            aes(color = "DT n=5"), linewidth = 1.2) +
  geom_point(data = all_mad %>% filter(method == "DT", n_reps_used == 5) %>%
               mutate(config = factor(config, levels = config_levels)),
             aes(color = "DT n=5"), size = 2.5) +
  facet_wrap(~config, ncol = 3, scales = "free_y") +
  scale_color_manual(
    values = c("DT n=5" = "steelblue", "DIC" = "coral", "WAIC" = "darkgoldenrod2"),
    breaks = c("DT n=5", "DIC", "WAIC")
  ) +
  labs(
    title = "Method Performance: Mean Absolute Deviation from Oracle Optimal",
    subtitle = "DT with n=5 averaging (solid line) vs DIC/WAIC (dashed lines). Lower is better. Ordered by sample size.",
    x = "Epsilon (ε)",
    y = "MAD (nbasis units)",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("analysis/figure/epsilon_curves_mad.png", p1, width = 12, height = 10, dpi = 300)
cat("✓ Saved: analysis/figure/epsilon_curves_mad.png\n")

# === Plot 2: Oracle MSE Penalty - DT n=5 only ===
p2 <- all_penalties %>%
  filter(method == "DT" & n_reps_used == 5 | method %in% c("DIC", "WAIC")) %>%
  mutate(config = factor(config, levels = config_levels)) %>%
  ggplot(aes(x = epsilon, y = mean_pct_increase)) +
  geom_hline(data = all_penalties %>% filter(method == "DIC") %>%
               mutate(config = factor(config, levels = config_levels)),
             aes(yintercept = mean_pct_increase, color = "DIC"),
             linetype = "dashed", linewidth = 0.8) +
  geom_hline(data = all_penalties %>% filter(method == "WAIC") %>%
               mutate(config = factor(config, levels = config_levels)),
             aes(yintercept = mean_pct_increase, color = "WAIC"),
             linetype = "dashed", linewidth = 0.8) +
  geom_line(data = all_penalties %>% filter(method == "DT", n_reps_used == 5) %>%
              mutate(config = factor(config, levels = config_levels)),
            aes(color = "DT n=5"), linewidth = 1.2) +
  geom_point(data = all_penalties %>% filter(method == "DT", n_reps_used == 5) %>%
               mutate(config = factor(config, levels = config_levels)),
             aes(color = "DT n=5"), size = 2.5) +
  facet_wrap(~config, ncol = 3, scales = "free_y") +
  scale_color_manual(
    values = c("DT n=5" = "steelblue", "DIC" = "coral", "WAIC" = "darkgoldenrod2"),
    breaks = c("DT n=5", "DIC", "WAIC")
  ) +
  labs(
    title = "Method Performance: Percentage Increase in MSE from Oracle Model",
    subtitle = "DT with n=5 averaging (solid line) vs DIC/WAIC (dashed lines). Lower is better. Ordered by sample size.",
    x = "Epsilon (ε)",
    y = "Percentage Increase in MSE (%)",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("analysis/figure/epsilon_curves_pct_mse.png", p2, width = 12, height = 10, dpi = 300)
cat("✓ Saved: analysis/figure/epsilon_curves_pct_mse.png\n")

# === Summary Statistics ===
cat("\n=== DT n=5: Best epsilon per config ===\n")
best_dt_n5 <- all_penalties %>%
  filter(method == "DT", n_reps_used == 5) %>%
  group_by(config) %>%
  slice_min(mean_pct_increase, n = 1, with_ties = FALSE) %>%
  select(config, best_epsilon = epsilon, best_pct = mean_pct_increase)

print(best_dt_n5, n = 20)

cat("\n=== DT n=5 vs DIC: Mean penalty across configs ===\n")
comparison <- all_penalties %>%
  filter((method == "DT" & n_reps_used == 5) | method == "DIC") %>%
  group_by(method, epsilon) %>%
  summarise(
    mean_pct = mean(mean_pct_increase),
    .groups = "drop"
  ) %>%
  arrange(mean_pct)

print(comparison, n = 20)
