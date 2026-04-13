#!/usr/bin/env Rscript
# ==============================================================================
# oracle_margin_analysis.R — Oracle margin and plateau analysis
# ==============================================================================
# Seeds 1-50, 3 PA designs (trio: 0.75% / 1.25% / 1.75%)
# Questions:
#   1. How "decisive" is the oracle? (margin between 1st and 2nd best model)
#   2. How wide is the near-optimal plateau? (range of models within X% of best)
#   3. Do method rankings change when filtering out ambiguous samples?
# ==============================================================================

library(tidyverse)

SEEDS   <- 1:50
EPS     <- 0.7
REPS    <- 5

# ==============================================================================
# PART 0: Load oracle data for trio designs
# ==============================================================================

cat("Loading trio from methodcomp_results.RDS...\n")
mc <- readRDS("results_summary/methodcomp_results.RDS")

# Trio oracle: all 20 per-model MSE values per comp_no × config
oracle_all <- mc$oracle %>% select(comp_no, config, model, nbasis, mse_true, DIC, WAIC)
dt_all     <- mc$dt
esim_all   <- mc$esim

cat(sprintf("\nOracle data: %d rows across %d designs\n",
            nrow(oracle_all), n_distinct(oracle_all$config)))

# ==============================================================================
# PART 1: Compute oracle margins and plateau width
# ==============================================================================

cat("\n=== PART 1: Oracle margins and plateaus ===\n")

margin_df <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  arrange(mse_true) %>%
  mutate(rank = row_number()) %>%
  summarise(
    oracle_nbasis    = nbasis[rank == 1],
    best_mse         = mse_true[rank == 1],
    second_mse       = mse_true[rank == 2],
    margin_1v2       = mse_true[rank == 2] - mse_true[rank == 1],
    # Plateau: how many models within 5% of best MSE?
    n_within_5pct    = sum(mse_true <= best_mse * 1.05),
    n_within_10pct   = sum(mse_true <= best_mse * 1.10),
    # Plateau width in nbasis units: range of models within 5% of best
    plateau_min_5pct = min(nbasis[mse_true <= best_mse * 1.05]),
    plateau_max_5pct = max(nbasis[mse_true <= best_mse * 1.05]),
    plateau_width_5pct  = plateau_max_5pct - plateau_min_5pct,
    plateau_min_10pct = min(nbasis[mse_true <= best_mse * 1.10]),
    plateau_max_10pct = max(nbasis[mse_true <= best_mse * 1.10]),
    plateau_width_10pct = plateau_max_10pct - plateau_min_10pct,
    .groups = "drop"
  )

# Summary stats per design
cat("\n--- Margin summary (margin_1v2 = 2nd best - best MSE) ---\n")
margin_df %>%
  group_by(config) %>%
  summarise(
    mean_margin   = round(mean(margin_1v2), 6),
    median_margin = round(median(margin_1v2), 6),
    sd_margin     = round(sd(margin_1v2), 6),
    q10           = round(quantile(margin_1v2, 0.10), 6),
    q25           = round(quantile(margin_1v2, 0.25), 6),
    q75           = round(quantile(margin_1v2, 0.75), 6),
    .groups = "drop"
  ) %>% print(width = Inf)

cat("\n--- Oracle basis selection ---\n")
margin_df %>%
  group_by(config) %>%
  summarise(
    mean_oracle = round(mean(oracle_nbasis), 1),
    median_oracle = median(oracle_nbasis),
    sd_oracle = round(sd(oracle_nbasis), 1),
    .groups = "drop"
  ) %>% print()

cat("\n--- Grand mean oracle nbasis (across all 5 designs) ---\n")
cat(sprintf("  %.1f\n", mean(margin_df$oracle_nbasis)))

cat("\n--- Plateau width (5% tolerance) ---\n")
margin_df %>%
  group_by(config) %>%
  summarise(
    mean_n_models   = round(mean(n_within_5pct), 1),
    mean_width_nbasis = round(mean(plateau_width_5pct), 1),
    median_width    = median(plateau_width_5pct),
    .groups = "drop"
  ) %>% print()

cat("\n--- Plateau width (10% tolerance) ---\n")
margin_df %>%
  group_by(config) %>%
  summarise(
    mean_n_models   = round(mean(n_within_10pct), 1),
    mean_width_nbasis = round(mean(plateau_width_10pct), 1),
    median_width    = median(plateau_width_10pct),
    .groups = "drop"
  ) %>% print()

# ==============================================================================
# PART 2: Distribution plots
# ==============================================================================

cat("\n=== PART 2: Distribution plots ===\n")

# Order configs by sample size
config_order <- c("prop_0.75pct", "prop_1p25pct", "prop_1p5pct", "prop_1p75pct", "prop_2p25pct")
config_labels <- c("0.75%", "1.25%", "1.5%", "1.75%", "2.25%")
margin_df$config_f <- factor(margin_df$config, levels = config_order, labels = config_labels)

# 2a: Histogram of margin_1v2
p1 <- ggplot(margin_df, aes(x = margin_1v2)) +
  geom_histogram(bins = 25, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~config_f, scales = "free_y", nrow = 1) +
  labs(x = "Oracle margin (2nd best - best MSE)", y = "Count",
       title = "Distribution of oracle margins by design") +
  theme_minimal()
ggsave("p1_analysis/plot_oracle_margin_hist.png", p1, width = 12, height = 3.5)
cat("  Saved: analysis/plot_oracle_margin_hist.png\n")

# 2b: Log-scale version
p2 <- p1 + scale_x_log10() +
  labs(x = "Oracle margin (log scale)", title = "Oracle margins (log scale)")
ggsave("p1_analysis/plot_oracle_margin_hist_log.png", p2, width = 12, height = 3.5)
cat("  Saved: analysis/plot_oracle_margin_hist_log.png\n")

# 2c: Plateau width distribution
p3 <- ggplot(margin_df, aes(x = plateau_width_5pct)) +
  geom_histogram(bins = 20, fill = "darkorange", alpha = 0.7) +
  facet_wrap(~config_f, scales = "free_y", nrow = 1) +
  labs(x = "Plateau width (nbasis range within 5% of best)",
       y = "Count",
       title = "Width of near-optimal plateau by design") +
  theme_minimal()
ggsave("p1_analysis/plot_plateau_width.png", p3, width = 12, height = 3.5)
cat("  Saved: analysis/plot_plateau_width.png\n")

# 2d: Scatter of margin vs plateau width
p4 <- ggplot(margin_df, aes(x = margin_1v2, y = plateau_width_5pct)) +
  geom_point(alpha = 0.4) +
  facet_wrap(~config_f, nrow = 1, scales = "free") +
  labs(x = "Oracle margin (1v2)", y = "Plateau width (5%)",
       title = "Margin vs plateau width") +
  theme_minimal()
ggsave("p1_analysis/plot_margin_vs_plateau.png", p4, width = 12, height = 3.5)
cat("  Saved: analysis/plot_margin_vs_plateau.png\n")

# ==============================================================================
# PART 3: Method rankings under filtering
# ==============================================================================

cat("\n=== PART 3: Method rankings under margin filtering ===\n")

# Build method selections for all 5 designs
oracle_sel <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, oracle_nbasis = nbasis)

dic_sel <- oracle_all %>% filter(!is.na(DIC)) %>%
  group_by(config, comp_no) %>% slice_min(DIC, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, selected_nbasis = nbasis) %>% mutate(method = "DIC")

waic_sel <- oracle_all %>% filter(!is.na(WAIC)) %>%
  group_by(config, comp_no) %>% slice_min(WAIC, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, selected_nbasis = nbasis) %>% mutate(method = "WAIC")

nll_sel <- dt_all %>%
  filter(loss_function == "plugin_NLL", abs(epsilon - EPS) < 0.001) %>%
  group_by(config, comp_no) %>% slice_min(metric, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, selected_nbasis = nbasis) %>% mutate(method = "DT-NLL")

mse_sel <- dt_all %>%
  filter(loss_function == "MSE", abs(epsilon - EPS) < 0.001) %>%
  group_by(config, comp_no) %>% slice_min(metric, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, selected_nbasis = nbasis) %>% mutate(method = "DT-MSE")

esim_sel <- esim_all %>%
  filter(metric_type == "MSE", !is.na(nbasis)) %>%
  group_by(config, comp_no) %>% slice_min(metric, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, selected_nbasis = nbasis) %>% mutate(method = "ESIM")

all_sel <- bind_rows(dic_sel, waic_sel, nll_sel, mse_sel, esim_sel) %>%
  left_join(oracle_sel, by = c("config", "comp_no")) %>%
  mutate(deviation = selected_nbasis - oracle_nbasis) %>%
  left_join(margin_df %>% select(config, comp_no, margin_1v2, plateau_width_5pct),
            by = c("config", "comp_no"))

# Compute filtered MAE at various margin quantile thresholds
# Use within-design quantiles
threshold_results <- list()

for (q in c(0, 0.10, 0.25, 0.50)) {
  filtered <- all_sel %>%
    group_by(config) %>%
    filter(margin_1v2 >= quantile(margin_1v2, q)) %>%
    ungroup()

  res <- filtered %>%
    group_by(config, method) %>%
    summarise(
      MAE  = round(mean(abs(deviation)), 2),
      bias = round(mean(deviation), 1),
      n    = n(),
      .groups = "drop"
    ) %>%
    mutate(threshold = sprintf("q%.0f", q * 100))

  threshold_results[[as.character(q)]] <- res
}

threshold_all <- bind_rows(threshold_results)
threshold_all$config_f <- factor(threshold_all$config, levels = config_order, labels = config_labels)

# Print tables
for (q_label in unique(threshold_all$threshold)) {
  cat(sprintf("\n--- %s threshold (samples retained per design) ---\n", q_label))
  threshold_all %>%
    filter(threshold == q_label) %>%
    select(config_f, method, MAE, bias, n) %>%
    pivot_wider(names_from = method, values_from = c(MAE, bias),
                names_glue = "{method}_{.value}") %>%
    print(width = Inf)
}

# Also filter by plateau width
cat("\n--- Filtering by plateau width (keep only wide-plateau samples) ---\n")
for (q in c(0, 0.25, 0.50)) {
  filtered <- all_sel %>%
    group_by(config) %>%
    filter(plateau_width_5pct >= quantile(plateau_width_5pct, q)) %>%
    ungroup()

  cat(sprintf("\nPlateau width >= q%.0f:\n", q * 100))
  filtered %>%
    group_by(config, method) %>%
    summarise(MAE = round(mean(abs(deviation)), 2),
              bias = round(mean(deviation), 1),
              n = n(), .groups = "drop") %>%
    mutate(config_f = factor(config, levels = config_order, labels = config_labels)) %>%
    select(config_f, method, MAE, n) %>%
    pivot_wider(names_from = method, values_from = c(MAE),
                names_glue = "{method}") %>%
    print(width = Inf)
}

# ==============================================================================
# PART 4: Summary
# ==============================================================================

cat("\n\n=== SUMMARY ===\n")
cat("Key numbers to check:\n")
cat(sprintf("  Total samples: %d (5 designs x 50 seeds)\n", nrow(margin_df)))
cat(sprintf("  Grand mean oracle nbasis: %.1f\n", mean(margin_df$oracle_nbasis)))
cat(sprintf("  Grand median margin (1v2): %.6f\n", median(margin_df$margin_1v2)))
cat(sprintf("  Grand mean plateau width (5%%): %.1f nbasis\n", mean(margin_df$plateau_width_5pct)))

cat("\nDone.\n")
