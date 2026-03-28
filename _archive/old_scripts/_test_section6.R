## Temporary test script for Section 6 plots
## Tests each plot independently, saves to analysis/figure/

library(tidyverse)

# ── Config ──────────────────────────────────────────────────────────────────
# Paper focuses on 4 candidate models; full grid available for robustness
PAPER_NBASIS <- c(6, 18, 30, 42)
USE_PAPER_MODELS <- TRUE  # toggle to FALSE for full 20-model grid

CONFIGS <- c("equal_30", "equal_40", "equal_50",
             "equal_75", "equal_100", "equal_125")

dir.create("analysis/figure", showWarnings = FALSE, recursive = TRUE)

# ── Load data ───────────────────────────────────────────────────────────────
dt_all <- readRDS("results_multi_config/dt_epsilon_8configs.RDS") %>%
  filter(config %in% CONFIGS, loss_function == "MSE")

oracle_all <- readRDS("results_multi_config/oracle_epsilon_8configs.RDS") %>%
  filter(config %in% CONFIGS)

thinned_oracle <- readRDS("results_multi_config/thinned_oracle_6configs.RDS") %>%
  filter(config %in% CONFIGS)

if (USE_PAPER_MODELS) {
  dt_all <- filter(dt_all, nbasis %in% PAPER_NBASIS)
  oracle_all <- filter(oracle_all, nbasis %in% PAPER_NBASIS)
  thinned_oracle <- filter(thinned_oracle, nbasis %in% PAPER_NBASIS)
}

cat(sprintf("DT rows: %d | Oracle rows: %d | Thinned oracle rows: %d\n",
            nrow(dt_all), nrow(oracle_all), nrow(thinned_oracle)))

# ── Helper: DT and oracle selections ───────────────────────────────────────
oracle_sel <- oracle_all %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(config, comp_no, oracle_nbasis = nbasis, oracle_mse = mse_true)

dt_sel <- dt_all %>%
  group_by(config, comp_no, epsilon, n_reps_used) %>%
  slice_min(metric, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(config, comp_no, epsilon, n_reps_used, dt_nbasis = nbasis)

# Join DT selection to oracle selection + oracle MSE for each selected model
sel <- dt_sel %>%
  left_join(oracle_sel, by = c("config", "comp_no")) %>%
  left_join(
    oracle_all %>% select(config, comp_no, nbasis, mse_true),
    by = c("config", "comp_no", "dt_nbasis" = "nbasis")
  ) %>%
  rename(dt_mse = mse_true) %>%
  mutate(
    deviation = dt_nbasis - oracle_nbasis,
    abs_deviation = abs(deviation),
    mse_penalty = dt_mse - oracle_mse
  )

# ═══════════════════════════════════════════════════════════════════════════
# PLOT 1: Fundamental tension (gap vs variance)
# Uses thinned oracle for gap; DT metric variance for variance
# ═══════════════════════════════════════════════════════════════════════════
cat("\n── Plot 1: Fundamental Tension ──\n")

gap_data <- thinned_oracle %>%
  filter(n_reps_used == 5) %>%
  left_join(oracle_all %>% select(comp_no, config, nbasis, mse_true),
            by = c("comp_no", "config", "nbasis")) %>%
  mutate(gap_pct = 100 * (thinned_mse - mse_true) / mse_true)

# Average gap across configs/comps/models at each epsilon
gap_by_eps <- gap_data %>%
  group_by(epsilon) %>%
  summarise(mean_gap_pct = mean(gap_pct, na.rm = TRUE), .groups = "drop")

# Variance of DT metric (across comps within each epsilon), averaged over configs
var_by_eps <- dt_all %>%
  filter(n_reps_used == 5) %>%
  group_by(config, epsilon, nbasis) %>%
  summarise(v = var(metric), .groups = "drop") %>%
  group_by(epsilon) %>%
  summarise(mean_variance = mean(v), .groups = "drop")

tension <- gap_by_eps %>%
  left_join(var_by_eps, by = "epsilon") %>%
  mutate(
    gap_norm = (mean_gap_pct - min(mean_gap_pct)) /
               (max(mean_gap_pct) - min(mean_gap_pct)),
    var_norm = (mean_variance - min(mean_variance)) /
               (max(mean_variance) - min(mean_variance))
  )

p1 <- tension %>%
  pivot_longer(c(gap_norm, var_norm), names_to = "metric_type") %>%
  ggplot(aes(epsilon, value, color = metric_type)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(gap_norm = "coral", var_norm = "steelblue"),
    labels = c("Oracle gap", "Estimation variance")
  ) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Fundamental Tension: Gap vs Variance",
       subtitle = sprintf("%d equal allocation configs, R=5", length(CONFIGS)),
       x = expression(epsilon), y = "Normalized [0, 1]", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("analysis/figure/test_s6_tension.png", p1, width = 8, height = 5, dpi = 150)
cat("  Saved test_s6_tension.png\n")
print(tension %>% select(epsilon, mean_gap_pct, mean_variance))


# ═══════════════════════════════════════════════════════════════════════════
# PLOT 2: MAD basis vs epsilon (NEW — primary metric per outline)
#   Separate lines for R = 1, 3, 5
# ═══════════════════════════════════════════════════════════════════════════
cat("\n── Plot 2: MAD Basis vs Epsilon ──\n")

mad_data <- sel %>%
  group_by(epsilon, n_reps_used) %>%
  summarise(
    mad = mean(abs_deviation),
    se = sd(abs_deviation) / sqrt(n()),
    .groups = "drop"
  )

p2 <- ggplot(mad_data, aes(epsilon, mad,
                            color = factor(n_reps_used),
                            fill = factor(n_reps_used))) +
  geom_ribbon(aes(ymin = mad - se, ymax = mad + se), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "MAD Basis vs Epsilon",
       subtitle = "Average |selected - oracle| across configs and samples",
       x = expression(epsilon), y = "MAD (# basis functions)",
       color = "R (repeats)", fill = "R (repeats)") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("analysis/figure/test_s6_mad.png", p2, width = 8, height = 5, dpi = 150)
cat("  Saved test_s6_mad.png\n")
print(mad_data %>% arrange(n_reps_used, epsilon))


# ═══════════════════════════════════════════════════════════════════════════
# PLOT 3: Directional bias vs epsilon (NEW)
#   Positive = over-selection (too complex), Negative = under-selection
# ═══════════════════════════════════════════════════════════════════════════
cat("\n── Plot 3: Directional Bias ──\n")

bias_data <- sel %>%
  group_by(epsilon, n_reps_used) %>%
  summarise(
    mean_bias = mean(deviation),
    se = sd(deviation) / sqrt(n()),
    .groups = "drop"
  )

p3 <- ggplot(bias_data, aes(epsilon, mean_bias,
                             color = factor(n_reps_used),
                             fill = factor(n_reps_used))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_ribbon(aes(ymin = mean_bias - se, ymax = mean_bias + se),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Directional Bias vs Epsilon",
       subtitle = "Mean(selected - oracle); negative = under-selection (favors simpler)",
       x = expression(epsilon), y = "Directional bias (# basis functions)",
       color = "R (repeats)", fill = "R (repeats)") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("analysis/figure/test_s6_directional_bias.png", p3, width = 8, height = 5, dpi = 150)
cat("  Saved test_s6_directional_bias.png\n")
print(bias_data %>% arrange(n_reps_used, epsilon))


# ═══════════════════════════════════════════════════════════════════════════
# PLOT 4: Differential gap (complex - simple) vs epsilon
#   Existing plot from Rmd — connects to Corollary 3.4
# ═══════════════════════════════════════════════════════════════════════════
cat("\n── Plot 4: Differential Gap ──\n")

if (USE_PAPER_MODELS) {
  simple_p <- min(PAPER_NBASIS)
  complex_p <- max(PAPER_NBASIS)
} else {
  simple_p <- 6; complex_p <- 42
}

diff_gap <- gap_data %>%
  filter(nbasis %in% c(simple_p, complex_p)) %>%
  select(config, comp_no, epsilon, nbasis, gap_pct) %>%
  pivot_wider(names_from = nbasis, values_from = gap_pct,
              names_prefix = "p") %>%
  mutate(diff = .data[[paste0("p", complex_p)]] - .data[[paste0("p", simple_p)]]) %>%
  group_by(epsilon) %>%
  summarise(mean_diff = mean(diff, na.rm = TRUE),
            se = sd(diff, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

p4 <- ggplot(diff_gap, aes(epsilon, mean_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_ribbon(aes(ymin = mean_diff - se, ymax = mean_diff + se),
              alpha = 0.2, fill = "steelblue") +
  geom_line(linewidth = 1.2, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
  labs(title = sprintf("Differential Gap: p=%d (complex) minus p=%d (simple)",
                       complex_p, simple_p),
       subtitle = "Crosses zero where model rankings stabilize",
       x = expression(epsilon), y = "Gap difference (% of oracle MSE)") +
  theme_minimal()

ggsave("analysis/figure/test_s6_diff_gap.png", p4, width = 8, height = 5, dpi = 150)
cat("  Saved test_s6_diff_gap.png\n")
print(diff_gap)


# ═══════════════════════════════════════════════════════════════════════════
# PLOT 5: Selection accuracy vs epsilon, by R
# ═══════════════════════════════════════════════════════════════════════════
cat("\n── Plot 5: Selection Accuracy ──\n")

acc_data <- sel %>%
  mutate(exact = dt_nbasis == oracle_nbasis) %>%
  group_by(epsilon, n_reps_used) %>%
  summarise(pct_exact = 100 * mean(exact), .groups = "drop")

p5 <- ggplot(acc_data, aes(epsilon, pct_exact,
                            color = factor(n_reps_used))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Exact Selection Accuracy vs Epsilon",
       subtitle = "% of samples where DT selects oracle-optimal model",
       x = expression(epsilon), y = "Exact match (%)",
       color = "R (repeats)") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("analysis/figure/test_s6_accuracy.png", p5, width = 8, height = 5, dpi = 150)
cat("  Saved test_s6_accuracy.png\n")
print(acc_data %>% pivot_wider(names_from = n_reps_used, values_from = pct_exact,
                                names_prefix = "R="))


# ═══════════════════════════════════════════════════════════════════════════
# PLOT 6: Oracle MSE penalty of DT selection vs epsilon, by R
#   "How much worse is the model DT picks vs the oracle-best?"
# ═══════════════════════════════════════════════════════════════════════════
cat("\n── Plot 6: Oracle MSE Penalty ──\n")

penalty_data <- sel %>%
  group_by(epsilon, n_reps_used) %>%
  summarise(
    mean_penalty = mean(mse_penalty),
    se = sd(mse_penalty) / sqrt(n()),
    .groups = "drop"
  )

p6 <- ggplot(penalty_data, aes(epsilon, mean_penalty,
                                color = factor(n_reps_used),
                                fill = factor(n_reps_used))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_ribbon(aes(ymin = mean_penalty - se, ymax = mean_penalty + se),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Oracle MSE Penalty of DT Selection",
       subtitle = "MSE(DT pick) - MSE(oracle pick); lower is better",
       x = expression(epsilon), y = "MSE penalty",
       color = "R (repeats)", fill = "R (repeats)") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("analysis/figure/test_s6_mse_penalty.png", p6, width = 8, height = 5, dpi = 150)
cat("  Saved test_s6_mse_penalty.png\n")
print(penalty_data %>% arrange(n_reps_used, epsilon))

cat("\n── All plots saved to analysis/figure/test_s6_*.png ──\n")
