## Section 6 — averaged across 6 designs, full nbasis grid (no filtering)
library(tidyverse)
library(patchwork)
library(latex2exp)

CONFIGS <- c("equal_30", "equal_40", "equal_50",
             "equal_75", "equal_100", "equal_125")

dir.create("analysis/figure", showWarnings = FALSE, recursive = TRUE)

# ── Load (no nbasis filter) ─────────────────────────────────────────────────
dt_all <- readRDS("results_multi_config/dt_epsilon_8configs.RDS") %>%
  filter(config %in% CONFIGS, loss_function == "MSE")

oracle_all <- readRDS("results_multi_config/oracle_epsilon_8configs.RDS") %>%
  filter(config %in% CONFIGS)

# ── Selections ───────────────────────────────────────────────────────────────
oracle_sel <- oracle_all %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(config, comp_no, oracle_nbasis = nbasis)

dt_sel <- dt_all %>%
  group_by(config, comp_no, epsilon, n_reps_used) %>%
  slice_min(metric, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(config, comp_no, epsilon, n_reps_used, dt_nbasis = nbasis) %>%
  left_join(oracle_sel, by = c("config", "comp_no")) %>%
  mutate(deviation = dt_nbasis - oracle_nbasis)

# ── (a) MAD basis — averaged across designs ─────────────────────────────────
mad_data <- dt_sel %>%
  group_by(epsilon, n_reps_used) %>%
  summarise(mad = mean(abs(deviation)),
            se = sd(abs(deviation)) / sqrt(n()), .groups = "drop")

p_a <- ggplot(mad_data, aes(epsilon, mad,
                             color = factor(n_reps_used),
                             fill = factor(n_reps_used))) +
  geom_ribbon(aes(ymin = mad - se, ymax = mad + se), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_color_brewer(palette = "Set1", name = "R") +
  scale_fill_brewer(palette = "Set1", name = "R") +
  labs(title = "(a) MAD basis",
       x = expression(epsilon),
       y = TeX("Mean $|p - p^*|$")) +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(size = 11, face = "bold"))

# ── (b) Directional bias — averaged across designs ──────────────────────────
bias_data <- dt_sel %>%
  group_by(epsilon, n_reps_used) %>%
  summarise(mean_bias = mean(deviation),
            se = sd(deviation) / sqrt(n()), .groups = "drop")

all_nbasis <- sort(unique(oracle_all$nbasis))
random_bias <- mean(all_nbasis) - mean(oracle_sel$oracle_nbasis)

p_b <- ggplot(bias_data, aes(epsilon, mean_bias,
                              color = factor(n_reps_used),
                              fill = factor(n_reps_used))) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "gray50") +
  geom_hline(yintercept = random_bias, linetype = "dashed", linewidth = 0.5, color = "gray35") +
  annotate("text", x = 0.55, y = random_bias - 1, label = sprintf("random selection (+%.1f)", random_bias),
           size = 3, color = "gray35") +
  geom_ribbon(aes(ymin = mean_bias - se, ymax = mean_bias + se),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "(b) Directional bias",
       x = expression(epsilon),
       y = TeX("Mean $p - p^*$")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"))

# ── Combine ──────────────────────────────────────────────────────────────────
fig <- p_a + p_b + plot_layout(widths = c(1.15, 1))
ggsave("analysis/figure/test_s6_full_grid_avg.png", fig, width = 10, height = 4.5, dpi = 150)
cat("Saved test_s6_full_grid_avg.png\n")

# ── Print summaries ──────────────────────────────────────────────────────────
cat(sprintf("\nRandom bias: %.1f\n", random_bias))
cat("\nMAD (R=5):\n")
print(filter(mad_data, n_reps_used == 5))
cat("\nDirectional bias (R=5):\n")
print(filter(bias_data, n_reps_used == 5))
