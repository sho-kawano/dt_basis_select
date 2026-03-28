## Section 6.2 — one combined figure: (a) MAD basis, (b) Directional bias
library(tidyverse)
library(patchwork)

PAPER_NBASIS <- c(6, 18, 30, 42)
CONFIGS <- c("equal_30", "equal_40", "equal_50",
             "equal_75", "equal_100", "equal_125")

dir.create("analysis/figure", showWarnings = FALSE, recursive = TRUE)

# ── Load & filter ───────────────────────────────────────────────────────────
dt_all <- readRDS("results_multi_config/dt_epsilon_8configs.RDS") %>%
  filter(config %in% CONFIGS, loss_function == "MSE", nbasis %in% PAPER_NBASIS)

oracle_all <- readRDS("results_multi_config/oracle_epsilon_8configs.RDS") %>%
  filter(config %in% CONFIGS, nbasis %in% PAPER_NBASIS)

# ── Selections ──────────────────────────────────────────────────────────────
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

# ── (a) MAD basis — all R values ───────────────────────────────────────────
mad_data <- dt_sel %>%
  group_by(epsilon, n_reps_used) %>%
  summarise(mad = mean(abs(deviation)),
            se = sd(abs(deviation)) / sqrt(n()), .groups = "drop")

shared_theme <- theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"))

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
       y = "Mean |selected − oracle| basis") +
  shared_theme +
  theme(legend.position = "right")

# ── (b) Directional bias — all R values ───────────────────────────────────
bias_data <- dt_sel %>%
  group_by(epsilon, n_reps_used) %>%
  summarise(mean_bias = mean(deviation),
            se = sd(deviation) / sqrt(n()), .groups = "drop")

# Random selection baseline: E[Uniform{6,18,30,42}] - E[oracle]
random_bias <- mean(PAPER_NBASIS) - mean(
  (oracle_sel %>% group_by(config, comp_no) %>% slice(1) %>% ungroup())$oracle_nbasis
)

p_b <- ggplot(bias_data, aes(epsilon, mean_bias,
                              color = factor(n_reps_used),
                              fill = factor(n_reps_used))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = random_bias, linetype = "dashed", color = "gray40") +
  annotate("text", x = 0.55, y = random_bias - 0.7,
           label = "random selection", size = 3, color = "gray30") +
  geom_ribbon(aes(ymin = mean_bias - se, ymax = mean_bias + se),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "(b) Directional bias",
       x = expression(epsilon),
       y = "Mean(selected − oracle) basis") +
  shared_theme

# ── Combine ────────────────────────────────────────────────────────────────
fig <- p_a + p_b + plot_layout(widths = c(1.15, 1))
ggsave("analysis/figure/test_s6v3_fig1.png", fig, width = 10, height = 4.5, dpi = 150)
cat("Saved test_s6v3_fig1.png\n")

# ── Print ──────────────────────────────────────────────────────────────────
cat("\nMAD:\n")
print(mad_data %>% pivot_wider(names_from = n_reps_used, values_from = c(mad, se), names_prefix = "R="))
cat("\nDirectional bias (R=5):\n")
print(bias_data)
