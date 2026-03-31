## Section 6 — all 8 designs, full nbasis grid (no filtering)
library(tidyverse)
library(patchwork)
library(latex2exp)

CONFIGS <- c("equal_50", "equal_100")

dir.create("analysis/figure", showWarnings = FALSE, recursive = TRUE)

# ── Load (no nbasis filter) ─────────────────────────────────────────────────
dt_all <- readRDS("results_multi_config/dt_epsilon_8configs.RDS") %>%
  filter(config %in% CONFIGS, loss_function == "MSE")

oracle_all <- readRDS("results_multi_config/oracle_epsilon_8configs.RDS") %>%
  filter(config %in% CONFIGS)

# Ordered parseable facet label
label_map <- setNames(
  paste0("n[i]==", gsub("equal_", "", CONFIGS)),
  CONFIGS
)
add_design <- function(df) mutate(df, design = factor(label_map[config], levels = label_map))

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
  mutate(deviation = dt_nbasis - oracle_nbasis) %>%
  add_design()

oracle_sel <- add_design(oracle_sel)

# Quick summary
cat("nbasis values in data:", sort(unique(oracle_all$nbasis)), "\n")
cat("Oracle selections by design:\n")
print(oracle_sel %>% count(design, oracle_nbasis) %>% pivot_wider(names_from = oracle_nbasis, values_from = n, values_fill = 0))

# ── (a) MAD basis by design ─────────────────────────────────────────────────
mad_data <- dt_sel %>%
  group_by(design, epsilon, n_reps_used) %>%
  summarise(mad = mean(abs(deviation)),
            se = sd(abs(deviation)) / sqrt(n()), .groups = "drop")

p_a <- ggplot(mad_data, aes(epsilon, mad,
                             color = factor(n_reps_used),
                             fill = factor(n_reps_used))) +
  geom_ribbon(aes(ymin = mad - se, ymax = mad + se), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ design, labeller = label_parsed, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = c(0.3, 0.6, 0.9)) +
  scale_color_brewer(palette = "Set1", name = "R") +
  scale_fill_brewer(palette = "Set1", name = "R") +
  labs(title = "(a) Mean absolute deviation of selected number of basis functions",
       x = expression(epsilon),
       y = TeX("Mean $|p - p^*|$")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 10))

# ── (b) Directional bias by design ──────────────────────────────────────────
bias_data <- dt_sel %>%
  group_by(design, epsilon, n_reps_used) %>%
  summarise(mean_bias = mean(deviation),
            se = sd(deviation) / sqrt(n()), .groups = "drop")

# Per-design random selection baseline
all_nbasis <- sort(unique(oracle_all$nbasis))
random_bias <- oracle_sel %>%
  group_by(design) %>%
  summarise(random_bias = mean(all_nbasis) - mean(oracle_nbasis), .groups = "drop")

random_label <- random_bias %>%
  mutate(label = sprintf("random (+%.1f)", random_bias))

p_b <- ggplot(bias_data, aes(epsilon, mean_bias,
                              color = factor(n_reps_used),
                              fill = factor(n_reps_used))) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "gray50") +
  geom_ribbon(aes(ymin = mean_bias - se, ymax = mean_bias + se),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  facet_wrap(~ design, labeller = label_parsed, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = c(0.3, 0.6, 0.9)) +
  scale_color_brewer(palette = "Set1", name = "R") +
  scale_fill_brewer(palette = "Set1", name = "R") +
  labs(title = "(b) Directional bias of selected number of basis functions",
       x = expression(epsilon),
       y = TeX("Mean $p - p^*$")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 10))

# ── Combine ──────────────────────────────────────────────────────────────────
fig <- p_a / p_b + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("analysis/figure/test_s6_full_grid.png", fig, width = 12, height = 10, dpi = 150)
cat("Saved test_s6_full_grid.png\n")
