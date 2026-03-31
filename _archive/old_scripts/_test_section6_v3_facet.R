## Section 6 — faceted by design (3 selected designs)
library(tidyverse)
library(patchwork)

PAPER_NBASIS <- c(6, 18, 30, 42)
CONFIGS <- c("equal_40", "equal_75", "equal_125")

dir.create("analysis/figure", showWarnings = FALSE, recursive = TRUE)

# ── Load & filter ───────────────────────────────────────────────────────────
dt_all <- readRDS("results_multi_config/dt_epsilon_8configs.RDS") %>%
  filter(config %in% CONFIGS, loss_function == "MSE", nbasis %in% PAPER_NBASIS)

oracle_all <- readRDS("results_multi_config/oracle_epsilon_8configs.RDS") %>%
  filter(config %in% CONFIGS, nbasis %in% PAPER_NBASIS)

# Ordered parseable facet label
label_map <- c(equal_40 = "n[i]==40", equal_75 = "n[i]==75", equal_125 = "n[i]==125")
add_design <- function(df) mutate(df, design = factor(label_map[config], levels = label_map))

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
  mutate(deviation = dt_nbasis - oracle_nbasis) %>%
  add_design()

oracle_sel <- add_design(oracle_sel)

# ── (a) MAD basis by design ────────────────────────────────────────────────
mad_data <- dt_sel %>%
  group_by(design, epsilon, n_reps_used) %>%
  summarise(mad = mean(abs(deviation)),
            se = sd(abs(deviation)) / sqrt(n()), .groups = "drop")

p_a <- ggplot(mad_data, aes(epsilon, mad,
                             color = factor(n_reps_used),
                             fill = factor(n_reps_used))) +
  geom_ribbon(aes(ymin = mad - se, ymax = mad + se), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ design, labeller = label_parsed, nrow = 1) +
  scale_x_continuous(breaks = c(0.3, 0.6, 0.9)) +
  scale_color_brewer(palette = "Set1", name = "R") +
  scale_fill_brewer(palette = "Set1", name = "R") +
  labs(title = "(a) MAD basis",
       x = expression(epsilon),
       y = "Mean |selected − oracle| basis") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 10))

# ── (b) Directional bias by design ────────────────────────────────────────
bias_data <- dt_sel %>%
  group_by(design, epsilon, n_reps_used) %>%
  summarise(mean_bias = mean(deviation),
            se = sd(deviation) / sqrt(n()), .groups = "drop")

# Per-design random selection baseline
random_bias <- oracle_sel %>%
  group_by(design) %>%
  summarise(random_bias = mean(PAPER_NBASIS) - mean(oracle_nbasis), .groups = "drop")

# Label for random baseline, positioned inside each facet
random_label <- random_bias %>%
  mutate(label = sprintf("random (%.1f)", random_bias))

p_b <- ggplot(bias_data, aes(epsilon, mean_bias,
                              color = factor(n_reps_used),
                              fill = factor(n_reps_used))) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "gray50") +
  geom_hline(data = random_bias, aes(yintercept = random_bias),
             linetype = "dashed", linewidth = 0.5, color = "gray35") +
  geom_text(data = random_label, aes(x = 0.6, y = random_bias - 0.8, label = label),
            inherit.aes = FALSE, size = 2.8, color = "gray35") +
  geom_ribbon(aes(ymin = mean_bias - se, ymax = mean_bias + se),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ design, labeller = label_parsed, nrow = 1) +
  scale_x_continuous(breaks = c(0.3, 0.6, 0.9)) +
  scale_color_brewer(palette = "Set1", name = "R") +
  scale_fill_brewer(palette = "Set1", name = "R") +
  labs(title = "(b) Directional bias",
       x = expression(epsilon),
       y = "Mean(selected − oracle) basis") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 10))

# ── Stack vertically ──────────────────────────────────────────────────────
fig <- p_a / p_b + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("analysis/figure/test_s6v3_facet.png", fig, width = 10, height = 7, dpi = 150)
cat("Saved test_s6v3_facet.png\n")
