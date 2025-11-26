library(dplyr)
library(ggplot2)
source("sim_functions/summary_oracle.R")

# Load oracle results
oracle_results <- lapply(1:10, function(comp_no) {
  summary_oracle(comp_no, "_results_wagp")
}) %>% bind_rows()

oracle <- oracle_results %>%
  mutate(nbasis = as.numeric(gsub("nbasis_", "", model))) %>%
  filter(!is.na(nbasis))

# Calculate average
oracle_avg <- oracle %>%
  group_by(nbasis) %>%
  summarise(mean_mse = mean(mse_true), .groups = "drop")

# Normalize each comparison to [0, 1] for easier comparison
oracle_norm <- oracle %>%
  group_by(comp_no) %>%
  mutate(
    mse_norm = (mse_true - min(mse_true)) / (max(mse_true) - min(mse_true))
  ) %>%
  ungroup()

# Create faceted plot
p1 <- ggplot(oracle_norm, aes(x = nbasis, y = mse_norm)) +
  geom_line(color = "blue", linewidth = 0.8) +
  geom_point(color = "blue", size = 2) +
  facet_wrap(~comp_no, ncol = 5, labeller = labeller(comp_no = function(x) paste("Comp", x))) +
  labs(
    title = "OD-Oracle MSE by Comparison (WAGP - Random Spatial Effects)",
    subtitle = "Normalized to [0,1] per comparison to show curve shape",
    x = "Number of basis functions",
    y = "Normalized MSE (0=best, 1=worst)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.text = element_text(size = 9),
    axis.text = element_text(size = 7)
  )

ggsave("/tmp/wagp_oracle_faceted.png", p1, width = 12, height = 6, dpi = 150)

# Also create overlay plot
p2 <- ggplot() +
  geom_line(data = oracle_norm, aes(x = nbasis, y = mse_norm, group = comp_no),
            alpha = 0.3, color = "blue") +
  stat_summary(data = oracle_norm, aes(x = nbasis, y = mse_norm),
               fun = mean, geom = "line", color = "red", linewidth = 1.5) +
  labs(
    title = "OD-Oracle MSE: Individual vs Average (WAGP)",
    subtitle = "Blue = individual comparisons (normalized), Red = average pattern",
    x = "Number of basis functions",
    y = "Normalized MSE"
  ) +
  theme_minimal()

ggsave("/tmp/wagp_oracle_overlay.png", p2, width = 10, height = 6, dpi = 150)

cat("✓ Plots saved to:\n")
cat("  /tmp/wagp_oracle_faceted.png\n")
cat("  /tmp/wagp_oracle_overlay.png\n")
