#!/usr/bin/env Rscript
# Plot faceted oracle MSE curves
# Usage: Rscript plot_oracle_faceted.R <results_dir>
# Example: Rscript plot_oracle_faceted.R _results_ca_full_comparison

library(dplyr)
library(ggplot2)
source("sim_functions/summary_oracle.R")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("Usage: Rscript plot_oracle_faceted.R <results_dir>\n")
  cat("Example: Rscript plot_oracle_faceted.R _results_ca_full_comparison\n")
  quit(status = 1)
}

results_dir <- args[1]

# Load config to get number of comparisons
config_file <- file.path(results_dir, "model_config.RDS")
if (!file.exists(config_file)) {
  stop(sprintf("Config file not found: %s", config_file))
}
model_config <- readRDS(config_file)
n_comp <- model_config$n_comp

cat(sprintf("Loading oracle results from %s (%d comparisons)...\n", results_dir, n_comp))

# Load oracle results
oracle_results <- lapply(1:n_comp, function(comp_no) {
  summary_oracle(comp_no, results_dir)
}) %>% bind_rows()

oracle <- oracle_results %>%
  mutate(nbasis = as.numeric(gsub("nbasis_", "", model))) %>%
  filter(!is.na(nbasis))

# Calculate average
oracle_avg <- oracle %>%
  group_by(nbasis) %>%
  summarise(mean_mse = mean(mse_true), .groups = "drop")

# Determine grid layout
ncols <- min(5, ceiling(sqrt(n_comp)))

# Create faceted plot (raw MSE values)
p1 <- ggplot(oracle, aes(x = nbasis, y = mse_true)) +
  geom_line(color = "blue", linewidth = 0.8) +
  geom_point(color = "blue", size = 2) +
  facet_wrap(~comp_no, ncol = ncols, scales = "free_y",
             labeller = labeller(comp_no = function(x) paste("Comp", x))) +
  labs(
    title = sprintf("OD-Oracle MSE by Comparison (%d comparisons)", n_comp),
    subtitle = "Raw MSE values (free y-axis per comparison)",
    x = "Number of basis functions",
    y = "MSE"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.text = element_text(size = 9),
    axis.text = element_text(size = 7)
  )

output_file <- file.path(results_dir, "oracle_faceted_raw.png")
ggsave(output_file, p1, width = 12, height = ceiling(n_comp/ncols) * 2, dpi = 150)

# Also create overlay plot (raw values)
p2 <- ggplot() +
  geom_line(data = oracle, aes(x = nbasis, y = mse_true, group = comp_no),
            alpha = 0.3, color = "blue") +
  stat_summary(data = oracle, aes(x = nbasis, y = mse_true),
               fun = mean, geom = "line", color = "red", linewidth = 1.5) +
  labs(
    title = sprintf("OD-Oracle MSE: Individual vs Average (%d comparisons)", n_comp),
    subtitle = "Blue = individual comparisons, Red = average pattern",
    x = "Number of basis functions",
    y = "MSE"
  ) +
  theme_minimal()

output_file2 <- file.path(results_dir, "oracle_overlay_raw.png")
ggsave(output_file2, p2, width = 10, height = 6, dpi = 150)

cat("\n✓ Plots saved to:\n")
cat(sprintf("  %s\n", output_file))
cat(sprintf("  %s\n", output_file2))
