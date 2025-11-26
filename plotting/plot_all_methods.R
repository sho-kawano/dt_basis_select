#!/usr/bin/env Rscript
# Plot all methods comparison with faceted view
# Usage: Rscript plot_all_methods.R <results_dir>
# Example: Rscript plot_all_methods.R _results_ca_full_comparison

library(dplyr)
library(ggplot2)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("Usage: Rscript plot_all_methods.R <results_dir>\n")
  cat("Example: Rscript plot_all_methods.R _results_ca_full_comparison\n")
  quit(status = 1)
}

results_dir <- args[1]
summaries_file <- file.path(results_dir, "results.RDS")

if (!file.exists(summaries_file)) {
  stop(sprintf("Summaries file not found: %s\nRun 'Rscript run_summary.R %s' first",
               summaries_file, results_dir))
}

cat(sprintf("Loading results from %s...\n", summaries_file))

# Load results
results <- readRDS(summaries_file)
results <- results %>%
  mutate(nbasis = as.numeric(gsub("nbasis_", "", model))) %>%
  filter(!is.na(nbasis))

# Prepare data for each method
# 1. Oracle (OD-Oracle MSE)
oracle <- results %>%
  filter(method_name == "oracle") %>%
  select(comp_no, nbasis, value = metric) %>%
  mutate(method = "OD-Oracle MSE")

# 2. WAIC
waic <- results %>%
  filter(method_name == "waic") %>%
  select(comp_no, nbasis, value = metric) %>%
  mutate(method = "WAIC")

# 3. DIC
dic <- results %>%
  filter(method_name == "dic") %>%
  select(comp_no, nbasis, value = metric) %>%
  mutate(method = "DIC")

# 4. DT 5-fold (use MSE loss)
dt5 <- results %>%
  filter(method_name == "dt_5fold", loss_function == "MSE") %>%
  select(comp_no, nbasis, value = metric) %>%
  mutate(method = "DT 5-fold")

# 5. DT 1-fold (eps=0.3, 3 reps, MSE loss)
dt1_03 <- results %>%
  filter(method_name == "dt_1fold", epsilon == 0.3, n_reps_used == 3, loss_function == "MSE") %>%
  select(comp_no, nbasis, value = metric) %>%
  mutate(method = "DT 1-fold (eps=0.3)")

# 6. DT 1-fold (eps=0.5, 3 reps, MSE loss)
dt1_05 <- results %>%
  filter(method_name == "dt_1fold", epsilon == 0.5, n_reps_used == 3, loss_function == "MSE") %>%
  select(comp_no, nbasis, value = metric) %>%
  mutate(method = "DT 1-fold (eps=0.5)")

# 7. ESIM Standard
esim_std <- results %>%
  filter(method_name == "esim_standard") %>%
  select(comp_no, nbasis, value = metric) %>%
  mutate(method = "ESIM Standard")

# 8. ESIM Data Fission
esim_fis <- results %>%
  filter(method_name == "esim_fission") %>%
  select(comp_no, nbasis, value = metric) %>%
  mutate(method = "ESIM Data Fission")

# Combine all methods
all_methods <- bind_rows(oracle, waic, dic, dt5, dt1_03, dt1_05, esim_std, esim_fis)

# Get number of comparisons
n_comp <- length(unique(all_methods$comp_no))
ncols <- min(5, ceiling(sqrt(n_comp)))

# Create faceted plot (raw values, free y-axis)
p1 <- ggplot(all_methods, aes(x = nbasis, y = value, color = method, group = method)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~comp_no, ncol = ncols, scales = "free_y",
             labeller = labeller(comp_no = function(x) paste("Comp", x))) +
  labs(
    title = sprintf("All Methods Comparison (%d comparisons)", n_comp),
    subtitle = "Raw values (free y-axis per comparison, lower = better)",
    x = "Number of basis functions",
    y = "Score",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.text = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2))

output_file1 <- file.path(results_dir, "all_methods_faceted_raw.png")
ggsave(output_file1, p1, width = 14, height = ceiling(n_comp/ncols) * 2 + 2, dpi = 150)

# Create overlay plot showing average patterns (raw values)
all_methods_avg <- all_methods %>%
  group_by(nbasis, method) %>%
  summarise(mean_value = mean(value), .groups = "drop")

p2 <- ggplot() +
  geom_line(data = all_methods,
            aes(x = nbasis, y = value, group = interaction(comp_no, method), color = method),
            alpha = 0.2, linewidth = 0.5) +
  geom_line(data = all_methods_avg,
            aes(x = nbasis, y = mean_value, color = method),
            linewidth = 1.5) +
  labs(
    title = sprintf("All Methods: Individual vs Average (%d comparisons)", n_comp),
    subtitle = "Thin lines = individual comparisons, Thick lines = average pattern",
    x = "Number of basis functions",
    y = "Score",
    color = "Method"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2))

output_file2 <- file.path(results_dir, "all_methods_overlay_raw.png")
ggsave(output_file2, p2, width = 12, height = 8, dpi = 150)

cat("\n✓ Plots saved to:\n")
cat(sprintf("  %s\n", output_file1))
cat(sprintf("  %s\n", output_file2))

# Print summary statistics
cat("\n=== OPTIMAL NBASIS BY METHOD ===\n")
optimal_summary <- all_methods %>%
  group_by(comp_no, method) %>%
  summarise(optimal_nbasis = nbasis[which.min(value)], .groups = "drop") %>%
  group_by(method) %>%
  summarise(
    mean_optimal = mean(optimal_nbasis),
    sd_optimal = sd(optimal_nbasis),
    min_optimal = min(optimal_nbasis),
    max_optimal = max(optimal_nbasis),
    .groups = "drop"
  ) %>%
  arrange(mean_optimal)

print(optimal_summary, n = Inf)
