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

# Create faceted plot by METHOD (each method in its own panel)
# Individual comparisons = light gray lines, Average = bold blue line
p1 <- ggplot(all_methods, aes(x = nbasis, y = value)) +
  geom_line(aes(group = comp_no), alpha = 0.25, color = "gray60", linewidth = 0.5) +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               color = "blue", linewidth = 1.2) +
  facet_wrap(~method, ncol = 3, scales = "free_y") +
  labs(
    title = sprintf("All Methods: Individual vs Average (%d comparisons)", n_comp),
    subtitle = "Gray = individual comparisons, Blue = average. Raw values (free y-axis per method)",
    x = "Number of basis functions",
    y = "Score"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.text = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 8)
  )

output_file1 <- file.path(results_dir, "all_methods_faceted_raw.png")
ggsave(output_file1, p1, width = 12, height = 8, dpi = 150)

cat("\n✓ Plot saved to:\n")
cat(sprintf("  %s\n", output_file1))

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
