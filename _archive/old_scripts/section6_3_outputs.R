library(tidyverse)
library(xtable)

# Load data
oracle_all <- readRDS("results_multi_config/methods_comparison/oracle_results.RDS")
selections <- readRDS("results_multi_config/methods_comparison/method_selections.RDS")

# Get oracle selections
oracle_selections <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, oracle_model = model, oracle_nbasis = nbasis) %>%
  ungroup()

# Compute deviations
deviations_all <- selections %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(
    deviation = selected_nbasis - oracle_nbasis,
    abs_deviation = abs(deviation),
    method_type = case_when(
      str_detect(method, "DT-MSE") ~ "DT-MSE",
      str_detect(method, "DT-plugin_NLL") ~ "DT-NLL",
      TRUE ~ method
    ),
    method_label = case_when(
      str_detect(method, "DT-MSE") ~ sprintf("DT-MSE ε=%.1f", epsilon),
      str_detect(method, "DT-plugin_NLL") ~ sprintf("DT-NLL ε=%.1f", epsilon),
      TRUE ~ method
    )
  )

cat("Loaded data:\n")
cat(sprintf("  - Oracle selections: %d samples × %d designs = %d\n",
            n_distinct(oracle_selections$comp_no),
            n_distinct(oracle_selections$config),
            nrow(oracle_selections)))
cat(sprintf("  - Deviations: %d rows\n", nrow(deviations_all)))

# ============================================================================
# OUTPUT 1: MAD Summary Table
# ============================================================================

cat("\n\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("OUTPUT 1: MAD Summary Table\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

# Filter to methods for table (exclude ε=0.5)
table_methods <- c(
  "DT-MSE ε=0.7",
  "DT-NLL ε=0.7",
  "DIC",
  "WAIC",
  "ESIM"
)

# Compute MAD and SD by method and design
mad_summary <- deviations_all %>%
  filter(method_label %in% table_methods) %>%
  group_by(method_label, config) %>%
  summarise(
    MAD = mean(abs_deviation),
    SD = sd(abs_deviation),
    n = n(),
    .groups = "drop"
  )

# Compute overall average (across all 60 samples)
mad_overall <- deviations_all %>%
  filter(method_label %in% table_methods) %>%
  group_by(method_label) %>%
  summarise(
    MAD = mean(abs_deviation),
    SD = sd(abs_deviation),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(config = "Overall")

# Combine
mad_combined <- bind_rows(mad_summary, mad_overall)

# Pivot to wide format
mad_wide <- mad_combined %>%
  mutate(
    cell_value = sprintf("%.2f (%.2f)", MAD, SD)
  ) %>%
  select(method_label, config, cell_value) %>%
  pivot_wider(names_from = config, values_from = cell_value)

# Get MAD values for sorting
mad_for_sorting <- mad_combined %>%
  filter(config == "Overall") %>%
  select(method_label, MAD)

# Reorder columns and sort by Overall MAD
mad_wide <- mad_wide %>%
  left_join(mad_for_sorting, by = "method_label") %>%
  arrange(MAD) %>%
  select(method_label, prop_0.75pct, prop_1p25pct, prop_1p75pct, Overall)

cat("MAD Summary (before formatting):\n")
print(mad_wide)

# Find minimum MAD in each column for bolding
min_vals <- mad_combined %>%
  group_by(config) %>%
  slice_min(MAD, n = 1, with_ties = FALSE) %>%
  select(config, method_label)

cat("\n\nMinimum MAD by column:\n")
print(min_vals)

# Create formatted table with bold minimums
mad_wide_bold <- mad_wide
for (i in 1:nrow(min_vals)) {
  cfg <- min_vals$config[i]
  method <- min_vals$method_label[i]

  # Find row and column
  row_idx <- which(mad_wide_bold$method_label == method)
  col_name <- cfg

  if (col_name == "prop_0.75pct") col_name <- "prop_0.75pct"
  if (col_name == "prop_1p25pct") col_name <- "prop_1p25pct"
  if (col_name == "prop_1p75pct") col_name <- "prop_1p75pct"

  # Bold only MAD (not SD in parentheses)
  current_val <- mad_wide_bold[[col_name]][row_idx]
  # Extract MAD and SD parts
  parts <- strsplit(current_val, " ")[[1]]
  mad_part <- parts[1]  # e.g., "7.20"
  sd_part <- parts[2]   # e.g., "(3.56)"
  # Bold only MAD
  mad_wide_bold[[col_name]][row_idx] <- paste0("\\textbf{", mad_part, "} ", sd_part)
}

cat("\n\nFormatted table with bold minimums:\n")
print(mad_wide_bold)

# Get oracle basis means
oracle_means <- oracle_selections %>%
  group_by(config) %>%
  summarise(mean_oracle = mean(oracle_nbasis))

cat("\n\nOracle basis means by design:\n")
print(oracle_means)

# Create LaTeX table
latex_table <- mad_wide_bold %>%
  rename(
    Method = method_label,
    `0.75\\% PA` = prop_0.75pct,
    `1.25\\% PA` = prop_1p25pct,
    `1.75\\% PA` = prop_1p75pct
  )

# Generate LaTeX
latex_code <- capture.output({
  print(xtable(latex_table,
               caption = "Mean and standard deviation of the absolute deviation (in parenthesis) from oracle basis by method and design across 20 simulation replicates. Bold indicates best performance in each column.",
               label = "tab:mad_comparison",
               align = c("l", "l", "r", "r", "r", "r")),
        include.rownames = FALSE,
        sanitize.text.function = identity,
        booktabs = TRUE)
})

# Write to file
writeLines(latex_code, "analysis/table_mad_comparison.tex")
cat("\n\nLaTeX table written to: analysis/table_mad_comparison.tex\n")

# ============================================================================
# OUTPUT 2: Boxplot of Deviations
# ============================================================================

cat("\n\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("OUTPUT 2: Boxplot of Deviations\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

# Filter to keep only ε=0.7 for DT methods
deviations_plot <- deviations_all %>%
  filter(is.na(epsilon) | abs(epsilon - 0.7) < 0.001) %>%
  mutate(
    config_label = case_when(
      config == "prop_0.75pct" ~ "0.75% PA",
      config == "prop_1p25pct" ~ "1.25% PA",
      config == "prop_1p75pct" ~ "1.75% PA"
    ),
    # Set factor order for x-axis
    method_label = factor(method_label, levels = c(
      "DT-NLL ε=0.7",
      "DT-MSE ε=0.7",
      "DIC",
      "WAIC",
      "ESIM"
    )),
    # Categorize as DT vs non-DT
    method_category = ifelse(method_type %in% c("DT-MSE", "DT-NLL"), "DT", "Non-DT")
  )

cat("Boxplot data summary:\n")
deviations_plot %>%
  group_by(method_label) %>%
  summarise(
    n = n(),
    mean_dev = mean(deviation),
    median_dev = median(deviation)
  ) %>%
  print()

# Create boxplot
p <- ggplot(deviations_plot, aes(x = method_label, y = deviation, fill = method_category)) +
  geom_boxplot(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~config_label, ncol = 3) +
  labs(
    x = "Method",
    y = "Deviation from oracle basis"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("DT" = "lightblue", "Non-DT" = "lightgray"))

# Save plot
ggsave("analysis/plot5_comparison.png", p, width = 10, height = 4, dpi = 300)
cat("\n\nBoxplot saved to: analysis/plot5_comparison.png\n")

cat("\n\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("DONE! All outputs generated successfully.\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
