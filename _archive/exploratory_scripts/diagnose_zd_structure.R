#!/usr/bin/env Rscript
# ==============================================================================
# Diagnose z/d structure across configs
# ==============================================================================
# Understand what makes prop_2pct different in terms of direct estimates

library(tidyverse)

cat("\n=== z/d STRUCTURE DIAGNOSTIC ===\n\n")

configs <- c("equal_40", "equal_50", "equal_75",
             "prop_0.5pct", "prop_1pct", "prop_2pct")

# Map config names to directory names
config_to_dir <- c(
  "equal_40" = "_results_equal40_comparison",
  "equal_50" = "_results_equal50_comparison",
  "equal_75" = "_results_equal75_comparison",
  "prop_0.5pct" = "_results_prop0.5pct_comparison",
  "prop_1pct" = "_results_prop1pct_comparison",
  "prop_2pct" = "_results_prop2pct_comparison"
)

# Analyze z/d for all comparisons in each config
zd_stats <- map_df(configs, function(cfg) {
  cat(sprintf("Processing %s...\n", cfg))

  # Get all comparison folders
  results_dir <- config_to_dir[cfg]
  if (!dir.exists(results_dir)) {
    cat(sprintf("  Directory %s not found, skipping...\n", results_dir))
    return(NULL)
  }

  comp_dirs <- list.files(results_dir,
                         pattern = "^comparison_",
                         full.names = TRUE)

  if (length(comp_dirs) == 0) {
    return(NULL)
  }

  # Collect z/d from all comparisons
  comp_data <- map_df(comp_dirs, function(comp_dir) {
    comp_no <- as.integer(str_extract(basename(comp_dir), "[0-9]+"))

    z_list <- readRDS(file.path(comp_dir, "z.RDS"))
    d_list <- readRDS(file.path(comp_dir, "d.RDS"))

    # Extract values (d and z are lists with $puma and $values)
    z <- z_list$values
    d <- d_list$values

    tibble(
      config = cfg,
      comp_no = comp_no,
      puma_id = 1:length(z),
      z = z,
      d = d,
      se = sqrt(d),  # Standard error
      cv = sqrt(d) / abs(z)  # Coefficient of variation
    )
  })

  # Summary stats per config
  comp_data %>%
    group_by(config) %>%
    summarise(
      n_comparisons = length(unique(comp_no)),
      n_pumas = n() / n_comparisons,
      # z statistics
      z_mean = mean(z),
      z_sd = sd(z),
      z_min = min(z),
      z_max = max(z),
      z_range = z_max - z_min,
      # d statistics
      d_mean = mean(d),
      d_sd = sd(d),
      d_min = min(d),
      d_max = max(d),
      d_ratio = d_max / d_min,
      # Effective sample size
      eff_n_mean = mean(1/d),
      eff_n_sd = sd(1/d),
      eff_n_min = min(1/d),
      eff_n_max = max(1/d),
      # Signal-to-noise
      mean_cv = mean(cv, na.rm = TRUE),  # Average coefficient of variation
      median_cv = median(cv, na.rm = TRUE),
      # Cross-comparison variability (z variability for same PUMA across comps)
      .groups = "drop"
    )
})

cat("\n=== SUMMARY STATISTICS ===\n\n")
print(zd_stats %>%
        select(config, n_comparisons, n_pumas) %>%
        knitr::kable())

cat("\n\n=== z (Direct Estimate) Characteristics ===\n\n")
print(zd_stats %>%
        select(config, z_mean, z_sd, z_range) %>%
        knitr::kable(digits = 4))

cat("\n\n=== d (Design Variance) Characteristics ===\n\n")
print(zd_stats %>%
        select(config, d_mean, d_sd, d_min, d_max, d_ratio) %>%
        knitr::kable(digits = 6))

cat("\n\n=== Effective Sample Size (1/d) ===\n\n")
print(zd_stats %>%
        select(config, eff_n_mean, eff_n_sd, eff_n_min, eff_n_max) %>%
        knitr::kable(digits = 1))

cat("\n\n=== Signal-to-Noise (Coefficient of Variation) ===\n\n")
print(zd_stats %>%
        select(config, mean_cv, median_cv) %>%
        knitr::kable(digits = 3))

# Deep dive: Cross-comparison variability
cat("\n\n=== CROSS-COMPARISON VARIABILITY ===\n")
cat("(How much does z vary for the same PUMA across comparisons?)\n\n")

cross_comp_var <- map_df(configs, function(cfg) {
  results_dir <- config_to_dir[cfg]
  if (!dir.exists(results_dir)) {
    return(NULL)
  }

  comp_dirs <- list.files(results_dir,
                         pattern = "^comparison_",
                         full.names = TRUE)

  comp_data <- map_df(comp_dirs, function(comp_dir) {
    comp_no <- as.integer(str_extract(basename(comp_dir), "[0-9]+"))
    z_list <- readRDS(file.path(comp_dir, "z.RDS"))
    z <- z_list$values
    tibble(
      config = cfg,
      comp_no = comp_no,
      puma_id = 1:length(z),
      z = z
    )
  })

  # For each PUMA, compute SD of z across comparisons
  puma_var <- comp_data %>%
    group_by(config, puma_id) %>%
    summarise(
      z_mean_across_comps = mean(z),
      z_sd_across_comps = sd(z),
      .groups = "drop"
    )

  # Summary across all PUMAs
  puma_var %>%
    group_by(config) %>%
    summarise(
      mean_z_sd = mean(z_sd_across_comps),
      median_z_sd = median(z_sd_across_comps),
      .groups = "drop"
    )
})

print(cross_comp_var %>% knitr::kable(digits = 5))

# Combine with oracle stats
oracle_all <- readRDS("results_multi_config/oracle_all.RDS")
oracle_summary <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  group_by(config) %>%
  summarise(
    oracle_mean_nbasis = mean(nbasis),
    oracle_sd_nbasis = sd(nbasis),
    .groups = "drop"
  )

# Actually compute MSE variation
mse_var <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  summarise(
    mse_range = max(mse_true) - min(mse_true),
    mse_min = min(mse_true),
    .groups = "drop"
  ) %>%
  group_by(config) %>%
  summarise(
    mse_variation_pct = mean(mse_range / mse_min) * 100,
    .groups = "drop"
  )

# Combine everything
final_summary <- zd_stats %>%
  left_join(cross_comp_var, by = "config") %>%
  left_join(mse_var, by = "config") %>%
  left_join(oracle_summary %>% select(config, oracle_mean_nbasis, oracle_sd_nbasis),
            by = "config")

cat("\n\n=== COMBINED SUMMARY (sorted by MSE variation) ===\n\n")
print(final_summary %>%
        select(config, d_mean, eff_n_mean, mean_cv, mean_z_sd,
               mse_variation_pct, oracle_mean_nbasis, oracle_sd_nbasis) %>%
        arrange(mse_variation_pct) %>%
        knitr::kable(digits = 3))

cat("\n\n=== KEY CORRELATIONS ===\n\n")
cat("Does low d (high eff_n) → low MSE variation → harder selection?\n\n")

correlations <- final_summary %>%
  summarise(
    cor_d_msevar = cor(d_mean, mse_variation_pct),
    cor_effn_msevar = cor(eff_n_mean, mse_variation_pct),
    cor_cv_msevar = cor(mean_cv, mse_variation_pct),
    cor_msevar_oracle_sd = cor(mse_variation_pct, oracle_sd_nbasis)
  )

print(correlations %>% knitr::kable(digits = 3))

cat("\n\nInterpretation:\n")
cat("- cor_d_msevar: Higher d (less data) → higher MSE variation (easier to differentiate models)\n")
cat("- cor_effn_msevar: Higher eff_n (more data) → lower MSE variation (harder to differentiate)\n")
cat("- cor_cv_msevar: Higher CV (more noise) → ? MSE variation\n")
cat("- cor_msevar_oracle_sd: Higher MSE variation → ? oracle stability\n")

cat("\n")
