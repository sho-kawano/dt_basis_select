#!/usr/bin/env Rscript
# ==============================================================================
# Multi-Config Summary Script
# ==============================================================================
# Purpose: Analyze and compare 6 sampling configurations
# ==============================================================================

library(tidyverse)
library(parallel)
library(doParallel)

cat("\n")
cat("================================================================================\n")
cat("  MULTI-CONFIG SUMMARY\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Configs analyzed (5 new + 1 baseline)
configs <- list(
  list(name = "equal_40", results_dir = "_results_equal40_comparison", n_comp = 20),
  list(name = "equal_50", results_dir = "_results_equal50_comparison", n_comp = 20),
  list(name = "equal_75", results_dir = "_results_equal75_rerun", n_comp = 20),
  list(name = "prop_0.5pct", results_dir = "_results_prop0.5pct_comparison", n_comp = 20),
  list(name = "prop_1pct", results_dir = "_results_prop1pct_comparison", n_comp = 20),
  list(name = "prop_2pct", results_dir = "_results_prop2pct_comparison", n_comp = 20)
)

cat("Configs to analyze:\n")
for (cfg in configs) {
  cat(sprintf("  - %s (%d comparisons) -> %s\n", cfg$name, cfg$n_comp, cfg$results_dir))
}
cat("\n")

n_cores <- 8

# ==============================================================================
# LOAD AND PROCESS RESULTS
# ==============================================================================

cat("=== LOADING RESULTS ===\n\n")

source("sim_functions/summary_oracle.R")
source("sim_functions/summary_dt.R")

all_results <- list()

for (cfg in configs) {
  cat(sprintf("Processing %s...\n", cfg$name))

  # Oracle summaries
  oracle_list <- lapply(1:cfg$n_comp, function(comp_no) {
    summary_oracle(comp_no, cfg$results_dir) %>%
      mutate(config = cfg$name)
  })
  oracle_df <- bind_rows(oracle_list)

  # DT 3-fold summaries
  dt3_list <- lapply(1:cfg$n_comp, function(comp_no) {
    summary_dt(comp_no, n_folds = 3, loss_function = "MSE",
               results_dir = cfg$results_dir, n_reps_to_use = 1) %>%
      mutate(config = cfg$name, method_name = "dt_3fold")
  })
  dt3_df <- bind_rows(dt3_list)

  # DT 5-fold summaries
  dt5_list <- lapply(1:cfg$n_comp, function(comp_no) {
    summary_dt(comp_no, n_folds = 5, loss_function = "MSE",
               results_dir = cfg$results_dir, n_reps_to_use = 1) %>%
      mutate(config = cfg$name, method_name = "dt_5fold")
  })
  dt5_df <- bind_rows(dt5_list)

  # DT 1-fold summaries (epsilon 0.3, 0.5, 0.7 with 5 reps)
  dt1_list <- list()
  for (eps in c(0.3, 0.5, 0.7)) {
    for (n_reps in c(1, 3, 5)) {
      dt1_eps <- lapply(1:cfg$n_comp, function(comp_no) {
        summary_dt(comp_no, n_folds = 1, loss_function = "MSE",
                   results_dir = cfg$results_dir,
                   n_reps_to_use = n_reps, eps = eps) %>%
          mutate(config = cfg$name, method_name = "dt_1fold",
                 epsilon = eps, n_reps_used = n_reps)
      })
      dt1_list <- c(dt1_list, dt1_eps)
    }
  }
  dt1_df <- bind_rows(dt1_list)

  all_results[[cfg$name]] <- list(
    oracle = oracle_df,
    dt3 = dt3_df,
    dt5 = dt5_df,
    dt1 = dt1_df
  )
}

cat("✓ All results loaded\n\n")

# Combine all oracle results
oracle_all <- bind_rows(lapply(all_results, function(x) x$oracle))

# Combine all DT results
dt_all <- bind_rows(
  bind_rows(lapply(all_results, function(x) x$dt3)),
  bind_rows(lapply(all_results, function(x) x$dt5)),
  bind_rows(lapply(all_results, function(x) x$dt1))
)

# ==============================================================================
# ORACLE ANALYSIS
# ==============================================================================

cat("=== ORACLE ANALYSIS ===\n\n")

# Compute oracle properties for each config
compute_oracle_properties <- function(oracle_df, config_name) {
  oracle_models <- oracle_df %>%
    filter(model != "Direct") %>%
    mutate(nbasis = as.numeric(gsub("nbasis_", "", model)))

  oracle_optimal <- oracle_models %>%
    group_by(comp_no) %>%
    summarise(
      optimal_nbasis = nbasis[which.min(mse_true)],
      min_mse = min(mse_true),
      max_mse = max(mse_true),
      .groups = "drop"
    )

  # Steepness metrics with asymmetry
  steepness_data <- oracle_models %>%
    left_join(oracle_optimal %>% select(comp_no, optimal_nbasis, min_mse), by = "comp_no") %>%
    mutate(
      dist_from_optimal = nbasis - optimal_nbasis,
      pct_from_min = (mse_true / min_mse - 1) * 100
    )

  # Left vs right asymmetry
  penalty_left_6 <- steepness_data %>%
    filter(dist_from_optimal == -6) %>%
    summarise(avg = mean(pct_from_min, na.rm = TRUE)) %>%
    pull(avg)

  penalty_right_6 <- steepness_data %>%
    filter(dist_from_optimal == 6) %>%
    summarise(avg = mean(pct_from_min, na.rm = TRUE)) %>%
    pull(avg)

  penalty_at_10 <- steepness_data %>%
    filter(abs(abs(dist_from_optimal) - 9) <= 1.5) %>%
    summarise(avg = mean(pct_from_min, na.rm = TRUE)) %>%
    pull(avg)

  data.frame(
    config = config_name,
    mean_optimal = round(mean(oracle_optimal$optimal_nbasis), 1),
    sd_optimal = round(sd(oracle_optimal$optimal_nbasis), 1),
    range_min = min(oracle_optimal$optimal_nbasis),
    range_max = max(oracle_optimal$optimal_nbasis),
    boundary_pct = round(100 * mean(oracle_optimal$optimal_nbasis %in% c(3, 60)), 1),
    mse_variation = round(100 * mean((oracle_optimal$max_mse - oracle_optimal$min_mse) / oracle_optimal$min_mse), 1),
    penalty_left_6 = round(penalty_left_6, 1),
    penalty_right_6 = round(penalty_right_6, 1),
    asymmetry_ratio = round(penalty_left_6 / penalty_right_6, 2),
    penalty_at_10 = round(penalty_at_10, 1)
  )
}

oracle_properties <- bind_rows(lapply(configs, function(cfg) {
  oracle_df <- all_results[[cfg$name]]$oracle
  compute_oracle_properties(oracle_df, cfg$name)
}))

print(oracle_properties)
cat("\n")

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

cat("=== SAVING RESULTS ===\n\n")

output_dir <- "results_multi_config"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(oracle_all, file.path(output_dir, "oracle_all.RDS"))
saveRDS(dt_all, file.path(output_dir, "dt_all.RDS"))
saveRDS(oracle_properties, file.path(output_dir, "oracle_properties.RDS"))

write.csv(oracle_properties, file.path(output_dir, "oracle_properties.csv"), row.names = FALSE)

cat(sprintf("✓ Results saved to %s/\n\n", output_dir))

cat("================================================================================\n")
cat("  SUMMARY COMPLETE\n")
cat("================================================================================\n\n")

cat("Next steps:\n")
cat("  - Review oracle_properties.csv for config comparison\n")
cat("  - Run plotting script (to be created)\n\n")
