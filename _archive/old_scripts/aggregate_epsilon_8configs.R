#!/usr/bin/env Rscript
# ==============================================================================
# Aggregate Epsilon Sensitivity Analysis Results (7 Configs)
# ==============================================================================
# Purpose: Aggregate all DT 1-fold results with full epsilon grid (0.1-0.9)
# Configs: equal_20, 30, 40, 50, 75, 100, 125
# DT configs: 9 eps × 3 n_reps × 2 loss = 54 per sample size
# ==============================================================================

library(parallel)
library(doParallel)
library(tidyverse)

cat("\n")
cat("================================================================================\n")
cat("  AGGREGATE EPSILON SENSITIVITY ANALYSIS (8 CONFIGS)\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# All 8 equal allocation configs
all_configs <- list(
  list(name = "equal_20", dir = "_results_equal20_epsilon", n_comp = 20, n_cores = 11),
  list(name = "equal_30", dir = "_results_equal30_epsilon", n_comp = 20, n_cores = 11),
  list(name = "equal_40", dir = "_results_equal40_comparison", n_comp = 20, n_cores = 11),
  list(name = "equal_50", dir = "_results_equal50_comparison", n_comp = 20, n_cores = 11),
  list(name = "equal_75", dir = "_results_equal75_rerun", n_comp = 20, n_cores = 11),
  list(name = "equal_100", dir = "_results_equal100_comparison", n_comp = 20, n_cores = 11),
  list(name = "equal_125", dir = "_results_equal125_epsilon", n_comp = 20, n_cores = 11),
  list(name = "equal_150", dir = "_results_equal150_epsilon", n_comp = 20, n_cores = 11)
)

eps_values <- seq(0.1, 0.9, by = 0.1)  # 9 values
n_reps_values <- c(1, 3, 5)           # 3 values
loss_functions <- c("MSE", "plugin_NLL")  # 2 values
# Total: 9 × 3 × 2 = 54 DT configurations per config

cat("Configuration:\n")
cat(sprintf("  Configs: %d\n", length(all_configs)))
cat(sprintf("  Comparisons per config: %d\n", all_configs[[1]]$n_comp))
cat(sprintf("  Epsilon values: %d (%s)\n", length(eps_values),
            paste(eps_values, collapse=", ")))
cat(sprintf("  n_reps values: %d (%s)\n", length(n_reps_values),
            paste(n_reps_values, collapse=", ")))
cat(sprintf("  Loss functions: %d (%s)\n", length(loss_functions),
            paste(loss_functions, collapse=", ")))
cat(sprintf("  Total DT configs per sample size: %d\n\n",
            length(eps_values) * length(n_reps_values) * length(loss_functions)))

# ==============================================================================
# AGGREGATE RESULTS
# ==============================================================================

dt_all_list <- list()
oracle_all_list <- list()

for (cfg_info in all_configs) {
  cat("\n")
  cat("================================================================================\n")
  cat(sprintf("  %s\n", toupper(cfg_info$name)))
  cat("================================================================================\n\n")

  if (!dir.exists(cfg_info$dir)) {
    cat(sprintf("⚠️  Directory not found: %s\n", cfg_info$dir))
    cat("Skipping this config...\n")
    next
  }

  # Check that model_config.RDS exists
  config_file <- file.path(cfg_info$dir, "model_config.RDS")
  if (!file.exists(config_file)) {
    cat(sprintf("⚠️  model_config.RDS not found in: %s\n", cfg_info$dir))
    cat("Skipping this config...\n")
    next
  }

  cat(sprintf("Processing: %s\n", cfg_info$dir))
  cat(sprintf("Comparisons: %d\n", cfg_info$n_comp))
  cat(sprintf("DT configs: %d\n\n",
              length(eps_values) * length(n_reps_values) * length(loss_functions)))

  # Create FORK cluster
  cl <- makeForkCluster(cfg_info$n_cores)
  registerDoParallel(cl)

  # --------------------------------------------------------------------------
  # DT 1-FOLD (ALL EPSILON × N_REPS × LOSS COMBINATIONS)
  # --------------------------------------------------------------------------

  cat("Aggregating DT 1-fold results...\n")

  # Create grid of all DT configurations
  dt_configs <- expand.grid(
    eps = eps_values,
    n_reps = n_reps_values,
    loss = loss_functions,
    stringsAsFactors = FALSE
  )

  cat(sprintf("  Total DT configs: %d\n", nrow(dt_configs)))

  dt_list <- list()
  dt_start <- proc.time()

  for (i in seq_len(nrow(dt_configs))) {
    dcfg <- dt_configs[i, ]

    # Progress indicator every 10 configs
    if (i %% 10 == 0) {
      cat(sprintf("    Progress: %d/%d (eps=%.1f, n_reps=%d, loss=%s)\n",
                  i, nrow(dt_configs), dcfg$eps, dcfg$n_reps, dcfg$loss))
    }

    result <- parLapply(cl, 1:cfg_info$n_comp, function(x, config_obj, res_dir) {
      source("sim_functions/summary_dt.R")
      summary_dt(x,
        n_folds = 1, loss_function = config_obj$loss, results_dir = res_dir,
        n_reps_to_use = config_obj$n_reps, eps = config_obj$eps
      ) %>%
        mutate(
          method_name = "dt_1fold",
          epsilon = config_obj$eps,
          n_reps_used = config_obj$n_reps
        )
    }, config_obj = dcfg, res_dir = cfg_info$dir)

    dt_list[[i]] <- bind_rows(result)
  }

  dt_elapsed <- proc.time() - dt_start
  cat(sprintf("✓ DT aggregation complete (%.1f sec)\n\n", dt_elapsed[3]))

  dt_1fold <- bind_rows(dt_list) %>% mutate(config = cfg_info$name)

  # --------------------------------------------------------------------------
  # ORACLE + FULL-DATA (DIC/WAIC)
  # --------------------------------------------------------------------------

  cat("Aggregating Oracle/DIC/WAIC results...\n")
  oracle_start <- proc.time()

  fulldata_oracle <- parLapply(cl, 1:cfg_info$n_comp, function(x, res_dir) {
    source("sim_functions/summary_oracle.R")
    summary_oracle(x, res_dir)
  }, res_dir = cfg_info$dir)

  fulldata_oracle <- bind_rows(fulldata_oracle) %>% mutate(config = cfg_info$name)

  oracle_elapsed <- proc.time() - oracle_start
  cat(sprintf("✓ Oracle/DIC/WAIC aggregation complete (%.1f sec)\n\n", oracle_elapsed[3]))

  stopCluster(cl)

  # Store results
  dt_all_list[[cfg_info$name]] <- dt_1fold
  oracle_all_list[[cfg_info$name]] <- fulldata_oracle

  # Summary
  cat("Summary for this config:\n")
  cat(sprintf("  DT rows: %d\n", nrow(dt_1fold)))
  cat(sprintf("  Oracle rows: %d\n", nrow(fulldata_oracle)))
  cat(sprintf("  Unique epsilon values: %d\n", n_distinct(dt_1fold$epsilon)))
  cat(sprintf("  Unique n_reps: %d\n", n_distinct(dt_1fold$n_reps_used)))
  cat(sprintf("  Unique loss functions: %d\n", n_distinct(dt_1fold$loss_function)))
  cat(sprintf("✓ %s aggregated successfully\n", cfg_info$name))
}

# ==============================================================================
# COMBINE AND SAVE
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  COMBINING ALL CONFIGS\n")
cat("================================================================================\n\n")

dt_all_combined <- bind_rows(dt_all_list)
oracle_all_combined <- bind_rows(oracle_all_list)

cat("Combined dataset summary:\n")
cat(sprintf("  Total configs: %d\n", n_distinct(dt_all_combined$config)))
cat(sprintf("  DT total rows: %d\n", nrow(dt_all_combined)))
cat(sprintf("  Oracle total rows: %d\n", nrow(oracle_all_combined)))
cat("\n")

# Create output directory
dir.create("results_multi_config", showWarnings = FALSE)

# Save
cat("Saving results...\n")
saveRDS(dt_all_combined, "results_multi_config/dt_epsilon_8configs.RDS")
saveRDS(oracle_all_combined, "results_multi_config/oracle_epsilon_8configs.RDS")

cat("\n")
cat("================================================================================\n")
cat("  COMPLETE!\n")
cat("================================================================================\n\n")

cat("✓ Saved:\n")
cat("  - results_multi_config/dt_epsilon_8configs.RDS\n")
cat("  - results_multi_config/oracle_epsilon_8configs.RDS\n\n")

cat("Dataset details:\n")
cat("\nConfigs included:\n")
print(sort(unique(dt_all_combined$config)))

cat("\nEpsilon values:\n")
print(sort(unique(dt_all_combined$epsilon)))

cat("\nn_reps values:\n")
print(sort(unique(dt_all_combined$n_reps_used)))

cat("\nLoss functions:\n")
print(sort(unique(dt_all_combined$loss_function)))

cat("\n")
cat("Next steps:\n")
cat("  1. Verify row counts match expectations\n")
cat("  2. Run Phase 4: Create analysis/epsilon_sensitivity.Rmd\n")
cat("  3. Analyze optimal epsilon by sample size\n\n")
