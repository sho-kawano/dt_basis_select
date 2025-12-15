#!/usr/bin/env Rscript
# ==============================================================================
# 10-Config Sampling Design Comparison
# ==============================================================================
# Purpose: Run all 10 sampling designs for comprehensive DT vs DIC/WAIC analysis
# Configs: 4 equal allocation + 6 proportional PPS
# ==============================================================================

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

cat("\n")
cat("================================================================================\n")
cat("  10-CONFIG SAMPLING DESIGN COMPARISON\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Base config (shared across all configs)
base_config <- list(
  # Data sources
  population_file = "data/ca_pums_population.rds",
  adjacency_file = "data/ca_puma_adjacency.RDA",

  # Response variable
  response_var = "employed",
  response_type = "binary",
  response_filter = NULL,

  # Covariates
  X_covariates = NULL, # Intercept-only model
  X_approach = "population",

  # Spatial basis function settings
  nbasis_values = seq(3, 60, by = 3),

  # Model specification
  spatial_type = "fixed",

  # Priors
  hyp = list(
    c = 0.001,
    d = 0.001
  ),

  # MCMC settings
  ndesired = 2000,
  nburn = 1500,
  nthin = 1,

  # Computational settings
  n_cores = 8
)

n_comp <- 20

# ALL 10 CONFIGS
configs_to_run <- list(
  # Equal allocation (4 configs)
  list(
    name = "equal_40",
    equal_allocation = TRUE,
    equal_n = 40,
    results_dir = "_results_equal40_comparison"
  ),
  list(
    name = "equal_50",
    equal_allocation = TRUE,
    equal_n = 50,
    results_dir = "_results_equal50_comparison"
  ),
  list(
    name = "equal_75",
    equal_allocation = TRUE,
    equal_n = 75,
    results_dir = "_results_equal75_rerun"
  ),
  list(
    name = "equal_100",
    equal_allocation = TRUE,
    equal_n = 100,
    results_dir = "_results_equal100_comparison"
  ),
  # Proportional PPS (6 configs)
  list(
    name = "prop_0.5pct",
    equal_allocation = FALSE,
    samp_frac = 0.005,
    min_sample_size = 10,
    results_dir = "_results_prop0.5pct_comparison"
  ),
  list(
    name = "prop_1pct",
    equal_allocation = FALSE,
    samp_frac = 0.01,
    min_sample_size = 10,
    results_dir = "_results_prop1pct_comparison"
  ),
  list(
    name = "prop_1p25pct",
    equal_allocation = FALSE,
    samp_frac = 0.0125,
    min_sample_size = 10,
    results_dir = "_results_prop1p25pct_comparison"
  ),
  list(
    name = "prop_1p5pct",
    equal_allocation = FALSE,
    samp_frac = 0.015,
    min_sample_size = 10,
    results_dir = "_results_prop1p5pct_comparison"
  ),
  list(
    name = "prop_1p75pct",
    equal_allocation = FALSE,
    samp_frac = 0.0175,
    min_sample_size = 10,
    results_dir = "_results_prop1p75pct_comparison"
  ),
  list(
    name = "prop_2pct",
    equal_allocation = FALSE,
    samp_frac = 0.02,
    min_sample_size = 10,
    results_dir = "_results_prop2pct_comparison"
  )
)

cat("Configurations to run:\n")
cat("\nEqual Allocation:\n")
for (i in 1:4) {
  cfg <- configs_to_run[[i]]
  cat(sprintf("  %d. %-15s n=%d -> %s\n", i, cfg$name, cfg$equal_n, cfg$results_dir))
}
cat("\nProportional PPS:\n")
for (i in 5:10) {
  cfg <- configs_to_run[[i]]
  cat(sprintf("  %d. %-15s %.2f%% -> %s\n", i, cfg$name, cfg$samp_frac * 100, cfg$results_dir))
}
cat(sprintf("\nComparisons per config: %d\n", n_comp))
cat(sprintf("Total comparisons: %d\n", n_comp * length(configs_to_run)))
cat(sprintf("Cores: %d\n\n", base_config$n_cores))

# ==============================================================================
# SOURCE FUNCTIONS
# ==============================================================================

source("sim_functions/sampling_and_setup.R")
source("sim_functions/full_data_fit.R")
source("sim_functions/run_dt.R")

# ==============================================================================
# LOOP THROUGH CONFIGS
# ==============================================================================

for (config_idx in seq_along(configs_to_run)) {
  config_spec <- configs_to_run[[config_idx]]

  cat("\n")
  cat("================================================================================\n")
  cat(sprintf("  CONFIG %d/%d: %s\n", config_idx, length(configs_to_run), config_spec$name))
  cat("================================================================================\n\n")

  # Merge base config with specific config
  config <- c(base_config, config_spec)
  config$n_comp <- n_comp  # Explicitly set n_comp

  # Print config details
  if (config$equal_allocation) {
    cat(sprintf("  Sampling: Equal allocation, n=%d per PUMA\n", config$equal_n))
  } else {
    cat(sprintf("  Sampling: Proportional PPS, samp_frac=%.3f%%, min=%d\n",
                config$samp_frac * 100, config$min_sample_size))
  }
  cat(sprintf("  Comparisons: %d\n", n_comp))
  cat(sprintf("  Results: %s\n\n", config$results_dir))

  # --------------------------------------------------------------------------
  # PART 1: SETUP COMPARISONS
  # --------------------------------------------------------------------------

  cat("=== PART 1: SETUP COMPARISONS ===\n\n")
  cat("Setting up comparisons (sampling + direct estimates)...\n")
  t1 <- proc.time()
  setup_comp(ncomps = n_comp, results_dir = config$results_dir, model_config = config)
  elapsed <- proc.time() - t1
  cat(sprintf("✓ %d comparisons set up (%.1f sec)\n\n", n_comp, elapsed[3]))

  # --------------------------------------------------------------------------
  # PART 2: FULL DATA FITS
  # --------------------------------------------------------------------------

  cat("=== PART 2: FULL DATA FITS (OD-ORACLE) ===\n\n")
  cat("Running full data fits...\n")
  t1 <- proc.time()

  cl <- makeCluster(config$n_cores, type = "FORK")
  registerDoParallel(cl)

  foreach(comp_no = 1:n_comp) %dopar% {
    full_data_fit(comp_no, config$results_dir, config)
  }

  stopCluster(cl)
  elapsed <- proc.time() - t1
  cat(sprintf("✓ Full data fits complete (%.1f sec, ~%.1f sec/comp)\n\n",
              elapsed[3], elapsed[3] / n_comp))

  # --------------------------------------------------------------------------
  # PART 3: DATA THINNING
  # --------------------------------------------------------------------------

  cat("=== PART 3: DATA THINNING ===\n\n")

  # DT 3-fold
  cat("Running DT 3-fold...\n")
  t1 <- proc.time()

  cl <- makeCluster(config$n_cores, type = "FORK")
  registerDoParallel(cl)

  foreach(comp_no = 1:n_comp) %dopar% {
    run_dt(comp_no, k = 3, results_dir = config$results_dir,
           model_config = config, eps = NA, n_reps = 1)
  }

  stopCluster(cl)
  elapsed <- proc.time() - t1
  cat(sprintf("✓ DT 3-fold complete (%.1f sec)\n\n", elapsed[3]))

  # DT 5-fold
  cat("Running DT 5-fold...\n")
  t1 <- proc.time()

  cl <- makeCluster(config$n_cores, type = "FORK")
  registerDoParallel(cl)

  foreach(comp_no = 1:n_comp) %dopar% {
    run_dt(comp_no, k = 5, results_dir = config$results_dir,
           model_config = config, eps = NA, n_reps = 1)
  }

  stopCluster(cl)
  elapsed <- proc.time() - t1
  cat(sprintf("✓ DT 5-fold complete (%.1f sec)\n\n", elapsed[3]))

  # DT 1-fold (epsilon 0.1-0.9 × 5 reps)
  cat("Running DT 1-fold...\n")
  t1 <- proc.time()

  eps_values <- seq(0.1, 0.9, by = 0.1)  # 9 epsilon values
  n_reps <- 5

  cl <- makeCluster(config$n_cores, type = "FORK")
  registerDoParallel(cl)

  for (eps in eps_values) {
    cat(sprintf("  eps=%.1f, n_reps=%d\n", eps, n_reps))

    foreach(comp_no = 1:n_comp) %dopar% {
      run_dt(comp_no, k = 1, results_dir = config$results_dir,
             model_config = config, eps = eps, n_reps = n_reps)
    }
  }

  stopCluster(cl)
  elapsed <- proc.time() - t1
  cat(sprintf("✓ DT 1-fold complete (%.1f sec)\n\n", elapsed[3]))

  # Save configuration
  cat("Saving model configuration...\n")
  saveRDS(config, file.path(config$results_dir, "model_config.RDS"))
  cat(sprintf("✓ Configuration saved\n\n"))

  cat(sprintf("✓ Config %d/%d (%s) COMPLETE\n", config_idx, length(configs_to_run), config_spec$name))
}

# ==============================================================================
# COMPLETE
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  ALL 10 CONFIGS COMPLETE\n")
cat("================================================================================\n\n")

cat("Results saved to:\n")
for (cfg in configs_to_run) {
  cat(sprintf("  - %s\n", cfg$results_dir))
}

cat("\nNext steps:\n")
cat("  1. Aggregate results: ./aggregate_additional_configs.R\n")
cat("  2. Run analysis: ./analyze_all_dt_configs.R\n\n")
