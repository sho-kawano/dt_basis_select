#!/usr/bin/env Rscript
# ==============================================================================
# Re-run equal_75 with DT 3-fold
# ==============================================================================
# Purpose: Re-run equal_75 config with 20 comparisons to match other configs
#          and include DT 3-fold method
# ==============================================================================

library(tidyverse)
library(survey)
library(parallel)
library(doParallel)
library(foreach)

cat("\n")
cat("================================================================================\n")
cat("  RE-RUN EQUAL_75 WITH DT 3-FOLD\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Model configuration (inline)
model_config <- list(
  # Data sources
  population_file = "data/ca_pums_population.rds",
  adjacency_file = "data/ca_puma_adjacency.RDA",

  # Response variable
  response_var = "employed",
  response_type = "binary",
  response_filter = NULL,

  # Covariates
  X_covariates = NULL,  # Intercept-only model
  X_approach = "population",

  # Sampling design
  equal_allocation = TRUE,
  equal_n = 75,

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
results_dir <- "_results_equal75_rerun"

cat("Configuration:\n")
cat(sprintf("  Sampling: Equal allocation, n=%d per PUMA\n", model_config$equal_n))
cat(sprintf("  Comparisons: %d\n", n_comp))
cat(sprintf("  Results: %s\n\n", results_dir))

# ==============================================================================
# SETUP
# ==============================================================================

source("sim_functions/sampling_and_setup.R")
source("sim_functions/full_data_fit.R")
source("sim_functions/run_dt.R")

# Create new results directory (DO NOT DELETE EXISTING)
if (dir.exists(results_dir)) {
  stop(sprintf("ERROR: Directory %s already exists. Choose a different name.", results_dir))
}
dir.create(results_dir, recursive = TRUE)

# Save model config
saveRDS(model_config, file.path(results_dir, "model_config.RDS"))
cat("✓ Model config saved\n\n")

# ==============================================================================
# PART 1: SETUP COMPARISONS
# ==============================================================================

cat("=== PART 1: SETUP COMPARISONS ===\n\n")
cat("Setting up comparisons (sampling + direct estimates)...\n")

t_start <- Sys.time()

setup_comp(ncomps = n_comp, results_dir = results_dir, model_config = model_config)

t_end <- Sys.time()
cat(sprintf("✓ %d comparisons set up (%.1f sec)\n\n",
            n_comp, as.numeric(difftime(t_end, t_start, units = "secs"))))

# ==============================================================================
# PART 2: FULL DATA FITS (OD-ORACLE)
# ==============================================================================

cat("=== PART 2: FULL DATA FITS (OD-ORACLE) ===\n\n")
cat("Running full data fits...\n")

t_start <- Sys.time()

cl <- makeCluster(model_config$n_cores, type = "FORK")
registerDoParallel(cl)

foreach(comp_no = 1:n_comp) %dopar% {
  full_data_fit(comp_no, results_dir, model_config)
}

stopCluster(cl)

t_end <- Sys.time()
elapsed <- as.numeric(difftime(t_end, t_start, units = "secs"))
cat(sprintf("✓ Full data fits complete (%.1f sec, ~%.1f sec/comp)\n\n",
            elapsed, elapsed / n_comp))

# ==============================================================================
# PART 3: DATA THINNING
# ==============================================================================

cat("=== PART 3: DATA THINNING ===\n\n")

# DT 3-fold
cat("Running DT 3-fold...\n")
t_start <- Sys.time()

cl <- makeCluster(model_config$n_cores, type = "FORK")
registerDoParallel(cl)

foreach(comp_no = 1:n_comp) %dopar% {
  run_dt(comp_no, k = 3, results_dir = results_dir,
         model_config = model_config, eps = NA, n_reps = 1)
}

stopCluster(cl)

t_end <- Sys.time()
cat(sprintf("✓ DT 3-fold complete (%.1f sec)\n\n",
            as.numeric(difftime(t_end, t_start, units = "secs"))))

# DT 5-fold
cat("Running DT 5-fold...\n")
t_start <- Sys.time()

cl <- makeCluster(model_config$n_cores, type = "FORK")
registerDoParallel(cl)

foreach(comp_no = 1:n_comp) %dopar% {
  run_dt(comp_no, k = 5, results_dir = results_dir,
         model_config = model_config, eps = NA, n_reps = 1)
}

stopCluster(cl)

t_end <- Sys.time()
cat(sprintf("✓ DT 5-fold complete (%.1f sec)\n\n",
            as.numeric(difftime(t_end, t_start, units = "secs"))))

# DT 1-fold (epsilon 0.3, 0.5, 0.7 with 5 reps each)
cat("Running DT 1-fold...\n")
t_start <- Sys.time()

epsilon_values <- c(0.3, 0.5, 0.7)
n_reps <- 5

cl <- makeCluster(model_config$n_cores, type = "FORK")
registerDoParallel(cl)

for (eps in epsilon_values) {
  cat(sprintf("  eps=%.1f, n_reps=%d\n", eps, n_reps))

  foreach(comp_no = 1:n_comp) %dopar% {
    run_dt(comp_no, k = 1, results_dir = results_dir,
           model_config = model_config, eps = eps, n_reps = n_reps)
  }
}

stopCluster(cl)

t_end <- Sys.time()
cat(sprintf("✓ DT 1-fold complete (%.1f sec)\n\n",
            as.numeric(difftime(t_end, t_start, units = "secs"))))

# ==============================================================================
# DONE
# ==============================================================================

cat("================================================================================\n")
cat("  EQUAL_75 RE-RUN COMPLETE\n")
cat("================================================================================\n\n")

cat("Results saved to: ", results_dir, "\n")
cat("Next steps:\n")
cat("  Rscript run_multi_summary.R\n\n")
