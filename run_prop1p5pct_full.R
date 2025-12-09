#!/usr/bin/env Rscript
# ==============================================================================
# PRODUCTION RUN: prop_1p5pct (1.5% sampling)
# ==============================================================================
# Purpose: Full comparison study for prop_1p5pct configuration
# Sample size: ~94 per PUMA (1.5% sampling fraction)
# Comparisons: 50 (can increase to 60-70 if time permits)
# Methods: DT (multiple eps/n_reps), DIC, WAIC
# ==============================================================================

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

cat("\n")
cat("================================================================================\n")
cat("  PRODUCTION RUN: prop_1p5pct (1.5% SAMPLING)\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

model_config <- list(
  # Study design
  n_comp = 50, # Total comparisons (can increase to 60-70)
  results_dir = "_results_prop1p5pct_production",

  # Data sources
  population_file = "data/ca_pums_population.rds",
  adjacency_file = "data/ca_puma_adjacency.RDA",

  # Sampling design (proportional PPS)
  equal_allocation = FALSE,
  samp_frac = 0.015, # 1.5% sampling fraction
  min_sample_size = 10,

  # Response variable
  response_var = "employed",
  response_type = "binary",
  response_filter = NULL, # NULL = employment-to-population ratio

  # Covariates
  X_covariates = NULL, # Intercept-only model
  X_approach = "population",

  # Spatial basis function settings
  nbasis_values = seq(3, 60, by = 3), # Grid: 3, 6, 9, ..., 60

  # Model specification
  spatial_type = "fixed",

  # Priors
  hyp = list(
    c = 0.001,
    d = 0.001
  ),

  # MCMC settings
  ndesired = 1500,
  nburn = 1500,
  nthin = 1,

  # Data thinning settings
  eps_values = rep(0.4, 0.7, by = 0.1), # Epsilon values to test
  n_reps_dt = c(1, 3, 5), # Number of repetitions
  k_folds = 5, # For k-fold DT

  # Empirical simulation settings
  n_iter_esim = 100,

  # Computational settings
  n_cores = 11  # Use 11 cores (one less than 12 total, keep system responsive)
)

# Quick summary (keep terse)
cat(sprintf(
  "Sampling: %.1f%% PPS; ~94/PUMA; comps=%d\n",
  model_config$samp_frac * 100, model_config$n_comp
))
cat(sprintf(
  "Response: %s; nbasis=%d values (%d-%d by %d)\n",
  model_config$response_var,
  length(model_config$nbasis_values),
  min(model_config$nbasis_values),
  max(model_config$nbasis_values),
  diff(model_config$nbasis_values)[1]
))
cat(sprintf(
  "DT eps: %s; reps: %s; k-fold: %d; cores: %d; out=%s/\n\n",
  paste(model_config$eps_values, collapse = ", "),
  paste(model_config$n_reps_dt, collapse = ", "),
  model_config$k_folds,
  model_config$n_cores,
  model_config$results_dir
))
cat(sprintf("ESIM iters: %d\n\n", model_config$n_iter_esim))

# Workload snapshot
cat(sprintf(
  "Fits: full=%d, dt1fold=%d, dtkfold=%d, esim_iters=%d; across %d nbasis\n\n",
  model_config$n_comp,
  model_config$n_comp * length(model_config$eps_values) *
    max(model_config$n_reps_dt),
  model_config$n_comp,
  model_config$n_comp * model_config$n_iter_esim,
  length(model_config$nbasis_values)
))

# ==============================================================================
# SOURCE FUNCTIONS
# ==============================================================================

source("sim_functions/sampling_and_setup.R")
source("sim_functions/full_data_fit.R")
source("sim_functions/run_dt.R")
source("sim_functions/run_esim.R")

# ==============================================================================
# PART 1: SETUP COMPARISONS
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  PART 1: SETUP COMPARISONS\n")
cat("================================================================================\n\n")

cat(sprintf(
  "Setting up %d comparisons (sampling + direct estimates)...\n",
  model_config$n_comp
))
t_start <- proc.time()

setup_comp(
  ncomps = model_config$n_comp,
  results_dir = model_config$results_dir,
  model_config = model_config
)

elapsed <- proc.time() - t_start
cat(sprintf("\n✓ Setup complete (%.1f sec)\n", elapsed[3]))

# ==============================================================================
# PART 2: SETUP EMPIRICAL SIMULATION
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  PART 2: SETUP EMPIRICAL SIMULATION\n")
cat("================================================================================\n\n")

cat(sprintf(
  "Setting up empirical simulation (%d iterations per comparison)...\n",
  model_config$n_iter_esim
))
t_start <- proc.time()
setup_esim(
  ncomps = model_config$n_comp,
  results_dir = model_config$results_dir,
  niters = model_config$n_iter_esim
)
elapsed <- proc.time() - t_start
cat(sprintf("\n✓ ESIM setup complete (%.1f sec)\n", elapsed[3]))

# ==============================================================================
# PART 3: FULL DATA FITS (Oracle + DIC + WAIC)
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  PART 3: FULL DATA FITS (Oracle + DIC + WAIC)\n")
cat("================================================================================\n\n")

cat(sprintf("Running full data fits for %d comparisons...\n", model_config$n_comp))
cat(sprintf("  - %d nbasis values per comparison\n", length(model_config$nbasis_values)))
cat(sprintf("  - Using %d cores\n\n", model_config$n_cores))

t_start <- proc.time()

cl <- makeCluster(model_config$n_cores, type = "FORK")
registerDoParallel(cl)

foreach(comp_no = 1:model_config$n_comp) %dopar% {
  full_data_fit(comp_no, model_config$results_dir, model_config)
}

stopCluster(cl)

elapsed <- proc.time() - t_start
cat(sprintf(
  "\n✓ Full data fits complete (%.1f min, ~%.1f sec/comp)\n",
  elapsed[3] / 60, elapsed[3] / model_config$n_comp
))

# ==============================================================================
# PART 4: DATA THINNING
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  PART 4: DATA THINNING\n")
cat("================================================================================\n\n")

# --------------------------------------------------------------------------
# DT 1-fold (multiple eps × n_reps combinations)
# --------------------------------------------------------------------------

cat("Running DT 1-fold...\n")
cat(sprintf("  Epsilon values: %s\n", paste(model_config$eps_values, collapse = ", ")))
cat(sprintf("  Repetitions requested: %s\n", paste(model_config$n_reps_dt, collapse = ", ")))

# Run once per eps using the max rep count; downstream summaries can reuse the first r reps
dt_reps_run <- max(model_config$n_reps_dt)
cat(sprintf(
  "  Repetitions run: %d (reuses first r reps for smaller settings)\n\n",
  dt_reps_run
))

t_start <- proc.time()

cl <- makeCluster(model_config$n_cores, type = "FORK")
registerDoParallel(cl)

for (eps in model_config$eps_values) {
  cat(sprintf("  eps=%.1f, n_reps=%d... ", eps, dt_reps_run))
  t_config <- proc.time()

  foreach(comp_no = 1:model_config$n_comp) %dopar% {
    run_dt(comp_no,
      k = 1,
      results_dir = model_config$results_dir,
      model_config = model_config,
      eps = eps, n_reps = dt_reps_run
    )
  }

  elapsed_config <- proc.time() - t_config
  cat(sprintf("done (%.1f min)\n", elapsed_config[3] / 60))
}

stopCluster(cl)

elapsed <- proc.time() - t_start
cat(sprintf("\n✓ DT 1-fold complete (%.1f min)\n", elapsed[3] / 60))

# --------------------------------------------------------------------------
# DT k-fold
# --------------------------------------------------------------------------

cat(sprintf("\nRunning DT %d-fold...\n", model_config$k_folds))
t_start <- proc.time()

cl <- makeCluster(model_config$n_cores, type = "FORK")
registerDoParallel(cl)

foreach(comp_no = 1:model_config$n_comp) %dopar% {
  run_dt(comp_no,
    k = model_config$k_folds,
    results_dir = model_config$results_dir,
    model_config = model_config,
    eps = NA, n_reps = 1
  )
}

stopCluster(cl)

elapsed <- proc.time() - t_start
cat(sprintf(
  "\n✓ DT %d-fold complete (%.1f min)\n",
  model_config$k_folds, elapsed[3] / 60
))

# ==============================================================================
# PART 5: EMPIRICAL SIMULATION
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  PART 5: EMPIRICAL SIMULATION\n")
cat("================================================================================\n\n")

cat(sprintf(
  "Running empirical simulation (%d iterations × %d comparisons)...\n",
  model_config$n_iter_esim, model_config$n_comp
))
t_start <- proc.time()

cl <- makeCluster(model_config$n_cores, type = "FORK")
registerDoParallel(cl)

foreach(comp_no = 1:model_config$n_comp) %dopar% {
  run_esim(
    comp_no = comp_no,
    n_cores = 1, # avoid nested parallel inside run_esim
    results_dir = model_config$results_dir,
    model_config = model_config,
    n_iters = model_config$n_iter_esim
  )
}

stopCluster(cl)

elapsed <- proc.time() - t_start
cat(sprintf(
  "\n✓ Empirical simulation complete (%.1f min, ~%.1f sec/comp)\n",
  elapsed[3] / 60, elapsed[3] / model_config$n_comp
))

# ==============================================================================
# SAVE CONFIGURATION
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  SAVING CONFIGURATION\n")
cat("================================================================================\n\n")

config_file <- file.path(model_config$results_dir, "model_config.RDS")
saveRDS(model_config, config_file)
cat(sprintf("✓ Configuration saved to: %s\n", config_file))

# ==============================================================================
# COMPLETE
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  PRODUCTION RUN COMPLETE!\n")
cat("================================================================================\n\n")

cat(sprintf("Results saved to: %s/\n\n", model_config$results_dir))

cat("Next steps:\n")
cat("  1. Run summary analysis:\n")
cat(sprintf("     Rscript run_summary.R %s\n\n", model_config$results_dir))

cat("  2. Generate plots and tables:\n")
cat("     # Open analysis/prop1p5pct_results.Rmd\n\n")

cat("Files to check:\n")
cat(sprintf("  - %s/results.RDS (aggregated results)\n", model_config$results_dir))
cat(sprintf("  - %s/model_config.RDS (configuration)\n", model_config$results_dir))
cat(sprintf(
  "  - %s/comparison_*/chains.RDS (MCMC chains)\n\n",
  model_config$results_dir
))
