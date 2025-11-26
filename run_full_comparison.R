#!/usr/bin/env Rscript
# ==============================================================================
# CA PRODUCTION RUN - 50 Comparisons
# ==============================================================================
# Purpose: Full production run for method comparison study
# Based on: CA employed, equal_75 (validated in test pipeline)
# ==============================================================================

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

cat("\n")
cat("================================================================================\n")
cat("  CA EMPLOYED: PRODUCTION RUN (50 COMPARISONS)\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION (Inline - explicit and auditable)
# ==============================================================================

model_config <- list(
  # Data sources
  population_file = "data/ca_pums_population.rds",
  adjacency_file = "data/ca_puma_adjacency.RDA",

  # Response variable
  response_var = "employed",
  response_type = "binary",
  response_filter = NULL, # NULL = employment-to-population ratio (includes children 0-15)

  # Covariates
  X_covariates = NULL, # Intercept-only model

  # Spatial basis function settings
  nbasis_values = seq(3, 60, by = 3), # Fine grid around expected optimal ~15

  # Model specification
  spatial_type = "fixed", # Fixed spatial basis (no random effects)

  # Priors
  hyp = list(
    c = 0.001, # sigma^2 ~ IG(c, d) for IID variance
    d = 0.001 # Weak priors for U-shaped curves
  ),

  # MCMC settings
  ndesired = 2000,
  nburn = 1500, # Validated convergence for CA PUMAs
  nthin = 1,

  # Sampling settings
  equal_allocation = TRUE,
  equal_n = 75, # Winner from sample size comparison
  samp_frac = 0.01, # Not used with equal allocation
  min_sample_size = 30,

  # Data thinning settings
  eps_values = seq(0.2, 0.8, by = 0.1), # many epsilon values
  n_reps_dt = 5, # Run 5 repetitions (summary will compare using 1, 3, 5)
  k_folds = 5, # Multi-fold DT

  # Empirical simulation settings
  n_iter_esim = 100, # Production setting (was 50 in test)

  # Computational settings
  n_cores = 11,

  # Output
  results_dir = "_results_ca_full_comparison",

  # Study design
  n_comp = 50 # PRODUCTION: 50 comparisons
)

cat("Configuration:\n")
cat(sprintf("  State: California (281 PUMAs)\n"))
cat(sprintf("  Response: %s (employment-to-population ratio)\n", model_config$response_var))
cat(sprintf("  Model: spatial_basis_fh (spatial_type=%s)\n", model_config$spatial_type))
cat(sprintf("  Comparisons: %d\n", model_config$n_comp))
cat(sprintf("  Sampling: equal_allocation, n=%d per PUMA\n", model_config$equal_n))
cat(sprintf(
  "  nbasis grid: %d-%d by %d (%d values)\n",
  min(model_config$nbasis_values), max(model_config$nbasis_values),
  diff(model_config$nbasis_values)[1], length(model_config$nbasis_values)
))
cat(sprintf("  Priors: c=%.0e, d=%.0e\n", model_config$hyp$c, model_config$hyp$d))
cat(sprintf("  MCMC: nburn=%d, ndesired=%d\n", model_config$nburn, model_config$ndesired))
cat(sprintf("  DT eps: %s\n", paste(model_config$eps_values, collapse = ", ")))
cat(sprintf("  DT reps: %d (summary will compare 1, 3, 5)\n", model_config$n_reps_dt))
cat(sprintf("  ESIM iters: %d\n", model_config$n_iter_esim))
cat(sprintf("  Cores: %d\n", model_config$n_cores))
cat(sprintf("  Results: %s\n\n", model_config$results_dir))

# ==============================================================================
# PART 1: SETUP COMPARISONS
# ==============================================================================

cat("=== PART 1: SETUP COMPARISONS ===\n\n")

source("sim_functions/sampling_and_setup.R")

cat("Setting up comparisons (sampling + direct estimates)...\n")
t1 <- proc.time()
setup_comp(ncomps = model_config$n_comp, results_dir = model_config$results_dir, model_config = model_config)
elapsed <- proc.time() - t1
cat(sprintf("✓ %d comparisons set up (%.1f sec)\n\n", model_config$n_comp, elapsed[3]))

# ==============================================================================
# PART 2: SETUP EMPIRICAL SIMULATION
# ==============================================================================

cat("=== PART 2: SETUP EMPIRICAL SIMULATION ===\n\n")

cat("Setting up empirical simulation (synthetic data)...\n")
t1 <- proc.time()
setup_esim(ncomps = model_config$n_comp, results_dir = model_config$results_dir, niters = model_config$n_iter_esim)
elapsed <- proc.time() - t1
cat(sprintf(
  "✓ Empirical simulation set up (%d iterations × %d comparisons, %.1f sec)\n\n",
  model_config$n_iter_esim, model_config$n_comp, elapsed[3]
))

# ==============================================================================
# PART 3: RUN METHODS
# ==============================================================================

cat("=== PART 3: RUN METHODS ===\n\n")

source("sim_functions/full_data_fit.R")
source("sim_functions/run_dt.R")
source("sim_functions/run_esim.R")

# --- Full Data Fit (OD-Oracle) ---
cat("Running full data fits (OD-Oracle)...\n")
t1 <- proc.time()

cl <- makeCluster(model_config$n_cores, type = "FORK")
registerDoParallel(cl)

foreach(comp_no = 1:model_config$n_comp) %dopar% {
  full_data_fit(comp_no, model_config$results_dir, model_config)
}

stopCluster(cl)
elapsed <- proc.time() - t1
cat(sprintf(
  "✓ Full data fits complete (%.1f sec, ~%.1f sec/comp)\n\n",
  elapsed[3], elapsed[3] / model_config$n_comp
))

# --- Data Thinning 1-fold ---
cat("Running data thinning (1-fold)...\n")

for (eps in model_config$eps_values) {
  cat(sprintf("  ε=%.2f, %d reps...\n", eps, model_config$n_reps_dt))
  t1 <- proc.time()

  cl <- makeCluster(model_config$n_cores, type = "FORK")
  registerDoParallel(cl)

  foreach(comp_no = 1:model_config$n_comp) %dopar% {
    run_dt(comp_no,
      k = 1, results_dir = model_config$results_dir,
      model_config = model_config, eps = eps, n_reps = model_config$n_reps_dt
    )
  }

  stopCluster(cl)
  elapsed <- proc.time() - t1
  cat(sprintf("    ✓ Done (%.1f sec)\n", elapsed[3]))
}
cat("\n")

# --- Data Thinning 5-fold ---
cat("Running data thinning (5-fold)...\n")
t1 <- proc.time()

cl <- makeCluster(model_config$n_cores, type = "FORK")
registerDoParallel(cl)

foreach(comp_no = 1:model_config$n_comp) %dopar% {
  run_dt(comp_no,
    k = model_config$k_folds, results_dir = model_config$results_dir,
    model_config = model_config, eps = NULL, n_reps = 1
  )
}

stopCluster(cl)
elapsed <- proc.time() - t1
cat(sprintf("✓ 5-fold DT complete (%.1f sec)\n\n", elapsed[3]))

# --- Empirical Simulation ---
cat("Running empirical simulation...\n")
t1 <- proc.time()

cl <- makeCluster(model_config$n_cores, type = "FORK")
registerDoParallel(cl)

foreach(comp_no = 1:model_config$n_comp) %dopar% {
  run_esim(comp_no, 1, model_config$results_dir,
    model_config,
    n_iters = model_config$n_iter_esim
  )
}

stopCluster(cl)
elapsed <- proc.time() - t1
cat(sprintf(
  "✓ Empirical simulation complete (%.1f sec, ~%.1f sec/comp)\n\n",
  elapsed[3], elapsed[3] / model_config$n_comp
))

# ==============================================================================
# SAVE CONFIGURATION AND COMPLETE
# ==============================================================================

# Save config
saveRDS(model_config, file.path(model_config$results_dir, "model_config.RDS"))
cat(sprintf("✓ Configuration saved: %s/model_config.RDS\n", model_config$results_dir))

cat("\n")
cat("================================================================================\n")
cat("  PRODUCTION RUN COMPLETE!\n")
cat("================================================================================\n\n")
cat("Next steps:\n")
cat(sprintf("  1. Run summaries: Rscript run_summary.R %s\n", model_config$results_dir))
cat("  2. Analyze results in results.RDS and results_list.RDS\n\n")
