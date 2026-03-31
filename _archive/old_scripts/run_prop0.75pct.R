#!/usr/bin/env Rscript
# ==============================================================================
# Run prop_0.75pct (0.75% PA) Design
# ==============================================================================
# Purpose: Generate simulations for lowest sample size PA design
# For: Section 6.3 methods comparison
# ==============================================================================

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

cat("\n")
cat("================================================================================\n")
cat("  PROP 0.75% PA DESIGN (Section 6.3 Methods Comparison)\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

model_config <- list(
  # Data sources
  population_file = "data/ca_pums_population.rds",
  adjacency_file = "data/ca_puma_adjacency.RDA",

  # Response variable
  response_var = "employed",
  response_type = "binary",
  response_filter = NULL, # Employment-to-population ratio

  # Covariates
  X_covariates = NULL, # Intercept-only model

  # Spatial basis function settings
  nbasis_values = seq(3, 60, by = 3), # 20 models

  # Model specification
  spatial_type = "fixed", # Fixed spatial basis

  # Priors
  hyp = list(
    c = 0.001, # sigma^2 ~ IG(c, d)
    d = 0.001  # Weak priors
  ),

  # MCMC settings
  ndesired = 2000,
  nburn = 1500,
  nthin = 1,

  # Sampling settings
  equal_allocation = FALSE, # Proportional allocation
  samp_frac = 0.0075,      # 0.75% sampling rate
  min_sample_size = 10,
  equal_n = NULL,          # Not used for PA

  # Data thinning settings
  eps_values = c(0.5, 0.6, 0.7), # Section 6 recommended range
  n_reps_dt = 5,                  # R=5 repetitions

  # Empirical simulation settings
  n_iter_esim = 100, # ESIM iterations

  # Computational settings
  n_cores = 11, # Use 11 cores for overnight run

  # Output
  results_dir = "_results_prop0.75pct_comparison",

  # Study design
  n_comp = 20 # 20 comparisons
)

cat("Configuration:\n")
cat(sprintf("  State: California (281 PUMAs)\n"))
cat(sprintf("  Response: %s (employment-to-population ratio)\n", model_config$response_var))
cat(sprintf("  Model: spatial_basis_fh (spatial_type=%s)\n", model_config$spatial_type))
cat(sprintf("  Comparisons: %d\n", model_config$n_comp))
cat(sprintf("  Sampling: proportional allocation, 0.75%% sampling rate\n"))
cat(sprintf("  Expected mean n/PUMA: ~47 (range: ~25-82)\n"))
cat(sprintf(
  "  nbasis grid: %d-%d by %d (%d values)\n",
  min(model_config$nbasis_values), max(model_config$nbasis_values),
  diff(model_config$nbasis_values)[1], length(model_config$nbasis_values)
))
cat(sprintf("  Priors: c=%.0e, d=%.0e\n", model_config$hyp$c, model_config$hyp$d))
cat(sprintf("  MCMC: nburn=%d, ndesired=%d\n", model_config$nburn, model_config$ndesired))
cat(sprintf("  DT eps: %s\n", paste(model_config$eps_values, collapse = ", ")))
cat(sprintf("  DT reps: %d\n", model_config$n_reps_dt))
cat(sprintf("  Cores: %d (leaving 2 free)\n", model_config$n_cores))
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
# PART 2: FULL DATA FITS (Oracle + DIC/WAIC)
# ==============================================================================

cat("=== PART 2: FULL DATA FITS ===\n\n")

source("sim_functions/full_data_fit.R")

cat("Running full data fits (Oracle MSE, DIC, WAIC)...\n")
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

# ==============================================================================
# PART 3: DATA THINNING (1-fold)
# ==============================================================================

cat("=== PART 3: DATA THINNING ===\n\n")

source("sim_functions/run_dt.R")

cat("Running data thinning (1-fold, R=5)...\n")

for (eps in model_config$eps_values) {
  cat(sprintf("  ε=%.1f, R=%d...\n", eps, model_config$n_reps_dt))
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
  cat(sprintf("    ✓ Done (%.1f sec, ~%.1f sec/comp)\n", elapsed[3], elapsed[3] / model_config$n_comp))
}
cat("\n")

# ==============================================================================
# PART 4: SETUP EMPIRICAL SIMULATION
# ==============================================================================

cat("=== PART 4: SETUP EMPIRICAL SIMULATION ===\n\n")

cat("Generating synthetic direct estimates...\n")
t1 <- proc.time()
setup_esim(ncomps = model_config$n_comp, results_dir = model_config$results_dir, niters = model_config$n_iter_esim)
elapsed <- proc.time() - t1
cat(sprintf(
  "✓ ESIM setup complete (%d iterations × %d comparisons, %.1f sec)\n\n",
  model_config$n_iter_esim, model_config$n_comp, elapsed[3]
))

# ==============================================================================
# PART 5: RUN EMPIRICAL SIMULATION
# ==============================================================================

cat("=== PART 5: RUN EMPIRICAL SIMULATION ===\n\n")

source("sim_functions/run_esim.R")

cat("Fitting models on synthetic datasets...\n")
cat(sprintf(
  "  (%d iterations × %d models × %d comparisons = %d total fits)\n",
  model_config$n_iter_esim, length(model_config$nbasis_values), model_config$n_comp,
  model_config$n_iter_esim * length(model_config$nbasis_values) * model_config$n_comp
))
cat("  (This will take ~3-4 hours)\n\n")

t1 <- proc.time()

cl <- makeCluster(model_config$n_cores, type = "FORK")
registerDoParallel(cl)

foreach(comp_no = 1:model_config$n_comp) %dopar% {
  run_esim(comp_no,
           n_cores = 1,  # Each worker processes all iterations sequentially
           results_dir = model_config$results_dir,
           model_config = model_config,
           n_iters = model_config$n_iter_esim)
}

stopCluster(cl)
elapsed <- proc.time() - t1
cat(sprintf(
  "✓ ESIM complete (%.1f sec = %.1f hours, ~%.1f sec/comp)\n\n",
  elapsed[3], elapsed[3] / 3600, elapsed[3] / model_config$n_comp
))

# ==============================================================================
# SAVE CONFIGURATION AND COMPLETE
# ==============================================================================

# Save config
saveRDS(model_config, file.path(model_config$results_dir, "model_config.RDS"))
cat(sprintf("✓ Configuration saved: %s/model_config.RDS\n", model_config$results_dir))

cat("\n")
cat("================================================================================\n")
cat("  PROP 0.75% PA DESIGN COMPLETE (INCLUDING ESIM)\n")
cat("================================================================================\n")
cat(sprintf("Results saved to: %s\n", model_config$results_dir))
cat("\n")
cat("Next step: Update aggregation to include all 3 PA designs\n")
cat("\n")
