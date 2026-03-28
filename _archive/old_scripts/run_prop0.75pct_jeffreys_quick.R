#!/usr/bin/env Rscript
# ==============================================================================
# Quick test: prop_0.75pct design with Jeffreys prior on (beta, sigma^2)
#   Prior: pi(beta, sigma^2) ∝ 1/sigma^2  ->  sigma^2 ~ IG(c=0, d=0); beta flat
#   Spatial basis fixed; no ESIM; 5 comparisons; epsilon = 0.7; R = 5
#   Results saved to a new directory: _results_prop0.75pct_jeffreys_test
# ==============================================================================

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

cat("\n")
cat("================================================================================\n")
cat("  QUICK TEST: PROP 0.75% PA DESIGN WITH JEFFREYS PRIOR (NO ESIM)\n")
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
  response_filter = NULL,

  # Covariates
  X_covariates = NULL, # Intercept-only model

  # Spatial basis function settings
  nbasis_values = seq(3, 60, by = 3), # 20 models (all < 114 positive eigenvalues)

  # Model specification
  spatial_type = "fixed",

  # Priors: Jeffreys on (beta, sigma^2)
  # sigma^2 ~ IG(c=0, d=0); beta flat (implicit in spatial_basis_fh)
  hyp = list(c = 0, d = 0),

  # MCMC settings (kept modest for quick test)
  ndesired = 1500,
  nburn = 1000,
  nthin = 1,

  # Sampling settings
  equal_allocation = FALSE,  # proportional allocation
  samp_frac = 0.0075,        # 0.75% sampling rate
  min_sample_size = 10,
  equal_n = NULL,

  # Data thinning settings
  eps_values = c(0.7), # single epsilon for quick test
  n_reps_dt = 5,       # R = 5 repetitions

  # Computational settings
  n_cores = max(1, detectCores() - 2),

  # Output
  results_dir = "_results_prop0.75pct_jeffreys_test",

  # Study design
  n_comp = 20 # Use all 20 comparisons for confirmation
)

cat("Configuration:\n")
cat(sprintf("  Comparisons: %d\n", model_config$n_comp))
cat("  Design: proportional allocation, 0.75%\n")
cat("  Prior: sigma^2 ~ IG(0,0), beta flat (Jeffreys)\n")
cat(sprintf("  nbasis grid: %d-%d by %d (%d values)\n",
            min(model_config$nbasis_values), max(model_config$nbasis_values),
            diff(model_config$nbasis_values)[1], length(model_config$nbasis_values)))
cat(sprintf("  MCMC: nburn=%d, ndesired=%d\n", model_config$nburn, model_config$ndesired))
cat(sprintf("  DT eps: %s | DT reps: %d\n", paste(model_config$eps_values, collapse = ", "), model_config$n_reps_dt))
cat(sprintf("  Cores: %d\n", model_config$n_cores))
cat(sprintf("  Results dir: %s\n\n", model_config$results_dir))

# ==============================================================================
# PART 1: SETUP COMPARISONS
# ==============================================================================

cat("=== PART 1: SETUP COMPARISONS ===\n\n")
source("sim_functions/sampling_and_setup.R")
t1 <- proc.time()
setup_comp(ncomps = model_config$n_comp, results_dir = model_config$results_dir, model_config = model_config)
elapsed <- proc.time() - t1
cat(sprintf("✓ Setup done (%.1f sec)\n\n", elapsed[3]))

# ==============================================================================
# PART 2: FULL DATA FITS (Oracle + DIC/WAIC)
# ==============================================================================

cat("=== PART 2: FULL DATA FITS ===\n\n")
source("sim_functions/full_data_fit.R")
t1 <- proc.time()
cl <- makeCluster(model_config$n_cores, type = "FORK")
registerDoParallel(cl)
foreach(comp_no = 1:model_config$n_comp) %dopar% {
  full_data_fit(comp_no, model_config$results_dir, model_config)
}
stopCluster(cl)
elapsed <- proc.time() - t1
cat(sprintf("✓ Full data fits done (%.1f sec, ~%.1f sec/comp)\n\n", elapsed[3], elapsed[3] / model_config$n_comp))

# ==============================================================================
# PART 3: DATA THINNING (1-fold, R=5)
# ==============================================================================

cat("=== PART 3: DATA THINNING (1-fold, R=5) ===\n\n")
source("sim_functions/run_dt.R")
for (eps in model_config$eps_values) {
  cat(sprintf("  ε=%.1f ...\n", eps))
  t1 <- proc.time()
  cl <- makeCluster(model_config$n_cores, type = "FORK")
  registerDoParallel(cl)
  foreach(comp_no = 1:model_config$n_comp) %dopar% {
    run_dt(comp_no,
           k = 1,
           results_dir = model_config$results_dir,
           model_config = model_config,
           eps = eps,
           n_reps = model_config$n_reps_dt)
  }
  stopCluster(cl)
  elapsed <- proc.time() - t1
  cat(sprintf("    ✓ Done (%.1f sec, ~%.1f sec/comp)\n", elapsed[3], elapsed[3] / model_config$n_comp))
}
cat("\n")

cat("=== COMPLETE ===\n")
cat(sprintf("Outputs in: %s\n", model_config$results_dir))
cat("Next: aggregate/compare to default-prior runs as needed.\n")
