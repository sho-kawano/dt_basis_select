#!/usr/bin/env Rscript
# ==============================================================================
# PA designs S=300 run — Section 6.3
# ==============================================================================
# Designs + seeds:
#   prop_0.75pct  : seeds 51-300 (extends existing _results_prop0.75pct,  seeds 1-50)
#   prop_1p25pct  : seeds 51-300 (extends existing _results_prop1p25pct,  seeds 1-50)
#   prop_1p5pct   : seeds  1-300 (new clean dir   _results_prop1p5pct)
#   prop_2p25pct  : seeds 51-300 (extends existing _results_prop2p25pct,  seeds 1-50)
#
# eps=0.7, R=5 only. ESIM off.
# Run fastest -> slowest (0.75, 1.25, 1.5, 2.25) to front-load.
# ==============================================================================

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

N_CORES   <- 11
EPSILON   <- 0.7
N_REPS_DT <- 5

base_config <- list(
  population_file  = "data/ca_pums_population.rds",
  adjacency_file   = "data/ca_puma_adjacency.RDA",
  response_var     = "employed",
  response_type    = "binary",
  response_filter  = NULL,
  X_covariates     = NULL,
  nbasis_values    = seq(3, 60, by = 3),
  spatial_type     = "fixed",
  hyp              = list(c = 0.001, d = 0.001),
  ndesired         = 2000,
  nburn            = 1500,
  nthin            = 1,
  n_cores          = N_CORES,
  equal_allocation = FALSE,
  min_sample_size  = 10
)

DESIGNS <- list(
  list(name = "prop_0.75pct",  samp_frac = 0.0075, results_dir = "_results_prop0.75pct",
       seeds = 51:300),
  list(name = "prop_1p25pct",  samp_frac = 0.0125, results_dir = "_results_prop1p25pct",
       seeds = 51:300),
  list(name = "prop_1p5pct",   samp_frac = 0.0150, results_dir = "_results_prop1p5pct",
       seeds = 1:300),
  list(name = "prop_2p25pct",  samp_frac = 0.0225, results_dir = "_results_prop2p25pct",
       seeds = 51:300)
)

source("sim_functions/sampling_and_setup.R")
source("sim_functions/full_data_fit.R")
source("sim_functions/run_dt.R")

for (d in DESIGNS) {
  cat(sprintf("\n================================================================================\n"))
  cat(sprintf("  %s  (seeds %d-%d)\n", toupper(d$name), min(d$seeds), max(d$seeds)))
  cat(sprintf("================================================================================\n\n"))

  model_config             <- base_config
  model_config$samp_frac   <- d$samp_frac
  model_config$n_comp      <- max(d$seeds)

  dir.create(d$results_dir, showWarnings = FALSE)

  # --- Setup ------------------------------------------------------------------
  cat("Setup...\n")
  t0 <- proc.time()
  setup_comp(ncomps      = max(d$seeds),
             results_dir = d$results_dir,
             model_config = model_config,
             comp_nos    = d$seeds)
  cat(sprintf("  Done (%.1f sec)\n\n", (proc.time() - t0)[3]))

  # --- Full data fits ----------------------------------------------------------
  cat("Full data fits...\n")
  t0 <- proc.time()
  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)
  foreach(comp_no = d$seeds) %dopar% {
    full_data_fit(comp_no, d$results_dir, model_config)
  }
  stopCluster(cl)
  cat(sprintf("  Done (%.1f sec)\n\n", (proc.time() - t0)[3]))

  # --- Data thinning (eps=0.7, R=5) -------------------------------------------
  cat(sprintf("DT eps=%.1f R=%d...\n", EPSILON, N_REPS_DT))
  t0 <- proc.time()
  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)
  foreach(comp_no = d$seeds) %dopar% {
    run_dt(comp_no, k = 1, results_dir = d$results_dir,
           model_config = model_config, eps = EPSILON, n_reps = N_REPS_DT)
  }
  stopCluster(cl)
  cat(sprintf("  Done (%.1f sec)\n\n", (proc.time() - t0)[3]))

  cat(sprintf("%s COMPLETE\n", toupper(d$name)))
}

cat("\n================================================================================\n")
cat("  ALL DESIGNS COMPLETE\n")
cat("  Next: run aggregate_section63.R\n")
cat("================================================================================\n\n")
