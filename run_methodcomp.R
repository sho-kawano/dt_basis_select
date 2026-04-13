#!/usr/bin/env Rscript
# ==============================================================================
# run_methodcomp.R — Section 6.3 Method Comparison
# ==============================================================================
# Reproduces all results for the method comparison section.
# Trio: prop_0.75pct / prop_1p25pct / prop_1p75pct
# S=50 seeds per design, eps=0.6 R=5
#
# Pipeline per design:
#   1. setup_comp   — draw samples, compute z/d
#   2. full_data_fit — fit 20 candidate FH models on each sample
#   3. run_dt       — data thinning (eps=0.6, R=5)
#   4. setup_esim   — generate synthetic direct estimates
#   5. run_esim     — fit models on synthetic data (L=100)
#
# Runtime: ~2-3 days total (ESIM dominates at ~6h/design)
# Output dirs: _results_prop0.75pct, _results_prop1p25pct, _results_prop1p75pct
# Next: run aggregate_methodcomp.R
# ==============================================================================

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

N_CORES    <- 11
SEEDS      <- 1:50
EPSILON    <- 0.6
N_REPS_DT  <- 5
ESIM_ITERS <- 100

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
  list(name = "prop_0.75pct",  samp_frac = 0.0075, results_dir = "_results_prop0.75pct"),
  list(name = "prop_1p25pct",  samp_frac = 0.0125, results_dir = "_results_prop1p25pct"),
  list(name = "prop_1p75pct",  samp_frac = 0.0175, results_dir = "_results_prop1p75pct")
)

source("sim_functions/sampling_and_setup.R")
source("sim_functions/full_data_fit.R")
source("sim_functions/run_dt.R")
source("sim_functions/run_esim.R")

for (d in DESIGNS) {
  cat(sprintf("\n================================================================================\n"))
  cat(sprintf("  %s  (seeds %d-%d)\n", toupper(d$name), min(SEEDS), max(SEEDS)))
  cat(sprintf("================================================================================\n\n"))

  model_config             <- base_config
  model_config$samp_frac   <- d$samp_frac
  model_config$n_comp      <- max(SEEDS)
  model_config$n_iter_esim <- ESIM_ITERS

  dir.create(d$results_dir, showWarnings = FALSE)

  # --- 1. Setup ----------------------------------------------------------------
  cat("Setup...\n")
  t0 <- proc.time()
  setup_comp(ncomps = max(SEEDS), results_dir = d$results_dir,
             model_config = model_config, comp_nos = SEEDS)
  cat(sprintf("  Done (%.1f sec)\n\n", (proc.time() - t0)[3]))

  # --- 2. Full data fits -------------------------------------------------------
  cat("Full data fits...\n")
  t0 <- proc.time()
  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)
  foreach(comp_no = SEEDS) %dopar% {
    full_data_fit(comp_no, d$results_dir, model_config)
  }
  stopCluster(cl)
  cat(sprintf("  Done (%.1f sec)\n\n", (proc.time() - t0)[3]))

  # --- 3. Data thinning (eps=0.6, R=5) ----------------------------------------
  cat(sprintf("DT eps=%.1f R=%d...\n", EPSILON, N_REPS_DT))
  t0 <- proc.time()
  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)
  foreach(comp_no = SEEDS) %dopar% {
    run_dt(comp_no, k = 1, results_dir = d$results_dir,
           model_config = model_config, eps = EPSILON, n_reps = N_REPS_DT)
  }
  stopCluster(cl)
  cat(sprintf("  Done (%.1f sec)\n\n", (proc.time() - t0)[3]))

  # --- 4. Setup ESIM -----------------------------------------------------------
  cat("Setup ESIM...\n")
  t0 <- proc.time()
  setup_esim(ncomps = max(SEEDS), results_dir = d$results_dir,
             niters = ESIM_ITERS, comp_nos = SEEDS)
  cat(sprintf("  Done (%.1f sec)\n\n", (proc.time() - t0)[3]))

  # --- 5. Run ESIM -------------------------------------------------------------
  cat(sprintf("ESIM (L=%d)...\n", ESIM_ITERS))
  cat(sprintf("  %d iters x %d models x %d seeds = %d fits\n",
              ESIM_ITERS, length(base_config$nbasis_values), length(SEEDS),
              ESIM_ITERS * length(base_config$nbasis_values) * length(SEEDS)))
  t0 <- proc.time()
  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)
  foreach(comp_no = SEEDS) %dopar% {
    run_esim(comp_no, n_cores = 1, results_dir = d$results_dir,
             model_config = model_config, n_iters = ESIM_ITERS)
  }
  stopCluster(cl)
  cat(sprintf("  Done (%.1f sec = %.1f hours)\n\n", (proc.time() - t0)[3],
              (proc.time() - t0)[3] / 3600))

  cat(sprintf("%s COMPLETE\n", toupper(d$name)))
}

cat("\n================================================================================\n")
cat("  ALL DESIGNS COMPLETE\n")
cat("  Next: run aggregate_methodcomp.R\n")
cat("================================================================================\n\n")
