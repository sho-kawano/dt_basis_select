#!/usr/bin/env Rscript
# ==============================================================================
# backfill_esim_methodcomp.R — One-off ESIM backfill for Section 6.3
# ==============================================================================
# Backfills ESIM to cover seeds 1-75 for all three Opt5 trio designs.
#
# After the initial backfill (seeds 1-50 for prop_1p5pct & prop_2p25pct),
# all three designs have ESIM for seeds 1-50. This script fills 51-75.
#
# Steps 1-3 (setup, full_data_fit, DT) already ran for all seeds 1-300.
#
# Runtime: ~10 hours (25 seeds × 3 designs)
# After: run aggregate_methodcomp.R
# ==============================================================================

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

N_CORES    <- 11
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
  list(name = "prop_0.75pct", samp_frac = 0.0075, results_dir = "_results_prop0.75pct",
       seeds = 51:75),
  list(name = "prop_1p5pct",  samp_frac = 0.0150, results_dir = "_results_prop1p5pct",
       seeds = 51:75),
  list(name = "prop_2p25pct", samp_frac = 0.0225, results_dir = "_results_prop2p25pct",
       seeds = 51:75)
)

source("sim_functions/sampling_and_setup.R")
source("sim_functions/run_esim.R")

for (d in DESIGNS) {
  cat(sprintf("\n================================================================================\n"))
  cat(sprintf("  ESIM BACKFILL: %s  (seeds %d-%d)\n", toupper(d$name),
              min(d$seeds), max(d$seeds)))
  cat(sprintf("================================================================================\n\n"))

  model_config             <- base_config
  model_config$samp_frac   <- d$samp_frac
  model_config$n_comp      <- max(d$seeds)
  model_config$n_iter_esim <- ESIM_ITERS

  # --- Setup ESIM (generate w.RDS files) --------------------------------------
  cat("Setup ESIM...\n")
  t0 <- proc.time()
  setup_esim(ncomps = max(d$seeds), results_dir = d$results_dir,
             niters = ESIM_ITERS, comp_nos = d$seeds)
  cat(sprintf("  Done (%.1f sec)\n\n", (proc.time() - t0)[3]))

  # --- Run ESIM ---------------------------------------------------------------
  cat(sprintf("ESIM (L=%d)...\n", ESIM_ITERS))
  cat(sprintf("  %d iters x %d models x %d seeds = %d fits\n",
              ESIM_ITERS, length(base_config$nbasis_values), length(d$seeds),
              ESIM_ITERS * length(base_config$nbasis_values) * length(d$seeds)))
  t0 <- proc.time()
  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)
  foreach(comp_no = d$seeds) %dopar% {
    run_esim(comp_no, n_cores = 1, results_dir = d$results_dir,
             model_config = model_config, n_iters = ESIM_ITERS)
  }
  stopCluster(cl)
  cat(sprintf("  Done (%.1f sec = %.1f hours)\n\n", (proc.time() - t0)[3],
              (proc.time() - t0)[3] / 3600))

  cat(sprintf("%s COMPLETE\n", toupper(d$name)))
}

cat("\n================================================================================\n")
cat("  ESIM BACKFILL COMPLETE\n")
cat("  Next: run aggregate_methodcomp.R\n")
cat("================================================================================\n\n")
