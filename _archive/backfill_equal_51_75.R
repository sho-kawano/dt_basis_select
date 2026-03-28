#!/usr/bin/env Rscript
# ==============================================================================
# backfill_equal_51_75.R — Extend equal allocation from S=50 to S=75
# ==============================================================================
# All 6 equal designs currently have seeds 1-50. This adds seeds 51-75.
# Runs: setup + full_data_fit + DT (full epsilon grid, R=5) + dt_5fold
# dt_5fold for equal_50, equal_75, equal_100 (Section 4)
#
# Runtime: ~2-3 hours (25 seeds × 6 designs, fast relative to PA)
# After: re-run aggregate_equal_allocation.R
# ==============================================================================

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

source("sim_functions/sampling_and_setup.R")
source("sim_functions/full_data_fit.R")
source("sim_functions/run_dt.R")

N_CORES      <- 11
SEEDS        <- 51:75
EPSILON_GRID <- seq(0.1, 0.9, 0.1)
N_REPS_DT    <- 5

DT_5FOLD_DESIGNS <- c("equal_50", "equal_75", "equal_100")

base_config <- list(
  population_file = "data/ca_pums_population.rds",
  adjacency_file  = "data/ca_puma_adjacency.RDA",
  response_var    = "employed",
  response_type   = "binary",
  response_filter = NULL,
  X_covariates    = NULL,
  nbasis_values   = seq(3, 60, by = 3),
  spatial_type    = "fixed",
  hyp             = list(c = 0.001, d = 0.001),
  ndesired        = 2000,
  nburn           = 1500,
  nthin           = 1,
  n_cores         = N_CORES,
  equal_allocation = TRUE
)

designs <- list(
  list(name = "equal_30",  equal_n = 30,  results_dir = "_results_equal30"),
  list(name = "equal_40",  equal_n = 40,  results_dir = "_results_equal40"),
  list(name = "equal_50",  equal_n = 50,  results_dir = "_results_equal50"),
  list(name = "equal_75",  equal_n = 75,  results_dir = "_results_equal75"),
  list(name = "equal_100", equal_n = 100, results_dir = "_results_equal100"),
  list(name = "equal_125", equal_n = 125, results_dir = "_results_equal125")
)

for (des in designs) {
  cat(sprintf("\n================================================================================\n"))
  cat(sprintf("  %s  (seeds %d-%d)\n", toupper(des$name), min(SEEDS), max(SEEDS)))
  cat(sprintf("================================================================================\n\n"))

  config <- c(base_config, list(equal_n = des$equal_n, results_dir = des$results_dir))
  config$n_comp <- max(SEEDS)

  # --- Setup ------------------------------------------------------------------
  cat("Setup...\n")
  t0 <- proc.time()
  setup_comp(ncomps = max(SEEDS), results_dir = des$results_dir,
             model_config = config, comp_nos = SEEDS)
  cat(sprintf("  Done (%.1f sec)\n\n", (proc.time() - t0)[3]))

  # --- Full data fits ---------------------------------------------------------
  cat("Full data fits...\n")
  t0 <- proc.time()
  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)
  foreach(comp_no = SEEDS) %dopar% {
    full_data_fit(comp_no, des$results_dir, config)
  }
  stopCluster(cl)
  cat(sprintf("  Done (%.1f sec)\n\n", (proc.time() - t0)[3]))

  # --- DT 1-fold (full epsilon grid) -----------------------------------------
  cat("DT 1-fold...\n")
  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)
  for (eps in EPSILON_GRID) {
    cat(sprintf("  eps=%.1f R=%d...\n", eps, N_REPS_DT))
    t0 <- proc.time()
    foreach(comp_no = SEEDS) %dopar% {
      run_dt(comp_no, k = 1, results_dir = des$results_dir,
             model_config = config, eps = eps, n_reps = N_REPS_DT)
    }
    cat(sprintf("    Done (%.1f sec)\n", (proc.time() - t0)[3]))
  }
  stopCluster(cl)

  # --- DT 5-fold (Section 4 designs only) ------------------------------------
  if (des$name %in% DT_5FOLD_DESIGNS) {
    cat("\nDT 5-fold...\n")
    t0 <- proc.time()
    cl <- makeCluster(N_CORES, type = "FORK")
    registerDoParallel(cl)
    foreach(comp_no = SEEDS) %dopar% {
      run_dt(comp_no, k = 5, results_dir = des$results_dir,
             model_config = config, n_reps = 1)
    }
    stopCluster(cl)
    cat(sprintf("  Done (%.1f sec)\n", (proc.time() - t0)[3]))
  }

  cat(sprintf("\n%s COMPLETE\n", toupper(des$name)))
}

cat("\n================================================================================\n")
cat("  EQUAL ALLOCATION BACKFILL COMPLETE (seeds 51-75)\n")
cat("  Next: re-run aggregate_equal_allocation.R\n")
cat("================================================================================\n\n")
