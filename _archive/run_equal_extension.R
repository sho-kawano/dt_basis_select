#!/usr/bin/env Rscript
# ==============================================================================
# Extend Equal Allocation to S=100 — n=30, 75, 125 only
# ==============================================================================
# Extends comparisons 51-100 for equal_30, equal_75, equal_125
# Runs: setup + full_data_fit + DT (eps=0.7, R=5 only — fast verification)
# Fill in remaining epsilons + ESIM separately after results confirmed
# ==============================================================================

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

cat("\n")
cat("================================================================================\n")
cat("  EXTEND EQUAL ALLOCATION TO S=100 (eps=0.7 only)\n")
cat("================================================================================\n\n")

DESIGNS <- list(
  list(name = "equal_30",  equal_n = 30,  results_dir = "_results_equal30"),
  list(name = "equal_75",  equal_n = 75,  results_dir = "_results_equal75"),
  list(name = "equal_125", equal_n = 125, results_dir = "_results_equal125")
)

START_FROM   <- 51
NCOMPS_TOTAL <- 100
EPSILON      <- 0.7
N_REPS_DT    <- 5
N_CORES      <- 11

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
  equal_allocation = TRUE,
  min_sample_size  = 10
)

cat(sprintf("Extending comparisons %d-%d for: %s\n",
  START_FROM, NCOMPS_TOTAL,
  paste(sapply(DESIGNS, `[[`, "name"), collapse=", ")))
cat(sprintf("DT: eps=%.1f, R=%d only\n\n", EPSILON, N_REPS_DT))

source("sim_functions/sampling_and_setup.R")
source("sim_functions/full_data_fit.R")
source("sim_functions/run_dt.R")

for (design in DESIGNS) {
  cat("\n")
  cat("================================================================================\n")
  cat(sprintf("  %s\n", toupper(design$name)))
  cat("================================================================================\n\n")

  model_config <- c(base_config, list(equal_n = design$equal_n,
                                       results_dir = design$results_dir,
                                       n_comp = NCOMPS_TOTAL))
  results_dir <- design$results_dir

  # --- Setup ---
  cat("Setup comparisons 51-100...\n")
  t1 <- proc.time()
  setup_comp(ncomps = NCOMPS_TOTAL, results_dir = results_dir,
             model_config = model_config, start_from = START_FROM)
  cat(sprintf("Done (%.1f sec)\n\n", (proc.time()-t1)[3]))

  # --- Full data fits ---
  cat("Full data fits 51-100...\n")
  t1 <- proc.time()
  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)
  foreach(comp_no = START_FROM:NCOMPS_TOTAL) %dopar% {
    full_data_fit(comp_no, results_dir, model_config)
  }
  stopCluster(cl)
  cat(sprintf("Done (%.1f sec)\n\n", (proc.time()-t1)[3]))

  # --- DT eps=0.7 only ---
  cat(sprintf("DT eps=%.1f, R=%d...\n", EPSILON, N_REPS_DT))
  t1 <- proc.time()
  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)
  foreach(comp_no = START_FROM:NCOMPS_TOTAL) %dopar% {
    run_dt(comp_no, k = 1, results_dir = results_dir,
           model_config = model_config, eps = EPSILON, n_reps = N_REPS_DT)
  }
  stopCluster(cl)
  cat(sprintf("Done (%.1f sec)\n\n", (proc.time()-t1)[3]))

  cat(sprintf("%s COMPLETE\n", toupper(design$name)))
}

cat("\n")
cat("================================================================================\n")
cat("  ALL DESIGNS COMPLETE\n")
cat("================================================================================\n")
cat("\nNext: run aggregate_equal_extension.R to check results.\n\n")
