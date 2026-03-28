#!/usr/bin/env Rscript
# ==============================================================================
# Extend PA Allocation Results to S=50
# ==============================================================================
# Purpose: Extend PA allocation simulations from S=20 to S=50
# For: Section 6.3 (method comparison)
#
# Designs: prop_0.75pct, prop_1p25pct, prop_1p75pct
# Extends comparisons 21-50 (keeps existing 1-20)
# Includes ESIM (100 iterations per comparison)
#
# WARNING: ESIM is computationally expensive (~48+ hours)
# ==============================================================================

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

cat("\n")
cat("================================================================================\n")
cat("  EXTEND PA ALLOCATION TO S=50 (WITH ESIM)\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

PA_DESIGNS <- c("prop_0.75pct", "prop_1p25pct", "prop_1p75pct")
NCOMPS_TOTAL <- 50
NCOMPS_EXISTING <- 20
START_FROM <- NCOMPS_EXISTING + 1
EPSILON_VALUES <- c(0.5, 0.6, 0.7)
N_REPS_DT <- 5
ESIM_ITERS <- 100
N_CORES <- 11

# Base model config
base_config <- list(
  population_file = "data/ca_pums_population.rds",
  adjacency_file = "data/ca_puma_adjacency.RDA",
  response_var = "employed",
  response_type = "binary",
  response_filter = NULL,
  X_covariates = NULL,
  nbasis_values = seq(3, 60, by = 3),
  spatial_type = "fixed",
  hyp = list(c = 0.001, d = 0.001),
  ndesired = 2000,
  nburn = 1500,
  nthin = 1,
  n_cores = N_CORES,
  equal_allocation = FALSE,
  min_sample_size = 10
)

# Design-specific parameters
design_params <- list(
  `prop_0.75pct` = list(samp_frac = 0.0075, results_dir = "_results_prop0.75pct"),
  `prop_1p25pct` = list(samp_frac = 0.0125, results_dir = "_results_prop1p25pct"),
  `prop_1p75pct` = list(samp_frac = 0.0175, results_dir = "_results_prop1p75pct_comparison")
)

cat("Configuration:\n")
cat(sprintf("  Designs: %s\n", paste(PA_DESIGNS, collapse = ", ")))
cat(sprintf("  Total comparisons: %d (extending from %d)\n", NCOMPS_TOTAL, NCOMPS_EXISTING))
cat(sprintf("  New comparisons: %d-%d\n", START_FROM, NCOMPS_TOTAL))
cat(sprintf("  Epsilon values: %s\n", paste(EPSILON_VALUES, collapse = ", ")))
cat(sprintf("  DT reps: %d\n", N_REPS_DT))
cat(sprintf("  ESIM iterations: %d\n", ESIM_ITERS))
cat(sprintf("  Cores: %d\n\n", N_CORES))

cat("WARNING: ESIM is computationally expensive.\n")
cat("Estimated time: 48+ hours for all PA designs.\n\n")

# ==============================================================================
# LOAD FUNCTIONS
# ==============================================================================

source("sim_functions/sampling_and_setup.R")
source("sim_functions/full_data_fit.R")
source("sim_functions/run_dt.R")
source("sim_functions/run_esim.R")

# ==============================================================================
# MAIN LOOP: PROCESS EACH DESIGN
# ==============================================================================

for (design in PA_DESIGNS) {
  cat("\n")
  cat("================================================================================\n")
  cat(sprintf("  PROCESSING: %s\n", toupper(design)))
  cat("================================================================================\n\n")

  # Build config for this design
  params <- design_params[[design]]
  model_config <- c(base_config, params)
  model_config$n_comp <- NCOMPS_TOTAL
  model_config$n_iter_esim <- ESIM_ITERS

  results_dir <- params$results_dir

  # --------------------------------------------------------------------------
  # PART 1: SETUP NEW COMPARISONS (21-100)
  # --------------------------------------------------------------------------

  cat("=== PART 1: SETUP COMPARISONS ===\n\n")
  cat(sprintf("Setting up comparisons %d-%d...\n", START_FROM, NCOMPS_TOTAL))

  t1 <- proc.time()
  setup_comp(
    ncomps = NCOMPS_TOTAL,
    results_dir = results_dir,
    model_config = model_config,
    start_from = START_FROM
  )
  elapsed <- proc.time() - t1
  cat(sprintf("Done (%.1f sec)\n\n", elapsed[3]))

  # --------------------------------------------------------------------------
  # PART 2: FULL DATA FITS FOR NEW COMPARISONS
  # --------------------------------------------------------------------------

  cat("=== PART 2: FULL DATA FITS ===\n\n")
  cat(sprintf("Running full data fits for comparisons %d-%d...\n", START_FROM, NCOMPS_TOTAL))

  t1 <- proc.time()
  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)

  foreach(comp_no = START_FROM:NCOMPS_TOTAL) %dopar% {
    full_data_fit(comp_no, results_dir, model_config)
  }

  stopCluster(cl)
  elapsed <- proc.time() - t1
  cat(sprintf("Done (%.1f sec, ~%.1f sec/comp)\n\n", elapsed[3], elapsed[3] / (NCOMPS_TOTAL - NCOMPS_EXISTING)))

  # --------------------------------------------------------------------------
  # PART 3: DATA THINNING (1-fold) FOR NEW COMPARISONS
  # --------------------------------------------------------------------------

  cat("=== PART 3: DATA THINNING ===\n\n")

  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)

  for (eps in EPSILON_VALUES) {
    cat(sprintf("  epsilon=%.1f, R=%d, comparisons %d-%d...\n", eps, N_REPS_DT, START_FROM, NCOMPS_TOTAL))
    t1 <- proc.time()

    foreach(comp_no = START_FROM:NCOMPS_TOTAL) %dopar% {
      run_dt(comp_no,
        k = 1, results_dir = results_dir,
        model_config = model_config, eps = eps, n_reps = N_REPS_DT
      )
    }

    elapsed <- proc.time() - t1
    cat(sprintf("    Done (%.1f sec)\n", elapsed[3]))
  }

  stopCluster(cl)
  cat("\n")

  # --------------------------------------------------------------------------
  # PART 4: SETUP ESIM FOR NEW COMPARISONS
  # --------------------------------------------------------------------------

  cat("=== PART 4: SETUP ESIM ===\n\n")
  cat(sprintf("Generating synthetic direct estimates for comparisons %d-%d...\n", START_FROM, NCOMPS_TOTAL))

  t1 <- proc.time()
  setup_esim(
    ncomps = NCOMPS_TOTAL,
    results_dir = results_dir,
    niters = ESIM_ITERS,
    start_from = START_FROM
  )
  elapsed <- proc.time() - t1
  cat(sprintf("Done (%.1f sec)\n\n", elapsed[3]))

  # --------------------------------------------------------------------------
  # PART 5: RUN ESIM FOR NEW COMPARISONS
  # --------------------------------------------------------------------------

  cat("=== PART 5: RUN ESIM ===\n\n")
  cat(sprintf("Running ESIM for comparisons %d-%d...\n", START_FROM, NCOMPS_TOTAL))
  cat(sprintf("  (%d iterations x %d models x %d comparisons = %d total fits)\n",
    ESIM_ITERS, length(model_config$nbasis_values), NCOMPS_TOTAL - NCOMPS_EXISTING,
    ESIM_ITERS * length(model_config$nbasis_values) * (NCOMPS_TOTAL - NCOMPS_EXISTING)
  ))
  cat("  (This will take many hours...)\n\n")

  t1 <- proc.time()
  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)

  foreach(comp_no = START_FROM:NCOMPS_TOTAL) %dopar% {
    run_esim(comp_no,
      n_cores = 1,
      results_dir = results_dir,
      model_config = model_config,
      n_iters = ESIM_ITERS
    )
  }

  stopCluster(cl)
  elapsed <- proc.time() - t1
  cat(sprintf("Done (%.1f sec = %.1f hours)\n\n", elapsed[3], elapsed[3] / 3600))

  # --------------------------------------------------------------------------
  # UPDATE CONFIG
  # --------------------------------------------------------------------------

  model_config$eps_values <- EPSILON_VALUES
  model_config$n_reps_dt <- N_REPS_DT
  saveRDS(model_config, file.path(results_dir, "model_config.RDS"))
  cat(sprintf("Config saved: %s/model_config.RDS\n", results_dir))

  cat(sprintf("\n%s COMPLETE\n", toupper(design)))
}

# ==============================================================================
# VERIFICATION
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  VERIFICATION\n")
cat("================================================================================\n\n")

for (design in PA_DESIGNS) {
  results_dir <- design_params[[design]]$results_dir
  n_comps <- length(list.dirs(results_dir, recursive = FALSE))
  cat(sprintf("%s: %d comparisons\n", design, n_comps))
}

cat("\n")
cat("================================================================================\n")
cat("  PA ALLOCATION EXTENSION COMPLETE\n")
cat("================================================================================\n")
cat("\nNext steps:\n")
cat("1. Run aggregate_pa_allocation.R to aggregate results\n")
cat("2. Re-render section6_methods_comparison.Rmd\n\n")
