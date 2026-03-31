#!/usr/bin/env Rscript
# ==============================================================================
# Equal Allocation Results (S=50)
# ==============================================================================
# Purpose: Run equal allocation simulations for S=50
# For: Sections 3 (theory figures), 4 (repeated thinning), 6.2 (epsilon analysis)
#
# Designs: equal_30, equal_40, equal_50, equal_75, equal_100, equal_125
# Seeds 1-50 per design
# NO ESIM needed for equal allocation
#
# Section 4 note: dt_5fold runs only for designs in DT_5FOLD_DESIGNS
# ==============================================================================

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

cat("\n")
cat("================================================================================\n")
cat("  EQUAL ALLOCATION S=50\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

EQUAL_DESIGNS <- c("equal_30", "equal_40", "equal_50", "equal_75", "equal_100", "equal_125")
NCOMPS_TOTAL <- 50
NCOMPS_EXISTING <- 0
START_FROM <- 1
EPSILON_GRID <- seq(0.1, 0.9, 0.1)
N_REPS_DT <- 5
N_CORES <- 11

# Designs that need dt_5fold (for Section 4)
DT_5FOLD_DESIGNS <- c("equal_50", "equal_75", "equal_100")

# Base model config (will be modified per design)
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
  equal_allocation = TRUE
)

# Design-specific equal_n values
design_params <- list(
  equal_30 = list(equal_n = 30, results_dir = "_results_equal30"),
  equal_40 = list(equal_n = 40, results_dir = "_results_equal40"),
  equal_50 = list(equal_n = 50, results_dir = "_results_equal50"),
  equal_75 = list(equal_n = 75, results_dir = "_results_equal75"),
  equal_100 = list(equal_n = 100, results_dir = "_results_equal100"),
  equal_125 = list(equal_n = 125, results_dir = "_results_equal125")
)

cat("Configuration:\n")
cat(sprintf("  Designs: %s\n", paste(EQUAL_DESIGNS, collapse = ", ")))
cat(sprintf("  Total comparisons: %d (extending from %d)\n", NCOMPS_TOTAL, NCOMPS_EXISTING))
cat(sprintf("  New comparisons: %d-%d\n", START_FROM, NCOMPS_TOTAL))
cat(sprintf("  Epsilon grid: %s\n", paste(EPSILON_GRID, collapse = ", ")))
cat(sprintf("  DT reps: %d\n", N_REPS_DT))
cat(sprintf("  dt_5fold designs: %s\n", paste(DT_5FOLD_DESIGNS, collapse = ", ")))
cat(sprintf("  Cores: %d\n\n", N_CORES))

# ==============================================================================
# LOAD FUNCTIONS
# ==============================================================================

source("sim_functions/sampling_and_setup.R")
source("sim_functions/full_data_fit.R")
source("sim_functions/run_dt.R")

# ==============================================================================
# MAIN LOOP: PROCESS EACH DESIGN
# ==============================================================================

for (design in EQUAL_DESIGNS) {
  cat("\n")
  cat("================================================================================\n")
  cat(sprintf("  PROCESSING: %s\n", toupper(design)))
  cat("================================================================================\n\n")

  # Build config for this design
  params <- design_params[[design]]
  model_config <- c(base_config, params)
  model_config$n_comp <- NCOMPS_TOTAL

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
    start_from = START_FROM,
    seed_offset = 11
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

  cat("=== PART 3: DATA THINNING (1-fold) ===\n\n")

  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)

  for (eps in EPSILON_GRID) {
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
  # PART 4: DATA THINNING (5-fold) FOR SECTION 4 DESIGNS ONLY
  # --------------------------------------------------------------------------

  if (design %in% DT_5FOLD_DESIGNS) {
    cat("=== PART 4: DATA THINNING (5-fold) [Section 4] ===\n\n")
    cat(sprintf("Running dt_5fold for comparisons %d-%d...\n", START_FROM, NCOMPS_TOTAL))

    t1 <- proc.time()
    cl <- makeCluster(N_CORES, type = "FORK")
    registerDoParallel(cl)

    foreach(comp_no = START_FROM:NCOMPS_TOTAL) %dopar% {
      run_dt(comp_no,
        k = 5, results_dir = results_dir,
        model_config = model_config, n_reps = 1
      )
    }

    stopCluster(cl)
    elapsed <- proc.time() - t1
    cat(sprintf("Done (%.1f sec)\n\n", elapsed[3]))
  }

  # --------------------------------------------------------------------------
  # UPDATE CONFIG
  # --------------------------------------------------------------------------

  # Update and save config with new n_comp
  model_config$eps_values <- EPSILON_GRID
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

for (design in EQUAL_DESIGNS) {
  results_dir <- design_params[[design]]$results_dir
  n_comps <- length(list.dirs(results_dir, recursive = FALSE))
  cat(sprintf("%s: %d comparisons\n", design, n_comps))
}

cat("\n")
cat("================================================================================\n")
cat("  EQUAL ALLOCATION COMPLETE\n")
cat("================================================================================\n")
cat("\nNext steps:\n")
cat("1. Run aggregate_equal_allocation.R to aggregate results\n")
cat("2. Re-render analysis notebooks to update figures\n\n")
