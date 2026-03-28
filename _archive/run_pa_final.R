#!/usr/bin/env Rscript
# ==============================================================================
# PA Allocation Final Run — Section 6.3
# ==============================================================================
# Designs: prop_0.75pct, prop_1p25pct, prop_1p5pct
# Seeds: seq(4, 400, by=4) → S=100 per design
# New clean directories (fresh run from scratch)
# ==============================================================================

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

cat("\n")
cat("================================================================================\n")
cat("  PA ALLOCATION FINAL RUN (S=100)\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

PA_DESIGNS   <- c(#"prop_0.75pct",  # DONE
                  #"prop_1p25pct",  # DONE (DT complete; ESIM run separately)
                  "prop_1p5pct")
SEEDS        <- seq(4, 200, by = 4)   # first 50 seeds (run seq(204, 400, by=4) for second batch)
EPSILON_VALUES <- c(0.5, 0.6, 0.7)
N_REPS_DT    <- 5
ESIM_ITERS   <- 100
N_CORES      <- 11

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
  equal_allocation = FALSE,
  min_sample_size = 10
)

design_params <- list(
  `prop_0.75pct` = list(samp_frac = 0.0075, results_dir = "_results_prop0.75pct_final"),
  `prop_1p25pct` = list(samp_frac = 0.0125, results_dir = "_results_prop1p25pct_final"),
  `prop_1p5pct`  = list(samp_frac = 0.0150, results_dir = "_results_prop1p5pct_final")
)

cat("Configuration:\n")
cat(sprintf("  Designs: %s\n", paste(PA_DESIGNS, collapse = ", ")))
cat(sprintf("  Seeds: seq(4, 400, by=4) → %d per design\n", length(SEEDS)))
cat(sprintf("  Epsilon values: %s\n", paste(EPSILON_VALUES, collapse = ", ")))
cat(sprintf("  DT reps: %d\n", N_REPS_DT))
cat(sprintf("  ESIM iterations: %d\n", ESIM_ITERS))
cat(sprintf("  Cores: %d\n\n", N_CORES))

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

  params       <- design_params[[design]]
  model_config <- c(base_config, params)
  model_config$n_comp      <- max(SEEDS)
  model_config$n_iter_esim <- ESIM_ITERS

  results_dir <- params$results_dir

  # --------------------------------------------------------------------------
  # PART 1: SETUP (create comparison directories + data)
  # --------------------------------------------------------------------------

  cat("=== PART 1: SETUP ===\n\n")
  cat(sprintf("Creating comparison directories for %d seeds...\n", length(SEEDS)))

  t1 <- proc.time()
  setup_comp(
    ncomps      = max(SEEDS),
    results_dir = results_dir,
    model_config = model_config,
    comp_nos    = SEEDS
  )
  elapsed <- proc.time() - t1
  cat(sprintf("Done (%.1f sec)\n\n", elapsed[3]))

  # --------------------------------------------------------------------------
  # PART 2: FULL DATA FITS
  # --------------------------------------------------------------------------

  cat("=== PART 2: FULL DATA FITS ===\n\n")

  t1 <- proc.time()
  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)

  foreach(comp_no = SEEDS) %dopar% {
    full_data_fit(comp_no, results_dir, model_config)
  }

  stopCluster(cl)
  elapsed <- proc.time() - t1
  cat(sprintf("Done (%.1f sec, ~%.1f sec/comp)\n\n", elapsed[3], elapsed[3] / length(SEEDS)))

  # --------------------------------------------------------------------------
  # PART 3: DATA THINNING (1-fold)
  # --------------------------------------------------------------------------

  cat("=== PART 3: DATA THINNING ===\n\n")

  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)

  for (eps in EPSILON_VALUES) {
    cat(sprintf("  epsilon=%.1f, R=%d...\n", eps, N_REPS_DT))
    t1 <- proc.time()

    foreach(comp_no = SEEDS) %dopar% {
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
  # PART 4: SETUP ESIM  [SKIPPED — run separately later]
  # --------------------------------------------------------------------------

  # setup_esim(
  #   ncomps      = max(SEEDS),
  #   results_dir = results_dir,
  #   niters      = ESIM_ITERS,
  #   comp_nos    = SEEDS
  # )

  # --------------------------------------------------------------------------
  # PART 5: RUN ESIM  [SKIPPED — run separately later]
  # --------------------------------------------------------------------------

  # cl <- makeCluster(N_CORES, type = "FORK")
  # registerDoParallel(cl)
  # foreach(comp_no = SEEDS) %dopar% {
  #   run_esim(comp_no, n_cores = 1, results_dir = results_dir,
  #            model_config = model_config, n_iters = ESIM_ITERS)
  # }
  # stopCluster(cl)

  # Save config
  model_config$seeds       <- SEEDS
  model_config$eps_values  <- EPSILON_VALUES
  model_config$n_reps_dt   <- N_REPS_DT
  saveRDS(model_config, file.path(results_dir, "model_config.RDS"))
  cat(sprintf("Config saved: %s/model_config.RDS\n", results_dir))

  cat(sprintf("\n%s COMPLETE\n", toupper(design)))
}

cat("\n")
cat("================================================================================\n")
cat("  ALL DESIGNS COMPLETE\n")
cat("================================================================================\n")
cat("\nNext: run aggregate_pa_final.R to aggregate results.\n\n")
