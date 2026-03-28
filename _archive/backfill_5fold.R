#!/usr/bin/env Rscript
# ==============================================================================
# Backfill 5-fold DT for equal_40 and equal_75 (seeds 21-50)
# ==============================================================================
# These designs only had 5-fold data for seeds 1-20.
# Needed for Section 4 variance ratio table.
# Setup + full-data fits already exist for all 50 seeds.
# ==============================================================================

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

source("sim_functions/sampling_and_setup.R")
source("sim_functions/full_data_fit.R")
source("sim_functions/run_dt.R")

N_CORES <- 11
START <- 21
END <- 50

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

designs <- list(
  list(equal_n = 40, results_dir = "_results_equal40"),
  list(equal_n = 75, results_dir = "_results_equal75")
)

for (des in designs) {
  cat(sprintf("\n=== %s: 5-fold backfill, seeds %d-%d ===\n", des$results_dir, START, END))

  config <- c(base_config, des)
  config$n_comp <- END

  t1 <- proc.time()
  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)

  foreach(comp_no = START:END) %dopar% {
    run_dt(comp_no, k = 5, results_dir = des$results_dir,
           model_config = config, n_reps = 1)
  }

  stopCluster(cl)
  elapsed <- proc.time() - t1
  cat(sprintf("Done (%.1f sec)\n", elapsed[3]))
}

# Verify
cat("\n=== Verification ===\n")
for (des in designs) {
  count <- sum(file.exists(
    sprintf("%s/comparison_%03d/dt_5fold/dt_split_5fold_rep1.RDA", des$results_dir, 1:50)
  ))
  cat(sprintf("%s: %d/50 seeds with 5-fold data\n", des$results_dir, count))
}
