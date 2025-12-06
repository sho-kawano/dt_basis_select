#!/usr/bin/env Rscript
# Run DT 1-fold with epsilon=0.3 for 4 configs missing this value

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

configs <- list(
  list(name = "equal_100", results_dir = "_results_equal100_comparison"),
  list(name = "prop_1p25pct", results_dir = "_results_prop1p25pct_comparison"),
  list(name = "prop_1p5pct", results_dir = "_results_prop1p5pct_comparison"),
  list(name = "prop_1p75pct", results_dir = "_results_prop1p75pct_comparison")
)

n_comp <- 20
eps_value <- 0.3
n_reps <- 5

source("sim_functions/run_dt.R")

for (cfg in configs) {
  cat(sprintf("\n=== %s ===\n", cfg$name))

  model_config <- readRDS(file.path(cfg$results_dir, "model_config.RDS"))

  cl <- makeForkCluster(model_config$n_cores)
  registerDoParallel(cl)

  foreach(comp_no = 1:n_comp) %dopar% {
    run_dt(comp_no, k = 1, results_dir = cfg$results_dir,
           model_config = model_config, eps = eps_value, n_reps = n_reps)
  }

  stopCluster(cl)
}

cat("\nNext: ./aggregate_additional_configs.R\n")
