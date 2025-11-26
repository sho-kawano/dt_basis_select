#!/usr/bin/env Rscript
# Quick test of screening logic on a single response variable

library(dplyr)

source("sim_functions/sampling_and_setup.R")
source("sim_functions/full_data_fit.R")
source("sim_functions/summary_oracle.R")

# Test with "employed" (highest Moran's I)
test_response <- data.frame(
  variable = "employed",
  label = "Employed",
  type = "binary",
  morans_i = 0.701,
  stringsAsFactors = FALSE
)

# Config
nbasis_grid <- c(10, 30, 50, 70, 90)
model_config <- list(
  equal_allocation = TRUE,
  equal_n = 50,
  min_sample_size = 30,
  samp_frac = 0.01,
  response_var = "employed",
  response_type = "binary",
  response_filter = NULL,
  X_covariates = NULL,
  c_prior = 0.001,
  d_prior = 0.001
)

cat("Testing response:", test_response$label, "\n\n")

# Test dataset 1 only
comp_no <- 1
results_dir <- "oracle_consistency_analysis/response_search/test_employed"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cat("Setting up comparison 1...\n")
setup_comp(comp_no, results_dir, model_config)

# Check for zero variance
d <- readRDS(file.path(results_dir, sprintf("comparison_%03d", comp_no), "d.RDS"))
if (any(d$values == 0)) {
  cat("ERROR: Zero variance detected!\n")
  quit(status = 1)
}
cat("  No zero variance - good!\n")

# Run oracle fits
cat("\nRunning oracle fits for nbasis =", paste(nbasis_grid, collapse = ", "), "...\n")

fit_config <- model_config
fit_config$nbasis_values <- nbasis_grid
fit_config$spatial_type <- "fixed"
fit_config$ndesired <- 2000
fit_config$nburn <- 1500
fit_config$nthin <- 1
fit_config$hyp <- list(a = 1, b = 1, c = model_config$c_prior, d = model_config$d_prior)

full_data_fit(comp_no, results_dir, fit_config)

# Summarize
cat("\nSummarizing oracle results...\n")
oracle_summary <- summary_oracle(comp_no, results_dir)
print(oracle_summary)

# Find optimal
optimal_nbasis <- oracle_summary$nbasis[which.min(oracle_summary$mse)]
cat(sprintf("\nOptimal nbasis: %d (MSE = %.6f)\n", optimal_nbasis, min(oracle_summary$mse)))

cat("\nTest completed successfully!\n")
