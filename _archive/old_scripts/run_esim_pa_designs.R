#!/usr/bin/env Rscript
# ==============================================================================
# Run ESIM for PA Designs
# ==============================================================================
# Purpose: Empirical simulation for methods comparison (Section 6.3)
# Designs: prop_1p25pct, prop_2pct, (prop_0.75pct when ready)
# ==============================================================================

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

cat("\n")
cat("================================================================================\n")
cat("  ESIM FOR PA DESIGNS (Section 6.3 Methods Comparison)\n")
cat("================================================================================\n\n")

# Source required functions
source("sim_functions/sampling_and_setup.R")
source("sim_functions/run_esim.R")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Shared ESIM settings
n_iter_esim <- 100  # 100 iterations per comparison
n_cores <- 10       # Leave 2 cores free for other work
n_comp <- 20        # 20 comparisons per design

# Model config for ESIM (needs spatial setup)
model_config <- list(
  adjacency_file = "data/ca_puma_adjacency.RDA",
  nbasis_values = seq(3, 60, by = 3),
  spatial_type = "fixed",
  ndesired = 2000,
  nburn = 1500,
  nthin = 1,
  hyp = list(c = 0.001, d = 0.001)
)

# PA designs to process
# Note: prop_0.75pct commented out until simulations complete
# Using 0.75%, 1.25%, 1.75% for evenly spaced 0.5% increments
pa_designs <- list(
  # list(name = "prop_0.75pct", dir = "_results_prop0.75pct_comparison"),
  list(name = "prop_1p25pct", dir = "_results_prop1p25pct_comparison"),
  list(name = "prop_1p75pct", dir = "_results_prop1p75pct_comparison")
)

cat("Configuration:\n")
cat(sprintf("  ESIM iterations: %d per comparison\n", n_iter_esim))
cat(sprintf("  Comparisons per design: %d\n", n_comp))
cat(sprintf("  Cores: %d (leaving 2 free)\n", n_cores))
cat(sprintf("  Models per iteration: %d\n", length(model_config$nbasis_values)))
cat(sprintf("  Total fits per design: %d × %d × %d = %d\n",
            n_comp, n_iter_esim, length(model_config$nbasis_values),
            n_comp * n_iter_esim * length(model_config$nbasis_values)))
cat("\n")

cat("Designs to process:\n")
for (design in pa_designs) {
  cat(sprintf("  - %s\n", design$name))
}
cat("\n")

# ==============================================================================
# PROCESS EACH DESIGN
# ==============================================================================

for (design in pa_designs) {
  cat("================================================================================\n")
  cat(sprintf("  PROCESSING: %s\n", design$name))
  cat("================================================================================\n\n")

  results_dir <- design$dir

  # Check if directory exists
  if (!dir.exists(results_dir)) {
    cat(sprintf("⚠️  Directory not found: %s (skipping)\n\n", results_dir))
    next
  }

  # --------------------------------------------------------------------------
  # STEP 1: Setup ESIM (generate synthetic data)
  # --------------------------------------------------------------------------

  cat(sprintf("=== STEP 1: Setup ESIM for %s ===\n\n", design$name))

  cat("Generating synthetic direct estimates...\n")
  t1 <- proc.time()

  setup_esim(ncomps = n_comp, results_dir = results_dir, niters = n_iter_esim)

  elapsed <- proc.time() - t1
  cat(sprintf(
    "✓ ESIM setup complete (%d iterations × %d comparisons, %.1f sec)\n\n",
    n_iter_esim, n_comp, elapsed[3]
  ))

  # --------------------------------------------------------------------------
  # STEP 2: Run ESIM (fit models on synthetic data)
  # --------------------------------------------------------------------------

  cat(sprintf("=== STEP 2: Run ESIM for %s ===\n\n", design$name))

  cat("Fitting models on synthetic datasets...\n")
  cat(sprintf("  (This will take ~6-12 hours for %d total fits)\n\n",
              n_comp * n_iter_esim * length(model_config$nbasis_values)))

  t1 <- proc.time()

  cl <- makeCluster(n_cores, type = "FORK")
  registerDoParallel(cl)

  foreach(comp_no = 1:n_comp) %dopar% {
    run_esim(comp_no,
             n_cores = 1,  # Each worker processes all iterations sequentially
             results_dir = results_dir,
             model_config = model_config,
             n_iters = n_iter_esim)
  }

  stopCluster(cl)
  elapsed <- proc.time() - t1

  cat(sprintf(
    "✓ ESIM complete for %s (%.1f sec = %.1f hours, ~%.1f sec/comp)\n\n",
    design$name, elapsed[3], elapsed[3] / 3600, elapsed[3] / n_comp
  ))
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  ESIM COMPLETE FOR ALL PA DESIGNS\n")
cat("================================================================================\n\n")

cat("Completed designs:\n")
for (design in pa_designs) {
  if (dir.exists(design$dir)) {
    cat(sprintf("  ✓ %s\n", design$name))
  }
}

cat("\n")
cat("Next steps:\n")
cat("  1. Run ESIM for prop_0.75pct once simulations complete\n")
cat("  2. Update aggregate_methods_comparison.R to include ESIM results\n")
cat("  3. Create Section 6.3 analysis notebook\n")
cat("\n")
