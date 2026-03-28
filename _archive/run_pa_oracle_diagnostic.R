#!/usr/bin/env Rscript
# ==============================================================================
# PA Oracle Diagnostic — Fresh Seeds 51-100
# ==============================================================================
# Purpose: Before committing to a full PA re-run (with DT + ESIM), check whether
# fresh seeds (51-100) produce oracle basis spikes above p=30.
# Only runs: setup_comp -> full_data_fit -> summary_oracle
# No DT, no ESIM.
#
# Designs: prop_0.75pct, prop_1p25pct, prop_1p75pct
# Seeds: 51-100 (via set.seed(sim) in setup_comp)
# Output: results_multi_config/pa_oracle_diagnostic.RDS
# ==============================================================================

library(parallel)
library(doParallel)
library(foreach)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

cat("\n")
cat("================================================================================\n")
cat("  PA ORACLE DIAGNOSTIC S=200 (no DT/ESIM)\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

PA_DESIGNS   <- c("prop_1p5pct")
START_FROM   <- 1
NCOMPS_TOTAL <- 400
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
  equal_allocation = FALSE,
  min_sample_size  = 10
)

design_params <- list(
  `prop_0.75pct` = list(samp_frac = 0.0075, results_dir = "_results_prop0.75pct_s200"),
  `prop_1p25pct` = list(samp_frac = 0.0125, results_dir = "_results_prop1p25pct_s200"),
  `prop_1p5pct`  = list(samp_frac = 0.015,  results_dir = "_results_prop1p5pct_s400"),
  `prop_1p75pct` = list(samp_frac = 0.0175, results_dir = "_results_prop1p75pct_s200")
)

cat(sprintf("Designs: %s\n", paste(PA_DESIGNS, collapse = ", ")))
cat(sprintf("Seeds:   %d-%d (comparisons %d-%d)\n", START_FROM, NCOMPS_TOTAL, START_FROM, NCOMPS_TOTAL))
cat(sprintf("Cores:   %d\n\n", N_CORES))

# ==============================================================================
# LOAD FUNCTIONS
# ==============================================================================

source("sim_functions/sampling_and_setup.R")
source("sim_functions/full_data_fit.R")
source("sim_functions/summary_oracle.R")

# ==============================================================================
# MAIN LOOP
# ==============================================================================

oracle_list <- list()

for (design in PA_DESIGNS) {
  cat("\n")
  cat("================================================================================\n")
  cat(sprintf("  %s\n", toupper(design)))
  cat("================================================================================\n\n")

  params      <- design_params[[design]]
  results_dir <- params$results_dir
  model_config <- c(base_config, list(samp_frac = params$samp_frac))

  # --------------------------------------------------------------------------
  # PART 1: SETUP (seeds 51-100)
  # --------------------------------------------------------------------------

  cat("--- Setup comparisons", START_FROM, "-", NCOMPS_TOTAL, "---\n")
  t1 <- proc.time()
  setup_comp(
    ncomps      = NCOMPS_TOTAL,
    results_dir = results_dir,
    model_config = model_config,
    start_from  = START_FROM
  )
  cat(sprintf("Done (%.1f sec)\n\n", (proc.time() - t1)[3]))

  # --------------------------------------------------------------------------
  # PART 2: FULL DATA FITS
  # --------------------------------------------------------------------------

  cat("--- Full data fits", START_FROM, "-", NCOMPS_TOTAL, "---\n")
  t1 <- proc.time()

  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)

  foreach(comp_no = START_FROM:NCOMPS_TOTAL) %dopar% {
    full_data_fit(comp_no, results_dir, model_config)
  }

  stopCluster(cl)
  cat(sprintf("Done (%.1f sec, ~%.1f sec/comp)\n\n",
              (proc.time() - t1)[3],
              (proc.time() - t1)[3] / (NCOMPS_TOTAL - START_FROM + 1)))

  # --------------------------------------------------------------------------
  # PART 3: ORACLE SUMMARY
  # --------------------------------------------------------------------------

  cat("--- Summary oracle", START_FROM, "-", NCOMPS_TOTAL, "---\n")
  t1 <- proc.time()

  cl <- makeCluster(N_CORES, type = "FORK")
  registerDoParallel(cl)

  oracle_raw <- foreach(comp_no = START_FROM:NCOMPS_TOTAL,
                        .combine  = bind_rows) %dopar% {
    tryCatch(
      summary_oracle(comp_no, results_dir) %>% mutate(config = design),
      error = function(e) {
        cat(sprintf("  [WARN] comp %d failed: %s\n", comp_no, conditionMessage(e)))
        NULL
      }
    )
  }

  stopCluster(cl)
  cat(sprintf("Done (%.1f sec)\n\n", (proc.time() - t1)[3]))

  oracle_list[[design]] <- oracle_raw
}

# ==============================================================================
# AGGREGATE AND SAVE
# ==============================================================================

oracle_diag <- bind_rows(oracle_list)
dir.create("results_multi_config", showWarnings = FALSE)
saveRDS(oracle_diag, "results_multi_config/pa_oracle_diagnostic.RDS")
cat("Saved: results_multi_config/pa_oracle_diagnostic.RDS\n\n")

# ==============================================================================
# INLINE ANALYSIS
# ==============================================================================

cat("================================================================================\n")
cat("  ORACLE BASIS DISTRIBUTION (S=400, prop_1p5pct)\n")
cat("================================================================================\n\n")

# Per-sample oracle and flatness metrics
oracle_stats <- oracle_diag %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  mutate(min_mse = min(mse_true)) %>%
  summarise(
    oracle_nbasis    = nbasis[which.min(mse_true)],
    # Models within 2% of oracle MSE: how wide is that band in basis count?
    flat_range       = {
      within2 <- nbasis[mse_true <= 1.02 * min(mse_true)]
      max(within2) - min(within2)
    },
    n_within_2pct    = sum(mse_true <= 1.02 * min(mse_true)),
    .groups = "drop"
  ) %>%
  mutate(
    flag_high_oracle = oracle_nbasis > 30,
    flag_flat        = flat_range > 15,   # 2%+ band spans >15 basis functions
    flag_any         = flag_high_oracle | flag_flat
  )

cat("--- Overall summary ---\n")
oracle_stats %>%
  group_by(config) %>%
  summarise(
    n              = n(),
    mean_oracle    = mean(oracle_nbasis),
    sd_oracle      = sd(oracle_nbasis),
    max_oracle     = max(oracle_nbasis),
    pct_gt30       = mean(flag_high_oracle),
    mean_flat_range = mean(flat_range),
    pct_flat       = mean(flag_flat),
    pct_flagged    = mean(flag_any),
    .groups = "drop"
  ) %>% print()

cat("\n--- By 50-seed block ---\n")
oracle_stats %>%
  mutate(block = paste0(((comp_no - 1) %/% 50) * 50 + 1, "-",
                        ((comp_no - 1) %/% 50 + 1) * 50)) %>%
  group_by(config, block) %>%
  summarise(
    n           = n(),
    pct_gt30    = mean(flag_high_oracle),
    pct_flat    = mean(flag_flat),
    pct_flagged = mean(flag_any),
    max_oracle  = max(oracle_nbasis),
    max_flat_range = max(flat_range),
    .groups = "drop"
  ) %>% arrange(config, pct_flagged, block) %>%
  print(n = Inf)

cat("\n--- Flagged seeds (oracle>30 OR flat_range>15) ---\n")
oracle_stats %>%
  filter(flag_any) %>%
  select(config, comp_no, oracle_nbasis, flat_range, n_within_2pct,
         flag_high_oracle, flag_flat) %>%
  arrange(config, comp_no) %>%
  print(n = Inf)

cat("\n")
cat("================================================================================\n")
cat("  DIAGNOSTIC COMPLETE\n")
cat("================================================================================\n\n")
