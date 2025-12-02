#!/usr/bin/env Rscript
# ==============================================================================
# Fine-Grid Proportional Sampling Analysis
# ==============================================================================
# Analyze prop_1p25pct, prop_1p5pct, prop_1p75pct configs
# ==============================================================================

library(tidyverse)

cat("\n")
cat("================================================================================\n")
cat("  PROPORTIONAL FINE-GRID ANALYSIS\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

configs <- list(
  list(name = "prop_1p25pct", results_dir = "_results_prop1p25pct_comparison", n_comp = 20),
  list(name = "prop_1p5pct", results_dir = "_results_prop1p5pct_comparison", n_comp = 20),
  list(name = "prop_1p75pct", results_dir = "_results_prop1p75pct_comparison", n_comp = 20)
)

cat("Configs to analyze:\n")
for (cfg in configs) {
  cat(sprintf("  - %s (%d comparisons) -> %s\n", cfg$name, cfg$n_comp, cfg$results_dir))
}
cat("\n")

# ==============================================================================
# LOAD RESULTS
# ==============================================================================

cat("=== LOADING RESULTS ===\n\n")

source("sim_functions/summary_oracle.R")
source("sim_functions/summary_dt.R")

all_results <- list()

for (cfg in configs) {
  cat(sprintf("Processing %s...\n", cfg$name))

  # Oracle summaries
  oracle_list <- lapply(1:cfg$n_comp, function(comp_no) {
    summary_oracle(comp_no, cfg$results_dir) %>%
      mutate(config = cfg$name)
  })
  oracle_df <- bind_rows(oracle_list)

  # DT 1-fold summaries (epsilon 0.5, 0.7 with reps 1, 5)
  dt1_list <- list()
  for (eps in c(0.5, 0.7)) {
    for (n_reps in c(1, 5)) {
      dt1_eps <- lapply(1:cfg$n_comp, function(comp_no) {
        summary_dt(comp_no, n_folds = 1, loss_function = "MSE",
                   results_dir = cfg$results_dir,
                   n_reps_to_use = n_reps, eps = eps) %>%
          mutate(config = cfg$name, method_name = "dt_1fold",
                 epsilon = eps, n_reps_used = n_reps)
      })
      dt1_list <- c(dt1_list, dt1_eps)
    }
  }
  dt1_df <- bind_rows(dt1_list)

  # DT 5-fold summaries
  dt5_list <- lapply(1:cfg$n_comp, function(comp_no) {
    summary_dt(comp_no, n_folds = 5, loss_function = "MSE",
               results_dir = cfg$results_dir, n_reps_to_use = 1) %>%
      mutate(config = cfg$name, method_name = "dt_5fold")
  })
  dt5_df <- bind_rows(dt5_list)

  all_results[[cfg$name]] <- list(
    oracle = oracle_df,
    dt1 = dt1_df,
    dt5 = dt5_df
  )
}

cat("✓ All results loaded\n\n")

# Combine all results
oracle_all <- bind_rows(lapply(all_results, function(x) x$oracle))
dt_all <- bind_rows(
  bind_rows(lapply(all_results, function(x) x$dt1)),
  bind_rows(lapply(all_results, function(x) x$dt5))
)

# ==============================================================================
# ORACLE ANALYSIS
# ==============================================================================

cat("=== ORACLE ANALYSIS ===\n\n")

compute_oracle_properties <- function(oracle_df, config_name) {
  oracle_models <- oracle_df %>%
    filter(model != "Direct" & model != "D.Est") %>%
    mutate(nbasis = as.numeric(gsub("nbasis_", "", model)))

  oracle_optimal <- oracle_models %>%
    group_by(comp_no) %>%
    summarise(
      optimal_nbasis = nbasis[which.min(mse_true)],
      min_mse = min(mse_true),
      max_mse = max(mse_true),
      .groups = "drop"
    )

  # Steepness metrics
  steepness_data <- oracle_models %>%
    left_join(oracle_optimal %>% select(comp_no, optimal_nbasis, min_mse), by = "comp_no") %>%
    mutate(
      dist_from_optimal = nbasis - optimal_nbasis,
      pct_from_min = (mse_true / min_mse - 1) * 100
    )

  # Left vs right asymmetry
  penalty_left_6 <- steepness_data %>%
    filter(dist_from_optimal == -6) %>%
    summarise(avg = mean(pct_from_min, na.rm = TRUE)) %>%
    pull(avg)

  penalty_right_6 <- steepness_data %>%
    filter(dist_from_optimal == 6) %>%
    summarise(avg = mean(pct_from_min, na.rm = TRUE)) %>%
    pull(avg)

  data.frame(
    config = config_name,
    mean_optimal = round(mean(oracle_optimal$optimal_nbasis), 1),
    sd_optimal = round(sd(oracle_optimal$optimal_nbasis), 1),
    range_min = min(oracle_optimal$optimal_nbasis),
    range_max = max(oracle_optimal$optimal_nbasis),
    boundary_pct = round(100 * mean(oracle_optimal$optimal_nbasis %in% c(3, 60)), 1),
    mse_variation = round(100 * mean((oracle_optimal$max_mse - oracle_optimal$min_mse) / oracle_optimal$min_mse), 1),
    penalty_left_6 = round(penalty_left_6, 1),
    penalty_right_6 = round(penalty_right_6, 1),
    asymmetry_ratio = round(penalty_left_6 / penalty_right_6, 2)
  )
}

oracle_properties <- bind_rows(lapply(configs, function(cfg) {
  oracle_df <- all_results[[cfg$name]]$oracle
  compute_oracle_properties(oracle_df, cfg$name)
}))

cat("Oracle Properties:\n")
print(oracle_properties)
cat("\n")

# ==============================================================================
# METHOD PERFORMANCE ANALYSIS
# ==============================================================================

cat("=== METHOD PERFORMANCE ===\n\n")

# Get oracle selections per comparison
oracle_selections <- oracle_all %>%
  filter(model != "Direct" & model != "D.Est") %>%
  mutate(nbasis = as.numeric(gsub("nbasis_", "", model))) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, oracle_nbasis = nbasis, oracle_mse = mse_true)

# Analyze DT method performance
dt_performance <- dt_all %>%
  group_by(config, comp_no, method_name, epsilon, n_reps_used) %>%
  slice_min(metric, n = 1, with_ties = FALSE) %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(deviation = nbasis - oracle_nbasis)

# Summarize by method and config
method_summary <- dt_performance %>%
  group_by(config, method_name, epsilon, n_reps_used) %>%
  summarise(
    mean_nbasis = round(mean(nbasis), 1),
    mean_dev = round(mean(deviation), 1),
    mad = round(mean(abs(deviation)), 1),
    exact_match_pct = round(100 * mean(deviation == 0), 1),
    within_3_pct = round(100 * mean(abs(deviation) <= 3), 1),
    .groups = "drop"
  )

cat("Method Performance Summary:\n")
print(method_summary)
cat("\n")

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

cat("=== SAVING RESULTS ===\n\n")

output_dir <- "results_prop_fine_grid"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(oracle_all, file.path(output_dir, "oracle_all.RDS"))
saveRDS(dt_all, file.path(output_dir, "dt_all.RDS"))
saveRDS(oracle_properties, file.path(output_dir, "oracle_properties.RDS"))
saveRDS(method_summary, file.path(output_dir, "method_summary.RDS"))

write.csv(oracle_properties, file.path(output_dir, "oracle_properties.csv"), row.names = FALSE)
write.csv(method_summary, file.path(output_dir, "method_summary.csv"), row.names = FALSE)

cat(sprintf("✓ Results saved to %s/\n\n", output_dir))

cat("================================================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")
