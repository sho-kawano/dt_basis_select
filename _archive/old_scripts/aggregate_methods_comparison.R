#!/usr/bin/env Rscript
# Aggregate results for Section 7 methods comparison
# Optimized aggregation: PA designs, both DT loss functions, pre-computed metrics

library(parallel)
library(doParallel)
library(tidyverse)

# Source helper functions
source("sim_functions/summary_dt.R")
source("sim_functions/summary_oracle.R")
source("sim_functions/compute_metrics.R")

#-----------------------------------------------------------------------------
# Configuration
#-----------------------------------------------------------------------------

# PA designs for Section 7
# Using 0.75%, 1.25%, 1.75% for evenly spaced 0.5% increments
pa_configs <- list(
  list(name = "prop_0.75pct", dir = "_results_prop0.75pct_comparison", n_comp = 20),
  list(name = "prop_1p25pct", dir = "_results_prop1p25pct_comparison", n_comp = 20),
  list(name = "prop_1p75pct", dir = "_results_prop1p75pct_comparison", n_comp = 20)
)

# DT parameters
epsilon_values <- c(0.5, 0.6, 0.7)  # Section 6 recommended range
n_reps <- 5                          # Repeated thinning
loss_functions <- c("MSE", "plugin_NLL")  # Both loss functions

# Parallel processing
n_cores <- 11

# Output directory
output_dir <- "results_multi_config/methods_comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#-----------------------------------------------------------------------------
# Aggregate DT results
#-----------------------------------------------------------------------------

cat("\n=== AGGREGATING DATA THINNING RESULTS ===\n")

dt_all_list <- list()

for (cfg_info in pa_configs) {
  cat(sprintf("\nProcessing %s...\n", cfg_info$name))

  if (!dir.exists(cfg_info$dir)) {
    cat(sprintf("⚠️  Directory not found: %s (skipping)\n", cfg_info$dir))
    next
  }

  cl <- makeForkCluster(n_cores)
  registerDoParallel(cl)

  # Generate all DT configurations
  dt_configs <- expand.grid(
    eps = epsilon_values,
    loss = loss_functions,
    stringsAsFactors = FALSE
  )

  dt_list <- list()
  for (i in seq_len(nrow(dt_configs))) {
    dcfg <- dt_configs[i, ]
    cat(sprintf("  DT: ε=%.1f, loss=%s, R=%d\n", dcfg$eps, dcfg$loss, n_reps))

    result <- parLapply(cl, 1:cfg_info$n_comp, function(x, config_obj, res_dir, n_r) {
      source("sim_functions/summary_dt.R")
      summary_dt(
        comp_no = x,
        n_folds = 1,
        loss_function = config_obj$loss,
        results_dir = res_dir,
        n_reps_to_use = n_r,
        eps = config_obj$eps
      ) %>%
        mutate(
          method_name = "dt_1fold",
          epsilon = config_obj$eps,
          n_reps_used = n_r
        )
    }, config_obj = dcfg, res_dir = cfg_info$dir, n_r = n_reps)

    dt_list[[i]] <- bind_rows(result)
  }

  stopCluster(cl)

  dt_all_list[[cfg_info$name]] <- bind_rows(dt_list) %>%
    mutate(config = cfg_info$name)

  cat(sprintf("✓ Aggregated DT for %s\n", cfg_info$name))
}

dt_all_combined <- bind_rows(dt_all_list)

#-----------------------------------------------------------------------------
# Aggregate Oracle/DIC/WAIC results
#-----------------------------------------------------------------------------

cat("\n=== AGGREGATING ORACLE/DIC/WAIC RESULTS ===\n")

oracle_all_list <- list()

for (cfg_info in pa_configs) {
  cat(sprintf("\nProcessing %s...\n", cfg_info$name))

  if (!dir.exists(cfg_info$dir)) {
    next
  }

  cl <- makeForkCluster(n_cores)
  registerDoParallel(cl)

  fulldata_oracle <- parLapply(cl, 1:cfg_info$n_comp, function(x, res_dir) {
    source("sim_functions/summary_oracle.R")
    summary_oracle(x, res_dir)
  }, res_dir = cfg_info$dir)

  stopCluster(cl)

  oracle_all_list[[cfg_info$name]] <- bind_rows(fulldata_oracle) %>%
    mutate(config = cfg_info$name)

  cat(sprintf("✓ Aggregated oracle for %s\n", cfg_info$name))
}

oracle_all_combined <- bind_rows(oracle_all_list)

#-----------------------------------------------------------------------------
# Aggregate ESIM results
#-----------------------------------------------------------------------------

cat("\n=== AGGREGATING ESIM RESULTS ===\n")

source("sim_functions/summary_esim.R")

esim_all_list <- list()

for (cfg_info in pa_configs) {
  cat(sprintf("\nProcessing ESIM for %s...\n", cfg_info$name))

  if (!dir.exists(cfg_info$dir)) {
    next
  }

  cl <- makeForkCluster(n_cores)
  registerDoParallel(cl)

  esim_results <- parLapply(cl, 1:cfg_info$n_comp, function(x, res_dir) {
    source("sim_functions/summary_esim.R")
    summary_esim(x, res_dir, validation = "standard", n_iters = 100)
  }, res_dir = cfg_info$dir)

  stopCluster(cl)

  esim_all_list[[cfg_info$name]] <- bind_rows(esim_results) %>%
    mutate(config = cfg_info$name)

  cat(sprintf("✓ Aggregated ESIM for %s\n", cfg_info$name))
}

esim_all_combined <- bind_rows(esim_all_list)

#-----------------------------------------------------------------------------
# Compute evaluation metrics
#-----------------------------------------------------------------------------

cat("\n=== COMPUTING EVALUATION METRICS ===\n")

# Get DT selections for both loss functions
dt_mse_selections <- get_dt_selections(
  dt_all_combined,
  loss_function = "MSE",
  epsilon_values = epsilon_values,
  n_reps_values = n_reps
)

dt_nll_selections <- get_dt_selections(
  dt_all_combined,
  loss_function = "plugin_NLL",
  epsilon_values = epsilon_values,
  n_reps_values = n_reps
)

# Get DIC/WAIC selections
dic_selections <- get_ic_selections(oracle_all_combined, criterion = "DIC")
waic_selections <- get_ic_selections(oracle_all_combined, criterion = "WAIC")

# Get ESIM selections (model with lowest average MSE)
esim_selections <- esim_all_combined %>%
  group_by(config, comp_no) %>%
  slice_min(metric, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(config, comp_no, selected_model = model, selected_nbasis = nbasis) %>%
  mutate(method = "ESIM")

# Combine all selections
all_selections <- bind_rows(
  dt_mse_selections,
  dt_nll_selections,
  dic_selections,
  waic_selections,
  esim_selections
)

# Compute metrics
cat("Computing MAD, directional bias, and oracle penalty...\n")

metrics_summary <- compute_metrics(
  oracle_data = oracle_all_combined,
  method_data = all_selections,
  method_col = "method",
  group_cols = c("epsilon", "n_reps_used")
)

# Add cleaner method labels for plotting
metrics_summary <- metrics_summary %>%
  mutate(
    method_type = case_when(
      grepl("DT-MSE", method) ~ "DT-MSE",
      grepl("DT-plugin_NLL", method) ~ "DT-Likelihood",
      method == "DIC" ~ "DIC",
      method == "WAIC" ~ "WAIC",
      method == "ESIM" ~ "ESIM",
      TRUE ~ "Other"
    )
  )

cat("✓ Metrics computed\n")

#-----------------------------------------------------------------------------
# Save results
#-----------------------------------------------------------------------------

cat("\n=== SAVING RESULTS ===\n")

# Save raw results
saveRDS(dt_all_combined, file.path(output_dir, "dt_results.RDS"))
saveRDS(oracle_all_combined, file.path(output_dir, "oracle_results.RDS"))
saveRDS(esim_all_combined, file.path(output_dir, "esim_results.RDS"))
saveRDS(metrics_summary, file.path(output_dir, "metrics_summary.RDS"))

# Also save selections for reference
saveRDS(all_selections, file.path(output_dir, "method_selections.RDS"))

cat("✓ Saved files:\n")
cat(sprintf("  - %s/dt_results.RDS (%d rows)\n",
            output_dir, nrow(dt_all_combined)))
cat(sprintf("  - %s/oracle_results.RDS (%d rows)\n",
            output_dir, nrow(oracle_all_combined)))
cat(sprintf("  - %s/esim_results.RDS (%d rows)\n",
            output_dir, nrow(esim_all_combined)))
cat(sprintf("  - %s/metrics_summary.RDS (%d rows)\n",
            output_dir, nrow(metrics_summary)))
cat(sprintf("  - %s/method_selections.RDS (%d rows)\n",
            output_dir, nrow(all_selections)))

#-----------------------------------------------------------------------------
# Summary statistics
#-----------------------------------------------------------------------------

cat("\n=== SUMMARY ===\n")

cat(sprintf("Configs aggregated: %d\n", length(unique(dt_all_combined$config))))
cat("  ", paste(sort(unique(dt_all_combined$config)), collapse = ", "), "\n")

cat(sprintf("\nDT configurations:\n"))
cat(sprintf("  Epsilon values: %s\n", paste(epsilon_values, collapse = ", ")))
cat(sprintf("  Repeats (R): %d\n", n_reps))
cat(sprintf("  Loss functions: %s\n", paste(loss_functions, collapse = ", ")))

cat(sprintf("\nMethods in metrics summary:\n"))
print(metrics_summary %>%
  group_by(method_type) %>%
  summarise(n_configs = n(), .groups = "drop"))

cat("\n✓ Aggregation complete!\n")
