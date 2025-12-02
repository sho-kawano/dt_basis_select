#!/usr/bin/env Rscript
# Aggregate 4 additional configs (prop1p25pct, prop1p5pct, prop1p75pct, equal_100)
# Skips ESIM since none have it

library(parallel)
library(doParallel)
library(tidyverse)

# Config mappings
new_configs <- list(
  list(name = "prop_1p25pct", dir = "_results_prop1p25pct_comparison", n_comp = 20, n_cores = 8),
  list(name = "prop_1p5pct", dir = "_results_prop1p5pct_comparison", n_comp = 20, n_cores = 8),
  list(name = "prop_1p75pct", dir = "_results_prop1p75pct_comparison", n_comp = 20, n_cores = 8),
  list(name = "equal_100", dir = "_results_equal100_comparison", n_comp = 20, n_cores = 8)
)

dt_all_new <- list()
oracle_all_new <- list()

for (cfg_info in new_configs) {
  cat(sprintf("\n=== Processing %s ===\n", cfg_info$name))

  if (!dir.exists(cfg_info$dir)) {
    cat(sprintf("⚠️  Directory not found: %s\n", cfg_info$dir))
    next
  }

  cl <- makeForkCluster(cfg_info$n_cores)
  registerDoParallel(cl)

  # DT 1-fold (only eps=0.5, 0.7 available)
  dt_configs <- expand.grid(
    eps = c(0.5, 0.7),
    n_reps = c(1, 3, 5),
    loss = "MSE",  # Only MSE for simplicity
    stringsAsFactors = FALSE
  )

  dt_list <- list()
  for (i in seq_len(nrow(dt_configs))) {
    dcfg <- dt_configs[i, ]
    result <- parLapply(cl, 1:cfg_info$n_comp, function(x, config_obj, res_dir) {
      source("sim_functions/summary_dt.R")
      summary_dt(x,
        n_folds = 1, loss_function = config_obj$loss, results_dir = res_dir,
        n_reps_to_use = config_obj$n_reps, eps = config_obj$eps
      ) %>%
        mutate(
          method_name = "dt_1fold",
          epsilon = config_obj$eps, n_reps_used = config_obj$n_reps
        )
    }, config_obj = dcfg, res_dir = cfg_info$dir)

    dt_list[[i]] <- bind_rows(result)
  }

  dt_1fold <- bind_rows(dt_list) %>% mutate(config = cfg_info$name)

  # Oracle + full-data
  fulldata_oracle <- parLapply(cl, 1:cfg_info$n_comp, function(x, res_dir) {
    source("sim_functions/summary_oracle.R")
    summary_oracle(x, res_dir)
  }, res_dir = cfg_info$dir)

  fulldata_oracle <- bind_rows(fulldata_oracle) %>% mutate(config = cfg_info$name)

  stopCluster(cl)

  dt_all_new[[cfg_info$name]] <- dt_1fold
  oracle_all_new[[cfg_info$name]] <- fulldata_oracle

  cat(sprintf("✓ Aggregated %s\n", cfg_info$name))
}

# Load existing data
dt_all_old <- readRDS("results_multi_config/dt_all.RDS")
oracle_all_old <- readRDS("results_multi_config/oracle_all.RDS")

# Combine with new data
dt_all_combined <- bind_rows(dt_all_old, bind_rows(dt_all_new))
oracle_all_combined <- bind_rows(oracle_all_old, bind_rows(oracle_all_new))

# Save
saveRDS(dt_all_combined, "results_multi_config/dt_all_10configs.RDS")
saveRDS(oracle_all_combined, "results_multi_config/oracle_all_10configs.RDS")

cat("\n=== SUMMARY ===\n")
cat(sprintf("Total configs: %d\n", length(unique(dt_all_combined$config))))
cat("\nConfigs:\n")
print(sort(unique(dt_all_combined$config)))

cat("\n✓ Saved:\n")
cat("  - results_multi_config/dt_all_10configs.RDS\n")
cat("  - results_multi_config/oracle_all_10configs.RDS\n")
