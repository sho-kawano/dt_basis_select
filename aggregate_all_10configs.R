#!/usr/bin/env Rscript
# Aggregate all 10 configs from scratch
# No dependency on old results_multi_config files

library(parallel)
library(doParallel)
library(tidyverse)

# All 10 config mappings
all_configs <- list(
  list(name = "equal_40", dir = "_results_equal40_comparison", n_comp = 20, n_cores = 11),
  list(name = "equal_50", dir = "_results_equal50_comparison", n_comp = 20, n_cores = 11),
  list(name = "equal_75", dir = "_results_equal75_rerun", n_comp = 20, n_cores = 11),
  list(name = "equal_100", dir = "_results_equal100_comparison", n_comp = 20, n_cores = 11),
  list(name = "prop_0.5pct", dir = "_results_prop0.5pct_comparison", n_comp = 20, n_cores = 11),
  list(name = "prop_1pct", dir = "_results_prop1pct_comparison", n_comp = 20, n_cores = 11),
  list(name = "prop_1p25pct", dir = "_results_prop1p25pct_comparison", n_comp = 20, n_cores = 11),
  list(name = "prop_1p5pct", dir = "_results_prop1p5pct_comparison", n_comp = 20, n_cores = 11),
  list(name = "prop_1p75pct", dir = "_results_prop1p75pct_comparison", n_comp = 20, n_cores = 11),
  list(name = "prop_2pct", dir = "_results_prop2pct_comparison", n_comp = 20, n_cores = 11)
)

dt_all_list <- list()
oracle_all_list <- list()

for (cfg_info in all_configs) {
  cat(sprintf("\n=== Processing %s ===\n", cfg_info$name))

  if (!dir.exists(cfg_info$dir)) {
    cat(sprintf("⚠️  Directory not found: %s\n", cfg_info$dir))
    next
  }

  cl <- makeForkCluster(cfg_info$n_cores)
  registerDoParallel(cl)

  # DT 1-fold (all epsilon values now available)
  dt_configs <- expand.grid(
    eps = c(0.3, 0.5, 0.7),
    n_reps = c(1, 3, 5),
    loss = "MSE",
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

  # Oracle + full-data (DIC/WAIC)
  fulldata_oracle <- parLapply(cl, 1:cfg_info$n_comp, function(x, res_dir) {
    source("sim_functions/summary_oracle.R")
    summary_oracle(x, res_dir)
  }, res_dir = cfg_info$dir)

  fulldata_oracle <- bind_rows(fulldata_oracle) %>% mutate(config = cfg_info$name)

  stopCluster(cl)

  dt_all_list[[cfg_info$name]] <- dt_1fold
  oracle_all_list[[cfg_info$name]] <- fulldata_oracle

  cat(sprintf("✓ Aggregated %s\n", cfg_info$name))
}

# Combine all
dt_all_combined <- bind_rows(dt_all_list)
oracle_all_combined <- bind_rows(oracle_all_list)

# Create output directory
dir.create("results_multi_config", showWarnings = FALSE)

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
