#!/usr/bin/env Rscript
# ==============================================================================
# Aggregate Equal Allocation Results (S=50)
# ==============================================================================
# Purpose: Aggregate all equal allocation results for Sections 3, 4, 6.2
# Designs: equal_30, equal_40, equal_50, equal_75, equal_100, equal_125
# S=50 comparisons per design
# ==============================================================================

library(parallel)
library(doParallel)
library(tidyverse)

cat("\n")
cat("================================================================================\n")
cat("  AGGREGATE EQUAL ALLOCATION RESULTS (S=50)\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

all_configs <- list(
  list(name = "equal_30", dir = "_results_equal30", n_comp = 50, n_cores = 11),
  list(name = "equal_40", dir = "_results_equal40", n_comp = 50, n_cores = 11),
  list(name = "equal_50", dir = "_results_equal50", n_comp = 50, n_cores = 11),
  list(name = "equal_75", dir = "_results_equal75", n_comp = 50, n_cores = 11),
  list(name = "equal_100", dir = "_results_equal100", n_comp = 50, n_cores = 11),
  list(name = "equal_125", dir = "_results_equal125", n_comp = 50, n_cores = 11)
)

eps_values <- seq(0.1, 0.9, by = 0.1)
n_reps_values <- c(1, 3, 5)
loss_functions <- c("MSE", "plugin_NLL")

# Designs that have dt_5fold (for Section 4)
DT_5FOLD_DESIGNS <- c("equal_50", "equal_75", "equal_100")

cat("Configuration:\n")
cat(sprintf("  Designs: %d\n", length(all_configs)))
cat(sprintf("  Comparisons per design: %d\n", all_configs[[1]]$n_comp))
cat(sprintf("  Epsilon values: %s\n", paste(eps_values, collapse = ", ")))
cat(sprintf("  n_reps values: %s\n", paste(n_reps_values, collapse = ", ")))
cat(sprintf("  Loss functions: %s\n", paste(loss_functions, collapse = ", ")))
cat(sprintf("  dt_5fold designs: %s\n\n", paste(DT_5FOLD_DESIGNS, collapse = ", ")))

# ==============================================================================
# AGGREGATE RESULTS
# ==============================================================================

dt_1fold_list <- list()
dt_5fold_list <- list()
oracle_list <- list()

for (cfg_info in all_configs) {
  cat("\n")
  cat("================================================================================\n")
  cat(sprintf("  %s\n", toupper(cfg_info$name)))
  cat("================================================================================\n\n")

  if (!dir.exists(cfg_info$dir)) {
    cat(sprintf("Directory not found: %s -- skipping\n", cfg_info$dir))
    next
  }

  comp <- 1:cfg_info$n_comp

  cl <- makeForkCluster(cfg_info$n_cores)
  registerDoParallel(cl)

  # --------------------------------------------------------------------------
  # DT 1-FOLD
  # --------------------------------------------------------------------------

  cat("Aggregating DT 1-fold...\n")

  dt_configs <- expand.grid(
    eps = eps_values,
    n_reps = n_reps_values,
    loss = loss_functions,
    stringsAsFactors = FALSE
  )

  dt_list_inner <- list()
  t1 <- proc.time()

  for (i in seq_len(nrow(dt_configs))) {
    dcfg <- dt_configs[i, ]

    result <- parLapply(cl, comp, function(x, config_obj, res_dir) {
      source("sim_functions/summary_dt.R")
      tryCatch({
        summary_dt(x,
          n_folds = 1, loss_function = config_obj$loss, results_dir = res_dir,
          n_reps_to_use = config_obj$n_reps, eps = config_obj$eps
        ) %>%
          mutate(
            method_name = "dt_1fold",
            epsilon = config_obj$eps,
            n_reps_used = config_obj$n_reps
          )
      }, error = function(e) NULL)
    }, config_obj = dcfg, res_dir = cfg_info$dir)

    dt_list_inner[[i]] <- bind_rows(result)
  }

  elapsed <- proc.time() - t1
  cat(sprintf("  Done (%.1f sec)\n", elapsed[3]))

  dt_1fold <- bind_rows(dt_list_inner) %>% mutate(config = cfg_info$name)
  dt_1fold_list[[cfg_info$name]] <- dt_1fold

  # --------------------------------------------------------------------------
  # DT 5-FOLD (for Section 4 designs only)
  # --------------------------------------------------------------------------

  if (cfg_info$name %in% DT_5FOLD_DESIGNS) {
    cat("Aggregating DT 5-fold (Section 4)...\n")
    t1 <- proc.time()

    dt5_list_inner <- list()
    for (loss in loss_functions) {
      result <- parLapply(cl, comp, function(x, res_dir, loss_fn) {
        source("sim_functions/summary_dt.R")
        tryCatch({
          summary_dt(x,
            n_folds = 5, loss_function = loss_fn, results_dir = res_dir,
            n_reps_to_use = 1
          ) %>%
            mutate(
              method_name = "dt_5fold",
              loss_function = loss_fn,
              epsilon = NA_real_,
              n_reps_used = NA_integer_
            )
        }, error = function(e) NULL)
      }, res_dir = cfg_info$dir, loss_fn = loss)

      dt5_list_inner[[loss]] <- bind_rows(result)
    }

    elapsed <- proc.time() - t1
    cat(sprintf("  Done (%.1f sec)\n", elapsed[3]))

    dt_5fold <- bind_rows(dt5_list_inner) %>% mutate(config = cfg_info$name)
    dt_5fold_list[[cfg_info$name]] <- dt_5fold
  }

  # --------------------------------------------------------------------------
  # ORACLE + DIC/WAIC
  # --------------------------------------------------------------------------

  cat("Aggregating Oracle/DIC/WAIC...\n")
  t1 <- proc.time()

  oracle_raw <- parLapply(cl, comp, function(x, res_dir) {
    source("sim_functions/summary_oracle.R")
    tryCatch({
      summary_oracle(x, res_dir)
    }, error = function(e) NULL)
  }, res_dir = cfg_info$dir)

  elapsed <- proc.time() - t1
  cat(sprintf("  Done (%.1f sec)\n", elapsed[3]))

  oracle <- bind_rows(oracle_raw) %>% mutate(config = cfg_info$name)
  oracle_list[[cfg_info$name]] <- oracle

  stopCluster(cl)

  cat(sprintf("  DT 1-fold rows: %d\n", nrow(dt_1fold)))
  if (cfg_info$name %in% DT_5FOLD_DESIGNS) {
    cat(sprintf("  DT 5-fold rows: %d\n", nrow(dt_5fold_list[[cfg_info$name]])))
  }
  cat(sprintf("  Oracle rows: %d\n", nrow(oracle)))
}

# ==============================================================================
# COMBINE AND SAVE
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  COMBINING AND SAVING\n")
cat("================================================================================\n\n")

dt_1fold_combined <- bind_rows(dt_1fold_list)
dt_5fold_combined <- bind_rows(dt_5fold_list)
oracle_combined <- bind_rows(oracle_list)

cat("Combined dataset summary:\n")
cat(sprintf("  Designs: %d\n", n_distinct(dt_1fold_combined$config)))
cat(sprintf("  DT 1-fold rows: %d\n", nrow(dt_1fold_combined)))
cat(sprintf("  DT 5-fold rows: %d\n", nrow(dt_5fold_combined)))
cat(sprintf("  Oracle rows: %d\n", nrow(oracle_combined)))
cat("\n")

# Create output directory
dir.create("results_summary", showWarnings = FALSE, recursive = TRUE)

# Save
results <- list(
  dt_1fold = dt_1fold_combined,
  dt_5fold = dt_5fold_combined,
  oracle = oracle_combined,
  metadata = list(
    designs = names(dt_1fold_list),
    eps_values = eps_values,
    n_reps_values = n_reps_values,
    loss_functions = loss_functions,
    n_comp = all_configs[[1]]$n_comp,
    created = Sys.time()
  )
)

saveRDS(results, "results_summary/equal_allocation_results.RDS")

cat("Saved: results_summary/equal_allocation_results.RDS\n")

cat("\n")
cat("================================================================================\n")
cat("  COMPLETE\n")
cat("================================================================================\n\n")
