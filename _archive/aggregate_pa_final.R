#!/usr/bin/env Rscript
# ==============================================================================
# Aggregate PA Final Results (S=100)
# ==============================================================================
# Designs: prop_0.75pct, prop_1p25pct, prop_1p5pct
# S=100 per design (seeds: seq(4, 400, by=4))
# Includes: DT, ESIM, DIC, WAIC, Oracle
# ==============================================================================

library(parallel)
library(doParallel)
library(tidyverse)

cat("\n")
cat("================================================================================\n")
cat("  AGGREGATE PA FINAL RESULTS (S=100)\n")
cat("================================================================================\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

all_configs <- list(
  list(name = "prop_0.75pct", dir = "_results_prop0.75pct_final", n_cores = 11),
  list(name = "prop_1p25pct", dir = "_results_prop1p25pct_final", n_cores = 11),
  list(name = "prop_1p5pct",  dir = "_results_prop1p5pct_final",  n_cores = 11)
)

SKIP_ESIM <- TRUE   # set FALSE when ESIM results are ready

eps_values    <- c(0.5, 0.6, 0.7)
n_reps_values <- c(1, 3, 5)
loss_functions <- c("MSE", "plugin_NLL")
n_iter_esim   <- 100

cat("Configuration:\n")
cat(sprintf("  Designs: %d\n", length(all_configs)))
cat(sprintf("  Epsilon values: %s\n", paste(eps_values, collapse = ", ")))
cat(sprintf("  n_reps values: %s\n", paste(n_reps_values, collapse = ", ")))
cat(sprintf("  Loss functions: %s\n", paste(loss_functions, collapse = ", ")))
cat(sprintf("  ESIM iterations: %d\n\n", n_iter_esim))

# ==============================================================================
# AGGREGATE RESULTS
# ==============================================================================

dt_list     <- list()
esim_list   <- list()
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

  comp_dirs  <- list.dirs(cfg_info$dir, recursive = FALSE)
  comp_dirs  <- comp_dirs[grepl("comparison_", comp_dirs)]
  n_comp     <- length(comp_dirs)
  comp_nos   <- as.integer(sub(".*comparison_0*", "", comp_dirs))

  cat(sprintf("Found %d comparisons in %s\n", n_comp, cfg_info$dir))

  if (n_comp == 0) {
    cat("No comparisons found -- skipping\n")
    next
  }

  cl <- makeForkCluster(cfg_info$n_cores)
  registerDoParallel(cl)

  # --------------------------------------------------------------------------
  # DT 1-FOLD
  # --------------------------------------------------------------------------

  cat("Aggregating DT 1-fold...\n")

  dt_configs <- expand.grid(
    eps    = eps_values,
    n_reps = n_reps_values,
    loss   = loss_functions,
    stringsAsFactors = FALSE
  )

  dt_list_inner <- list()
  t1 <- proc.time()

  for (i in seq_len(nrow(dt_configs))) {
    dcfg <- dt_configs[i, ]

    result <- parLapply(cl, comp_nos, function(x, config_obj, res_dir) {
      source("sim_functions/summary_dt.R")
      tryCatch({
        summary_dt(x,
          n_folds = 1, loss_function = config_obj$loss, results_dir = res_dir,
          n_reps_to_use = config_obj$n_reps, eps = config_obj$eps
        ) %>%
          mutate(
            method_name  = "dt_1fold",
            epsilon      = config_obj$eps,
            n_reps_used  = config_obj$n_reps
          )
      }, error = function(e) NULL)
    }, config_obj = dcfg, res_dir = cfg_info$dir)

    dt_list_inner[[i]] <- bind_rows(result)
  }

  elapsed <- proc.time() - t1
  cat(sprintf("  Done (%.1f sec)\n", elapsed[3]))

  dt_results <- bind_rows(dt_list_inner) %>% mutate(config = cfg_info$name)
  dt_list[[cfg_info$name]] <- dt_results

  # --------------------------------------------------------------------------
  # ESIM
  # --------------------------------------------------------------------------

  if (!SKIP_ESIM) {
    cat("Aggregating ESIM...\n")
    t1 <- proc.time()

    esim_list_inner <- list()
    for (version in c("standard", "data_fission")) {
      method <- sub("data_fission", "fission", paste0("esim_", version))

      result <- parLapply(cl, comp_nos, function(x, val_method, meth_name, res_dir, n_iters) {
        source("sim_functions/summary_esim.R")
        tryCatch({
          summary_esim(x, results_dir = res_dir, validation = val_method, n_iters = n_iters) %>%
            mutate(
              method_name = meth_name,
              method      = ifelse(val_method == "standard", "ESIM (standard)", "ESIM (data fission)"),
              epsilon     = NA_real_,
              n_reps_used = NA_integer_
            )
        }, error = function(e) NULL)
      }, val_method = version, meth_name = method, res_dir = cfg_info$dir, n_iters = n_iter_esim)

      esim_list_inner[[version]] <- bind_rows(result)
    }

    elapsed <- proc.time() - t1
    cat(sprintf("  Done (%.1f sec)\n", elapsed[3]))

    esim_results <- bind_rows(esim_list_inner) %>% mutate(config = cfg_info$name)
  } else {
    cat("Skipping ESIM (SKIP_ESIM = TRUE)\n")
    esim_results <- tibble()
  }
  esim_list[[cfg_info$name]] <- esim_results

  # --------------------------------------------------------------------------
  # ORACLE + DIC/WAIC
  # --------------------------------------------------------------------------

  cat("Aggregating Oracle/DIC/WAIC...\n")
  t1 <- proc.time()

  oracle_raw <- parLapply(cl, comp_nos, function(x, res_dir) {
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

  cat(sprintf("  DT rows: %d\n", nrow(dt_results)))
  cat(sprintf("  ESIM rows: %d\n", nrow(esim_results)))
  cat(sprintf("  Oracle rows: %d\n", nrow(oracle)))
}

# ==============================================================================
# COMBINE AND SAVE
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("  COMBINING AND SAVING\n")
cat("================================================================================\n\n")

dt_combined     <- bind_rows(dt_list)
esim_combined   <- bind_rows(esim_list)
oracle_combined <- bind_rows(oracle_list)

cat("Combined dataset summary:\n")
cat(sprintf("  Designs: %d\n",   n_distinct(dt_combined$config)))
cat(sprintf("  DT rows: %d\n",   nrow(dt_combined)))
cat(sprintf("  ESIM rows: %d\n", nrow(esim_combined)))
cat(sprintf("  Oracle rows: %d\n", nrow(oracle_combined)))
cat("\n")

dir.create("results_multi_config/paper_final", showWarnings = FALSE, recursive = TRUE)

results <- list(
  dt     = dt_combined,
  esim   = esim_combined,
  oracle = oracle_combined,
  metadata = list(
    designs        = names(dt_list),
    seeds          = seq(4, 400, by = 4),
    eps_values     = eps_values,
    n_reps_values  = n_reps_values,
    loss_functions = loss_functions,
    n_iter_esim    = n_iter_esim,
    n_comp         = 100,
    created        = Sys.time()
  )
)

saveRDS(results, "results_multi_config/paper_final/pa_method_comparison_final.RDS")
cat("Saved: results_multi_config/paper_final/pa_method_comparison_final.RDS\n")

cat("\n")
cat("================================================================================\n")
cat("  COMPLETE\n")
cat("================================================================================\n\n")
