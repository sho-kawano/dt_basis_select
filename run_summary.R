#!/usr/bin/env Rscript
# ==============================================================================
# Comparison Study Summary Script (run after the study script is done)
# ==============================================================================
# Usage: Rscript run_summary.R <results_directory>
# ==============================================================================

library(parallel)
library(doParallel)
library(tidyverse)

# ==============================================================================
# SETUP
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)

results_dir <- args[1]
if (!dir.exists(results_dir)) stop(sprintf("Directory '%s' not found", results_dir))

config_file <- file.path(results_dir, "model_config.RDS")
if (!file.exists(config_file)) stop(sprintf("Config file not found: %s", config_file))

model_config <- readRDS(config_file)
n_comp <- model_config$n_comp

cat("\n================================================================================\n")
cat("  CA PIPELINE: COMPUTING SUMMARIES\n")
cat("================================================================================\n\n")
cat(sprintf(
  "Results: %s\nComparisons: %d\nCores: %d\n\n",
  results_dir, n_comp, model_config$n_cores
))

# ==============================================================================
# COMPUTE SUMMARIES
# ==============================================================================

cl <- makeForkCluster(model_config$n_cores)
registerDoParallel(cl)

cat("=== DT SUMMARIES ===\n")

# ----------- DT 1-fold: all epsilon × reps × loss combinations -----------
# Note: run_dt was run with n_reps=5, we compare performance using first 1, 3, 5
dt_configs <- expand.grid(
  eps = model_config$eps_values,
  n_reps = c(1, 3, 5), # Compare using different numbers of reps
  loss = c("MSE", "plugin_NLL", "predictive_NLL"),
  stringsAsFactors = FALSE
)

cat(sprintf("Total DT 1-fold configurations: %d\n\n", nrow(dt_configs)))
dt_1fold_list <- list()

for (i in seq_len(nrow(dt_configs))) {
  cfg <- dt_configs[i, ]
  cat(sprintf("DT 1-fold (ε=%.2f, reps=%d, %s)...\n", cfg$eps, cfg$n_reps, cfg$loss))

  result <- parLapply(cl, 1:n_comp, function(x, config_obj, res_dir) {
    source("sim_functions/summary_dt.R")
    summary_dt(x,
      n_folds = 1, loss_function = config_obj$loss, results_dir = res_dir,
      n_reps_to_use = config_obj$n_reps, eps = config_obj$eps
    ) %>%
      group_by(comp_no) %>%
      mutate(
        rank = rank(metric), method_name = "dt_1fold",
        loss_function = config_obj$loss, epsilon = config_obj$eps, n_reps_used = config_obj$n_reps,
        method = sprintf("DT 1-fold (ε=%.1f, reps=%d, %s)", config_obj$eps, config_obj$n_reps, config_obj$loss)
      )
  }, config_obj = cfg, res_dir = results_dir)

  result <- bind_rows(result)
  dt_1fold_list[[i]] <- result
}

dt_1fold <- bind_rows(dt_1fold_list)

#  ----------- DT 5-fold: all loss functions  -----------
cat("\n=== DT 5-FOLD SUMMARIES ===\n")

dt_5fold_list <- list()

for (loss in c("MSE", "plugin_NLL", "predictive_NLL")) {
  cat(sprintf("DT 5-fold (%s)...\n", loss))

  result <- parLapply(cl, 1:n_comp, function(x, res_dir, loss_fn) {
    source("sim_functions/summary_dt.R")
    summary_dt(x,
      n_folds = 5, loss_function = loss_fn, results_dir = res_dir,
      n_reps_to_use = 1
    ) %>%
      group_by(comp_no) %>%
      mutate(
        rank = rank(metric), method_name = "dt_5fold", loss_function = loss_fn,
        method = sprintf("DT 5-fold (%s)", loss_fn),
        epsilon = NA_real_, n_reps_used = NA_integer_
      )
  }, res_dir = results_dir, loss_fn = loss)

  result <- bind_rows(result)
  dt_5fold_list[[loss]] <- result
}

dt_5fold <- bind_rows(dt_5fold_list)
dt_results <- bind_rows(dt_1fold, dt_5fold)

#  ----------- ESIM: both validation methods  -----------
cat("\n=== ESIM SUMMARIES ===\n")

esim_list <- list()

for (version in c("standard", "data_fission")) {
  method <- sub("data_fission", "fission", paste0("esim_", version))
  cat(sprintf("ESIM (%s)...\n", version))

  result <- parLapply(cl, 1:n_comp, function(x, val_method, meth_name, res_dir, n_iters) {
    source("sim_functions/summary_esim.R")
    summary_esim(x, results_dir = res_dir, validation = val_method, n_iters = n_iters) %>%
      group_by(comp_no) %>%
      mutate(
        rank = rank(metric), method_name = meth_name, loss_function = NA_character_,
        method = ifelse(val_method == "standard", "ESIM (standard)", "ESIM (data fission)"),
        epsilon = NA_real_, n_reps_used = NA_integer_
      )
  }, val_method = version, meth_name = method, res_dir = results_dir, n_iters = model_config$n_iter_esim)

  result <- bind_rows(result)
  esim_list[[version]] <- result
}

esim_results <- bind_rows(esim_list)

#  ----------- Full-Data Methods & Oracle  -----------
cat("\n=== FULL-DATA & ORACLE SUMMARIES ===\n")
cat("Computing OD-Oracle MSE (oracle), DIC, WAIC (full-data)...\n")

fulldata_oracle_raw <- parLapply(cl, 1:n_comp, function(x, res_dir) {
  source("sim_functions/summary_oracle.R")
  summary_oracle(x, res_dir) %>%
    group_by(comp_no) %>%
    mutate(mse_true_rank = rank(mse_true))
}, res_dir = results_dir)

fulldata_oracle_raw <- bind_rows(fulldata_oracle_raw)

# Split into separate results for oracle, DIC, and WAIC
oracle_results <- fulldata_oracle_raw %>%
  transmute(
    comp_no, model, nbasis,
    method_name = "oracle", loss_function = "mse_true",
    metric = mse_true, rank = mse_true_rank,
    method = "OD-Oracle MSE",
    epsilon = NA_real_, n_reps_used = NA_integer_
  )

dic_results <- fulldata_oracle_raw %>%
  group_by(comp_no) %>%
  transmute(
    comp_no, model, nbasis,
    method_name = "dic", loss_function = "DIC",
    metric = DIC, rank = rank(DIC),
    method = "DIC",
    epsilon = NA_real_, n_reps_used = NA_integer_
  ) %>%
  ungroup()

waic_results <- fulldata_oracle_raw %>%
  group_by(comp_no) %>%
  transmute(
    comp_no, model, nbasis,
    method_name = "waic", loss_function = "WAIC",
    metric = WAIC, rank = rank(WAIC),
    method = "WAIC",
    epsilon = NA_real_, n_reps_used = NA_integer_
  ) %>%
  ungroup()

stopCluster(cl)

# ==============================================================================
# COMBINE AND SAVE
# ==============================================================================

cat("\n=== SAVING RESULTS ===\n")

all_results <- bind_rows(
  dt_results,
  esim_results,
  oracle_results,
  dic_results,
  waic_results
)

saveRDS(all_results, file.path(results_dir, "results.RDS"))
cat(sprintf("✓ Saved: %s\n", file.path(results_dir, "results.RDS")))

# ==============================================================================
# PERFORMANCE ANALYSIS
# ==============================================================================

cat("\n=== MODEL SELECTION PERFORMANCE ===\n\n")

# Oracle best per comparison
oracle_best <- all_results %>%
  filter(method_name == "oracle", loss_function == "mse_true") %>%
  group_by(comp_no) %>%
  slice_min(metric, n = 1, with_ties = FALSE) %>%
  select(comp_no, best_model = model)

# Performance function
calc_perf <- function(data, oracle) {
  oracle_models <- all_results %>%
    filter(method_name == "oracle", loss_function == "mse_true") %>%
    pull(model) %>%
    unique()

  data_filt <- data %>% filter(model %in% oracle_models)

  # Accuracy
  selections <- data_filt %>%
    group_by(comp_no) %>%
    slice_min(metric, n = 1, with_ties = FALSE) %>%
    left_join(oracle, by = "comp_no")

  acc <- selections %>%
    ungroup() %>%
    summarise(correct = sum(model == best_model, na.rm = TRUE)) %>%
    pull(correct)

  # Median absolute difference from oracle-selected nbasis
  oracle_nbasis <- all_results %>%
    filter(method_name == "oracle", loss_function == "mse_true") %>%
    group_by(comp_no) %>%
    slice_min(metric, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(comp_no, oracle_nbasis = nbasis)

  sel_with_oracle <- selections %>% left_join(oracle_nbasis, by = "comp_no")
  nbasis_diff_med <- sel_with_oracle %>%
    transmute(diff = abs(nbasis - oracle_nbasis)) %>%
    pull(diff) %>%
    median(na.rm = TRUE)

  # Kendall's tau
  oracle_ranks <- all_results %>%
    filter(method_name == "oracle", loss_function == "mse_true") %>%
    select(comp_no, model, oracle_rank = rank)

  tau <- data_filt %>%
    select(comp_no, model, rank) %>%
    left_join(oracle_ranks, by = c("comp_no", "model")) %>%
    group_by(comp_no) %>%
    summarise(tau = cor(rank, oracle_rank, method = "kendall"), .groups = "drop") %>%
    pull(tau) %>%
    mean(na.rm = TRUE)

  tibble(accuracy = acc, tau = tau, nbasis_diff_med = nbasis_diff_med)
}

# Compute for all methods
perf_list <- list(
  dt1 = all_results %>%
    filter(method_name == "dt_1fold") %>%
    group_by(epsilon, n_reps_used, loss_function) %>%
    group_modify(~ calc_perf(.x, oracle_best)) %>%
    ungroup() %>%
    arrange(desc(accuracy), desc(tau)),
  dt5 = all_results %>%
    filter(method_name == "dt_5fold") %>%
    group_by(loss_function) %>%
    group_modify(~ calc_perf(.x, oracle_best)) %>%
    arrange(desc(accuracy), desc(tau)),
  esim = all_results %>%
    filter(method_name %in% c("esim_standard", "esim_fission")) %>%
    group_by(method_name) %>%
    group_modify(~ calc_perf(.x, oracle_best)) %>%
    arrange(desc(accuracy), desc(tau)),
  fulldata = all_results %>%
    filter(method_name %in% c("dic", "waic")) %>%
    group_by(method_name) %>%
    group_modify(~ calc_perf(.x, oracle_best)) %>%
    arrange(desc(accuracy), desc(tau))
)

# Print summaries
cat("DT 1-FOLD:\n")
print(perf_list$dt1, n = 20)
cat("\nDT 5-FOLD:\n")
print(perf_list$dt5)
cat("\nESIM:\n")
print(perf_list$esim)
cat("\nFULL-DATA (DIC, WAIC):\n")
print(perf_list$fulldata)

# Overall ranking
cat("\n=== OVERALL RANKING ===\n\n")

overall <- bind_rows(
  perf_list$dt1 %>%
    mutate(method = sprintf(
      "DT 1-fold (ε=%.1f, reps=%d, %s)",
      epsilon, n_reps_used, loss_function
    )) %>%
    select(method, accuracy, tau, nbasis_diff_med),
  perf_list$dt5 %>%
    mutate(method = sprintf("DT 5-fold (%s)", loss_function)) %>%
    select(method, accuracy, tau, nbasis_diff_med),
  perf_list$esim %>%
    mutate(method = ifelse(method_name == "esim_standard",
      "ESIM (standard)", "ESIM (data fission)"
    )) %>%
    select(method, accuracy, tau, nbasis_diff_med),
  perf_list$fulldata %>%
    mutate(method = toupper(method_name)) %>%
    select(method, accuracy, tau, nbasis_diff_med)
) %>% arrange(desc(accuracy), desc(tau), nbasis_diff_med)

print(overall, n = 100)

cat(sprintf(
  "\nBest: %s (accuracy=%d/%d, tau=%.3f, nbasis_diff_med=%.1f)\n",
  overall$method[1], overall$accuracy[1], n_comp, overall$tau[1], overall$nbasis_diff_med[1]
))

# Save performance
saveRDS(
  c(perf_list, list(overall = overall)),
  file.path(results_dir, "results_list.RDS")
)

cat("\n================================================================================\n")
cat("  COMPLETE!\n")
cat("================================================================================\n\n")
