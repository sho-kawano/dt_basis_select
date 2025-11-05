#!/usr/bin/env Rscript
# Run only the summary calculation phase (Part 4)
# Use this to complete a run that failed during summary calculation

# ==============================================================================
# CONFIGURATION: Set to save results with a specific name
# ==============================================================================
SAVE_AS <- NULL  # Set to NULL to skip, or "run_name" to save to saved_results/run_name/

library(parallel)
library(doParallel)
library(tidyverse)

cat("\n================================================================================\n")
cat("  COMPUTING SUMMARIES ONLY\n")
cat("================================================================================\n\n")

# Load configuration
config <- readRDS("sim_config.RDS")
sim_config <- config

# Load data
load("data/data_IL.RDA")
chosen_var <- config$chosen_var
n_comp <- config$n_comp
n_cores <- config$n_cores

# Define epsilon and repeat configurations (must match what was run)
# NOTE: Only include epsilon values that have been computed
epsilon_values <- c(0.01, 0.02, 0.03, 0.04, seq(from = 0.05, to = 0.8, by = 0.05))
repeat_counts <- c(1, 3, 5)
loss_functions <- c("MSE", "plugin_NLL", "predictive_NLL")

# Load sources
source("sim_functions/summary_esim.R")
source("sim_functions/summary_dt.R")
source("sim_functions/summary_zfit.R")

cat("Calculating metrics for all methods in parallel...\n\n")

# Set up cluster for summary calculations
cl <- makeForkCluster(n_cores)
registerDoParallel(cl)

#-------------------------------------------------------------------------------
# Compute summaries for all DT 1-fold configurations (epsilon × repeats × loss)
cat("Computing DT 1-fold summaries for all configurations...\n\n")

# Show configured epsilons and total configurations
cat(sprintf("Epsilon values: %s\n", paste(epsilon_values, collapse = ", ")))
cat(sprintf(
  "Total DT 1-fold configurations: %d (eps × reps = %d × %d)\n\n",
  length(epsilon_values) * length(repeat_counts) * length(loss_functions),
  length(epsilon_values), length(repeat_counts)
))

dt_results_list <- list()

for (eps in epsilon_values) {
  cat(sprintf("\n=== EPSILON = %.2f ===\n", eps))

  for (n_reps in repeat_counts) {
    for (loss in loss_functions) {
      cat(sprintf("Computing %d reps, %s loss...\n", n_reps, loss))
      timer <- proc.time()

      result <- parLapply(cl, 1:n_comp,
        function(x, cfg) {
          source("sim_functions/summary_dt.R")
          summary_dt(x,
            n_folds = 1, cfg$chosen_var, loss_function = cfg$loss,
            all_data = cfg$all_data, results_dir = cfg$results_dir,
            n_reps_to_use = cfg$n_reps, eps = cfg$eps
          ) %>%
            group_by(comp_no) %>%
            mutate(
              rank = rank(metric),
              method_name = "dt_1fold",
              loss_function = cfg$loss,
              epsilon = cfg$eps,
              n_reps_used = cfg$n_reps
            )
        },
        cfg = list(
          eps = eps, n_reps = n_reps, loss = loss,
          chosen_var = chosen_var, all_data = all_data,
          results_dir = sim_config$results_dir
        )
      ) %>% bind_rows()

      key <- sprintf("dt1_eps%.2f_rep%d_%s", eps, n_reps, loss)
      dt_results_list[[key]] <- result

      elapsed <- (proc.time() - timer)[3]
      cat(sprintf("  ✓ Completed in %.1f seconds\n", elapsed))
    }
  }
}

#-------------------------------------------------------------------------------
# Compute summaries for DT 5-fold (3 loss functions)
cat("\n\n=== DT 5-FOLD SUMMARIES ===\n")

for (loss in loss_functions) {
  cat(sprintf("Computing %s loss...\n", loss))
  timer <- proc.time()

  result <- parLapply(cl, 1:n_comp,
    function(x, cfg) {
      source("sim_functions/summary_dt.R")
      summary_dt(x,
        n_folds = 5, cfg$chosen_var, loss_function = cfg$loss,
        all_data = cfg$all_data, results_dir = cfg$results_dir,
        n_reps_to_use = 1
      ) %>%
        group_by(comp_no) %>%
        mutate(
          rank = rank(metric),
          method_name = "dt_5fold",
          loss_function = cfg$loss
        )
    },
    cfg = list(
      loss = loss, chosen_var = chosen_var,
      all_data = all_data, results_dir = sim_config$results_dir
    )
  ) %>%
    bind_rows()

  key <- sprintf("dt5_%s", loss)
  dt_results_list[[key]] <- result

  elapsed <- (proc.time() - timer)[3]
  cat(sprintf("  ✓ Completed in %.1f seconds\n\n", elapsed))
}

# Combine all DT results
dt_results <- bind_rows(dt_results_list)

#-------------------------------------------------------------------------------
cat("Computing ESIM metrics (standard validation)...\n")
timer <- proc.time()

esim_standard <- parLapply(
  cl, 1:n_comp,
  function(x) {
    summary_esim(x,
      all_data = all_data,
      results_dir = sim_config$results_dir,
      validation = "standard"
    ) %>%
      group_by(comp_no) %>%
      mutate(
        rank = rank(metric),
        method_name = "esim_standard",
        loss_function = NA
      )
  }
) %>% bind_rows()

elapsed <- (proc.time() - timer)[3]
cat(sprintf("  ✓ Completed in %.1f seconds\n\n", elapsed))

#-------------------------------------------------------------------------------
cat("Computing ESIM metrics (data fission validation)...\n")
timer <- proc.time()

esim_fission <- parLapply(
  cl, 1:n_comp,
  function(x) {
    summary_esim(x,
      all_data = all_data,
      results_dir = sim_config$results_dir,
      validation = "data_fission"
    ) %>%
      group_by(comp_no) %>%
      mutate(
        rank = rank(metric),
        method_name = "esim_fission",
        loss_function = NA
      )
  }
) %>% bind_rows()

elapsed <- (proc.time() - timer)[3]
cat(sprintf("  ✓ Completed in %.1f seconds\n\n", elapsed))

# Combine both ESIM approaches
esim_results <- bind_rows(esim_standard, esim_fission)

#-------------------------------------------------------------------------------
# Process ZFIT results
cat("Processing ZFIT results...\n")

truth <- config$truth
all_covs <- config$all_covs

zfit_res <- parLapply(cl, 1:n_comp,
  function(x, truth_val, var_val, covs_val, dir_val) {
    source("sim_functions/summary_zfit.R")
    summary_zfit(x, truth_val, var_val, covs_val, dir_val) %>%
      group_by(comp_no) %>%
      mutate(mse_true_rank = rank(mse_true))
  },
  truth_val = truth, var_val = chosen_var,
  covs_val = all_covs, dir_val = sim_config$results_dir
) %>% bind_rows()

zfit_results <- zfit_res %>%
  mutate(
    method_name = "zfit",
    loss_function = NA
  )

stopCluster(cl)

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

cat("\n================================================================================\n")
cat("=== SAVING RESULTS ===\n")
cat("================================================================================\n\n")

# Combine all results into single unified structure
cat("Creating unified results structure...\n")

all_results <- bind_rows(
  dt_results %>% select(comp_no, method_name, loss_function, model, metric, rank, method, metric_type, epsilon, n_reps_used),
  esim_results %>% select(comp_no, method_name, loss_function, model, metric, rank, method, metric_type),
  zfit_results %>% select(comp_no, method_name, loss_function, model, mse_true, DIC, WAIC, mse_true_rank) %>%
    pivot_longer(cols = c(mse_true, DIC, WAIC), names_to = "metric_type", values_to = "metric") %>%
    mutate(
      method = case_when(
        metric_type == "mse_true" ~ "Oracle MSE (per-comp)",
        metric_type == "DIC" ~ "DIC",
        metric_type == "WAIC" ~ "WAIC",
        TRUE ~ metric_type
      )
    ) %>%
    group_by(comp_no, metric_type) %>%
    mutate(rank = rank(metric)) %>%
    select(comp_no, method_name, loss_function, model, metric, rank, method, metric_type)
)

cat("Saving results to results.RDS...\n")
saveRDS(all_results, "results.RDS")

cat("\n✓ Results saved successfully!\n\n")

cat("Summary of unified results:\n")
cat(sprintf("  Total rows: %d\n", nrow(all_results)))
cat(sprintf("  Comparisons: %d\n", n_distinct(all_results$comp_no)))
cat(sprintf("  Methods: %s\n", paste(unique(all_results$method_name), collapse = ", ")))
cat(sprintf("  Models tested: %s\n\n", paste(unique(all_results$model), collapse = ", ")))

# Show breakdown by method
cat("Results by method:\n")
all_results %>%
  group_by(method_name) %>%
  summarise(n_rows = n(), .groups = "drop") %>%
  arrange(desc(n_rows)) %>%
  print()

cat("\n================================================================================\n")
cat("  SUMMARY COMPLETE!\n")
cat("================================================================================\n\n")

# ===============================================================================
# PREVIEW: Model Selection Performance
# ===============================================================================

cat("\n================================================================================\n")
cat("=== MODEL SELECTION PERFORMANCE PREVIEW ===\n")
cat("================================================================================\n\n")

# Get oracle best model for each comparison (lowest oracle MSE)
oracle_best <- all_results %>%
  filter(method_name == "zfit", metric_type == "mse_true") %>%
  group_by(comp_no) %>%
  slice_min(metric, n = 1, with_ties = FALSE) %>%
  select(comp_no, best_model = model)

# Function to compute selection accuracy and Kendall's Tau
compute_performance <- function(data, oracle) {
  # Get all models that exist in oracle (1_cov, 4_cov, 7_cov)
  # This excludes ESIM's "Direct" model
  oracle_ranks_temp <- all_results %>%
    filter(method_name == "zfit", metric_type == "mse_true") %>%
    select(comp_no, model)

  oracle_all_models <- unique(oracle_ranks_temp$model)

  # Filter to models that exist in oracle
  data_filtered <- data %>%
    filter(model %in% oracle_all_models)

  # ACCURACY: Count how many times method selects same model as per-comparison oracle
  selections <- data_filtered %>%
    group_by(comp_no) %>%
    slice_min(metric, n = 1, with_ties = FALSE) %>%
    ungroup()

  # Join with oracle to compare selections
  accuracy <- selections %>%
    left_join(oracle, by = "comp_no") %>%
    summarise(correct = sum(model == best_model, na.rm = TRUE)) %>%
    pull(correct)

  # Compute Kendall's Tau
  # Compare to per-comparison oracle ranks
  oracle_ranks <- all_results %>%
    filter(method_name == "zfit", metric_type == "mse_true") %>%
    select(comp_no, model, oracle_rank = rank) %>%
    arrange(comp_no, model)

  method_ranks <- data %>%
    filter(model %in% oracle_all_models) %>%
    select(comp_no, model, rank) %>%
    arrange(comp_no, model)

  combined <- method_ranks %>%
    left_join(oracle_ranks, by = c("comp_no", "model"))

  # Compute tau for each comparison separately, then average
  per_comp_tau <- combined %>%
    group_by(comp_no) %>%
    summarise(tau = cor(rank, oracle_rank, method = "kendall"), .groups = "drop")

  tau <- mean(per_comp_tau$tau, na.rm = TRUE)

  return(list(accuracy = accuracy, tau = tau))
}

# Compute for all DT 1-fold configurations
dt1_count <- length(epsilon_values) * length(repeat_counts) * length(loss_functions)
cat(sprintf("DT 1-FOLD RESULTS (%d configurations):\n", dt1_count))
cat(sprintf(
  "  %d epsilon values × %d repeat counts × %d loss functions\n\n",
  length(epsilon_values), length(repeat_counts), length(loss_functions)
))

dt1_performance <- all_results %>%
  filter(method_name == "dt_1fold") %>%
  group_by(epsilon, n_reps_used, loss_function) %>%
  group_modify(~ {
    perf <- compute_performance(.x, oracle_best)
    tibble(accuracy = perf$accuracy, tau = perf$tau)
  }) %>%
  ungroup() %>%
  arrange(desc(accuracy), desc(tau))

print(dt1_performance, n = 50)

# Compute for DT 5-fold
cat("\n\nDT 5-FOLD RESULTS:\n\n")

dt5_performance <- all_results %>%
  filter(method_name == "dt_5fold") %>%
  group_by(loss_function) %>%
  group_modify(~ {
    perf <- compute_performance(.x, oracle_best)
    tibble(accuracy = perf$accuracy, tau = perf$tau)
  }) %>%
  ungroup() %>%
  arrange(desc(accuracy), desc(tau))

print(dt5_performance)

# Compute for ESIM
cat("\n\nESIM RESULTS:\n\n")

esim_performance <- all_results %>%
  filter(method_name %in% c("esim_standard", "esim_fission")) %>%
  group_by(method_name) %>%
  group_modify(~ {
    perf <- compute_performance(.x, oracle_best)
    tibble(accuracy = perf$accuracy, tau = perf$tau)
  }) %>%
  ungroup() %>%
  arrange(desc(accuracy), desc(tau))

print(esim_performance)

# Compute for full-data methods (DIC, WAIC)
cat("\n\nFULL-DATA METHODS:\n\n")

zfit_performance <- all_results %>%
  filter(method_name == "zfit", metric_type %in% c("DIC", "WAIC")) %>%
  group_by(metric_type) %>%
  group_modify(~ {
    perf <- compute_performance(.x, oracle_best)
    tibble(accuracy = perf$accuracy, tau = perf$tau)
  }) %>%
  ungroup() %>%
  mutate(
    method_label = case_when(
      metric_type == "DIC" ~ "DIC",
      metric_type == "WAIC" ~ "WAIC",
      TRUE ~ metric_type
    )
  ) %>%
  select(method = method_label, accuracy, tau) %>%
  arrange(desc(accuracy), desc(tau))

print(zfit_performance)

# Combined summary sorted by performance
cat("\n\n================================================================================\n")
# Compute dynamic total methods
dt1_methods <- length(epsilon_values) * length(repeat_counts) * length(loss_functions)
dt5_methods <- length(loss_functions)
esim_methods <- 2
zfit_methods <- 2  # DIC, WAIC
total_methods <- dt1_methods + dt5_methods + esim_methods + zfit_methods
cat(sprintf("=== OVERALL RANKING (all %d methods, sorted by accuracy) ===\n", total_methods))
cat("================================================================================\n\n")

overall_summary <- bind_rows(
  dt1_performance %>%
    mutate(method = sprintf(
      "DT 1-fold (ε=%.2f, reps=%d, %s)",
      epsilon, n_reps_used, loss_function
    )) %>%
    select(method, accuracy, tau),
  dt5_performance %>%
    mutate(method = sprintf("DT 5-fold (%s)", loss_function)) %>%
    select(method, accuracy, tau),
  esim_performance %>%
    mutate(method = ifelse(method_name == "esim_standard",
      "ESIM (standard)", "ESIM (data fission)"
    )) %>%
    select(method, accuracy, tau),
  zfit_performance  # Already has method, accuracy, tau columns with proper names
) %>%
  arrange(desc(accuracy), desc(tau))

print(overall_summary, n = 100)

cat("\n")
cat(sprintf("Best performing method: %s\n", overall_summary$method[1]))
cat(sprintf(
  "  Accuracy: %d/%d (%.1f%%)\n",
  overall_summary$accuracy[1], n_comp,
  100 * overall_summary$accuracy[1] / n_comp
))
cat(sprintf("  Kendall's Tau: %.3f\n", overall_summary$tau[1]))

cat("\n================================================================================\n")
cat("  PREVIEW COMPLETE!\n")
cat("================================================================================\n\n")

results_list = list(dt1_performance, dt5_performance, esim_performance, zfit_performance)
saveRDS(results_list, file="results_list.RDS")

# ==============================================================================
# SAVE ARCHIVED COPY (if SAVE_AS is set)
# ==============================================================================
if (!is.null(SAVE_AS)) {
  save_dir <- file.path("saved_results", SAVE_AS)
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("\nSaving archived copy to: %s\n", save_dir))
  file.copy("results.RDS", file.path(save_dir, "results.RDS"), overwrite = TRUE)
  file.copy("results_list.RDS", file.path(save_dir, "results_list.RDS"), overwrite = TRUE)
  file.copy("sim_config.RDS", file.path(save_dir, "sim_config.RDS"), overwrite = TRUE)

  # Save timestamp
  writeLines(sprintf("Run: %s\nDate: %s\nComparisons: %d",
                     SAVE_AS, Sys.time(), n_comp),
             file.path(save_dir, "info.txt"))

  cat(sprintf("✓ Saved to %s\n", save_dir))
}
