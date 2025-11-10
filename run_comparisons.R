#!/usr/bin/env Rscript
# ==============================================================================
# CONSOLIDATED COMPARISON PIPELINE
# ==============================================================================
# This script runs a complete comparison study of model selection methods:
# - ZFIT: Oracle fits using true means
# - DT 1-fold: Data thinning with single fold (configurable reps)
# - DT 5-fold: Data thinning with 5 folds
# - ESIM: Empirical simulation (100 iterations)
#
# Results are saved to results.RDS for analysis
# ==============================================================================

# ==============================================================================
# CONFIGURATION: Set to save results with a specific name
# ==============================================================================
SAVE_AS <- "HICOV" # Set to NULL to skip, or "run_name" to save to saved_results/run_name/

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

cat("\n")
cat("================================================================================\n")
cat("  MODEL SELECTION COMPARISON STUDY\n")
cat("================================================================================\n\n")

# ==============================================================================
# PART 1: DATA SETUP & CONFIGURATION
# ==============================================================================

cat("=== PART 1: DATA SETUP & CONFIGURATION ===\n\n")

# ---- Load PUMS data ----
cat("Loading CA PUMS data...\n")

# Load individual-level population data for sampling
acs_pop <- readRDS("data/ca_pums_population.rds")
cat(sprintf("✓ Population data loaded: %d individuals\n", nrow(acs_pop)))

# Define predictors (13 variables - matches create_X in sampling_and_setup.R)
predictor_names <- c(
  # Demographics (8)
  "mean_age", "pct_male", "pct_white", "pct_black",
  "pct_asian", "pct_hispanic", "pct_married", "pct_citizen",
  # Economic (1)
  "employment_rate",
  # Education (1)
  "pct_bachelor",
  # Housing (1)
  "homeownership_rate",
  # Health (2)
  "disability_rate", "has_hicov"
)

cat(sprintf("✓ Using %d predictors\n", length(predictor_names)))

# ---- Simulation parameters ----
n_comp <- 50
n_cores <- 10 # Use 10 out of 12 cores (leave cores for OS)
thinning_param <- 0.5 # Data thinning parameter (only applies to k=1; for k>1 it is 1/k)
n_reps_dt <- 1 # Number of repetitions for data thinning (independent splits)

# Create configuration object with all shared parameters
sim_config <- list(
  # Simulation parameters
  n_comp = n_comp,
  n_cores = n_cores,
  thinning_param = thinning_param,
  n_reps_dt = n_reps_dt,

  # Data (NEW: design-based sampling from population)
  acs_pop = acs_pop,                    # Individual-level population for PPS sampling
  samp_frac = 0.0225,                   # Sampling fraction (2.25% - tested optimal)
  response_var = "medi_cal_qualified",  # Response variable
  X_approach = "population",            # Compute X from population (vs "estimated" from sample)

  # Paths
  results_dir = "_results",
  data_dir = "data"
)

# Save config object
saveRDS(sim_config, "sim_config.RDS")

cat(sprintf("\n✓ Configuration created:\n"))
cat(sprintf("  - %d comparisons\n", n_comp))
cat(sprintf("  - %d cores\n", n_cores))
cat(sprintf("  - Thinning parameter: %.2f (for k=1)\n", thinning_param))
cat(sprintf("  - DT repetitions: %d\n", n_reps_dt))
cat(sprintf("  - %d predictors\n", length(predictor_names)))
cat(sprintf("  - Response variable: %s\n", sim_config$response_var))
cat(sprintf("  - Sampling fraction: %.4f (%.2f%%)\n", sim_config$samp_frac, 100*sim_config$samp_frac))
cat(sprintf("  - X approach: %s\n\n", sim_config$X_approach))

# ==============================================================================
# PART 2: SETUP COMPARISONS
# ==============================================================================

cat("=== PART 2: SETUP COMPARISONS ===\n\n")
cat(sprintf("Setting up folders for %d comparisons (each with unique z & d)...\n", n_comp))

source("sim_functions/sampling_and_setup.R")
setup_comp(ncomps = n_comp, results_dir = sim_config$results_dir)
setup_esim(ncomps = n_comp, results_dir = sim_config$results_dir)

cat("✓ Comparison folders created with design-based variances\n\n")

# ==============================================================================
# PART 3: RUN METHODS
# ==============================================================================

cat("================================================================================\n")
cat("=== PART 3: RUN MODEL SELECTION METHODS ===\n")
cat("================================================================================\n\n")

# ==============================================================================
# METHOD 1: FULL DATA FIT (Oracle - fit on observed z, evaluate on true means)
# ==============================================================================
# COMMENTED OUT: Spatial model not yet implemented
# Will uncomment when models/spatial_basis_fit.R is ready
#
# #-------------------------------------------------------------------------------
# cat("--- METHOD 1: FULL DATA FIT (Oracle fits) ---\n")
# cat(sprintf("Expected time: ~10-12 minutes for %d comparisons (3 models each)\n\n", n_comp))
#
# source("sim_functions/full_data_fit.R")
# timer <- proc.time()
# options(digits.secs = 0)
# cat(sprintf("Start time: %s\n", Sys.time()))
#
# cl <- makeForkCluster(n_cores)
# registerDoParallel(cl)
# cat(sprintf("Using %d cores\n\n", length(cl)))
#
# foreach(j = 1:n_comp) %dopar% {
#   tryCatch(
#     {
#       full_data_fit(j, sim_config$results_dir)  # Updated signature
#       print(paste0("✓ Full data fit comparison #", j, " completed"))
#     },
#     error = function(e) {
#       message(sprintf("✗ Error in comparison %d: %s", j, e$message))
#       error_file <- file.path(sim_config$results_dir, sprintf("comparison_%03d_fit_ERROR.RDS", j))
#       saveRDS(list(comp_no = j, error = e$message, traceback = capture.output(traceback())),
#         file = error_file
#       )
#       return(NULL)
#     }
#   )
# }
#
# stopCluster(cl)
# elapsed <- (proc.time() - timer)[3]
# cat(sprintf("\n✓ Full data fit completed in %.1f minutes\n", elapsed / 60))
#
# # Preview results
# source("sim_functions/summary_oracle.R")
# cl <- makeForkCluster(n_cores)
# zfit_res <- parLapply(
#   cl, 1:n_comp,
#   function(x) {
#     summary_oracle(x, sim_config$results_dir) %>%  # Updated signature
#       group_by(comp_no) %>%
#       mutate(mse_true_rank = rank(mse_true))
#   }
# ) %>% bind_rows()
# stopCluster(cl)
#
# cat("\nTop model selections (by true MSE):\n")
# zfit_res %>%
#   filter(mse_true_rank == 1) %>%
#   group_by(model) %>%
#   reframe(count = n()) %>%
#   arrange(desc(count)) %>%
#   rename(`True No.1 count` = count) %>%
#   print()
#
# cat("\nAverage MSE by model:\n")
# zfit_res %>%
#   group_by(model) %>%
#   reframe(true_mse = mean(mse_true)) %>%
#   arrange(true_mse) %>%
#   rename(Av_MSE = true_mse) %>%
#   mutate(pct_diff = 100 * (Av_MSE - first(Av_MSE)) / first(Av_MSE)) %>%
#   print()
# cat("\n")
# cat(sprintf("\n========================================\n"))

cat("⚠ Skipping METHOD 1 (Full Data Fit) - spatial model not yet implemented\n")
cat(sprintf("========================================\n\n"))

#-------------------------------------------------------------------------------
cat("--- METHOD 2: DT 1-FOLD (Data Thinning - Single Fold) ---\n")
cat("Testing 15 epsilon values x 3 repeat counts = 45 configurations\n")
cat("Approx. time: ~2.5 hr for 50 comparisons\n\n")

source("sim_functions/run_dt.R")

# Define epsilon and repeat configurations
epsilon_values <- c(0.01, 0.02, 0.03, 0.04, seq(from = 0.05, to = 0.8, by = 0.05))
repeat_counts <- c(1, 3, 5)

# Sanity check epsilon values
cat(sprintf("Epsilon values to test: %s\n", paste(epsilon_values, collapse = ", ")))
cat(sprintf(
  "Total configurations: %d (epsilon) x %d (repeats) = %d\n\n",
  length(epsilon_values), length(repeat_counts),
  length(epsilon_values) * length(repeat_counts)
))

overall_timer <- proc.time()
cat(sprintf("Start time: %s\n\n", Sys.time()))

total_configs <- length(epsilon_values) * length(repeat_counts)
current_config <- 0

for (eps in epsilon_values) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("EPSILON = %.2f\n", eps))
  cat(sprintf("========================================\n\n"))

  for (n_reps in repeat_counts) {
    cat(sprintf("--- Running DT with eps=%.2f, %d rep(s) ---\n", eps, n_reps))
    timer <- proc.time()

    # RUN PHASE: Parallel over comparisons
    cl <- makeForkCluster(n_cores)
    registerDoParallel(cl)

    foreach(comp = 1:n_comp) %dopar% {
      tryCatch(
        {
          source("sim_functions/run_dt.R")
          run_dt(
            comp_no = comp, k = 1,
            results_dir = sim_config$results_dir,
            eps = eps, n_reps = n_reps
          )
        },
        error = function(e) {
          message(sprintf("✗ Error comp %d eps %.2f reps %d: %s", comp, eps, n_reps, e$message))
          return(NULL)
        }
      )
    }

    stopCluster(cl)
    elapsed <- (proc.time() - timer)[3]
    current_config <- current_config + 1
    cat(sprintf(
      "✓ Completed in %.1f minutes (configuration %d/%d)\n\n",
      elapsed / 60, current_config, total_configs
    ))
  }
}

overall_elapsed <- (proc.time() - overall_timer)[3]
cat(sprintf("\n✓ All DT 1-fold runs completed in %.1f minutes\n\n", overall_elapsed / 60))

#-------------------------------------------------------------------------------
cat("--- METHOD 3: DT 5-FOLD (Data Thinning - 5 Folds) ---\n")
cat(sprintf("Expected time: ~11 minutes for %d comparisons (3 models each)\n\n", n_comp))

source("sim_functions/run_dt.R")
timer <- proc.time()
cat(sprintf("Start time: %s\n", Sys.time()))

cl <- makeForkCluster(n_cores)
registerDoParallel(cl)
cat(sprintf("Using %d cores\n\n", length(cl)))

foreach(comp = 1:n_comp) %dopar% {
  tryCatch(
    {
      source("sim_functions/run_dt.R")
      run_dt(
        comp_no = comp, k = 5,
        results_dir = sim_config$results_dir
      )
      print(paste0("✓ DT 5-fold comparison #", comp, " completed"))
    },
    error = function(e) {
      error_file <- file.path(sim_config$results_dir, sprintf("comparison_%03d/dt5_ERROR.txt", comp))
      writeLines(c(
        sprintf("Error in DT 5-fold comparison %d:", comp),
        e$message,
        "",
        "Traceback:",
        capture.output(traceback())
      ), error_file)
      message(sprintf("✗ Error in DT 5-fold comparison %d: %s (see %s)", comp, e$message, error_file))
      return(NULL)
    }
  )
}

stopCluster(cl)
elapsed <- (proc.time() - timer)[3]
cat(sprintf("\n✓ DT 5-fold completed in %.1f minutes\n\n", elapsed / 60))

#-------------------------------------------------------------------------------
cat("--- METHOD 4: ESIM (Empirical Simulation) ---\n")
cat(sprintf("Using 100 iterations per comparison\n"))
cat(sprintf("Expected time: ~30 minutes for 50 comparisons (3 models each)\n\n", n_comp))

source("sim_functions/run_esim.R")
timer <- proc.time()
cat(sprintf("Start time: %s\n", Sys.time()))

cl <- makeForkCluster(n_cores)
registerDoParallel(cl)
cat(sprintf("Using %d cores\n\n", length(cl)))

foreach(comp = 1:n_comp) %dopar% {
  tryCatch(
    {
      source("sim_functions/run_esim.R")
      run_esim(comp,
        n_cores = 1,
        results_dir = sim_config$results_dir, n_iters = 100
      )
      print(paste0("✓ ESIM comparison #", comp, " completed"))
    },
    error = function(e) {
      message(sprintf("✗ Error in ESIM comparison %d: %s", comp, e$message))
      error_file <- file.path(sim_config$results_dir, sprintf("comparison_%03d_esim_ERROR.RDS", comp))
      saveRDS(list(comp_no = comp, error = e$message), file = error_file)
      return(NULL)
    }
  )
}

stopCluster(cl)
elapsed <- (proc.time() - timer)[3]
cat(sprintf("\n✓ ESIM completed in %.1f minutes\n\n", elapsed / 60))

# ==============================================================================
# PART 4: CALCULATE SUMMARY METRICS (FLEXIBLE)
# ==============================================================================

cat("================================================================================\n")
cat("=== PART 4: CALCULATE SUMMARY METRICS ===\n")
cat("================================================================================\n\n")

source("sim_functions/summary_esim.R")
source("sim_functions/summary_dt.R")
source("sim_functions/summary_zfit.R")

cat("Calculating metrics for all methods in parallel...\n\n")

# Set up cluster for summary calculations
cl <- makeForkCluster(n_cores)
registerDoParallel(cl)

#-------------------------------------------------------------------------------
# Compute summaries for all DT 1-fold configurations (epsilon × repeats × loss)
cat("Computing DT 1-fold summaries (135 configurations: 15 eps × 3 reps × 3 loss)...\n")

loss_functions <- c("MSE", "plugin_NLL", "predictive_NLL")
dt_results_list <- list()
summary_timer <- proc.time()

for (eps in epsilon_values) {
  for (n_reps in repeat_counts) {
    for (loss in loss_functions) {
      result <- parLapply(cl, 1:n_comp,
        function(x, cfg) {
          source("sim_functions/summary_dt.R")
          summary_dt(x,
            n_folds = 1, loss_function = cfg$loss,
            results_dir = cfg$results_dir,
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
          results_dir = sim_config$results_dir
        )
      ) %>% bind_rows()

      key <- sprintf("dt1_eps%.2f_rep%d_%s", eps, n_reps, loss)
      dt_results_list[[key]] <- result
    }
  }
}
cat(sprintf("✓ Completed in %.1f seconds\n\n", (proc.time() - summary_timer)[3]))

#-------------------------------------------------------------------------------
# Compute summaries for DT 5-fold (3 loss functions)
cat("Computing DT 5-fold summaries (3 loss functions)...\n")
summary_timer <- proc.time()

for (loss in loss_functions) {
  result <- parLapply(cl, 1:n_comp,
    function(x, cfg) {
      source("sim_functions/summary_dt.R")
      summary_dt(x,
        n_folds = 5, loss_function = cfg$loss,
        results_dir = cfg$results_dir,
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
      loss = loss, results_dir = sim_config$results_dir
    )
  ) %>%
    bind_rows()

  key <- sprintf("dt5_%s", loss)
  dt_results_list[[key]] <- result
}
cat(sprintf("✓ Completed in %.1f seconds\n\n", (proc.time() - summary_timer)[3]))

# Combine all DT results
dt_results <- bind_rows(dt_results_list)

#-------------------------------------------------------------------------------
cat("Computing ESIM summaries (2 methods: standard & data fission)...\n")
summary_timer <- proc.time()

esim_standard <- parLapply(
  cl, 1:n_comp,
  function(x, cfg) {
    summary_esim(x,
      results_dir = cfg$results_dir,
      validation = "standard"
    ) %>%
      group_by(comp_no) %>%
      mutate(
        rank = rank(metric),
        method_name = "esim_standard",
        loss_function = NA
      )
  },
  cfg = list(results_dir = sim_config$results_dir)
) %>% bind_rows()

esim_fission <- parLapply(
  cl, 1:n_comp,
  function(x, cfg) {
    summary_esim(x,
      results_dir = cfg$results_dir,
      validation = "data_fission"
    ) %>%
      group_by(comp_no) %>%
      mutate(
        rank = rank(metric),
        method_name = "esim_fission",
        loss_function = NA
      )
  },
  cfg = list(results_dir = sim_config$results_dir)
) %>% bind_rows()

cat(sprintf("✓ Completed in %.1f seconds\n\n", (proc.time() - summary_timer)[3]))

# Combine both ESIM approaches
esim_results <- bind_rows(esim_standard, esim_fission)

#-------------------------------------------------------------------------------
# ZFIT results and performance evaluation COMMENTED OUT (model not ready)
# Will uncomment when spatial model is implemented
#
# # Process ZFIT results computed earlier
# cat("Processing ZFIT results...\n")
# zfit_results <- zfit_res %>%
#   mutate(
#     method_name = "zfit",
#     loss_function = NA
#   )

stopCluster(cl)

cat("\n⚠ Skipping performance evaluation - depends on full data fit (not yet implemented)\n\n")

# ===============================================================================
# PERFORMANCE EVALUATION (COMMENTED OUT - depends on zfit_results)
# ===============================================================================
#
# cat("\n================================================================================\n")
# cat("=== MODEL SELECTION PERFORMANCE PREVIEW ===\n")
# cat("================================================================================\n\n")

# # Get oracle best model for each comparison (lowest oracle MSE)
# oracle_best <- zfit_results %>%
#   select(comp_no, model, mse_true) %>%
#   group_by(comp_no) %>%
#   slice_min(mse_true, n = 1, with_ties = FALSE) %>%
#   select(comp_no, best_model = model)
#
# # Function to compute selection accuracy and Kendall's Tau
# compute_performance <- function(data, oracle) {
#   # Get all models that exist in oracle (1_cov, 4_cov, 7_cov)
#   # This excludes ESIM's "Direct" model
#   oracle_all_models <- unique(zfit_results$model)
#
#   # Filter to models that exist in oracle
#   data_filtered <- data %>%
#     filter(model %in% oracle_all_models)
#
#   # ACCURACY: Count how many times method selects same model as per-comparison oracle
#   selections <- data_filtered %>%
#     group_by(comp_no) %>%
#     slice_min(metric, n = 1, with_ties = FALSE) %>%
#     ungroup()
#
#   # Join with oracle to compare selections
#   accuracy <- selections %>%
#     left_join(oracle, by = "comp_no") %>%
#     summarise(correct = sum(model == best_model, na.rm = TRUE)) %>%
#     pull(correct)
#
#   # Compute Kendall's Tau
#   # Compare to per-comparison oracle ranks
#   oracle_ranks <- zfit_results %>%
#     select(comp_no, model, oracle_rank = mse_true_rank) %>%
#     arrange(comp_no, model)
#
#   method_ranks <- data %>%
#     filter(model %in% oracle_all_models) %>%
#     select(comp_no, model, rank) %>%
#     arrange(comp_no, model)
#
#   combined <- method_ranks %>%
#     left_join(oracle_ranks, by = c("comp_no", "model"))
#
#   # Compute tau for each comparison separately, then average
#   per_comp_tau <- combined %>%
#     group_by(comp_no) %>%
#     summarise(tau = cor(rank, oracle_rank, method = "kendall"), .groups = "drop")
#
#   tau <- mean(per_comp_tau$tau, na.rm = TRUE)
#
#   return(list(accuracy = accuracy, tau = tau))
# }
#
# # Combine all results for performance evaluation
# all_results_temp <- bind_rows(
#   dt_results,
#   esim_results,
#   zfit_results %>%
#     select(comp_no, model, DIC, WAIC, mse_true_rank) %>%
#     pivot_longer(cols = c(DIC, WAIC), names_to = "metric_type", values_to = "metric") %>%
#     group_by(comp_no, metric_type) %>%
#     mutate(
#       rank = rank(metric),
#       method_name = "zfit",
#       loss_function = NA
#     ) %>%
#     select(comp_no, method_name, loss_function, model, metric, rank, metric_type)
# )
#
# # Compute for all DT 1-fold configurations
# dt1_performance <- all_results_temp %>%
#   filter(method_name == "dt_1fold") %>%
#   group_by(epsilon, n_reps_used, loss_function) %>%
#   group_modify(~ {
#     perf <- compute_performance(.x, oracle_best)
#     tibble(accuracy = perf$accuracy, tau = perf$tau)
#   }) %>%
#   ungroup() %>%
#   arrange(desc(accuracy), desc(tau))
#
# # Compute for DT 5-fold
# dt5_performance <- all_results_temp %>%
#   filter(method_name == "dt_5fold") %>%
#   group_by(loss_function) %>%
#   group_modify(~ {
#     perf <- compute_performance(.x, oracle_best)
#     tibble(accuracy = perf$accuracy, tau = perf$tau)
#   }) %>%
#   ungroup() %>%
#   arrange(desc(accuracy), desc(tau))
#
# # Compute for ESIM
# esim_performance <- all_results_temp %>%
#   filter(method_name %in% c("esim_standard", "esim_fission")) %>%
#   group_by(method_name) %>%
#   group_modify(~ {
#     perf <- compute_performance(.x, oracle_best)
#     tibble(accuracy = perf$accuracy, tau = perf$tau)
#   }) %>%
#   ungroup() %>%
#   arrange(desc(accuracy), desc(tau))
#
# # Compute for full-data methods (DIC, WAIC)
# zfit_performance <- all_results_temp %>%
#   filter(method_name == "zfit", metric_type %in% c("DIC", "WAIC")) %>%
#   group_by(metric_type) %>%
#   group_modify(~ {
#     perf <- compute_performance(.x, oracle_best)
#     tibble(accuracy = perf$accuracy, tau = perf$tau)
#   }) %>%
#   ungroup() %>%
#   mutate(
#     method_label = case_when(
#       metric_type == "DIC" ~ "DIC",
#       metric_type == "WAIC" ~ "WAIC",
#       TRUE ~ metric_type
#     )
#   ) %>%
#   select(method = method_label, accuracy, tau) %>%
#   arrange(desc(accuracy), desc(tau))
#
# # Combined summary sorted by performance
# cat("\n================================================================================\n")
# # Compute method counts dynamically
# dt1_methods <- length(epsilon_values) * length(repeat_counts) * length(loss_functions)
# dt5_methods <- length(loss_functions)
# esim_methods <- 2
# zfit_methods <- 2 # DIC, WAIC
# total_methods <- dt1_methods + dt5_methods + esim_methods + zfit_methods
# cat(sprintf(
#   "=== OVERALL RANKING (all %d methods, sorted by accuracy) ===\n",
#   total_methods
# ))
# cat("================================================================================\n\n")
#
# overall_summary <- bind_rows(
#   dt1_performance %>%
#     mutate(method = sprintf(
#       "DT 1-fold (ε=%.2f, reps=%d, %s)",
#       epsilon, n_reps_used, loss_function
#     )) %>%
#     select(method, accuracy, tau),
#   dt5_performance %>%
#     mutate(method = sprintf("DT 5-fold (%s)", loss_function)) %>%
#     select(method, accuracy, tau),
#   esim_performance %>%
#     mutate(method = ifelse(method_name == "esim_standard",
#       "ESIM (standard)", "ESIM (data fission)"
#     )) %>%
#     select(method, accuracy, tau),
#   zfit_performance # Already has method, accuracy, tau columns with proper names
# ) %>%
#   arrange(desc(accuracy), desc(tau))
#
# print(overall_summary, n = 100)
#
# cat("\n")
# cat(sprintf("Best performing method: %s\n", overall_summary$method[1]))
# cat(sprintf(
#   "  Accuracy: %d/%d (%.1f%%)\n",
#   overall_summary$accuracy[1], n_comp,
#   100 * overall_summary$accuracy[1] / n_comp
# ))
# cat(sprintf("  Kendall's Tau: %.3f\n", overall_summary$tau[1]))
# cat("\n")

# ==============================================================================
# PART 5: SAVE RESULTS
# ==============================================================================

cat("================================================================================\n")
cat("=== PART 5: SAVE RESULTS ===\n")
cat("================================================================================\n\n")

# Combine DT and ESIM results (zfit not available yet)
cat("Saving DT and ESIM results...\n")

all_results <- bind_rows(
  dt_results,
  esim_results
)

cat("Saving results to results.RDS...\n")
saveRDS(all_results, "results.RDS")

# Save individual method results for easier access
results_list <- list(
  dt_results = dt_results,
  esim_results = esim_results
)
saveRDS(results_list, "results_list.RDS")

cat("\n✓ Results saved successfully!\n")
cat("  Note: Full data fit results not included (model not yet implemented)\n\n")

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
  writeLines(
    sprintf(
      "Run: %s\nDate: %s\nComparisons: %d",
      SAVE_AS, Sys.time(), sim_config$n_comp
    ),
    file.path(save_dir, "info.txt")
  )

  cat(sprintf("✓ Saved to %s\n", save_dir))
}

cat("================================================================================\n")
cat("  PIPELINE COMPLETE!\n")
cat("================================================================================\n\n")

cat("Next steps:\n")
cat("  1. Review results: results <- readRDS('results.RDS')\n")
cat("  2. Run analysis: Open analysis.Rmd\n")
cat("  3. Check for errors: source('rerun_failed.R'); check_status()\n\n")
