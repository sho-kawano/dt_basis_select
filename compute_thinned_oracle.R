#!/usr/bin/env Rscript
# Compute thinned oracle MSE for S=50 comparisons
# Thinned oracle = MSE of DT fits (on z_train) vs true theta
# This measures the "gap" introduced by thinning

library(dplyr)
library(tidyr)
library(parallel)

n_cores <- min(10, detectCores() - 1)
cat(sprintf("Using %d cores\n\n", n_cores))

# Config mapping to S=50 directories
config_dirs <- c(
  equal_30 = "_results_equal30",
  equal_40 = "_results_equal40",
  equal_50 = "_results_equal50",
  equal_75 = "_results_equal75",
  equal_100 = "_results_equal100",
  equal_125 = "_results_equal125"
)

epsilon_grid <- seq(0.1, 0.9, by = 0.1)
n_comparisons <- 50

# Compute thinned MSE for one (comp, eps, rep)
compute_thinned_single <- function(results_dir, comp_no, epsilon, rep_no, theta_true) {
  comp_dir <- sprintf("%s/comparison_%03d", results_dir, comp_no)
  eps_dir <- sprintf("%s/dt_1fold/eps_%.2f", comp_dir, epsilon)
  fit_file <- sprintf("%s/1fold_rep%d.RDS", eps_dir, rep_no)

  if (!file.exists(fit_file)) return(NULL)

  fit_data <- readRDS(fit_file)

  fit_with_truth <- fit_data %>%
    mutate(theta_true = theta_true$values[match(domain, theta_true$puma)])

  if (any(is.na(fit_with_truth$theta_true))) return(NULL)

  # DT fits estimate eps*theta, scale by 1/eps
  fit_with_truth %>%
    group_by(nbasis) %>%
    summarise(thinned_mse = mean((mean/epsilon - theta_true)^2), .groups = "drop") %>%
    mutate(comp_no = comp_no, epsilon = epsilon, rep_no = rep_no)
}

# Compute for all (eps, rep) for one comparison
compute_thinned_comparison <- function(comp_no, results_dir, theta_true, epsilon_grid, n_reps = 5) {
  results_list <- list()
  idx <- 1
  for (eps in epsilon_grid) {
    for (rep_no in 1:n_reps) {
      result <- compute_thinned_single(results_dir, comp_no, eps, rep_no, theta_true)
      if (!is.null(result)) {
        results_list[[idx]] <- result
        idx <- idx + 1
      }
    }
  }
  bind_rows(results_list)
}

# Compute for one config
compute_thinned_config <- function(config_name, results_dir, n_comparisons, epsilon_grid, n_reps, n_cores) {
  cat(sprintf("\n=== %s ===\n", config_name))

  theta_true <- readRDS(file.path(results_dir, "theta_true.RDS"))
  cat(sprintf("Loaded theta_true: %d PUMAs\n", length(theta_true$values)))

  results_list <- mclapply(1:n_comparisons, function(comp_no) {
    compute_thinned_comparison(comp_no, results_dir, theta_true, epsilon_grid, n_reps)
  }, mc.cores = n_cores)

  bind_rows(results_list) %>% mutate(config = config_name)
}

# Aggregate by n_reps
aggregate_by_n_reps <- function(raw_data, n_reps_options = c(1, 3, 5)) {
  results_list <- lapply(n_reps_options, function(n_reps) {
    raw_data %>%
      filter(rep_no <= n_reps) %>%
      group_by(config, comp_no, epsilon, nbasis) %>%
      summarise(thinned_mse = mean(thinned_mse), .groups = "drop") %>%
      mutate(n_reps_used = n_reps)
  })
  bind_rows(results_list)
}

# Main
cat("Computing thinned oracle MSE for S=50\n")
cat("Configs:", paste(names(config_dirs), collapse=", "), "\n")

raw_results_all <- list()
for (cfg in names(config_dirs)) {
  raw_results_all[[cfg]] <- compute_thinned_config(
    config_name = cfg,
    results_dir = config_dirs[cfg],
    n_comparisons = n_comparisons,
    epsilon_grid = epsilon_grid,
    n_reps = 5,
    n_cores = n_cores
  )
}

all_raw_data <- bind_rows(raw_results_all)
aggregated_data <- aggregate_by_n_reps(all_raw_data, c(1, 3, 5))

# Verify
cat("\n=== Verification ===\n")
cat("comp_no range:", range(aggregated_data$comp_no), "\n")
cat("n_distinct comp_no:", n_distinct(aggregated_data$comp_no), "\n")
cat("Total rows:", nrow(aggregated_data), "\n")

# Save
dir.create("results_summary", showWarnings = FALSE, recursive = TRUE)
output_file <- "results_summary/thinned_oracle_s50.RDS"
saveRDS(aggregated_data, output_file)
cat(sprintf("\nSaved to: %s\n", output_file))
