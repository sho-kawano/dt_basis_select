#!/usr/bin/env Rscript
# Compute thinned oracle MSE for gap decomposition analysis
# For all equal allocation configs: 30, 40, 50, 75, 100, 125
# Computes MSE of DT fits (on z_train) vs true population theta

library(dplyr)
library(tidyr)
library(parallel)

# Detect cores for parallel processing
n_cores <- min(10, detectCores() - 1)
cat(sprintf("Using %d cores for parallel processing\n\n", n_cores))

# Config mapping
config_dirs <- c(
  equal_30 = "_results_equal30_epsilon",
  equal_40 = "_results_equal40_comparison",
  equal_50 = "_results_equal50_comparison",
  equal_75 = "_results_equal75_rerun",
  equal_100 = "_results_equal100_comparison",
  equal_125 = "_results_equal125_epsilon"
)

# nbasis grid (20 values)
nbasis_grid <- seq(3, 60, by = 3)

# epsilon grid (9 values for epsilon_8configs)
epsilon_grid <- seq(0.1, 0.9, by = 0.1)

# reps (1-5)
n_reps_options <- c(1, 3, 5)

# ============================================================================
# Function: Compute thinned oracle MSE for one (comp, eps, rep)
# ============================================================================
compute_thinned_oracle_single <- function(results_dir, comp_no, epsilon, rep_no, theta_true) {
  # Construct file path
  comp_dir <- sprintf("%s/comparison_%03d", results_dir, comp_no)
  eps_dir <- sprintf("%s/dt_1fold/eps_%.2f", comp_dir, epsilon)
  fit_file <- sprintf("%s/1fold_rep%d.RDS", eps_dir, rep_no)

  # Check file exists
  if (!file.exists(fit_file)) {
    warning(sprintf("File not found: %s", fit_file))
    return(NULL)
  }

  # Load DT fit (contains all nbasis values)
  fit_data <- readRDS(fit_file)

  # Join true theta by PUMA ID
  fit_with_truth <- fit_data %>%
    mutate(theta_true = theta_true$values[match(domain, theta_true$puma)])

  # Check alignment
  if (any(is.na(fit_with_truth$theta_true))) {
    warning(sprintf("PUMA alignment failed for comp=%d, eps=%.2f, rep=%d",
                    comp_no, epsilon, rep_no))
    return(NULL)
  }

  # Compute MSE for each nbasis
  # NOTE: DT fits estimate eps*theta, so scale by 1/eps to get theta estimates
  mse_results <- fit_with_truth %>%
    group_by(nbasis) %>%
    summarise(
      thinned_mse = mean((mean/epsilon - theta_true)^2),
      .groups = "drop"
    ) %>%
    mutate(
      comp_no = comp_no,
      epsilon = epsilon,
      rep_no = rep_no
    )

  return(mse_results)
}

# ============================================================================
# Function: Compute thinned oracle MSE for one comparison (all eps, reps)
# ============================================================================
compute_thinned_oracle_comparison <- function(comp_no, results_dir, theta_true,
                                               epsilon_grid, n_reps = 5) {
  cat(sprintf("  Comp %02d...", comp_no))

  results_list <- list()
  idx <- 1

  # Loop over epsilon values
  for (eps in epsilon_grid) {
    # Loop over reps
    for (rep_no in 1:n_reps) {
      result <- compute_thinned_oracle_single(results_dir, comp_no, eps, rep_no, theta_true)

      if (!is.null(result)) {
        results_list[[idx]] <- result
        idx <- idx + 1
      }
    }
  }

  cat(" done\n")
  return(bind_rows(results_list))
}

# ============================================================================
# Function: Compute thinned oracle MSE for one config
# ============================================================================
compute_thinned_oracle_config <- function(config_name, results_dir,
                                           n_comparisons = 20,
                                           epsilon_grid = seq(0.1, 0.9, 0.1),
                                           n_reps = 5,
                                           n_cores = 10) {
  cat(sprintf("\n=== Processing config: %s ===\n", config_name))
  cat(sprintf("Results dir: %s\n", results_dir))

  # Load theta_true once for this config
  theta_true_file <- file.path(results_dir, "theta_true.RDS")
  if (!file.exists(theta_true_file)) {
    stop(sprintf("theta_true.RDS not found in %s", results_dir))
  }

  theta_true <- readRDS(theta_true_file)
  cat(sprintf("Loaded theta_true: %d PUMAs\n", length(theta_true$values)))

  # Parallel over comparisons
  cat(sprintf("Computing over %d comparisons (parallel):\n", n_comparisons))

  results_list <- mclapply(1:n_comparisons, function(comp_no) {
    compute_thinned_oracle_comparison(
      comp_no = comp_no,
      results_dir = results_dir,
      theta_true = theta_true,
      epsilon_grid = epsilon_grid,
      n_reps = n_reps
    )
  }, mc.cores = n_cores)

  # Combine all comparisons
  all_data <- bind_rows(results_list) %>%
    mutate(config = config_name)

  cat(sprintf("Computed %d rows for config %s\n", nrow(all_data), config_name))

  return(all_data)
}

# ============================================================================
# Function: Aggregate by n_reps
# ============================================================================
aggregate_by_n_reps <- function(raw_data, n_reps_options = c(1, 3, 5)) {
  cat("\nAggregating by n_reps...\n")

  results_list <- list()

  for (n_reps in n_reps_options) {
    cat(sprintf("  n_reps = %d...", n_reps))

    agg <- raw_data %>%
      filter(rep_no <= n_reps) %>%
      group_by(config, comp_no, epsilon, nbasis) %>%
      summarise(
        thinned_mse = mean(thinned_mse),
        .groups = "drop"
      ) %>%
      mutate(n_reps_used = n_reps)

    results_list[[as.character(n_reps)]] <- agg
    cat(sprintf(" %d rows\n", nrow(agg)))
  }

  combined <- bind_rows(results_list)
  cat(sprintf("Total aggregated rows: %d\n", nrow(combined)))

  return(combined)
}

# ============================================================================
# Main execution
# ============================================================================
cat("\n")
cat(paste(rep("=", 70), collapse=""))
cat("\nTHINNED ORACLE MSE COMPUTATION\n")
cat(paste(rep("=", 70), collapse=""))
cat("\n\n")

cat("Configs to process:\n")
for (cfg in names(config_dirs)) {
  cat(sprintf("  %s: %s\n", cfg, config_dirs[cfg]))
}

# Store raw results (all reps separate)
raw_results_all <- list()

# Process each config
for (cfg in names(config_dirs)) {
  raw_results <- compute_thinned_oracle_config(
    config_name = cfg,
    results_dir = config_dirs[cfg],
    n_comparisons = 20,
    epsilon_grid = epsilon_grid,
    n_reps = 5,
    n_cores = n_cores
  )

  raw_results_all[[cfg]] <- raw_results
}

# Combine all configs
cat("\n=== Combining all configs ===\n")
all_raw_data <- bind_rows(raw_results_all)
cat(sprintf("Total raw rows: %d\n", nrow(all_raw_data)))
cat(sprintf("Expected: %d configs × 20 comps × 9 eps × 5 reps × 20 nbasis = %d\n",
            length(config_dirs), 6 * 20 * 9 * 5 * 20))

# Aggregate by n_reps
aggregated_data <- aggregate_by_n_reps(all_raw_data, n_reps_options)

# Verify structure
cat("\n=== Data verification ===\n")
cat("Columns:", paste(names(aggregated_data), collapse=", "), "\n")
cat("Configs:", paste(unique(aggregated_data$config), collapse=", "), "\n")
cat(sprintf("Epsilon range: [%.1f, %.1f]\n",
            min(aggregated_data$epsilon), max(aggregated_data$epsilon)))
cat(sprintf("nbasis range: [%d, %d] (%d values)\n",
            min(aggregated_data$nbasis), max(aggregated_data$nbasis),
            length(unique(aggregated_data$nbasis))))
cat(sprintf("n_reps values: %s\n",
            paste(sort(unique(aggregated_data$n_reps_used)), collapse=", ")))

# Save output
output_file <- "results_multi_config/thinned_oracle_6configs.RDS"
cat(sprintf("\nSaving to: %s\n", output_file))

# Ensure output directory exists
dir.create("results_multi_config", showWarnings = FALSE, recursive = TRUE)

saveRDS(aggregated_data, output_file)

output_size <- file.info(output_file)$size / 1024^2
cat(sprintf("File saved: %.2f MB\n", output_size))

# Summary statistics
cat("\n=== Summary statistics ===\n")

summary_stats <- aggregated_data %>%
  group_by(config, n_reps_used) %>%
  summarise(
    n_rows = n(),
    mean_mse = mean(thinned_mse),
    min_mse = min(thinned_mse),
    max_mse = max(thinned_mse),
    .groups = "drop"
  )

print(summary_stats)

# Quick validation: Check one case
cat("\n=== Quick validation ===\n")
test_case <- aggregated_data %>%
  filter(config == "equal_40", comp_no == 1, epsilon == 0.5,
         n_reps_used == 5, nbasis == 18)

if (nrow(test_case) == 1) {
  cat(sprintf("equal_40, comp=1, eps=0.5, n_reps=5, nbasis=18: MSE = %.6f\n",
              test_case$thinned_mse))
} else {
  warning(sprintf("Expected 1 row for test case, got %d", nrow(test_case)))
}

cat("\n")
cat(paste(rep("=", 70), collapse=""))
cat("\nCOMPUTATION COMPLETE ✓\n")
cat(paste(rep("=", 70), collapse=""))
cat("\n\n")
cat("Output saved to:", output_file, "\n")
cat("Next: Run analysis/gap_decomposition.Rmd to analyze results\n\n")
