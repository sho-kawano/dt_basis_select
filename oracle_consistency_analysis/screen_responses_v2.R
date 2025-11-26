#!/usr/bin/env Rscript
# Phase 1 Screening with Steepness Filter (Screening Criteria v2)
# Tests single dataset per response to filter candidates

library(dplyr)
library(doParallel)
library(foreach)

source("sim_functions/sampling_and_setup.R")
source("sim_functions/full_data_fit.R")
source("sim_functions/summary_oracle.R")

# === CONFIGURATION ===
# Set dataset and output directory (allow override from command line)
if (!exists("DATASET")) DATASET <- "ny"  # Options: "ca", "tx", "ny", or "il"
nbasis_grid <- seq(10, 100, by = 10)  # Full grid for Phase 1
N_DATASETS <- 6  # Number of datasets for Phase 1 screening
BASE_SEED <- 12345  # Base seed for reproducibility
SKIP_MORANS_I <- TRUE  # Skip Moran's I pre-screening (uninformative)

# Phase 1 Filter Criteria (Screening Criteria v2)
INTERIOR_MIN <- 20
INTERIOR_MAX <- 90
MSE_DIFF_MIN <- 15
MSE_DIFF_MAX <- 1000  # Effectively no maximum (was 30)
PENALTY_AT_10_MIN <- 5.0  # Strict steepness threshold
MODELS_WITHIN_5PCT_MAX <- 3  # Alternative steepness metric

# Sampling configuration (equal_50)
base_config <- list(
  population_data = sprintf("data/%s_pums_population.rds", DATASET),
  adjacency_data = sprintf("data/%s_puma_adjacency.RDA", DATASET),
  equal_allocation = TRUE,
  equal_n = 50,
  min_sample_size = 30,
  samp_frac = 0.01,
  response_filter = NULL,
  X_covariates = NULL,   # Intercept-only models
  c_prior = 0.001,       # Weak spatial variance prior
  d_prior = 0.001        # Weak IID variance prior
)

# Candidate response variables
candidates <- tribble(
  ~variable,      ~label,                    ~type,
  # Binary responses
  "PUBCOV",       "public_health_coverage",  "binary",
  "employed",     "employed",                "binary",
  "owns_home",    "owns_home",               "binary",
  "rents_home",   "rents_home",              "binary",
  "married",      "married",                 "binary",
  "citizen",      "citizen",                 "binary",
  "bachelor_plus", "bachelor_plus",          "binary",
  "disabled",     "disabled",                "binary",
  "has_vehicle",  "has_vehicle",             "binary",
  "crowded",      "crowded",                 "binary",
  # Continuous responses
  "WAGP",         "annual_wages",            "continuous",
  "POVPIP",       "poverty_ratio",           "continuous",
  "AGEP",         "age",                     "continuous",
  "JWMNP",        "commute_time",            "continuous",
  "BDSP",         "num_bedrooms",            "continuous"
)

cat(sprintf("\n=== PHASE 1 SCREENING: %s Response Variables ===\n\n", toupper(DATASET)))
cat(sprintf("Initial candidates: %d responses\n\n", nrow(candidates)))

# Create output directory
results_dir_base <- sprintf("oracle_consistency_analysis/response_search/%s_phase1", DATASET)
dir.create(results_dir_base, recursive = TRUE, showWarnings = FALSE)

if (!SKIP_MORANS_I) {
  # === STEP 0: Moran's I Pre-screening ===
  cat("STEP 0: Moran's I Spatial Autocorrelation Pre-screening\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")

  MORANS_I_MIN <- 0.5  # Minimum spatial autocorrelation threshold

  # Load population data and adjacency matrix
  cat(sprintf("Loading %s population and adjacency data...\n", toupper(DATASET)))
  pop_data <- readRDS(base_config$population_data)
  load(base_config$adjacency_data)  # Loads A (adjacency matrix)

  library(spdep)

# Calculate Moran's I for each response
morans_results <- data.frame()

for (i in 1:nrow(candidates)) {
  var_name <- candidates$variable[i]
  var_label <- candidates$label[i]
  var_type <- candidates$type[i]

  # Create derived variables (match sampling_and_setup.R logic)
  pop_data_derived <- pop_data

  if (var_name == "employed") {
    pop_data_derived$employed <- as.numeric(pop_data_derived$ESR %in% c('1', '2', '4', '5'))
  } else if (var_name == "owns_home") {
    pop_data_derived$owns_home <- as.numeric(pop_data_derived$TEN %in% c(1, 2))
  } else if (var_name == "rents_home") {
    pop_data_derived$rents_home <- as.numeric(pop_data_derived$TEN == 3)
  } else if (var_name == "married") {
    pop_data_derived$married <- as.numeric(pop_data_derived$MAR == 1)
  } else if (var_name == "citizen") {
    pop_data_derived$citizen <- as.numeric(pop_data_derived$CIT %in% c(1, 2, 3, 4))
  } else if (var_name == "bachelor_plus") {
    pop_data_derived$bachelor_plus <- as.numeric(as.numeric(pop_data_derived$SCHL) >= 21)
  } else if (var_name == "disabled") {
    pop_data_derived$disabled <- as.numeric(pop_data_derived$DIS == 1)
  } else if (var_name == "has_vehicle") {
    pop_data_derived$has_vehicle <- as.numeric(pop_data_derived$VEH > 0)
  } else if (var_name == "crowded") {
    # Crowded: more than 1 person per room (assuming BDSP is rooms)
    pop_data_derived <- pop_data_derived %>%
      group_by(SERIALNO) %>%
      mutate(hh_size = n(),
             crowded = as.numeric(hh_size / BDSP > 1)) %>%
      ungroup()
  } else if (var_name == "PUBCOV") {
    # PUBCOV: 1=covered, 2=not covered
    pop_data_derived$PUBCOV <- as.numeric(pop_data_derived$PUBCOV == 1)
  }

  # Compute population-level response means per PUMA
  puma_means <- pop_data_derived %>%
    group_by(PUMA) %>%
    summarise(value = mean(get(var_name), na.rm = TRUE), .groups = "drop") %>%
    arrange(PUMA)

  # Create spatial weights from adjacency matrix
  W_list <- mat2listw(A, style = "W")

  # Calculate Moran's I
  morans_test <- moran.test(puma_means$value, W_list)
  morans_i <- morans_test$estimate["Moran I statistic"]

  morans_results <- rbind(morans_results, data.frame(
    variable = var_name,
    label = var_label,
    type = var_type,
    morans_i = morans_i,
    pass_morans = morans_i >= MORANS_I_MIN,
    stringsAsFactors = FALSE
  ))

  cat(sprintf("%-25s %s: Moran's I = %.3f %s\n",
              var_label,
              sprintf("[%s]", var_type),
              morans_i,
              ifelse(morans_i >= MORANS_I_MIN, "✓", "✗")))
}

cat(paste(rep("-", 80), collapse = ""), "\n")

# Filter candidates by Moran's I
candidates_filtered <- morans_results %>%
  filter(pass_morans) %>%
  select(variable, label, type)

cat(sprintf("\nPassed Moran's I >= %.2f: %d / %d responses\n\n",
            MORANS_I_MIN, nrow(candidates_filtered), nrow(candidates)))

if (nrow(candidates_filtered) == 0) {
  cat("ERROR: No responses passed Moran's I pre-screening!\n")
  cat("Consider lowering threshold or testing all responses.\n")
  quit(status = 1)
}

# Save Moran's I results
write.csv(morans_results,
          file.path(results_dir_base, "morans_i_prescreening.csv"),
          row.names = FALSE)

# Update candidates to filtered set
candidates <- candidates_filtered

} else {
  cat("Skipping Moran's I pre-screening (SKIP_MORANS_I = TRUE)\n\n")
}

cat(sprintf("=== EMPIRICAL SCREENING: %d Candidates ===\n\n", nrow(candidates)))
cat(sprintf("nbasis grid: %s\n", paste(nbasis_grid, collapse = ", ")))
cat(sprintf("Datasets per response: %d\n\n", N_DATASETS))

cat("Screening Criteria v2:\n")
if (!SKIP_MORANS_I) {
  cat(sprintf("  Filter 0: Moran's I >= %.2f (pre-screen)\n", MORANS_I_MIN))
}
cat(sprintf("  Filter 1: Interior minimum [%d, %d]\n", INTERIOR_MIN, INTERIOR_MAX))
cat(sprintf("  Filter 2: MSE differentiation [%d%%, %d%%]\n", MSE_DIFF_MIN, MSE_DIFF_MAX))
cat(sprintf("  Filter 3a: Penalty at ±10 units >= %.1f%%\n", PENALTY_AT_10_MIN))
cat(sprintf("  Filter 3b: Models within 5%% of optimal <= %d\n\n", MODELS_WITHIN_5PCT_MAX))

# Setup parallel cluster
n_cores <- min(10, parallel::detectCores() - 1)
cl <- makeCluster(n_cores, type = "FORK")
registerDoParallel(cl)

# Create all tasks (response × dataset combinations)
tasks <- expand.grid(
  response_idx = 1:nrow(candidates),
  dataset_id = 1:N_DATASETS,
  stringsAsFactors = FALSE
)

cat(sprintf("Total tasks: %d (responses × datasets)\n", nrow(tasks)))
cat(sprintf("Running in parallel (%d cores)...\n\n", n_cores))

# Run all tasks in parallel
results_list <- foreach(i = 1:nrow(tasks), .errorhandling = "pass") %dopar% {
  task <- tasks[i, ]
  response_info <- candidates[task$response_idx, ]
  var_name <- response_info$variable
  var_label <- response_info$label
  var_type <- response_info$type
  dataset_id <- task$dataset_id

  comp_no <- dataset_id
  results_dir <- sprintf("%s/%s", results_dir_base, var_name)

  # Set response variable in config
  config <- base_config
  config$response_var <- var_name
  config$response_type <- var_type

  tryCatch({
    # Setup comparison
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    set.seed(BASE_SEED + dataset_id)
    setup_comp(comp_no, results_dir, config)

    # Check for zero variance PUMAs
    d_data <- readRDS(file.path(results_dir, sprintf("comparison_%03d", comp_no), "d.RDS"))
    if (any(d_data$values == 0)) {
      return(data.frame(
        variable = var_name,
        label = var_label,
        type = var_type,
        dataset_id = dataset_id,
        status = "zero_variance",
        stringsAsFactors = FALSE
      ))
    }

    # Run oracle fits
    fit_config <- config
    fit_config$nbasis_values <- nbasis_grid
    fit_config$spatial_type <- "fixed"
    fit_config$ndesired <- 2000
    fit_config$nburn <- 1500
    fit_config$nthin <- 1
    fit_config$hyp <- list(a = 1, b = 1, c = config$c_prior, d = config$d_prior)

    full_data_fit(comp_no, results_dir, fit_config)

    # Summarize oracle results
    oracle_summary <- summary_oracle(comp_no, results_dir)
    oracle_summary$variable <- var_name
    oracle_summary$label <- var_label
    oracle_summary$type <- var_type
    oracle_summary$dataset_id <- dataset_id
    oracle_summary$status <- "success"

    return(oracle_summary)

  }, error = function(e) {
    return(data.frame(
      variable = var_name,
      label = var_label,
      type = var_type,
      dataset_id = dataset_id,
      status = paste("error:", e$message),
      stringsAsFactors = FALSE
    ))
  })
}

stopCluster(cl)

cat("\n\nParallel execution complete!\n\n")

# Combine results
cat("Combining results...\n")
results <- bind_rows(results_list)

# Filter successful results
successful_results <- results %>%
  filter(status == "success" | is.na(status))  # NA status means it has oracle data

# Check if we have results
if (nrow(successful_results) == 0) {
  cat("ERROR: No successful runs!\n")
  cat("\nError summary:\n")
  print(table(results$status))
  quit(status = 1)
}

cat(sprintf("Successful runs: %d / %d\n\n",
    nrow(successful_results %>% filter(!is.na(nbasis))),
    nrow(tasks) * length(nbasis_grid)))

# Calculate Phase 1 screening metrics per response
screening_results <- data.frame()

for (var in unique(successful_results$variable)) {
  var_data <- successful_results %>% filter(variable == var, !is.na(nbasis))
  var_label <- unique(successful_results %>% filter(variable == var) %>% pull(label))[1]
  var_type <- unique(successful_results %>% filter(variable == var) %>% pull(type))[1]

  if (nrow(var_data) == 0) {
    cat(sprintf("WARNING: %s - No valid data\n\n", var_label))
    next
  }

  # Find optimal nbasis per dataset
  optimal_per_dataset <- var_data %>%
    group_by(dataset_id) %>%
    summarise(optimal_nbasis = nbasis[which.min(mse_true)], .groups = "drop")

  n_datasets <- nrow(optimal_per_dataset)
  opt_sd <- sd(optimal_per_dataset$optimal_nbasis)
  opt_median <- median(optimal_per_dataset$optimal_nbasis)

  # Boundary selections
  n_boundary <- sum(optimal_per_dataset$optimal_nbasis %in% c(min(nbasis_grid), max(nbasis_grid)))
  pct_boundary <- (n_boundary / n_datasets) * 100

  # MSE variation (average across datasets)
  mse_by_nbasis <- var_data %>%
    group_by(nbasis) %>%
    summarise(mean_mse = mean(mse_true), .groups = "drop")

  min_mse <- min(mse_by_nbasis$mean_mse)
  max_mse <- max(mse_by_nbasis$mean_mse)

  # Filter 1: Interior minimum (use median)
  interior_min <- opt_median >= INTERIOR_MIN && opt_median <= INTERIOR_MAX

  # Filter 2: MSE differentiation
  mse_diff <- ((max_mse - min_mse) / min_mse) * 100
  pass_mse_diff <- mse_diff >= MSE_DIFF_MIN && mse_diff <= MSE_DIFF_MAX

  # Filter 3a: Steepness (penalty at ±10 units) - average across datasets
  penalty_at_10_list <- var_data %>%
    group_by(dataset_id) %>%
    mutate(
      min_mse_ds = min(mse_true),
      optimal_ds = nbasis[which.min(mse_true)],
      dist_from_opt = abs(nbasis - optimal_ds)
    ) %>%
    filter(dist_from_opt == 10) %>%
    summarise(penalty = mean((mse_true / min_mse_ds - 1) * 100), .groups = "drop")

  if (nrow(penalty_at_10_list) > 0) {
    penalty_at_10 <- mean(penalty_at_10_list$penalty)
  } else {
    penalty_at_10 <- NA
  }
  pass_penalty <- !is.na(penalty_at_10) && penalty_at_10 >= PENALTY_AT_10_MIN

  # Filter 3b: Steepness (models within 5% of optimal) - average across datasets
  n_within_5pct_list <- var_data %>%
    group_by(dataset_id) %>%
    mutate(
      min_mse_ds = min(mse_true),
      pct_from_min = (mse_true / min_mse_ds - 1) * 100
    ) %>%
    summarise(n_within_5 = sum(pct_from_min <= 5), .groups = "drop")

  n_within_5pct <- mean(n_within_5pct_list$n_within_5)
  pass_within_5pct <- n_within_5pct <= MODELS_WITHIN_5PCT_MAX

  # Overall steepness (pass either 3a OR 3b)
  pass_steepness <- pass_penalty || pass_within_5pct

  # Pass all filters
  pass_all <- interior_min && pass_mse_diff && pass_steepness

  # Record results
  screening_results <- rbind(screening_results, data.frame(
    variable = var,
    label = var_label,
    type = var_type,
    n_datasets = n_datasets,
    opt_median = opt_median,
    opt_sd = opt_sd,
    pct_boundary = pct_boundary,
    mse_diff = mse_diff,
    penalty_at_10 = penalty_at_10,
    n_within_5pct = n_within_5pct,
    interior_min = interior_min,
    pass_mse_diff = pass_mse_diff,
    pass_penalty = pass_penalty,
    pass_within_5pct = pass_within_5pct,
    pass_steepness = pass_steepness,
    pass_all = pass_all,
    stringsAsFactors = FALSE
  ))
}

# Sort by pass_all, then steepness metrics
screening_results <- screening_results %>%
  arrange(desc(pass_all), desc(penalty_at_10), n_within_5pct)

# Summary
cat("=== PHASE 1 SCREENING RESULTS ===\n\n")
cat(sprintf("Tested: %d responses\n", nrow(screening_results)))

# Finalists
finalists <- screening_results %>% filter(pass_all)
cat(sprintf("PASS ALL FILTERS: %d responses\n\n", nrow(finalists)))

if (nrow(finalists) > 0) {
  cat("FINALISTS:\n")
  cat(paste(rep("-", 110), collapse = ""), "\n")
  cat(sprintf("%-25s %-12s %8s %10s %10s %12s %10s %10s %s\n",
              "Label", "Type", "Median", "SD", "Boundary%", "MSE Diff", "Penalty@10", "N≤5%", "Result"))
  cat(paste(rep("-", 110), collapse = ""), "\n")

  for (i in 1:nrow(finalists)) {
    f <- finalists[i, ]
    cat(sprintf("%-25s %-12s %8.0f %10.1f %9.0f%% %11.1f%% %9.1f%% %10.1f ✓ PASS\n",
                f$label, f$type, f$opt_median, f$opt_sd, f$pct_boundary,
                f$mse_diff, f$penalty_at_10, f$n_within_5pct))
  }
} else {
  cat("WARNING: No responses passed all filters!\n\n")
}

# Show all results
cat("\n\nALL RESULTS (sorted by steepness):\n")
cat(paste(rep("-", 120), collapse = ""), "\n")
cat(sprintf("%-25s %-12s %8s %10s %10s %12s %10s %10s %s\n",
            "Label", "Type", "Median", "SD", "Boundary%", "MSE Diff", "Penalty@10", "N≤5%", "Filters"))
cat(paste(rep("-", 120), collapse = ""), "\n")

for (i in 1:nrow(screening_results)) {
  r <- screening_results[i, ]
  filters_pass <- sum(c(r$interior_min, r$pass_mse_diff, r$pass_steepness))
  penalty_str <- ifelse(is.na(r$penalty_at_10), "NA", sprintf("%.1f%%", r$penalty_at_10))

  cat(sprintf("%-25s %-12s %8.0f %10.1f %9.0f%% %11.1f%% %10s %10.1f %d/3 %s\n",
              r$label, r$type, r$opt_median, r$opt_sd, r$pct_boundary,
              r$mse_diff, penalty_str, r$n_within_5pct, filters_pass,
              ifelse(r$pass_all, "✓", "✗")))
}

# Save results
write.csv(screening_results,
          file.path(results_dir_base, "screening_results.csv"),
          row.names = FALSE)

cat("\n\nResults saved to:", file.path(results_dir_base, "screening_results.csv"), "\n")

if (nrow(finalists) > 0) {
  write.csv(finalists,
            file.path(results_dir_base, "finalists.csv"),
            row.names = FALSE)
  cat("Finalists saved to:", file.path(results_dir_base, "finalists.csv"), "\n")

  cat("\n=== NEXT STEPS ===\n")
  cat("Phase 2: Run multi-dataset consistency test for finalists\n")
  cat("  - 10 independent datasets per finalist\n")
  cat("  - Check SD(optimal_nbasis) < 15\n")
  cat("  - Check 0% boundary selections\n")
} else {
  cat("\n=== RECOMMENDATIONS ===\n")
  cat("No responses passed strict criteria. Consider:\n")
  cat("  1. Try moderate thresholds (penalty >= 3%, N≤5%% <= 4)\n")
  cat("  2. Test different model configurations (add covariates, adjust priors)\n")
  cat("  3. Try different sampling schemes (equal_30, proportional)\n")
  cat("  4. Test responses from other states\n")
}

cat("\n")
