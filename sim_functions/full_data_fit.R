library(SUMMER)
library(parallel)
library(doParallel)
library(Matrix)
library(tidyverse)
# library(rstan)

# Fits models on 'z' (direct estimates) using full observed data
# Saves full posterior samples for calculating WAIC, DIC, oracle MSE, etc.
full_data_fit <- function(comp_no, results_dir) {
  # ----- loading data /samplers -----
  source("models/fh_fit.R")

  # specify folder/comparison
  comp_folder <- file.path(getwd(), results_dir, sprintf("comparison_%03d", comp_no))

  # Load comparison-specific data (X, z, d)
  X_data <- readRDS(file.path(comp_folder, "X.RDS"))
  z_data <- readRDS(file.path(comp_folder, "z.RDS"))
  d_data <- readRDS(file.path(comp_folder, "d.RDS"))

  # Extract values
  puma_ids <- X_data$puma
  X_df <- X_data$X  # Data frame with predictor columns
  z <- z_data$values
  d <- d_data$values
  predictor_names <- colnames(X_df)

  # Fit models with 1, 4, 7 covariates (matching main analysis)
  cov_counts <- c(1, 4, 7)
  all_results <- vector("list", length(cov_counts))

  for (i in seq_along(cov_counts)) {
    ncovs <- cov_counts[i]
    # Convert data frame to model matrix (adds intercept)
    X <- model.matrix(~., X_df[, predictor_names[1:ncovs], drop = F])

    # Set seed for reproducibility
    set.seed(comp_no * 100 + ncovs)

    # Fit model (use minimal iterations for testing - increase for production)
    fh_chain <- fh_fit(X, z, d, ndesired = 10, nburn = 10, nthin = 1, verbose = F)

    # Save results
    all_results[[i]] <- list(
      model = paste0(ncovs, "_cov"),
      chain = fh_chain,
      puma_ids = puma_ids
    )
  }

  # Save all results
  file_name <- file.path(comp_folder, "fit_on_z", "chains.RDS")
  saveRDS(all_results, file = file_name)
}
