library(SUMMER)
library(parallel)
library(doParallel)
library(Matrix)
library(tidyverse)

# Fits models on 'z' (direct estimates) using full observed data
# Saves full posterior samples for calculating WAIC, DIC, oracle MSE, etc.
#
# @param comp_no comparison number
# @param results_dir results directory
# @param model_config configuration list (from configs/model_configs.R)
full_data_fit <- function(comp_no, results_dir, model_config) {
  # ----- Load model function -----
  source("models/spatial_basis_fh.R")
  source("models/mcmc_helper.R")

  # Specify folder/comparison
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

  # Convert data frame to model matrix
  # If X_df has no columns (intercept-only), create intercept matrix
  if (ncol(X_df) == 0) {
    X <- matrix(1, nrow = length(z), ncol = 1)
    colnames(X) <- "(Intercept)"
  } else {
    X <- model.matrix(~., X_df)
  }

  # ----- Load adjacency matrix and compute spatial basis -----
  if (is.null(model_config$adjacency_file)) {
    stop("model_config$adjacency_file is required! Specify the adjacency matrix file.")
  }
  load(model_config$adjacency_file)  # Loads A (adjacency matrix)

  # Compute Moran eigenvectors once (for efficiency)
  eig_moran <- compute_moran_eigen(A, X)

  # ----- Fit models across nbasis values -----
  nbasis_values <- model_config$nbasis_values
  all_results <- vector("list", length(nbasis_values))

  for (i in seq_along(nbasis_values)) {
    nb <- nbasis_values[i]

    # Set seed for reproducibility
    set.seed(comp_no * 1000 + nb)

    cat(sprintf("  Fitting nbasis=%d...\n", nb))

    # Fit model
    fit <- spatial_basis_fh(
      X = X,
      y = z,
      d.var = d,
      nbasis = nb,
      eig_moran = eig_moran,
      spatial_type = model_config$spatial_type,
      ndesired = model_config$ndesired,
      nburn = model_config$nburn,
      nthin = model_config$nthin,
      hyp = model_config$hyp,
      verbose = FALSE
    )

    # Save results
    all_results[[i]] <- list(
      model = sprintf("nbasis_%03d", nb),
      nbasis = nb,
      chain = fit,
      puma_ids = puma_ids,
      spatial_type = model_config$spatial_type
    )
  }

  # Save all results
  file_name <- file.path(comp_folder, "fit_on_z", "chains.RDS")
  saveRDS(all_results, file = file_name)

  cat(sprintf("  Saved results to %s\n", file_name))
}
