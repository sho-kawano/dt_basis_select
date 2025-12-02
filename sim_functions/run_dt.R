library(parallel)
library(doParallel)
library(Matrix)
library(tidyverse)

# k is the number of 'folds' in data thinning
# eps is the thinning parameter for k=1 (proportion of signal reserved for training)
# n_reps is the number of repetitions (independent data thinning splits)
# model_config is the configuration list from configs/model_configs.R
run_dt <- function(comp_no, k, results_dir, model_config, eps = 0.5, n_reps = 1) {
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

  # Convert to model matrix
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

  # Construct DT folder path based on k and eps
  if (k == 1) {
    dt_folder <- file.path(comp_folder, "dt_1fold", sprintf("eps_%.2f", eps))
  } else {
    dt_folder <- file.path(comp_folder, sprintf("dt_%dfold", k))
  }
  # Create directory (recursive=TRUE creates eps subdirs for k=1)
  dir.create(dt_folder, showWarnings = FALSE, recursive = TRUE)

  n <- length(z)  # Number of areas (PUMAs)
  nbasis_values <- model_config$nbasis_values

  # ----- Loop over repetitions -----
  for (rep in 1:n_reps) {
    # Set seed for this repetition to ensure different splits
    set.seed(comp_no * 1000 + rep)

    # ----- Create train/test sets via data thinning -----
    if (k == 1) {
      # Single-fold: define z_train first (see dt paper)
      # eps parameter is passed as function argument (default 0.5)
      # eps = training proportion, (1-eps) = test proportion
      z_train <- rnorm(n = length(z), mean = eps * z, sd = sqrt(eps * (1 - eps) * d))
      z_test <- z - z_train
      # Correlation between test & validation set in each fold (for diagnostics)
      fold_cors <- cor(z_train, z_test) %>% round(4)
      # Convert to matrices
      z_train <- matrix(z_train, nrow = length(z_train), ncol = 1)
      z_test <- matrix(z_test, nrow = length(z_test), ncol = 1)
    } else {
      # Multi-fold: define z_test first (see dt paper)
      e0 <- 1 / k # test proportion (0.2 for k=5)
      eps_mf <- matrix(rep(e0, k), ncol = 1) # thinning parameter
      z_test <- matrix(nrow = n, ncol = k)
      for (i in 1:n) {
        Sig <- d[i] * (diag(e0, nrow = k) - eps_mf %*% t(eps_mf))
        z_test[i, ] <- mvtnorm::rmvnorm(n = 1, mean = eps_mf * z[i], sigma = Sig) %>%
          as.numeric()
      }
      # Create training set
      z_train <- z - z_test
      # Correlation between test & validation set in each fold (for diagnostics)
      fold_cors <- sapply(
        1:ncol(z_train),
        function(u) cor(z_train[, u], z_test[, u]) %>% round(4)
      )
    }

    # save the results from data thinning
    # Note: For k>1, eps is deterministic (1/k), no need to save it
    split_file <- sprintf("dt_split_%dfold_rep%d.RDA", k, rep)
    if (k == 1) {
      save(z_train, z_test, fold_cors, eps,
        file = file.path(dt_folder, split_file)
      )
    } else {
      save(z_train, z_test, fold_cors,
        file = file.path(dt_folder, split_file)
      )
    }

    # ----- Loop over folds -----
    # Run all folds sequentially (no nested parallelization)
    # Outer loop (comparisons) is already parallelized, which is more efficient
    for (fold in 1:k) {
      tryCatch(
        {
          # Loop over nbasis values
          all_results <- vector("list", length(nbasis_values))

          for (i in seq_along(nbasis_values)) {
            nb <- nbasis_values[i]

            # Set unique seed for each model fit: varies by comparison, rep, fold, and nbasis
            # This ensures independent MCMC chains across all dimensions
            set.seed(comp_no * 10000 + rep * 100 + fold * 10 + nb)

            # Variance for training: eps*d for k=1, (1-e0)*d for k>1
            train_d <- if (k == 1) eps * d else (1 - e0) * d

            # Fit model
            fit <- spatial_basis_fh(
              X = X,
              y = z_train[, fold],
              d.var = train_d,
              nbasis = nb,
              eig_moran = eig_moran,
              A = A,
              spatial_type = model_config$spatial_type,
              ndesired = model_config$ndesired,
              nburn = model_config$nburn,
              nthin = model_config$nthin,
              hyp = model_config$hyp,
              verbose = FALSE
            )

            # Save posterior mean and variance for theta
            all_results[[i]] <- data.frame(
              domain = puma_ids,
              mean = colMeans(fit$theta),
              var = apply(fit$theta, 2, var),
              method = sprintf("nbasis_%03d", nb),
              nbasis = nb
            )
          }
          all_results <- all_results %>% bind_rows()
          row.names(all_results) <- NULL

          # File naming: always use rep naming
          if (k == 1) {
            file_name <- sprintf("%dfold_rep%d.RDS", k, rep)
          } else {
            file_name <- sprintf("%dfold_rep%d_fold%d.RDS", k, rep, fold)
          }
          saveRDS(all_results, file = file.path(dt_folder, file_name))
        },
        error = function(e) {
          error_msg <- sprintf("Error in comp %d rep %d fold %d: %s", comp_no, rep, fold, e$message)
          message(error_msg)
          error_file <- file.path(
            dt_folder,
            sprintf("dt_rep%d_fold%d_ERROR.RDS", rep, fold)
          )
          saveRDS(list(comp_no = comp_no, rep = rep, fold = fold, error = e$message),
            file = error_file
          )
          return(NULL)
        }
      )
    }
  } # End of repetitions loop
}
