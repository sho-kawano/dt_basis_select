library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

# runs empirical simulation study for comparison # `comp_no`
# model_config is the configuration list from configs/model_configs.R
run_esim <- function(comp_no, n_cores, results_dir, model_config, n_iters=100){

  # ---- Set up ----
  # Folder where the comparison is stored
  comp_folder = file.path(getwd(), results_dir, sprintf("comparison_%03d", comp_no))

  # Load comparison-specific data (X, d)
  X_data <- readRDS(file.path(comp_folder, "X.RDS"))
  d_data <- readRDS(file.path(comp_folder, "d.RDS"))

  # Extract values
  puma_ids <- X_data$puma
  X_df <- X_data$X  # Data frame with predictor columns
  d_var <- d_data$values

  # Convert to model matrix
  # If X_df has no columns (intercept-only), create intercept matrix
  if (ncol(X_df) == 0) {
    X <- matrix(1, nrow = length(d_var), ncol = 1)
    colnames(X) <- "(Intercept)"
  } else {
    X <- model.matrix(~., X_df)
  }

  # Load model functions
  source("models/spatial_basis_fh.R")
  source("models/mcmc_helper.R")

  # Load adjacency matrix and compute spatial basis
  if (is.null(model_config$adjacency_file)) {
    stop("model_config$adjacency_file is required! Specify the adjacency matrix file.")
  }
  load(model_config$adjacency_file)  # Loads A (adjacency matrix)
  eig_moran <- compute_moran_eigen(A, X)

  # Get nbasis values from config
  nbasis_values <- model_config$nbasis_values

  # ---- decide parallelization strategy based on n_cores ----
  if(n_cores > 1) {
    # Parallel mode: set up cluster
    progress_txt = file.path(comp_folder, "esim_outfile.txt")
    if (file.exists(progress_txt)) {
      file.remove(progress_txt)
    }
    cl = makeForkCluster(n_cores, outfile=progress_txt)
    registerDoParallel(cl)
    if(length(cl)!=n_cores){warning("No. clusters less than desired.")}

    # Run iterations in parallel
    foreach(iter=1:n_iters) %dopar% {
    tryCatch({
      set.seed(iter*comp_no)
      # Each of the models + direct estimate - preallocate
      results = vector("list", length(nbasis_values)+1)

      # Load synthetic direct estimate w
      sub_folder = file.path(comp_folder, "emp_sim", sprintf("%03d", iter))
      w = readRDS(file.path(sub_folder, "w.RDS")) %>% as.numeric()

      # Save the direct estimate to results
      results[[1]] = data.frame(domain=puma_ids, mean=w,
                            median=NA, var=d_var, lower=NA, upper=NA, method="Direct")

      # Run each of the models - save results
      for(i in seq_along(nbasis_values)){
        nb <- nbasis_values[i]

        # Set unique seed for each model fit: varies by comparison, iteration, and nbasis
        set.seed(comp_no * 10000 + iter * 10 + nb)

        # Fit model
        fit <- spatial_basis_fh(
          X = X,
          y = w,
          d.var = d_var,
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

        results[[i+1]] = fit %>%
          mcmc_summary("theta", sprintf("nbasis_%03d", nb)) %>%
          mutate(domain=puma_ids, nbasis=nb) %>%
          relocate(domain)
      }

      # ---- Save results & wrap up  ----
      results = results %>% bind_rows()
      row.names(results) = NULL
      saveRDS(results, file=file.path(sub_folder, "results.RDS"))
      print("__________")
      print(paste0("Iteration #", iter, " is finished."))
      print(Sys.time())
    }, error = function(e) {
      error_msg <- sprintf("Error in comp %d iter %d: %s", comp_no, iter, e$message)
      message(error_msg)
      error_file <- file.path(comp_folder, sprintf("esim_iter%03d_ERROR.RDS", iter))
      saveRDS(list(comp_no = comp_no, iter = iter, error = e$message), file = error_file)
      return(NULL)
    })
    }
    stopCluster(cl)
  } else {
    # Sequential mode (n_cores=1): run iterations in a regular for loop
    for(iter in 1:n_iters) {
      tryCatch({
        set.seed(iter*comp_no)
        # Each of the models + direct estimate - preallocate
        results = vector("list", length(nbasis_values)+1)

        # Load synthetic direct estimate w
        sub_folder = file.path(comp_folder, "emp_sim", sprintf("%03d", iter))
        w = readRDS(file.path(sub_folder, "w.RDS")) %>% as.numeric()

        # Save the direct estimate to results
        results[[1]] = data.frame(domain=puma_ids, mean=w,
                              median=NA, var=d_var, lower=NA, upper=NA, method="Direct")

        # Run each of the models - save results
        for(i in seq_along(nbasis_values)){
          nb <- nbasis_values[i]

          # Set unique seed for each model fit: varies by comparison, iteration, and nbasis
          set.seed(comp_no * 10000 + iter * 10 + nb)

          # Fit model
          fit <- spatial_basis_fh(
            X = X,
            y = w,
            d.var = d_var,
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

          results[[i+1]] = fit %>%
            mcmc_summary("theta", sprintf("nbasis_%03d", nb)) %>%
            mutate(domain=puma_ids, nbasis=nb) %>%
            relocate(domain)
        }

        # ---- Save results & wrap up  ----
        results = results %>% bind_rows()
        row.names(results) = NULL
        saveRDS(results, file=file.path(sub_folder, "results.RDS"))
      }, error = function(e) {
        error_msg <- sprintf("Error in comp %d iter %d: %s", comp_no, iter, e$message)
        message(error_msg)
        error_file <- file.path(comp_folder, sprintf("esim_iter%03d_ERROR.RDS", iter))
        saveRDS(list(comp_no = comp_no, iter = iter, error = e$message), file = error_file)
      })
    }
  }
}
