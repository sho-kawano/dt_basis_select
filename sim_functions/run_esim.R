library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)
#library(rstan)

#(comp_no, k, s, cov_names)
# runs empirical simulation study for comparison # `comp_no`
run_esim <- function(comp_no, s, all_covs, n_cores, all_data, results_dir, n_iters=100){

  # ---- set up ----
  # folder where the comparison is stored
  comp_folder = file.path(getwd(), results_dir, sprintf("comparison_%03d", comp_no))

  # load sampler
  source("models/fh_fit.R")

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
      # Only fit models with 1, 4, 7 covariates
      cov_counts <- c(1, 4, 7)
      # each of the models + direct estimate - preallocate
      results = vector("list", length(cov_counts)+1)

      # load synthetic direct estimate w
      sub_folder = file.path(comp_folder, "emp_sim", sprintf("%03d", iter))
      w = readRDS(file.path(sub_folder, "w.RDS")) %>% as.numeric()

      # save the direct estimate to results
      results[[1]] = data.frame(domain=all_data$fips, mean=w,
                            median=NA, var=s, lower=NA, upper=NA, method="Direct")

      # run each of the models - save results
      for(i in seq_along(cov_counts)){
        ncovs <- cov_counts[i]
        X = model.matrix(~., all_data[, all_covs[1:ncovs], drop=F])

        # Set unique seed for each model fit: varies by comparison, iteration, and model
        set.seed(comp_no * 10000 + iter * 10 + ncovs)
        fh_chain = fh_fit(X, w, s, ndesired=1000, nburn=10000, nthin=1, verbose=F)

        results[[i+1]] = fh_chain %>%
          mcmc_summary("theta", paste0(ncovs, "_cov_model")) %>%
          mutate(domain=all_data$fips) %>%
          relocate(domain)
      }

      # ---- save results & wrap up  ----
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
        # Only fit models with 1, 4, 7 covariates
        cov_counts <- c(1, 4, 7)
        # each of the models + direct estimate - preallocate
        results = vector("list", length(cov_counts)+1)

        # load synthetic direct estimate w
        sub_folder = file.path(comp_folder, "emp_sim", sprintf("%03d", iter))
        w = readRDS(file.path(sub_folder, "w.RDS")) %>% as.numeric()

        # save the direct estimate to results
        results[[1]] = data.frame(domain=all_data$fips, mean=w,
                              median=NA, var=s, lower=NA, upper=NA, method="Direct")

        # run each of the models - save results
        for(i in seq_along(cov_counts)){
          ncovs <- cov_counts[i]
          X = model.matrix(~., all_data[, all_covs[1:ncovs], drop=F])

          # Set unique seed for each model fit: varies by comparison, iteration, and model
          set.seed(comp_no * 10000 + iter * 10 + ncovs)
          fh_chain = fh_fit(X, w, s, ndesired=1000, nburn=10000, nthin=1, verbose=F)

          results[[i+1]] = fh_chain %>%
            mcmc_summary("theta", paste0(ncovs, "_cov_model")) %>%
            mutate(domain=all_data$fips) %>%
            relocate(domain)
        }

        # ---- save results & wrap up  ----
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
