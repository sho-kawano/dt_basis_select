library(SUMMER)
library(parallel)
library(doParallel)
library(Matrix)
library(tidyverse)
#library(rstan)

# simply fit each of the models on 'z' the direct estimates
# saves the full sample (to use for CRPS, WAIC, DIC, etc.)
zfit <- function(comp_no, s, all_covs, all_data, results_dir){
  # ----- loading data /samplers -----
  source("samplers/fh_fit.R")

  # specify folder/comparison
  comp_folder = file.path(getwd(), results_dir, sprintf("comparison_%03d", comp_no))

  # load direct estimates
  z = readRDS(file.path(comp_folder, "z.RDS")) %>% as.numeric()

  # Only fit models with 1, 4, 7 covariates
  cov_counts <- c(1, 4, 7)
  results = vector("list", length(cov_counts))
  for(i in seq_along(cov_counts)){
    ncovs <- cov_counts[i]
    X = model.matrix(~., all_data[, all_covs[1:ncovs], drop=F])
    # Set unique seed for each model fit: varies by comparison and model
    set.seed(comp_no * 100 + ncovs)
    # fit the chain
    fh_chain = fh_fit(X, z, s, ndesired=1000, nburn=8000, nthin=1, verbose=F)

    results[[i]] = fh_chain$theta
  }
  # save results
  file_name = file.path(comp_folder, "fit_on_z", "chains.RDS")
  saveRDS(results, file=file_name)
}

