library(parallel)
library(doParallel)
library(Matrix)
library(tidyverse)

# k is the number of 'folds' in data thinning
# eps is the thinning parameter for k=1 (proportion of signal reserved for training)
# n_reps is the number of repetitions (independent data thinning splits)
run_dt <- function(comp_no, k, s, all_covs, all_data, results_dir, eps = 0.5, n_reps = 1) {
  # ----- loading data /samplers -----
  source("models/fh_fit.R")

  # specify folder/comparison
  comp_folder <- file.path(getwd(), results_dir, sprintf("comparison_%03d", comp_no))

  # Construct DT folder path based on k and eps
  if (k == 1) {
    dt_folder <- file.path(comp_folder, "dt_1fold", sprintf("eps_%.2f", eps))
  } else {
    dt_folder <- file.path(comp_folder, "dt_5fold")
  }
  # Create directory (recursive=TRUE creates eps subdirs for k=1)
  dir.create(dt_folder, showWarnings = FALSE, recursive = TRUE)

  n <- nrow(all_data)
  z <- readRDS(file.path(comp_folder, "z.RDS")) %>% as.numeric()

  # ----- Loop over repetitions -----
  for (rep in 1:n_reps) {
    # Set seed for this repetition to ensure different splits
    set.seed(comp_no * 1000 + rep)

    # ----- create train/test sets via data thinning -----
    if (k == 1) {
      # in single-fold it is easier to define z_train first (see dt paper)
      # eps parameter is passed as function argument (default 0.5)
      # eps = training proportion, (1-eps) = test proportion
      z_train <- rnorm(n = length(z), mean = eps * z, sd = sqrt(eps * (1 - eps) * s))
      z_test <- z - z_train
      # correlation between test & validation set in each fold (for diagnostics)
      fold_cors <- cor(z_train, z_test) %>% round(4)
      # convert to matrices
      z_train <- matrix(z_train, nrow = length(z_train), ncol = 1)
      z_test <- matrix(z_test, nrow = length(z_test), ncol = 1)
    } else {
      # in multi-fold it is easier to define z_test first (see dt paper)
      e0 <- 1 / k # test proportion (0.2 for k=5)
      eps <- matrix(rep(e0, k), ncol = 1) # thinning parameter
      z_test <- matrix(nrow = n, ncol = k)
      check <- c()
      for (i in 1:n) {
        Sig <- s[i] * (diag(e0, nrow = k) - eps %*% t(eps))
        z_test[i, ] <- mvtnorm::rmvnorm(n = 1, mean = eps * z[i], sigma = Sig) %>%
          as.numeric()
      }
      # create training set
      z_train <- z - z_test
      # correlation between test & validation set in each fold (for diagnostics)
      fold_cors <- sapply(
        1:ncol(z_train),
        function(u) cor(z_train[, u], z_test[, u]) %>% round(4)
      )
    }

    # save the results from data thinning (including eps for summary function)
    split_file <- sprintf("dt_split_%dfold_rep%d.RDA", k, rep)
    save(z_train, z_test, fold_cors, eps,
      file = file.path(dt_folder, split_file)
    )

    # ----- Loop over folds -----
    # Run all folds sequentially (no nested parallelization)
    # Outer loop (comparisons) is already parallelized, which is more efficient
    for (fold in 1:k) {
      tryCatch(
        {
          # Only fit models with 1, 4, 7 covariates (matching main.tex)
          cov_counts <- c(1, 4, 7)
          all_results <- vector("list", length(cov_counts))

          for (i in seq_along(cov_counts)) {
            ncovs <- cov_counts[i]
            X <- model.matrix(~., all_data[, all_covs[1:ncovs], drop = F])

            # Set unique seed for each model fit: varies by comparison, rep, fold, and model
            # This ensures independent MCMC chains across all dimensions
            set.seed(comp_no * 10000 + rep * 100 + fold * 10 + ncovs)
            # Variance for training: eps*s for k=1, (1-e0)*s for k>1
            train_s <- if (k == 1) eps * s else (1 - e0) * s
            fh_chain <- fh_fit(X, z_train[, fold], train_s,
              ndesired = 1000, nburn = 10000, nthin = 1, verbose = F
            )

            # Save posterior mean and variance for theta
            all_results[[i]] <- data.frame(
              domain = all_data$fips,
              mean = colMeans(fh_chain$theta),
              var = apply(fh_chain$theta, 2, var),
              method = paste0(ncovs, "_cov_model")
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
