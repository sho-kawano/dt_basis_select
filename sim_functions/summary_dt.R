library(dplyr)
library(foreach)

#-----------------------------------------------------------------------------
# Loss Function Factories
#-----------------------------------------------------------------------------
# These functions return loss functions configured with data thinning parameters
# NOTE: Formulas differ between k=1 and k>1 due to different data thinning setups

make_plugin_nll <- function(z_test_rep, eps_rep, d_var, n_folds) {
  function(theta.hat, j) {
    if (n_folds == 1) {
      # For k=1: z_train~N(eps*theta,...), z_test~N((1-eps)*theta,...)
      # theta.hat estimates eps*theta
      mean_pred <- ((1 - eps_rep) / eps_rep) * theta.hat
      sd_pred <- sqrt((1 - eps_rep) * d_var)
    } else {
      # For k>1: z_train~N((1-eps)*theta,...), z_test~N(eps*theta,...)
      # theta.hat estimates (1-eps)*theta
      mean_pred <- (eps_rep / (1 - eps_rep)) * theta.hat
      sd_pred <- sqrt(eps_rep * d_var)
    }
    -1 * sum(dnorm(z_test_rep[, j], mean_pred, sd_pred, log = TRUE))
  }
}

make_predictive_nll <- function(z_test_rep, eps_rep, d_var, n_folds) {
  function(theta.hat, theta.var, j) {
    if (n_folds == 1) {
      # For k=1: model estimates eps*theta, need to transform to theta
      mean_pred <- ((1 - eps_rep) / eps_rep) * theta.hat
      sd_pred <- sqrt((1 - eps_rep) * d_var + ((1 - eps_rep) / eps_rep)^2 * theta.var)
    } else {
      # For k>1: model estimates (1-eps)*theta, need to transform to theta
      mean_pred <- (eps_rep / (1 - eps_rep)) * theta.hat
      sd_pred <- sqrt(eps_rep * d_var + (eps_rep / (1 - eps_rep))^2 * theta.var)
    }
    -1 * sum(dnorm(z_test_rep[, j], mean_pred, sd_pred, log = TRUE))
  }
}

make_mse <- function(z_test_rep, eps_rep, d_var, n_folds) {
  function(theta.hat, j) {
    if (n_folds == 1) {
      # For k=1: theta.hat estimates eps*theta, z_test represents (1-eps)*theta
      de_thinned_hat <- theta.hat / eps_rep
      de_thinned_test <- z_test_rep[, j] / (1 - eps_rep)
      bias_correction <- mean(d_var / (1 - eps_rep))
    } else {
      # For k>1: theta.hat estimates (1-eps)*theta, z_test represents eps*theta
      de_thinned_hat <- theta.hat / (1 - eps_rep)
      de_thinned_test <- z_test_rep[, j] / eps_rep
      bias_correction <- mean(d_var / eps_rep)
    }
    naive_mse <- mean((de_thinned_hat - de_thinned_test)^2)
    return(naive_mse - bias_correction)
  }
}

make_naive_mse <- function(z_test_rep, eps_rep, d_var, n_folds) {
  function(theta.hat, j) {
    if (n_folds == 1) {
      de_thinned_hat <- theta.hat / eps_rep
      de_thinned_test <- z_test_rep[, j] / (1 - eps_rep)
    } else {
      de_thinned_hat <- theta.hat / (1 - eps_rep)
      de_thinned_test <- z_test_rep[, j] / eps_rep
    }
    naive_mse <- mean((de_thinned_hat - de_thinned_test)^2)
    return(naive_mse)
  }
}

#-----------------------------------------------------------------------------
# Main Summary Function
#-----------------------------------------------------------------------------

summary_dt <- function(comp_no, n_folds, loss_function, results_dir, n_reps_to_use, eps = NULL) {
  #-----------------------------------------------------------------------------
  # folder for comparison
  comp_folder <- file.path(getwd(), results_dir, sprintf("comparison_%03d", comp_no))

  # Load comparison-specific variance
  d_data <- readRDS(file.path(comp_folder, "d.RDS"))
  d_var <- d_data$values

  # Construct DT folder path based on n_folds and eps
  if (n_folds == 1) {
    if (is.null(eps)) stop("eps parameter required for n_folds=1")
    dt_folder <- file.path(comp_folder, "dt_1fold", sprintf("eps_%.2f", eps))
  } else {
    dt_folder <- file.path(comp_folder, sprintf("dt_%dfold", n_folds))
  }

  #-----------------------------------------------------------------------------
  # Loop over repetitions (outer loop)
  result_all_reps <- foreach(rep = 1:n_reps_to_use, .combine = rbind) %do% {
    # Load split file for this repetition
    split_file <- sprintf("dt_split_%dfold_rep%d.RDA", n_folds, rep)

    # Load into a new environment to avoid scope issues
    split_env <- new.env()
    load(file.path(dt_folder, split_file), envir = split_env)

    # Extract variables from the loaded environment
    z_train_rep <- split_env$z_train
    z_test_rep <- split_env$z_test
    fold_cors_rep <- split_env$fold_cors

    # Determine eps value
    # For k=1, eps is saved in the file (user-specified thinning parameter)
    # For k>1, eps is deterministic: 1/k (not saved in file)
    if (n_folds == 1) {
      eps_scalar <- split_env$eps
    } else {
      eps_scalar <- 1 / n_folds
    }

    #-----------------------------------------------------------------------------
    # iterates over folds and combines via rbind
    result_one_rep <- foreach(j = 1:n_folds, .combine = rbind) %do% {
      # load the summary for a given fold
      if (n_folds == 1) {
        result_file <- sprintf("%dfold_rep%d.RDS", n_folds, rep)
      } else {
        result_file <- sprintf("%dfold_rep%d_fold%d.RDS", n_folds, rep, j)
      }
      results <- readRDS(file.path(dt_folder, result_file))

      # create the results table for a given fold
      res <- results %>%
        select(domain, mean, var, method, nbasis) %>%
        rename(fips = domain, estim = mean, estim_var = var)

      # Create the appropriate loss function using factory
      # Use eps_scalar (not eps_rep) to ensure scalar value
      # Pass n_folds to use correct formula for k=1 vs k>1
      loss_fn <- switch(loss_function,
        "plugin_NLL" = make_plugin_nll(z_test_rep, eps_scalar, d_var, n_folds),
        "predictive_NLL" = make_predictive_nll(z_test_rep, eps_scalar, d_var, n_folds),
        "MSE" = make_mse(z_test_rep, eps_scalar, d_var, n_folds),
        "naive_MSE" = make_naive_mse(z_test_rep, eps_scalar, d_var, n_folds),
        stop("Unknown loss_function: ", loss_function)
      )

      # calculate the loss for each model
      res %>%
        group_by(method) %>%
        reframe(metric = if (loss_function == "predictive_NLL") {
          loss_fn(estim, estim_var, j)
        } else {
          loss_fn(estim, j)
        }) %>%
        rename(model = method) %>%
        mutate(
          comp_no = comp_no, fold_no = j, rep_no = rep,
          method = paste0("DT ", n_folds, " fold - ", loss_function),
          loss_function = loss_function
        ) %>%
        left_join(res %>% select(fips, model = method, nbasis) %>% distinct(model, nbasis), by = "model") %>%
        relocate(rep_no, fold_no, comp_no, method, model, metric, loss_function) %>%
        arrange(metric)
    } # End of folds loop

    # Average across folds for this repetition
    result_one_rep %>%
      group_by(comp_no, method, model, loss_function, rep_no, nbasis) %>%
      reframe(metric = mean(metric))
  } # End of repetitions loop

  #-----------------------------------------------------------------------------
  # Average across repetitions and return final result
  result_all_reps %>%
    group_by(comp_no, method, model, loss_function, nbasis) %>%
    reframe(metric = mean(metric)) %>%
    arrange(metric) %>%
    relocate(loss_function, .after = last_col())
}
