library(tidyverse)

# function to caculate DIC for one posterior sample 
calc_DIC <- function(theta_samples, z, d){
  # theta_samples: MCMC samples, each row is a sample
  # z: observed direct estimates
  # d: known variances (length n)
  # log-likelihood function
  log_lik <- function(theta, z, d) {
    sum(dnorm(z, mean = theta, sd = sqrt(d), log = TRUE))
  }
  # Compute deviance for each posterior sample
  deviances <- apply(theta_samples, 1, function(theta) {
    -2 * log_lik(theta, z, d)
  })
  # the components
  average_deviances = mean(deviances)
  deviance_at_post_mean = -2 * log_lik(colMeans(theta_samples), z, d)
  2 * average_deviances - deviance_at_post_mean
}

# function to caculate WAIC for one posterior sample 
calc_WAIC <- function(theta_samples, z, d){
  S <- nrow(theta_samples)
  n <- length(z)
  
  # Log-likelihood matrix: S rows (samples) × n columns (data points)
  log_lik_matrix <- matrix(NA, nrow = S, ncol = n)
  
  for (s in 1:S) {
    log_lik_matrix[s, ] <- dnorm(z, mean = theta_samples[s, ],
                                 sd = sqrt(d), log = TRUE)
  }
  
  # lppd: log pointwise predictive density
  lppd <- sum(log(colMeans(exp(log_lik_matrix))))
  
  # p_waic: effective number of parameters
  p_waic <- sum(apply(log_lik_matrix, 2, var))
  
  WAIC <- -2 * (lppd - p_waic)
  return(WAIC)
}

# wrapper function to help create the resulting table
# calculates true MSE, WAIC, DIC and puts them into one table
summary_oracle <- function(comp_no, results_dir){
  comp_folder = file.path(getwd(), results_dir, sprintf("comparison_%03d", comp_no))
  all_results = readRDS(file.path(comp_folder, "fit_on_z", "chains.RDS"))

  # Load comparison-specific estimates and variances (new list format)
  z_data = readRDS(file.path(comp_folder, "z.RDS"))
  d_data = readRDS(file.path(comp_folder, "d.RDS"))

  z = z_data$values
  d = d_data$values
  puma_ids = z_data$puma

  # Load true population values from results directory
  base_folder = file.path(getwd(), results_dir)
  theta_true_data = readRDS(file.path(base_folder, "theta_true.RDS"))

  # Ensure alignment with z_data (same PUMAs, same order)
  theta_true_df <- data.frame(
    PUMA = theta_true_data$puma,
    truth = theta_true_data$values
  ) %>%
    filter(PUMA %in% puma_ids) %>%
    arrange(PUMA)

  truth <- theta_true_df$truth

  # Extract chains and model info from new structure
  post_samples <- lapply(all_results, function(x) x$chain$theta)
  model_names <- sapply(all_results, function(x) x$model)
  nbasis_values <- sapply(all_results, function(x) x$nbasis)

  # Compute posterior means for each model
  estimates = sapply(post_samples, function(x){x %>% colMeans()})
  estimates = cbind(estimates, z)

  # Calculate MSE against truth (not z!)
  mse_df = data.frame(comp_no=comp_no,
                      model=c(model_names, "D.Est"),
                      nbasis = c(nbasis_values, NA),
                      mse_true=colMeans((estimates-truth)^2))
  waic_df = data.frame(comp_no=comp_no,
                      model=model_names,
                      nbasis = nbasis_values,
                      WAIC = sapply(post_samples, calc_WAIC, z=z, d=d))
  dic_df = data.frame(comp_no=comp_no,
                      model=model_names,
                      nbasis = nbasis_values,
                      DIC = sapply(post_samples, calc_DIC, z=z, d=d))

  summary_df = mse_df %>%
    left_join(waic_df, by = join_by(comp_no, model, nbasis)) %>%
    left_join(dic_df, by = join_by(comp_no, model, nbasis)) %>%
    arrange(mse_true)

  return(summary_df)
}

# Alias for backward compatibility
summary_zfit <- summary_oracle


