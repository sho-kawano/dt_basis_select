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
summary_zfit <- function(comp_no, truth, d, all_covs, results_dir){
  comp_folder = file.path(getwd(), results_dir, sprintf("comparison_%03d", comp_no))
  post_samples = readRDS(file.path(comp_folder, "fit_on_z", "chains.RDS"))
  z = readRDS(file.path(comp_folder, "z.RDS"))

  # Only 3 models: 1, 4, 7 covariates
  cov_counts <- c(1, 4, 7)

  # workaround for adding in direct estimate
  estimates = sapply(post_samples, function(x){x %>% colMeans()})
  estimates = cbind(estimates, z %>% t())

  mse_df = data.frame(comp_no=comp_no,
                      model=c(paste0(cov_counts, "_cov_model"), "D.Est"),
                      no_covs = c(cov_counts, NA),
                            mse_true=colMeans((estimates-truth)^2))
  waic_df = data.frame(comp_no=comp_no,
                      model=c(paste0(cov_counts, "_cov_model")),
                      no_covs = cov_counts,
                      WAIC = sapply(post_samples, calc_WAIC, z=z, d=d))
  dic_df = data.frame(comp_no=comp_no,
                      model=c(paste0(cov_counts, "_cov_model")),
                      no_covs = cov_counts,
                      DIC = sapply(post_samples, calc_DIC, z=z, d=d))

  summary_df = mse_df %>%
    left_join(waic_df, by = join_by(comp_no, model, no_covs)) %>%
    right_join(dic_df, by = join_by(comp_no, model, no_covs)) %>%
    arrange(mse_true)

  return(summary_df)
}


