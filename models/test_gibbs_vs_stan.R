# Test R Gibbs sampler with a=c=2, b=d=0.25
# Compare timing and results with Stan

library(Matrix)
library(rstan)
library(dplyr)
source("models/spatial_basis_fh.R")

rstan_options(auto_write = TRUE)
options(mc.cores = 1)  # Use 1 core for fair comparison

# Load test data
load("models/nc_testdata.RDA")

covs <- c("degree", "assistance")
X <- model.matrix(~., all_data[, covs, drop = FALSE])
y <- log(all_data$povPerc)
d_var <- all_data$povPercSE^2 /all_data$povPerc^2

cat("=== Gibbs vs Stan: Timing Comparison ===\n\n")
cat("Prior: a=b=c=d=1 (IG(1,1) for both variance components)\n")
cat("MCMC: sample=2000, nburn=5000 \n")
cat("STAN: sample=2000, warmup=1000")

# Test configurations
nbasis_values <- seq(5, 35, by=5)
results <- data.frame()

for (nbasis_val in nbasis_values) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("nbasis = %d\n", nbasis_val))
  cat(sprintf("========================================\n\n"))

  # ========================================
  # Test R Gibbs Sampler
  # ========================================
  cat("--- R Gibbs Sampler ---\n")

  gibbs_result <- tryCatch({
    time_start <- Sys.time()

    fit_gibbs <- spatial_basis_fh(
      X = X,
      y = y,
      d.var = d_var,
      A = A,
      nbasis = nbasis_val,
      ndesired = 2000,
      nburn = 5000,
      nthin = 1,
      verbose = FALSE
    )

    time_end <- Sys.time()
    elapsed <- as.numeric(difftime(time_end, time_start, units = "secs"))

    # Check results
    mean_tau_sq <- mean(fit_gibbs$tau.sq)
    mean_sigma_sq <- mean(fit_gibbs$sigma.sq)

    cat(sprintf("  ✓ SUCCESS: %.2f seconds\n", elapsed))
    cat(sprintf("    tau^2: %.6f, sigma^2: %.6f\n", mean_tau_sq, mean_sigma_sq))

    list(
      success = TRUE,
      time = elapsed,
      mean_tau_sq = mean_tau_sq,
      mean_sigma_sq = mean_sigma_sq
    )

  }, error = function(e) {
    cat(sprintf("  ✗ FAILED: %s\n", conditionMessage(e)))
    list(
      success = FALSE,
      error = conditionMessage(e),
      time = NA
    )
  })

  # ========================================
  # Test Stan
  # ========================================
  cat("\n--- Stan (1 chain for fair comparison) ---\n")

  stan_result <- tryCatch({
    # Pre-compute basis (same as Gibbs would do internally)
    eig_moran <- compute_moran_eigen(A, X)
    S <- select_basis_from_eigen(eig_moran, nbasis_val)
    r <- ncol(S)

    stan_data <- list(
      m = nrow(X),
      p = ncol(X),
      r = r,
      y = as.numeric(y),
      d_var = as.numeric(d_var),
      X = X,
      S = S,
      prior_sigma_shape = 1,
      prior_sigma_scale = 1,
      prior_tau_shape = 1,
      prior_tau_scale = 1,
      use_ncp = 1L
    )

    time_start <- Sys.time()

    fit_stan <- stan(
      file = "models/spatial_basis_fh.stan",
      data = stan_data,
      chains = 1,  # 1 chain for fair comparison
      warmup = 1000,
      iter = 3000,  # 1000 warmup + 2000 sampling
      cores = 1,
      refresh = 0,
      verbose = FALSE
    )

    time_end <- Sys.time()
    elapsed <- as.numeric(difftime(time_end, time_start, units = "secs"))

    # Extract results
    mean_tau_sq <- mean(rstan::extract(fit_stan, "tau_sq")$tau_sq)
    mean_sigma_sq <- mean(rstan::extract(fit_stan, "sigma_sq")$sigma_sq)

    # Diagnostics (Rhat omitted - not meaningful with 1 chain)
    summary_df <- summary(fit_stan)$summary %>% as.data.frame()
    min_ess <- min(summary_df$n_eff, na.rm = TRUE)

    sampler_params <- get_sampler_params(fit_stan, inc_warmup = FALSE)
    n_divergences <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))

    cat(sprintf("  ✓ SUCCESS: %.2f seconds\n", elapsed))
    cat(sprintf("    tau^2: %.6f, sigma^2: %.6f\n", mean_tau_sq, mean_sigma_sq))
    cat(sprintf("    ESS: %.0f, divergences: %d\n",
                min_ess, n_divergences))

    list(
      success = TRUE,
      time = elapsed,
      mean_tau_sq = mean_tau_sq,
      mean_sigma_sq = mean_sigma_sq,
      min_ess = min_ess,
      n_divergences = n_divergences
    )

  }, error = function(e) {
    cat(sprintf("  ✗ FAILED: %s\n", conditionMessage(e)))
    list(
      success = FALSE,
      error = conditionMessage(e),
      time = NA
    )
  })

  # Store results
  results <- rbind(results, data.frame(
    nbasis = nbasis_val,
    method = "R Gibbs",
    success = gibbs_result$success,
    time_sec = gibbs_result$time,
    mean_tau_sq = ifelse(gibbs_result$success, gibbs_result$mean_tau_sq, NA),
    mean_sigma_sq = ifelse(gibbs_result$success, gibbs_result$mean_sigma_sq, NA),
    min_ess = NA,
    n_divergences = NA,
    stringsAsFactors = FALSE
  ))

  results <- rbind(results, data.frame(
    nbasis = nbasis_val,
    method = "Stan",
    success = stan_result$success,
    time_sec = stan_result$time,
    mean_tau_sq = ifelse(stan_result$success, stan_result$mean_tau_sq, NA),
    mean_sigma_sq = ifelse(stan_result$success, stan_result$mean_sigma_sq, NA),
    min_ess = ifelse(stan_result$success, stan_result$min_ess, NA),
    n_divergences = ifelse(stan_result$success, stan_result$n_divergences, NA),
    stringsAsFactors = FALSE
  ))
}

# Summary
cat("\n\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n\n")

cat("Success rates:\n")
print(results %>%
  group_by(method) %>%
  summarize(
    n_total = n(),
    n_success = sum(success),
    pct_success = 100 * mean(success),
    .groups = "drop"
  ))

cat("\n\nTiming comparison (successful runs only):\n")
timing <- results %>%
  filter(success == TRUE) %>%
  select(nbasis, method, time_sec, mean_tau_sq, mean_sigma_sq)

if (nrow(timing) > 0) {
  print(timing %>% arrange(nbasis, method))

  # Calculate speedup
  cat("\n\nSpeedup (Gibbs time / Stan time):\n")
  timing_wide <- timing %>%
    select(nbasis, method, time_sec) %>%
    tidyr::pivot_wider(names_from = method, values_from = time_sec)

  if ("R Gibbs" %in% names(timing_wide) && "Stan" %in% names(timing_wide)) {
    timing_wide <- timing_wide %>%
      mutate(speedup = `R Gibbs` / Stan)
    print(timing_wide)
  }
}

# Save results
write.csv(results, "results/gibbs_vs_stan_timing.csv", row.names = FALSE)
cat("\n\nResults saved to results/gibbs_vs_stan_timing.csv\n")
