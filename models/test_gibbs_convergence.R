# Test Gibbs sampler convergence across multiple chains
# Assess Rhat and ESS for different nbasis values

library(Matrix)
library(coda)  # For Rhat and ESS calculations
library(dplyr)
source("models/spatial_basis_fh.R")

# Load test data
load("models/nc_testdata.RDA")

covs <- c("degree", "assistance")
X <- model.matrix(~., all_data[, covs, drop = FALSE])
y <- log(all_data$povPerc)
d_var <- all_data$povPercSE^2 /all_data$povPerc^2


cat("=== Gibbs Sampler: Multi-Chain Convergence Assessment ===\n\n")
cat("Configuration:\n")
cat("  Chains: 5\n")
cat("  Prior: a=b=c=d=1 (IG(1,1) for both variance components)\n")
cat("  MCMC: sample=2000, nburn=5000\n\n")

# Test configurations
# Note: NC test data has 37 positive eigenvalues, so nbasis capped at 35
nbasis_values <- c(5, 10, 20, 30, 35)
n_chains <- 5
results <- data.frame()

for (nbasis_val in nbasis_values) {
  cat(sprintf("========================================\n"))
  cat(sprintf("nbasis = %d\n", nbasis_val))
  cat(sprintf("========================================\n\n"))

  # Storage for chains
  tau_sq_chains <- list()
  sigma_sq_chains <- list()
  theta_chains <- list()

  # Run multiple chains with different seeds
  for (chain in 1:n_chains) {
    cat(sprintf("  Running chain %d/%d...", chain, n_chains))

    set.seed(1000 + chain)  # Different seed for each chain

    time_start <- Sys.time()

    fit <- spatial_basis_fh(
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

    # Store chains
    tau_sq_chains[[chain]] <- fit$tau.sq
    sigma_sq_chains[[chain]] <- fit$sigma.sq
    theta_chains[[chain]] <- fit$theta

    cat(sprintf(" %.1fs\n", elapsed))
  }

  # Convert to mcmc.list for coda diagnostics
  tau_sq_mcmc <- mcmc.list(lapply(tau_sq_chains, mcmc))
  sigma_sq_mcmc <- mcmc.list(lapply(sigma_sq_chains, mcmc))

  # Calculate Rhat for tau^2 and sigma^2
  tau_sq_rhat <- gelman.diag(tau_sq_mcmc)$psrf[1, "Point est."]
  sigma_sq_rhat <- gelman.diag(sigma_sq_mcmc)$psrf[1, "Point est."]

  # Calculate ESS for tau^2 and sigma^2
  tau_sq_ess <- effectiveSize(tau_sq_mcmc)
  sigma_sq_ess <- effectiveSize(sigma_sq_mcmc)

  # Calculate Rhat for theta parameters (all areas)
  theta_mcmc <- mcmc.list(lapply(theta_chains, mcmc))
  theta_gelman <- gelman.diag(theta_mcmc)
  theta_rhats <- theta_gelman$psrf[, "Point est."]

  max_theta_rhat <- max(theta_rhats)
  mean_theta_rhat <- mean(theta_rhats)

  # Posterior means
  mean_tau_sq <- mean(unlist(tau_sq_chains))
  mean_sigma_sq <- mean(unlist(sigma_sq_chains))

  cat("\n  Convergence diagnostics:\n")
  cat(sprintf("    tau^2   - Rhat: %.4f, ESS: %.0f\n", tau_sq_rhat, tau_sq_ess))
  cat(sprintf("    sigma^2 - Rhat: %.4f, ESS: %.0f\n", sigma_sq_rhat, sigma_sq_ess))
  cat(sprintf("    theta   - Rhat (max): %.4f, Rhat (mean): %.4f\n",
              max_theta_rhat, mean_theta_rhat))
  cat(sprintf("  Posterior means: tau^2 = %.6f, sigma^2 = %.6f\n\n",
              mean_tau_sq, mean_sigma_sq))

  # Store results
  results <- rbind(results, data.frame(
    nbasis = nbasis_val,
    tau_sq_rhat = tau_sq_rhat,
    sigma_sq_rhat = sigma_sq_rhat,
    max_theta_rhat = max_theta_rhat,
    mean_theta_rhat = mean_theta_rhat,
    tau_sq_ess = tau_sq_ess,
    sigma_sq_ess = sigma_sq_ess,
    mean_tau_sq = mean_tau_sq,
    mean_sigma_sq = mean_sigma_sq,
    stringsAsFactors = FALSE
  ))
}

# Summary
cat("\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n\n")

cat("Convergence across nbasis values:\n")
print(results %>%
  select(nbasis, tau_sq_rhat, sigma_sq_rhat, max_theta_rhat,
         tau_sq_ess, sigma_sq_ess))

cat("\n\nConvergence assessment:\n")
cat("  Good convergence: Rhat < 1.1\n")
cat("  Excellent convergence: Rhat < 1.05\n\n")

# Flag any issues
max_rhat <- max(c(results$tau_sq_rhat, results$sigma_sq_rhat, results$max_theta_rhat))
if (max_rhat < 1.05) {
  cat("  ✓ All parameters show excellent convergence (Rhat < 1.05)\n")
} else if (max_rhat < 1.1) {
  cat("  ✓ All parameters show good convergence (Rhat < 1.1)\n")
} else {
  cat("  ⚠ Some parameters show poor convergence (Rhat >= 1.1)\n")
  cat("    Consider increasing burn-in or number of iterations\n")
}

# Save results
write.csv(results, "results/gibbs_convergence_assessment.csv", row.names = FALSE)
cat("\nResults saved to results/gibbs_convergence_assessment.csv\n")
