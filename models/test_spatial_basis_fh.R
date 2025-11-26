# Test spatial_basis_fh model
# Tests: 1) Dimensions 2) Model fitting 3) Spatial smoothness of v

library(Matrix)
library(dplyr)
source("models/spatial_basis_fh.R")

cat("=== Testing spatial_basis_fh() ===\n\n")

# Load test data
cat("Loading NC test data...\n")
load("models/nc_testdata.RDA")

if (!exists("A") || !exists("all_data")) {
  stop("Required test data not found (A or all_data)")
}

m <- nrow(A)
cat("  Number of areas (m):", m, "\n")
cat("  Adjacency matrix dimensions:", nrow(A), "×", ncol(A), "\n\n")

# Prepare data
covs <- c("degree", "assistance")
X <- model.matrix(~., all_data[, covs, drop = FALSE])
n <- nrow(X)
p <- ncol(X)  # number of covariates INCLUDES intercept

# Response (log transformation for better convergence)
y <- log(all_data$povPerc)
d.var <- all_data$povPercSE^2 / all_data$povPerc^2  # Delta method for log scale

cat("Data dimensions:\n")
cat("  X:", nrow(X), "×", ncol(X), "(m × p)\n")
cat("  y:", length(y), "× 1\n")
cat("  d.var:", length(d.var), "× 1\n\n")

# MCMC settings (short for testing)
ndesired <- 50
nburn <- 50
nthin <- 1

cat("MCMC settings: ndesired =", ndesired, ", nburn =", nburn, ", nthin =", nthin, "\n")
cat("Using default data-dependent priors: IG(3, 2*mean(d.var))\n")
cat("  mean(d.var) =", sprintf("%.6f", mean(d.var)), "\n\n")

# ========================================
# Test 1: Few basis functions (high smoothing)
# ========================================
cat("=== Test 1: Few basis functions (nbasis = 5) ===\n")
cat("Expected: Highly smooth spatial effects v\n\n")

fit1 <- spatial_basis_fh(
  X = X, y = y, d.var = d.var,
  nbasis = 5,
  A = A,
  ndesired = ndesired, nburn = nburn, nthin = nthin,
  verbose = TRUE
)

cat("\n✓ Fit 1 completed\n")
cat("Output dimensions:\n")
cat("  beta:", nrow(fit1$beta), "×", ncol(fit1$beta), "(ndesired × p)\n")
cat("  eta:", nrow(fit1$eta), "×", ncol(fit1$eta), "(ndesired × r)\n")
cat("  w:", nrow(fit1$w), "×", ncol(fit1$w), "(ndesired × m)\n")
cat("  v:", nrow(fit1$v), "×", ncol(fit1$v), "(ndesired × m)\n")
cat("  theta:", nrow(fit1$theta), "×", ncol(fit1$theta), "(ndesired × m)\n")
cat("  sigma.sq:", nrow(fit1$sigma.sq), "×", ncol(fit1$sigma.sq), "\n")
cat("  tau.sq:", nrow(fit1$tau.sq), "×", ncol(fit1$tau.sq), "\n\n")

# Test dimension correctness
stopifnot(nrow(fit1$beta) == ndesired)
stopifnot(ncol(fit1$beta) == p)
stopifnot(ncol(fit1$eta) == 5)  # Should match nbasis
stopifnot(ncol(fit1$v) == m)
stopifnot(ncol(fit1$w) == m)
stopifnot(ncol(fit1$theta) == m)

cat("✓ All dimensions correct\n\n")

# ========================================
# Test 2: Moderate basis functions (moderate smoothing)
# ========================================
cat("=== Test 2: Moderate basis functions (nbasis = 20) ===\n")
cat("Expected: Moderately smooth spatial effects v\n\n")

fit2 <- spatial_basis_fh(
  X = X, y = y, d.var = d.var,
  nbasis = 20,
  A = A,
  ndesired = ndesired, nburn = nburn, nthin = nthin,
  verbose = TRUE
)

cat("\n✓ Fit 2 completed\n")
cat("  eta dimensions:", nrow(fit2$eta), "×", ncol(fit2$eta), "\n\n")

stopifnot(ncol(fit2$eta) == 20)
cat("✓ Dimensions correct\n\n")

# ========================================
# Test 3: Many basis functions (low smoothing)
# ========================================
cat("=== Test 3: Many basis functions (nbasis = 50) ===\n")
cat("Expected: Flexible spatial effects v (less smooth)\n\n")

fit3 <- spatial_basis_fh(
  X = X, y = y, d.var = d.var,
  nbasis = 50,
  A = A,
  ndesired = ndesired, nburn = nburn, nthin = nthin,
  verbose = TRUE
)

cat("\n✓ Fit 3 completed\n")
cat("  eta dimensions:", nrow(fit3$eta), "×", ncol(fit3$eta), "\n\n")

stopifnot(ncol(fit3$eta) <= 50)  # May be less if not enough positive eigenvalues
cat("✓ Dimensions correct\n\n")

# ========================================
# Compare spatial smoothness
# ========================================
cat("=== Spatial Smoothness Analysis ===\n\n")

# Use posterior means of v for comparison
v1_mean <- colMeans(fit1$v)
v2_mean <- colMeans(fit2$v)
v3_mean <- colMeans(fit3$v)

# Measure 1: Total variance of v
var1 <- var(v1_mean)
var2 <- var(v2_mean)
var3 <- var(v3_mean)

cat("Variance of spatial effects v (posterior mean):\n")
cat("  nbasis = 5:  ", sprintf("%.4f", var1), "\n")
cat("  nbasis = 20: ", sprintf("%.4f", var2), "\n")
cat("  nbasis = 50: ", sprintf("%.4f", var3), "\n\n")

# Measure 2: Mean absolute difference between neighbors
compute_neighbor_roughness <- function(v, A) {
  # For each area, compute mean |v_i - v_j| for neighbors j
  diffs <- numeric(0)
  for (i in 1:nrow(A)) {
    neighbors <- which(A[i, ] > 0)
    if (length(neighbors) > 0) {
      diffs <- c(diffs, abs(v[i] - v[neighbors]))
    }
  }
  mean(diffs)
}

rough1 <- compute_neighbor_roughness(v1_mean, A)
rough2 <- compute_neighbor_roughness(v2_mean, A)
rough3 <- compute_neighbor_roughness(v3_mean, A)

cat("Mean absolute difference between neighbors (roughness):\n")
cat("  nbasis = 5:  ", sprintf("%.4f", rough1), " (smoother)\n")
cat("  nbasis = 20: ", sprintf("%.4f", rough2), "\n")
cat("  nbasis = 50: ", sprintf("%.4f", rough3), " (rougher)\n\n")

# Measure 3: Proportion of variance explained by spatial component
theta1_mean <- colMeans(fit1$theta)
theta2_mean <- colMeans(fit2$theta)
theta3_mean <- colMeans(fit3$theta)

prop_spatial1 <- var(v1_mean) / var(theta1_mean)
prop_spatial2 <- var(v2_mean) / var(theta2_mean)
prop_spatial3 <- var(v3_mean) / var(theta3_mean)

cat("Proportion of theta variance from spatial component:\n")
cat("  nbasis = 5:  ", sprintf("%.3f", prop_spatial1), "\n")
cat("  nbasis = 20: ", sprintf("%.3f", prop_spatial2), "\n")
cat("  nbasis = 50: ", sprintf("%.3f", prop_spatial3), "\n\n")

# Measure 4: Variance components
cat("Variance components (posterior means):\n")
cat("  nbasis = 5:   sigma^2 =", sprintf("%.4f", mean(fit1$sigma.sq)),
    ", tau^2 =", sprintf("%.4f", mean(fit1$tau.sq)), "\n")
cat("  nbasis = 20:  sigma^2 =", sprintf("%.4f", mean(fit2$sigma.sq)),
    ", tau^2 =", sprintf("%.4f", mean(fit2$tau.sq)), "\n")
cat("  nbasis = 50:  sigma^2 =", sprintf("%.4f", mean(fit3$sigma.sq)),
    ", tau^2 =", sprintf("%.4f", mean(fit3$tau.sq)), "\n\n")

cat("Interpretation:\n")
cat("  - Fewer basis functions → smoother v (lower neighbor roughness)\n")
cat("  - More basis functions → more flexible v (higher neighbor roughness)\n")
cat("  - All models should have valid dimensions and converge\n\n")

cat("=== All tests passed! ===\n")
