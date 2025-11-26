# Test construct_moran_basis function

library(Matrix)
source("models/spatial_basis_fh.R")

cat("=== Testing construct_moran_basis() ===\n\n")

# Test 1: NC data
cat("Test 1: NC test data\n")
load("models/nc_testdata.RDA")
cat("Loaded objects:", ls(), "\n")

if (exists("A") && exists("all_data")) {
  cat("  Adjacency matrix A: ", nrow(A), "×", ncol(A), "\n")

  # Create X from all_data
  X <- model.matrix(~., all_data[, c("degree", "assistance")])
  cat("  Covariate matrix X: ", nrow(X), "×", ncol(X), "\n")

  cat("\n  Running construct_moran_basis(A, X, nbasis=10)...\n")
  basis_nc <- construct_moran_basis(A, X, nbasis = 10)

  cat("  ✓ Success!\n")
  cat("    - S dimensions:", nrow(basis_nc$S), "×", ncol(basis_nc$S), "\n")
  cat("    - Number of basis functions (r):", basis_nc$r, "\n")
  cat("    - Eigenvalues (first 5):", head(basis_nc$eigvals, 5), "\n")
} else {
  cat("  ✗ Error: A or all_data not found in nc_testdata.RDA\n")
  cat("    Available objects:", ls(), "\n")
}

cat("\n")

# Test 2: IL data
cat("Test 2: Illinois test data\n")
rm(list = setdiff(ls(), "construct_moran_basis"))  # Clear workspace but keep function
load("models/il_testdata.RDA")
cat("Loaded objects:", ls(), "\n")

if (exists("A") && exists("all_data")) {
  cat("  Adjacency matrix A: ", nrow(A), "×", ncol(A), "\n")

  # Create X from all_data
  X <- model.matrix(~., all_data[, c("degree", "assistance")])
  cat("  Covariate matrix X: ", nrow(X), "×", ncol(X), "\n")

  cat("\n  Running construct_moran_basis(A, X, nbasis=20)...\n")
  basis_il <- construct_moran_basis(A, X, nbasis = 20)

  cat("  ✓ Success!\n")
  cat("    - S dimensions:", nrow(basis_il$S), "×", ncol(basis_il$S), "\n")
  cat("    - Number of basis functions (r):", basis_il$r, "\n")
  cat("    - Eigenvalues (first 5):", head(basis_il$eigvals, 5), "\n")

  # Check orthogonality
  cat("\n  Checking basis properties:\n")
  cat("    - S'S is diagonal?", max(abs(crossprod(basis_il$S) - diag(basis_il$r))) < 1e-10, "\n")
  cat("    - All eigenvalues positive?", all(basis_il$eigvals > 0), "\n")

} else if (exists("A")) {
  cat("  ✗ Error: all_data not found in il_testdata.RDA\n")
  cat("    Available objects:", ls(), "\n")
} else {
  cat("  ✗ Error: A not found in il_testdata.RDA\n")
  cat("    Available objects:", ls(), "\n")
}

cat("\n=== Tests complete ===\n")
