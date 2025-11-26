library(Matrix)
library(LaplacesDemon)
library(dplyr)

# Source mcmc_helper.R (same directory)
if (file.exists("mcmc_helper.R")) {
  source("mcmc_helper.R")
} else if (file.exists("models/mcmc_helper.R")) {
  source("models/mcmc_helper.R")
}

#' Description: Spatial Basis Function Fay-Herriot Model with IID + Spatial Random Effects (BYM-style)
#' @param X the covariates (from model.matrix, include intercept)
#' @param y the response
#' @param d.var known sample *variances*
#' @param nbasis number of spatial basis functions to use (required)
#'   Note: Limited by number of positive eigenvalues from Moran's I decomposition
#'   Requesting nbasis > available eigenvalues will use all available and trigger warning
#' @param eig_moran (Optional) pre-computed eigendecomposition of Moran's I operator. Ouput from `compute_moran_eigen`
#' @param A adjacency matrix (required only if eig_moran not provided)
#' @param spatial_type (Optional) either "random" (default) or "fixed"
#'   - "random": spatial basis coefficients treated as random effects (eta ~ N(0, tau^2 I))
#'   - "fixed": spatial basis functions added as fixed covariates (no eta, no tau^2)
#' @param ndesired desired size of MCMC sample
#' @param nburn # of burn in iterations
#' @param nthin (Optional) how many iterations to thin by
#' @param hyp (Optional) list of the hyperparameters for IG priors (default is IG(1, 1))
#'   - c, d: hyperparameters for sigma^2 ~ IG(c, d) [IID variance]
#'   - a, b: hyperparameters for tau^2 ~ IG(a, b) [spatial variance] (only for spatial_type="random")
#' @param ini (Optional) list of the initial values for any of the parameters
#' @param verbose (Optional) Default T. Print statements with progress bar & computation time
#' @returns list of MCMC containers (matrices) for each set of parameters

spatial_basis_fh <- function(X, y, d.var, nbasis, eig_moran = NULL, A = NULL,
                             spatial_type = "random",
                             ndesired, nburn, nthin = 1,
                             hyp = list(), ini = list(), verbose = TRUE) {
  nsim <- nthin * ndesired
  m <- nrow(X) # number of areas
  p <- ncol(X) # number of covariates (including intercept)

  # ---- Construct spatial basis functions ----
  if (is.null(eig_moran)) {
    # Compute eigendecomposition on-the-fly
    if (is.null(A)) {
      stop("Either eig_moran or A must be provided")
    }
    eig_moran <- compute_moran_eigen(A, X)
  }

  # Select first nbasis positive eigenvectors
  S <- select_basis_from_eigen(eig_moran, nbasis)
  r <- ncol(S) # actual number of basis functions

  # ---- Spatial type determines model structure ----
  if (spatial_type == "fixed") {
    # Fixed effects model: spatial basis functions are covariates
    # theta = X_aug * beta + w where X_aug = [X, S]
    X_aug <- cbind(X, S)
    colnames(X_aug) <- c(colnames(X), paste0("basis_", 1:r))
    p_total <- ncol(X_aug) # total number of fixed effects
    use_random_spatial <- FALSE
  } else if (spatial_type == "random") {
    # Random effects model: spatial basis coefficients are random
    # theta = X * beta + S * eta + w where eta ~ N(0, tau^2 I)
    X_aug <- X
    p_total <- p
    use_random_spatial <- TRUE
    Z <- cbind(X, S) # Z = [X  S], dimension m × (p+r)
  } else {
    stop("spatial_type must be 'random' or 'fixed'")
  }

  # ---- Hyperparameters for variance priors ----
  # sigma^2 ~ IG(c, d) for IID effects w
  c <- ifelse(is.null(hyp$c), 5e-05, hyp$c)
  d <- ifelse(is.null(hyp$d), 5e-05, hyp$d)

  if (use_random_spatial) {
    # tau^2 ~ IG(a, b) for spatial basis coefficients eta
    a <- ifelse(is.null(hyp$a), 5e-05, hyp$a)
    b <- ifelse(is.null(hyp$b), 5e-05, hyp$b)
    details <- paste0(
      "Model: Random spatial effects | ",
      "Priors: sigma^2 ~ IG(", c, ", ", d, "); ",
      "tau^2 ~ IG(", a, ", ", b, ")"
    )
  } else {
    details <- paste0(
      "Model: Fixed spatial basis covariates | ",
      "Priors: sigma^2 ~ IG(", c, ", ", d, ")"
    )
  }

  if (verbose) {
    cat(details, "\n")
  }

  # ---- MCMC containers ----
  Res_beta <- mcmc_mat(nsim / nthin, p_total, colnames(X_aug))
  Res_w <- mcmc_mat(nsim / nthin, m, paste0("w_", 1:m))
  Res_sigma.sq <- mcmc_mat(nsim / nthin, 1, "sigma.sq")
  Res_theta <- mcmc_mat(nsim / nthin, m, paste0("theta_", 1:m))

  if (use_random_spatial) {
    Res_eta <- mcmc_mat(nsim / nthin, r, paste0("eta_", 1:r))
    Res_v <- mcmc_mat(nsim / nthin, m, paste0("v_", 1:m))
    Res_tau.sq <- mcmc_mat(nsim / nthin, 1, "tau.sq")
  }

  # ---- Initial values ----
  # fixed effects
  ls_fit <- lm(y ~ X_aug - 1)
  beta <- if (is.null(ini$beta)) ls_fit$coefficients else ini$beta

  if (use_random_spatial) {
    # spatial effects (random)
    eta <- if (is.null(ini$eta)) rnorm(r) else ini$eta
    v <- as.numeric(S %*% eta)
  }

  # IID random effects
  w <- if (is.null(ini$w)) rnorm(m) else ini$w

  # Variance components
  sigma.sq <- if (is.null(ini$sigma.sq)) mean(ls_fit$residuals^2) else ini$sigma.sq
  if (use_random_spatial) {
    tau.sq <- if (is.null(ini$tau.sq)) 1 else ini$tau.sq
  }

  # ---- Progress tracking ----
  if (verbose) {
    ptm <- start_chain(nsim, nthin, nburn)
    pb <- txtProgressBar(min = 0, max = nsim + nburn, style = 3)
  }

  # ---- MCMC chain ----
  for (index in 1:(nsim + nburn)) {
    if (use_random_spatial) {
      # ========================================
      # RANDOM SPATIAL MODEL
      # ========================================

      # Step 1: Update gamma = (beta, eta) jointly
      # Precision matrix: Omega_gamma = Z'D^{-1}Z + diag(0_p, tau^{-2}I_r)
      Omega_gamma <- crossprod(Z, (1 / d.var) * Z)
      # Add prior precision for eta (last r diagonal elements)
      diag(Omega_gamma)[(p + 1):(p + r)] <- diag(Omega_gamma)[(p + 1):(p + r)] + (1 / tau.sq)

      # Cholesky decomposition
      U <- chol(forceSymmetric(Omega_gamma))

      # Mean: mu_gamma = Omega_gamma^{-1} Z'D^{-1}(y - w)
      a_vec <- crossprod(Z, (y - w) / d.var)
      b_vec <- rnorm(p + r)

      # Sample gamma using https://www.stat.berkeley.edu/~paciorek/research/techVignettes/techVignette6.pdf
      gamma <- backsolve(U, backsolve(U, a_vec, transpose = TRUE) + b_vec)

      # Partition gamma into beta and eta
      beta <- gamma[1:p]
      eta <- gamma[(p + 1):(p + r)]

      # Step 2: Update w (IID random effects)
      var_w <- (d.var * sigma.sq) / (d.var + sigma.sq)
      mean_w <- var_w * (y - Z %*% gamma) / d.var
      w <- rnorm(m, mean = mean_w, sd = sqrt(var_w))

      # Step 3: Update sigma^2
      sigma.sq <- 1 / rgamma(1, shape = c + m / 2, rate = d + sum(w^2) / 2)

      # Step 4: Update tau^2
      tau.sq <- 1 / rgamma(1, shape = a + r / 2, rate = b + sum(eta^2) / 2)

      # Derived quantities
      v <- as.numeric(S %*% eta)
      theta <- as.numeric(X %*% beta + v + w)
    } else {
      # ========================================
      # FIXED SPATIAL MODEL (basis functions as covariates)
      # ========================================

      # Step 1: Update beta (now includes spatial basis coefficients)
      # Precision matrix: Omega_beta = X_aug'D^{-1}X_aug (no prior on fixed effects)
      Omega_beta <- crossprod(X_aug, (1 / d.var) * X_aug)

      # Cholesky decomposition
      U <- chol(forceSymmetric(Omega_beta))

      # Mean: mu_beta = Omega_beta^{-1} X_aug'D^{-1}(y - w)
      a_vec <- crossprod(X_aug, (y - w) / d.var)
      b_vec <- rnorm(p_total)

      # Sample beta
      beta <- backsolve(U, backsolve(U, a_vec, transpose = TRUE) + b_vec)

      # Step 2: Update w (IID random effects)
      var_w <- (d.var * sigma.sq) / (d.var + sigma.sq)
      mean_w <- var_w * (y - X_aug %*% beta) / d.var
      w <- rnorm(m, mean = mean_w, sd = sqrt(var_w))

      # Step 3: Update sigma^2
      sigma.sq <- 1 / rgamma(1, shape = c + m / 2, rate = d + sum(w^2) / 2)

      # Derived quantities
      theta <- as.numeric(X_aug %*% beta + w)
    }

    # ========================================
    # Save results (after burn-in, with thinning)
    # ========================================
    if (index > nburn && (index - nburn) %% nthin == 0) {
      idx <- (index - nburn) / nthin
      Res_beta[idx, ] <- beta
      Res_w[idx, ] <- w
      Res_sigma.sq[idx, ] <- sigma.sq
      Res_theta[idx, ] <- theta

      if (use_random_spatial) {
        Res_eta[idx, ] <- eta
        Res_v[idx, ] <- v
        Res_tau.sq[idx, ] <- tau.sq
      }
    }

    if (verbose) setTxtProgressBar(pb, index)
  }

  if (verbose) {
    writeLines("")
    print(proc.time() - ptm)
  }

  # Return results
  results <- list(
    beta = Res_beta,
    w = Res_w,
    sigma.sq = Res_sigma.sq,
    theta = Res_theta,
    spatial_type = spatial_type,
    prior_details = details
  )

  if (use_random_spatial) {
    results$eta <- Res_eta
    results$v <- Res_v
    results$tau.sq <- Res_tau.sq
  }

  return(results)
}

#' Compute eigendecomposition of Moran's I operator
#' @param A adjacency matrix
#' @param X covariate matrix
#' @returns eigendecomposition list with $values and $vectors
compute_moran_eigen <- function(A, X) {
  m <- nrow(A)

  # Projection matrix: P_X = X(X'X)^{-1}X' (using Cholesky for numerical stability)
  XtX_inv <- chol2inv(chol(crossprod(X)))
  P_X <- X %*% XtX_inv %*% t(X)

  # Moran's I operator: G = (I - P_X) A (I - P_X)
  I_m <- diag(m)
  G <- (I_m - P_X) %*% A %*% (I_m - P_X)

  # Eigendecomposition
  eig_out <- eigen(G, symmetric = TRUE)

  return(eig_out)
}

#' Select spatial basis functions from eigendecomposition
#' @param eig_moran eigendecomposition of Moran's I operator (list with $values and $vectors)
#' @param nbasis number of basis functions to select
#' @returns S the spatial basis matrix (m x r)
select_basis_from_eigen <- function(eig_moran, nbasis) {
  # Select positive eigenvalues
  eigvals <- eig_moran$values
  idx_pos <- which(eigvals > 1e-10) # Small tolerance for numerical stability

  # Select first nbasis eigenvectors from positive spectrum
  r <- min(nbasis, length(idx_pos)) # Actual number of basis functions
  if (r < nbasis) {
    warning(sprintf("Only %d positive eigenvalues available (requested %d)", r, nbasis))
  }

  idx_select <- idx_pos[1:r]
  S <- eig_moran$vectors[, idx_select] # m × r basis matrix

  return(S)
}
