library(Matrix)
library(LaplacesDemon)
library(dplyr)

source("models/mcmc_helper.R")

#' Description: Spatial Basis Function Fay-Herriot Model with IID + Spatial Random Effects (BYM-style)
#'
#' Model: theta_i = X_i'beta + v_i + w_i
#'        v = S * eta, eta ~ N(0, tau^2 I_r) [spatial via basis functions]
#'        w ~ N(0, sigma^2 I_m) [IID random effects]
#'
#' Dimensions:
#'   X: m × p (covariates, includes intercept)
#'   S: m × r (spatial basis functions, r < m)
#'   Z = [X S]: m × (p+r) (augmented design)
#'   gamma = (beta', eta')': (p+r) × 1 (joint parameter)
#'   w: m × 1 (IID random effects)
#'   v = S*eta: m × 1 (spatial random effects)
#'   theta = Z*gamma + w: m × 1 (area parameters)
#'
#' @param X the covariates (from model.matrix, include intercept)
#' @param y the response
#' @param d.var known sample *variances*
#' @param nbasis number of spatial basis functions to use (required)
#' @param eig_moran (Optional) pre-computed eigendecomposition of Moran's I operator.
#'   Should be a list with $values and $vectors from eigen(). If provided, A is not needed.
#' @param A adjacency matrix (required only if eig_moran not provided)
#' @param ndesired desired size of MCMC sample
#' @param nburn # of burn in iterations
#' @param nthin how many iterations to thin by
#' @param hyp (Optional) list of the hyperparameters
#'   - c, d: hyperparameters for sigma^2 ~ IG(c, d) [IID variance]
#'   - a, b: hyperparameters for tau^2 ~ IG(a, b) [spatial variance]
#' @param ini (Optional) list of the initial values for any of the parameters
#' @param verbose (Optional) Default T. Print statements with progress bar & computation time
#' @returns list of MCMC containers (matrices) for each set of parameters

spatial_basis_fh <- function(X, y, d.var, nbasis, eig_moran = NULL, A = NULL,
                             ndesired, nburn, nthin,
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

  # ---- Hyperparameters for variance priors ----
  # Data-dependent priors (standard SAE practice)
  # sigma^2 ~ IG(c, d) for IID effects w
  c <- ifelse(is.null(hyp$c), 3, hyp$c)
  d <- ifelse(is.null(hyp$d), 2 * mean(d.var), hyp$d)

  # tau^2 ~ IG(a, b) for spatial basis coefficients eta
  a <- ifelse(is.null(hyp$a), 3, hyp$a)
  b <- ifelse(is.null(hyp$b), 2 * mean(d.var), hyp$b)

  details <- paste0(
    "Priors: sigma^2 ~ IG(", c, ", ", d, "); ",
    "tau^2 ~ IG(", a, ", ", b, ")"
  )

  if (verbose) {
    cat(details, "\n")
  }

  # ---- MCMC containers ----
  Res_beta <- mcmc_mat(nsim / nthin, p, colnames(X))
  Res_eta <- mcmc_mat(nsim / nthin, r, paste0("eta_", 1:r))
  Res_w <- mcmc_mat(nsim / nthin, m, paste0("w_", 1:m))
  Res_v <- mcmc_mat(nsim / nthin, m, paste0("v_", 1:m))
  Res_sigma.sq <- mcmc_mat(nsim / nthin, 1, "sigma.sq")
  Res_tau.sq <- mcmc_mat(nsim / nthin, 1, "tau.sq")
  Res_theta <- mcmc_mat(nsim / nthin, m, paste0("theta_", 1:m))

  # ---- Initial values ----
  ls_fit <- lm(y ~ X - 1)

  # beta, eta (combined as gamma in joint sampling)
  beta <- if (is.null(ini$beta)) ls_fit$coefficients else ini$beta
  eta <- if (is.null(ini$eta)) rnorm(r) else ini$eta

  # IID random effects w
  w <- if (is.null(ini$w)) rnorm(m) else ini$w

  # Variance components
  sigma.sq <- if (is.null(ini$sigma.sq)) mean(ls_fit$residuals^2) else ini$sigma.sq
  tau.sq <- if (is.null(ini$tau.sq)) 1 else ini$tau.sq

  # Derived: spatial effects v = S * eta
  v <- as.numeric(S %*% eta)

  # ---- Construct augmented design matrix (fixed across iterations) ----
  Z <- cbind(X, S)  # Z = [X  S], dimension m × (p+r)

  # ---- Progress tracking ----
  if (verbose) {
    ptm <- start_chain(nsim, nthin, nburn)
    pb <- txtProgressBar(min = 0, max = nsim + nburn, style = 3)
  }

  # ---- MCMC chain ----
  for (index in 1:(nsim + nburn)) {
    # ========================================
    # Step 1: Update gamma = (beta, eta) jointly
    # ========================================

    # Precision matrix: Omega_gamma = Z'D^{-1}Z + diag(0_p, tau^{-2}I_r)
    Omega_gamma <- crossprod(Z, (1 / d.var) * Z)
    # Add prior precision for eta (last r diagonal elements)
    diag(Omega_gamma)[(p + 1):(p + r)] <- diag(Omega_gamma)[(p + 1):(p + r)] + (1 / tau.sq)

    # Cholesky decomposition
    U <- chol(forceSymmetric(Omega_gamma))

    # Mean: mu_gamma = Omega_gamma^{-1} Z'D^{-1}(y - w)
    a_vec <- crossprod(Z, (y - w) / d.var)
    b_vec <- rnorm(p + r)

    # Sample gamma using Paciorek's method
    gamma <- backsolve(U, backsolve(U, a_vec, transpose = TRUE) + b_vec)

    # Partition gamma into beta and eta
    beta <- gamma[1:p]
    eta <- gamma[(p + 1):(p + r)]


    # ========================================
    # Step 2: Update w (IID random effects)
    # ========================================
    # Compute Z*gamma (X*beta + S*eta)
    Z_gamma <- as.numeric(Z %*% gamma)

    # Elementwise variance: s_i^2 = (d_i^{-1} + sigma^{-2})^{-1}
    var_w <- (d.var * sigma.sq) / (d.var + sigma.sq)

    # Elementwise mean: m_i = s_i^2 * (y_i - z_i'gamma) / d_i
    mean_w <- var_w * (y - Z_gamma) / d.var

    # Sample w elementwise
    w <- rnorm(m, mean = mean_w, sd = sqrt(var_w))


    # ========================================
    # Step 3: Update sigma^2 (IID variance)
    # ========================================
    # sigma^2 ~ IG(c + m/2, d + w'w/2)
    sigma.sq <- 1 / rgamma(1, shape = c + m / 2, rate = d + sum(w^2) / 2)


    # ========================================
    # Step 4: Update tau^2 (spatial variance)
    # ========================================
    # tau^2 ~ IG(a + r/2, b + eta'eta/2)
    tau.sq <- 1 / rgamma(1, shape = a + r / 2, rate = b + sum(eta^2) / 2)


    # ========================================
    # Derived quantities
    # ========================================
    # Compute v = S * eta
    v <- as.numeric(S %*% eta)

    # Compute theta = X*beta + v + w
    theta <- as.numeric(X %*% beta + v + w)

    # ========================================
    # Save results (after burn-in, with thinning)
    # ========================================
    if (index > nburn && (index - nburn) %% nthin == 0) {
      idx <- (index - nburn) / nthin
      Res_beta[idx, ] <- beta
      Res_eta[idx, ] <- eta
      Res_w[idx, ] <- w
      Res_v[idx, ] <- v
      Res_sigma.sq[idx, ] <- sigma.sq
      Res_tau.sq[idx, ] <- tau.sq
      Res_theta[idx, ] <- theta
    }

    if (verbose) setTxtProgressBar(pb, index)
  }

  if (verbose) {
    writeLines("")
    print(proc.time() - ptm)
  }

  return(list(
    beta = Res_beta,
    eta = Res_eta,
    w = Res_w,
    v = Res_v,
    sigma.sq = Res_sigma.sq,
    tau.sq = Res_tau.sq,
    theta = Res_theta,
    prior_details = details
  ))
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
