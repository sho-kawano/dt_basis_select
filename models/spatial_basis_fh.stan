// spatial_basis_fh.stan
// Spatial Basis Function Fay-Herriot Model (Non-Centered Parameterization)
//
// Model:
//   y_i ~ N(theta_i, d_i)           [observation with known variance d_i]
//   theta_i = X_i'beta + v_i + w_i  [area-specific means]
//   v = S * eta                      [spatial effects via basis functions]
//   eta ~ N(0, tau^2 I_r)           [spatial basis coefficients]
//   w ~ N(0, sigma^2 I_m)           [IID random effects]
//   sigma^2 ~ IG(c, d)              [IID variance prior]
//   tau^2 ~ IG(a, b)                [spatial variance prior]
//
// Parameterization:
//   Uses non-centered parameterization (NCP) for eta and w to avoid funnel
//   geometries when variance parameters can be small. This improves HMC sampling
//   efficiency compared to centered parameterization.
//
//   NCP: eta = sqrt(tau^2) * eta_raw where eta_raw ~ N(0,1)
//        w = sqrt(sigma^2) * w_raw where w_raw ~ N(0,1)
//
//   Toggle with use_ncp flag: 1 = NCP (recommended), 0 = centered
//
// Dimensions:
//   m = number of areas
//   p = number of covariates (includes intercept)
//   r = number of spatial basis functions

data {
  int<lower=1> m;                  // number of areas
  int<lower=1> p;                  // number of covariates
  int<lower=1> r;                  // number of basis functions

  vector[m] y;                     // response variable
  vector<lower=0>[m] d_var;        // known sampling variances
  matrix[m, p] X;                  // covariate matrix
  matrix[m, r] S;                  // spatial basis functions

  // Prior hyperparameters
  real<lower=0> prior_sigma_shape; // shape for sigma^2 ~ IG(c, d)
  real<lower=0> prior_sigma_scale; // scale for sigma^2 ~ IG(c, d)
  real<lower=0> prior_tau_shape;   // shape for tau^2 ~ IG(a, b)
  real<lower=0> prior_tau_scale;   // scale for tau^2 ~ IG(a, b)

  // Control for parameterization
  int<lower=0, upper=1> use_ncp;   // 1 = non-centered, 0 = centered
}

transformed data {
  vector[m] d_sd = sqrt(d_var);    // sampling standard deviations

  // Compute QR decomposition of X for numerical stability
  matrix[m, p] Q_X = qr_thin_Q(X) * sqrt(m - 1);
  matrix[p, p] R_X = qr_thin_R(X) / sqrt(m - 1);
  matrix[p, p] R_X_inv = inverse(R_X);
}

parameters {
  vector[p] beta_raw;              // regression coefficients (QR scale)
  vector[r] eta_raw;               // spatial basis coefficients (raw scale if NCP)
  vector[m] w_raw;                 // IID random effects (raw scale if NCP)

  real<lower=0> sigma_sq;          // IID variance
  real<lower=0> tau_sq;            // spatial variance
}

transformed parameters {
  vector[p] beta = R_X_inv * beta_raw; // back-transform to original scale
  vector[r] eta;
  vector[m] w;
  vector[m] v;                     // spatial random effects v = S * eta
  vector[m] theta;                 // area parameters

  if (use_ncp == 1) {
    // Non-centered parameterization
    eta = sqrt(tau_sq) * eta_raw;
    w = sqrt(sigma_sq) * w_raw;
  } else {
    // Centered parameterization
    eta = eta_raw;
    w = w_raw;
  }

  v = S * eta;                     // compute spatial effects
  theta = X * beta + v + w;        // compute area parameters
}

model {
  // Priors on variance components (IG in rate parameterization)
  // Stan uses rate parameterization: p(x) ∝ x^{-α-1} exp(-β/x)
  // So IG(shape=c, rate=d) has mean d/(c-1) for c > 1
  target += -prior_sigma_shape * log(sigma_sq) - prior_sigma_scale / sigma_sq;
  target += -prior_tau_shape * log(tau_sq) - prior_tau_scale / tau_sq;

  // Flat prior on beta (implicit via beta_raw ~ uniform)

  // Prior on random effects
  if (use_ncp == 1) {
    eta_raw ~ std_normal();
    w_raw ~ std_normal();
  } else {
    eta ~ normal(0, sqrt(tau_sq));
    w ~ normal(0, sqrt(sigma_sq));
  }

  // Likelihood
  y ~ normal(theta, d_sd);
}

generated quantities {
  vector[m] log_lik;               // pointwise log-likelihood
  vector[m] y_rep;                 // posterior predictive draws

  for (i in 1:m) {
    log_lik[i] = normal_lpdf(y[i] | theta[i], d_sd[i]);
    y_rep[i] = normal_rng(theta[i], d_sd[i]);
  }
}
