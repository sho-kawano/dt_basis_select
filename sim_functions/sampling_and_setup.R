library(LaplacesDemon)
library(tidyverse)
library(sampling)

# ==================================================
# Sampling functions used to create the direct estimates
# ==================================================

# Performs probability proportional to size (PPS) sampling from a population dataframe
get_PS <- function(pop_df, type = "PPS", samp_frac = .01) { # type is PPS or some other size variable
  sample_size <- max(floor(nrow(pop_df) * samp_frac), 20) # need to ensure a minimum
  # create the size variable that will be sampled in proportion to
  # using PWGTP directly - people with higher weights represent more people
  size_var <- as.numeric(pop_df$PWGTP)
  inclusion_probs <- inclusionprobabilities(size_var, sample_size)
  inclusion_probs <- inclusion_probs / sum(inclusion_probs) * sample_size
  # survey weights are inverse probabilities of selection
  weights <- 1 / inclusion_probs
  # Draw Poisson sample: https://en.wikipedia.org/wiki/Poisson_sampling
  sample_idx <- which(UPpoisson(inclusion_probs) == 1)
  sample_size <- length(sample_idx)
  weights <- weights[sample_idx]

  return(list(weights = weights, idx = sample_idx))
}

# Performs stratified PPS sampling (using `get_PS`) across PUMA (Public Use Microdata Area)
get_strat_PS <- function(pop_df, samp_frac) {
  PUMAs <- unique(pop_df$PUMA)
  weights <- idx <- c()
  for (PUMA in PUMAs) {
    puma_ids <- which(pop_df$PUMA == PUMA)
    tmp <- get_PS(pop_df = pop_df[puma_ids, ], samp_frac = samp_frac)
    weights <- c(weights, tmp$weights)
    idx <- c(idx, tmp$idx)
  }
  return(list(weights = weights, idx = idx))
}


# ==================================================
# Core function to set up the simulation study
# ==================================================

setup_comp <- function(z_mean, z_var, ncomps, results_dir) {
  # set up the sim_files folder containing all the results
  base_folder <- file.path(getwd(), results_dir)
  dir.create(base_folder, showWarnings = FALSE)

  # main loop
  for (sim in 1:ncomps) {
    # create simulation folder
    sim_folder <- file.path(base_folder, sprintf("comparison_%03d", sim))
    dir.create(sim_folder, showWarnings = FALSE)

    # create dt_1fold, dt_5fold, emp_sim, & fit_on_z  within each simulation folder
    dir.create(file.path(sim_folder, "dt_1fold"), showWarnings = FALSE)
    dir.create(file.path(sim_folder, "dt_5fold"), showWarnings = FALSE)
    dir.create(file.path(sim_folder, "emp_sim"), showWarnings = FALSE)
    dir.create(file.path(sim_folder, "fit_on_z"), showWarnings = FALSE)

    # generates 'z' - the direct estimates
    set.seed(sim)
    # draw samples
    sample <- get_strat_PS(pop_df = acs_pop, samp_frac = .002)
    sample_df <- acs_pop[sample$idx, ]
    sample_df$design_weight <- sample$weights

    # get direct estimate & design-based variance
    sample_design <- svydesign(ids = ~1, weights = ~design_weight, data = sample_df)
    direst <- svyby(~HICOV, ~PUMA, sample_design, svymean, vartype = "var")
    z <- direst$HICOV
    d <- direst$var

    # save results
    saveRDS(z, file = file.path(sim_folder, "z.RDS"))
    saveRDS(d, file = file.path(sim_folder, "d.RDS"))
  }
}


# ==================================================
# Function to set up the Empirical Simulation Study
# ==================================================

# NOTE: niters is the number of iteration per empirical simulation study
esim_helper <- function(comp_no, w_var, niters, results_dir) {
  # set up the folder
  base_folder <- file.path(getwd(), results_dir)
  sim_folder <- file.path(base_folder, sprintf("comparison_%03d", comp_no))

  # load direct estimates / truth for a given emp.study
  z <- readRDS(file.path(sim_folder, "z.RDS")) %>% as.numeric()

  # generate synthetic data for one simulation study
  for (iter in 1:niters) {
    sub_folder <- sprintf("%03d", iter)
    sub_path <- file.path(sim_folder, "emp_sim", sub_folder)
    dir.create(sub_path, recursive = TRUE, showWarnings = FALSE)
    set.seed(iter * 2)
    w <- rmvn(mu = z, Sigma = diag(w_var))
    # save results
    saveRDS(w, file = file.path(sub_path, "w.RDS"))
  }
}

# iterates over each comparison
setup_esim <- function(w_var, ncomps, results_dir) {
  lapply(1:ncomps, function(k) {
    esim_helper(comp_no = k, w_var = w_var, niters = 100, results_dir = results_dir)
  })
  return("Setup complete!")
}
