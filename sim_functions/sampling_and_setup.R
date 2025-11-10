library(LaplacesDemon)
library(tidyverse)
library(sampling)
library(survey)

# ==================================================
# Core function to generate the data and set things up
# ==================================================

setup_comp <- function(ncomps, results_dir) {
  # Load config to get population data and sampling parameters
  sim_config <- readRDS("sim_config.RDS")
  acs_pop <- sim_config$acs_pop
  samp_frac <- sim_config$samp_frac
  response_var <- sim_config$response_var
  X_approach <- sim_config$X_approach

  # set up the sim_files folder containing all the results
  base_folder <- file.path(getwd(), results_dir)
  dir.create(base_folder, showWarnings = FALSE)

  if (X_approach == "population") {
    X_pop <- create_X(approach = "population", pop_data = acs_pop)
  }
  # main loop - each comparison has unique z & d via sampling
  for (sim in 1:ncomps) {
    # create simulation folder
    sim_folder <- file.path(base_folder, sprintf("comparison_%03d", sim))
    dir.create(sim_folder, showWarnings = FALSE)

    # create dt_1fold, dt_5fold, emp_sim, & fit_on_z  within each simulation folder
    dir.create(file.path(sim_folder, "dt_1fold"), showWarnings = FALSE)
    dir.create(file.path(sim_folder, "dt_5fold"), showWarnings = FALSE)
    dir.create(file.path(sim_folder, "emp_sim"), showWarnings = FALSE)
    dir.create(file.path(sim_folder, "fit_on_z"), showWarnings = FALSE)

    # Draw sample
    set.seed(sim)
    sample <- get_strat_PS(pop_df = acs_pop, samp_frac = samp_frac)
    sample_df <- acs_pop[sample$idx, ]
    sample_df$design_weight <- sample$weights

    # Compute medi_cal_qualified if using that as response
    if (response_var == "medi_cal_qualified") {
      sample_df$medi_cal_qualified <- ifelse(sample_df$POVPIP < 138, 1, 0)
    }

    # Get direct estimates & design-based variance for response
    sample_design <- svydesign(ids = ~1, weights = ~design_weight, data = sample_df)
    direst <- svyby(as.formula(paste0("~", response_var)), ~PUMA, sample_design, svymean, vartype = "var", na.rm = TRUE) %>%
      arrange(PUMA)  # Explicit sort to ensure alignment with X_pop

    puma_ids <- direst$PUMA
    z_vals <- direst[[response_var]]
    d_vals <- direst$var

    # Get X based on approach
    if (X_approach == "population") {
      X_to_save <- X_pop # computed once above
    } else { # "estimated"
      X_to_save <- create_X(approach = "estimated", sample_design = sample_design)
    }

    # Save results
    saveRDS(list(puma = puma_ids, values = z_vals), file.path(sim_folder, "z.RDS"))
    saveRDS(list(puma = puma_ids, values = d_vals), file.path(sim_folder, "d.RDS"))
    saveRDS(list(puma = X_to_save$PUMA, X = X_to_save %>% select(-PUMA)), file.path(sim_folder, "X.RDS"))
  }
}


# ==================================================
# Data Generation: Sampling functions used to create the direct estimates
# ==================================================

# Performs probability proportional to size (PPS) sampling from a population dataframe
get_PS <- function(pop_df, samp_frac = .0225, min_sample_size = 30) {
  sample_size <- max(floor(nrow(pop_df) * samp_frac), min_sample_size) # ensure minimum for stable variance estimation

  size_var <- pop_df$PWGTP
  # trying an informative sampling design
  # size_var <- as.numeric(exp(.2 * scale(-pop_df$WAGP) + .4 * scale(pop_df$PWGTP)))

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
    idx <- c(idx, puma_ids[tmp$idx]) # Map local indices back to full dataset
  }
  return(list(weights = weights, idx = idx))
}

# ==================================================
# Functions to create X (from population or direct estimates)
# ==================================================

# create X based on population data or sample
create_X <- function(approach, pop_data = NA, sample_design = NA) {
  if (approach == "population") {
    X <- pop_data %>%
      group_by(PUMA) %>%
      summarise(
        mean_age = mean(AGEP, na.rm = TRUE),
        pct_male = mean(SEX == 1, na.rm = TRUE),
        pct_white = mean(RAC1P == 1, na.rm = TRUE),
        pct_black = mean(RAC1P == 2, na.rm = TRUE),
        pct_asian = mean(RAC1P == 6, na.rm = TRUE),
        pct_hispanic = mean(HISP != "01", na.rm = TRUE),
        pct_married = mean(MAR == 1, na.rm = TRUE),
        pct_citizen = mean(CIT %in% c(1, 2, 3, 4), na.rm = TRUE),
        employment_rate = mean(ESR == 1, na.rm = TRUE),
        pct_bachelor = mean(SCHL >= 21, na.rm = TRUE),
        homeownership_rate = mean(TEN %in% c(1, 2), na.rm = TRUE),
        disability_rate = mean(DIS == 1, na.rm = TRUE),
        has_hicov = mean(HICOV == 1, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(PUMA)
  } else if (approach == "estimated") {
    # Get direct estimates for all 13 predictors using survey design
    mean_age <- svyby(~AGEP, ~PUMA, sample_design, svymean, na.rm = TRUE)
    pct_male <- svyby(~ I(SEX == 1), ~PUMA, sample_design, svymean, na.rm = TRUE)
    pct_white <- svyby(~ I(RAC1P == 1), ~PUMA, sample_design, svymean, na.rm = TRUE)
    pct_black <- svyby(~ I(RAC1P == 2), ~PUMA, sample_design, svymean, na.rm = TRUE)
    pct_asian <- svyby(~ I(RAC1P == 6), ~PUMA, sample_design, svymean, na.rm = TRUE)
    pct_hispanic <- svyby(~ I(HISP != "01"), ~PUMA, sample_design, svymean, na.rm = TRUE)
    pct_married <- svyby(~ I(MAR == 1), ~PUMA, sample_design, svymean, na.rm = TRUE)
    pct_citizen <- svyby(~ I(CIT %in% c(1, 2, 3, 4)), ~PUMA, sample_design, svymean, na.rm = TRUE)
    employment_rate <- svyby(~ I(ESR == 1), ~PUMA, sample_design, svymean, na.rm = TRUE)
    pct_bachelor <- svyby(~ I(SCHL >= 21), ~PUMA, sample_design, svymean, na.rm = TRUE)
    homeownership_rate <- svyby(~ I(TEN %in% c(1, 2)), ~PUMA, sample_design, svymean, na.rm = TRUE)
    disability_rate <- svyby(~ I(DIS == 1), ~PUMA, sample_design, svymean, na.rm = TRUE)
    has_hicov <- svyby(~ I(HICOV == 1), ~PUMA, sample_design, svymean, na.rm = TRUE)

    # Combine into data frame
    X <- data.frame(
      PUMA = mean_age$PUMA,
      mean_age = mean_age$AGEP,
      pct_male = pct_male$`I(SEX == 1)`,
      pct_white = pct_white$`I(RAC1P == 1)`,
      pct_black = pct_black$`I(RAC1P == 2)`,
      pct_asian = pct_asian$`I(RAC1P == 6)`,
      pct_hispanic = pct_hispanic$`I(HISP != "01")`,
      pct_married = pct_married$`I(MAR == 1)`,
      pct_citizen = pct_citizen$`I(CIT %in% c(1, 2, 3, 4))`,
      employment_rate = employment_rate$`I(ESR == 1)`,
      pct_bachelor = pct_bachelor$`I(SCHL >= 21)`,
      homeownership_rate = homeownership_rate$`I(TEN %in% c(1, 2))`,
      disability_rate = disability_rate$`I(DIS == 1)`,
      has_hicov = has_hicov$`I(HICOV == 1)`
    ) %>%
      arrange(PUMA)
  } else {
    stop("Pick one approach: `population` or `estimated`")
  }
  return(X)
}

# ==================================================
# Function to set up the Empirical Simulation Study
# ==================================================

# NOTE: niters is the number of iteration per empirical simulation study
esim_helper <- function(comp_no, niters, results_dir) {
  # set up the folder
  base_folder <- file.path(getwd(), results_dir)
  sim_folder <- file.path(base_folder, sprintf("comparison_%03d", comp_no))

  # load direct estimates and variance for this comparison (now saved as lists)
  z_data <- readRDS(file.path(sim_folder, "z.RDS"))
  d_data <- readRDS(file.path(sim_folder, "d.RDS"))

  z <- z_data$values
  d <- d_data$values

  # generate synthetic data for one simulation study
  for (iter in 1:niters) {
    sub_folder <- sprintf("%03d", iter)
    sub_path <- file.path(sim_folder, "emp_sim", sub_folder)
    dir.create(sub_path, recursive = TRUE, showWarnings = FALSE)
    set.seed(iter * 2)
    # Use comparison-specific variance d (not global w_var)
    w <- rmvn(mu = z, Sigma = diag(d))
    # save results
    saveRDS(w, file = file.path(sub_path, "w.RDS"))
  }
}

# iterates over each comparison
setup_esim <- function(ncomps, results_dir) {
  lapply(1:ncomps, function(k) {
    esim_helper(comp_no = k, niters = 100, results_dir = results_dir)
  })
  return("Setup complete!")
}
