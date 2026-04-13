library(LaplacesDemon)
library(tidyverse)
library(sampling)
library(survey)

# ==================================================
# Helper function to create response variables
# ==================================================
create_response_variable <- function(df, response_var) {
  # Returns data frame with response variable added
  # Handles both continuous (already in data) and binary (derived) variables

  if (response_var %in% names(df)) {
    # Continuous variable already in data (WAGP, POVPIP, AGEP)
    return(df)
  }

  # Binary derived variables
  if (response_var == "rents_home") {
    df$rents_home <- df$TEN == 3
  } else if (response_var == "below_poverty") {
    df$below_poverty <- as.numeric(df$POVPIP < 138)
  } else if (response_var == "hispanic") {
    df$hispanic <- as.numeric(df$HISP != "01")
  } else if (response_var == "employed") {
    df$employed <- as.numeric(df$ESR %in% c("1", "2", "4", "5"))
  } else if (response_var == "college_plus") {
    df$college_plus <- as.numeric(as.numeric(df$SCHL) >= 21)
  } else if (response_var == "citizen") {
    df$citizen <- as.numeric(df$CIT == "1")
  } else if (response_var == "disabled") {
    df$disabled <- as.numeric(df$DIS == "1")
  } else if (response_var == "elderly") {
    df$elderly <- as.numeric(df$AGEP >= 65)
  } else if (response_var == "young_adult") {
    df$young_adult <- as.numeric(df$AGEP >= 18 & df$AGEP < 35)
  } else if (response_var == "middle_age") {
    df$middle_age <- as.numeric(df$AGEP >= 35 & df$AGEP < 55)
  } else if (response_var == "female") {
    df$female <- as.numeric(df$SEX == "2")
  } else if (response_var == "asian") {
    df$asian <- as.numeric(df$RAC1P == "6")
  } else if (response_var == "black") {
    df$black <- as.numeric(df$RAC1P == "2")
  } else if (response_var == "white") {
    df$white <- as.numeric(df$RAC1P == "1")
  } else if (response_var == "married") {
    df$married <- as.numeric(df$MAR == "1")
  } else if (response_var == "long_commute") {
    df$long_commute <- as.numeric(df$JWMNP > 45)
  } else if (response_var == "pubcov") {
    df$pubcov <- as.numeric(df$PUBCOV == "1")
  } else {
    stop(paste("Unknown response variable:", response_var))
  }

  return(df)
}

# ==================================================
# Core function to generate the data and set things up
# ==================================================

setup_comp <- function(ncomps, results_dir, model_config, start_from = 1, comp_nos = NULL, seed_offset = 0) {
  # Load population data (REQUIRED - no defaults)
  if (is.null(model_config$population_file)) {
    stop("model_config$population_file is required! Specify the population data file.")
  }
  acs_pop <- readRDS(model_config$population_file)

  # Get parameters from model_config
  samp_frac <- model_config$samp_frac
  response_var <- model_config$response_var
  response_filter <- model_config$response_filter
  X_approach <- "population" # Use population X for both contenders

  # Sampling parameters (with defaults if not specified)
  min_sample_size <- if (!is.null(model_config$min_sample_size)) model_config$min_sample_size else 30
  equal_allocation <- if (!is.null(model_config$equal_allocation)) model_config$equal_allocation else FALSE
  equal_n <- if (!is.null(model_config$equal_n)) model_config$equal_n else 50

  # set up the sim_files folder containing all the results
  base_folder <- file.path(getwd(), results_dir)
  dir.create(base_folder, showWarnings = FALSE)

  # Determine which variables to exclude from X based on response
  exclude_from_X <- NULL
  if (response_var == "rents_home") {
    exclude_from_X <- c("homeownership_rate")
  }
  # WAGP: no exclusions needed

  # Get covariates from model_config
  X_covariates <- model_config$X_covariates

  if (X_approach == "population") {
    X_pop <- create_X(
      approach = "population", pop_data = acs_pop,
      covariates = X_covariates, exclude_vars = exclude_from_X
    )
  }

  # Compute and save theta_true if not already present
  if (!file.exists(file.path(base_folder, "theta_true.RDS"))) {
    # Compute true population means (theta_true) - ONCE for entire study
    # Apply same response variable computation and filter as in sampling
    acs_pop_for_truth <- create_response_variable(acs_pop, response_var)

    # Apply response filter if specified
    if (!is.null(response_filter) && response_filter != "") {
      acs_pop_for_truth <- acs_pop_for_truth %>% filter(eval(parse(text = response_filter)))
    }

    # Compute true population means by PUMA
    theta_true_df <- acs_pop_for_truth %>%
      group_by(PUMA) %>%
      summarise(theta_true = mean(.data[[response_var]], na.rm = TRUE), .groups = "drop") %>%
      arrange(PUMA)

    # Save theta_true to results directory (once, not per comparison)
    saveRDS(
      list(puma = theta_true_df$PUMA, values = theta_true_df$theta_true),
      file.path(base_folder, "theta_true.RDS")
    )
  }

  # main loop - each comparison has unique z & d via sampling
  indices <- if (!is.null(comp_nos)) comp_nos else start_from:ncomps
  for (sim in indices) {
    # create simulation folder
    sim_folder <- file.path(base_folder, sprintf("comparison_%03d", sim))
    dir.create(sim_folder, showWarnings = FALSE)

    # create dt_1fold, dt_5fold, emp_sim, & fit_on_z  within each simulation folder
    dir.create(file.path(sim_folder, "dt_1fold"), showWarnings = FALSE)
    dir.create(file.path(sim_folder, "dt_5fold"), showWarnings = FALSE)
    dir.create(file.path(sim_folder, "emp_sim"), showWarnings = FALSE)
    dir.create(file.path(sim_folder, "fit_on_z"), showWarnings = FALSE)

    # Draw sample
    set.seed(sim + seed_offset)
    sample <- get_strat_PS(
      pop_df = acs_pop,
      samp_frac = samp_frac,
      min_sample_size = min_sample_size,
      equal_allocation = equal_allocation,
      equal_n = equal_n
    )
    sample_df <- acs_pop[sample$idx, ]
    sample_df$design_weight <- sample$weights

    # Compute response variable (handles both continuous and binary)
    sample_df <- create_response_variable(sample_df, response_var)

    # Apply response filter if specified (e.g., "!is.na(WAGP) & WAGP >= 0")
    if (!is.null(response_filter) && response_filter != "") {
      sample_df <- sample_df %>% filter(eval(parse(text = response_filter)))
    }

    # Get direct estimates & design-based variance for response
    sample_design <- svydesign(ids = ~1, weights = ~design_weight, data = sample_df)
    direst <- svyby(as.formula(paste0("~", response_var)), ~PUMA, sample_design, svymean, vartype = "var", na.rm = TRUE) %>%
      arrange(PUMA) # Explicit sort to ensure alignment with X_pop

    # Extract y and d_var correctly for both binary and continuous variables
    # For logical/binary: svyby creates "varnameTRUE"/"varnameFALSE" and "var.varnameTRUE"/"var.varnameFALSE"
    # For continuous: svyby creates "varname" and "var"
    all_names <- names(direst)
    response_cols <- setdiff(all_names, c("PUMA", all_names[grepl("^var", all_names)]))

    if (length(response_cols) > 1) {
      # Binary: multiple columns, select TRUE
      response_col <- response_cols[grepl("TRUE$", response_cols)]
      var_col <- paste0("var.", response_col)
    } else {
      # Continuous: single column
      response_col <- response_cols[1]
      var_col <- "var"
    }

    puma_ids <- direst$PUMA
    z_vals <- direst[[response_col]]
    d_vals <- direst[[var_col]]

    # Get X based on approach
    if (X_approach == "population") {
      X_to_save <- X_pop # computed once above
    } else { # "estimated"
      X_to_save <- create_X(
        approach = "estimated", sample_design = sample_design,
        covariates = X_covariates, exclude_vars = exclude_from_X
      )
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
get_PS <- function(
    pop_df, samp_frac = .0225, min_sample_size = 30,
    equal_allocation = FALSE, equal_n = 50) {
  if (equal_allocation) {
    sample_size <- equal_n
  } else {
    sample_size <- max(floor(nrow(pop_df) * samp_frac), min_sample_size)
    # ensure minimum for stable variance estimation
  }

  size_var <- pop_df$PWGTP

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
# samp_frac: either a scalar (same rate for all PUMAs) or a named vector keyed by PUMA id
#            (per-PUMA rates for mixed precision designs).
get_strat_PS <- function(pop_df, samp_frac = 0.01, min_sample_size = 30,
                         equal_allocation = FALSE, equal_n = 50) {
  PUMAs <- unique(pop_df$PUMA)
  weights <- idx <- c()
  for (PUMA in PUMAs) {
    puma_ids <- which(pop_df$PUMA == PUMA)
    frac_i <- if (length(samp_frac) == 1) samp_frac else samp_frac[[as.character(PUMA)]]
    tmp <- get_PS(
      pop_df = pop_df[puma_ids, ],
      samp_frac = frac_i,
      min_sample_size = min_sample_size,
      equal_allocation = equal_allocation,
      equal_n = equal_n
    )
    weights <- c(weights, tmp$weights)
    idx <- c(idx, puma_ids[tmp$idx]) # Map local indices back to full dataset
  }
  return(list(weights = weights, idx = idx))
}

# ==================================================
# Functions to create X (from population or direct estimates)
# ==================================================

# create X based on population data or sample
# covariates: NULL for intercept-only, or character vector of covariate names to include
# exclude_vars: optional character vector of variable names to exclude from X
create_X <- function(approach, pop_data = NA, sample_design = NA, covariates = NULL, exclude_vars = NULL) {
  # If covariates = NULL, return intercept-only (just PUMA column, no predictors)
  if (is.null(covariates)) {
    if (approach == "population") {
      X <- pop_data %>%
        group_by(PUMA) %>%
        summarise(.groups = "drop") %>%
        arrange(PUMA)
    } else if (approach == "estimated") {
      # Get unique PUMAs from sample
      puma_df <- data.frame(PUMA = unique(sample_design$variables$PUMA)) %>%
        arrange(PUMA)
      X <- puma_df
    }
    return(X)
  }

  # Otherwise, compute all covariates then select requested ones
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
        # Housing characteristics (for pilot testing)
        mean_rooms = mean(BDSP, na.rm = TRUE),
        pct_crowded = mean(BDSP < 2, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(PUMA)
  } else if (approach == "estimated") {
    # Get direct estimates for all predictors using survey design
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
    # Housing characteristics
    mean_rooms <- svyby(~BDSP, ~PUMA, sample_design, svymean, na.rm = TRUE)
    pct_crowded <- svyby(~ I(BDSP < 2), ~PUMA, sample_design, svymean, na.rm = TRUE)

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
      has_hicov = has_hicov$`I(HICOV == 1)`,
      mean_rooms = mean_rooms$BDSP,
      pct_crowded = pct_crowded$`I(BDSP < 2)`
    ) %>%
      arrange(PUMA)
  } else {
    stop("Pick one approach: `population` or `estimated`")
  }

  # Select only requested covariates (if covariates specified)
  if (!is.null(covariates)) {
    missing_covs <- setdiff(covariates, names(X))
    if (length(missing_covs) > 0) {
      stop("Requested covariates not available: ", paste(missing_covs, collapse = ", "))
    }
    X <- X %>% select(PUMA, all_of(covariates))
  }

  # Exclude specified variables if any
  if (!is.null(exclude_vars)) {
    X <- X %>% select(-any_of(exclude_vars))
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
setup_esim <- function(ncomps, results_dir, niters = 100, start_from = 1, comp_nos = NULL) {
  indices <- if (!is.null(comp_nos)) comp_nos else start_from:ncomps
  lapply(indices, function(k) {
    esim_helper(comp_no = k, niters = niters, results_dir = results_dir)
  })
  return("Setup complete!")
}
