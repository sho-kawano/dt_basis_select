# Architecture Change Analysis: Fixed Data → Design-Based Sampling

## Summary of Change

**OLD APPROACH (from previous codebase):**
- Fixed `chosen_var` and `truth` passed via `sim_config.RDS` to all functions
- Fixed covariate matrix `X` loaded from external data file
- Same data used across all comparisons

**NEW APPROACH (this implementation):**
- Each comparison = independent sample from `ca_pums_population.rds`
- Each comparison generates its own `X` (predictors), `z` (direct estimates), `d` (variances)
- All three (X, z, d) come from the SAME sampled units per comparison
- Saved to `comparison_XXX/` folders
- **NO `create_puma_data.R`** - we don't pre-aggregate the population!

**Key Insight:** In a real survey, you don't know the "true" area characteristics - you estimate them from your sample. This approach mimics that: both response AND predictors are sample-based.

**Note on file naming:**
- **Full data fit files:** `sim_functions/full_data_fit.R` (fits models on full observed data) and `sim_functions/summary_oracle.R` (calculates MSE/WAIC/DIC)
- Renamed from `zfit.R` → `full_data_fit.R` for clarity

---

## PART 1: Data Setup in run_comparisons.R (Lines 36-111)

**KEY INSIGHT:** Data generation now happens in `sampling_and_setup.R`, not here. This section just loads population data and creates config.

### NEW CODE NEEDED:
```r
# Load population data (individual-level PUMS data)
acs_pop <- readRDS("data/ca_pums_population.rds")

# Define which predictors to aggregate from sampled units
predictor_names <- c(
  # Demographics (8)
  "mean_age", "pct_male", "pct_white", "pct_black",
  "pct_asian", "pct_hispanic", "pct_married", "pct_citizen",
  # Economic (1)
  "employment_rate",
  # Education (1)
  "pct_bachelor",
  # Housing (1)
  "homeownership_rate",
  # Health (2)
  "disability_rate", "has_hicov"
)

# sim_config NO LONGER contains:
# - truth, chosen_var (never did in new approach)
# - all_data (deleted - no pre-aggregated data!)
# - all_covs (replaced with predictor_names)

sim_config <- list(
  acs_pop = acs_pop,                  # Individual-level population for sampling
  samp_frac = 0.0225,                 # Sampling fraction (2.25% - tested optimal)
  response_var = "medi_cal_qualified", # Response variable name
  X_approach = "population",          # How to compute X (population or estimated)
  predictor_names = predictor_names,  # Which variables to aggregate per comparison
  n_comp = n_comp,
  n_cores = n_cores,
  results_dir = "_results"
)
```

**What happened to `all_data` and `create_puma_data.R`?**
- ❌ **DELETED** - We don't pre-aggregate the population!
- Each comparison generates its own X matrix from its own sample
- X varies across comparisons (just like z and d)
- This is more realistic: in a real survey, predictors are also sample-based estimates

---

## PART 2: Files That Need Changes

**IMPORTANT:** Implement in this order:
1. **Start with 2.2** (`sampling_and_setup.R`) - This is where actual data generation happens
2. Then 2.1 (`run_comparisons.R`) - Update function calls to remove variance parameters
3. Then 2.3-2.8 - Update all downstream functions

### 2.1 `run_comparisons.R`

**Lines 36-111: Data Generation & Config**
- ❌ REMOVE: `all_data <- readRDS("data/puma_aggregated.rds")` (file doesn't exist!)
- ❌ REMOVE: `all_covs` definition (replaced with `predictor_names`)
- ✅ ADD: Load `acs_pop`, define `predictor_names`
- ✅ MODIFY: `sim_config` structure (see Part 1 above)

**Line 135: setup_comp call**
```r
# OLD:
setup_comp(z_mean = truth, z_var = chosen_var, ncomps = n_comp, results_dir)

# NEW:
setup_comp(ncomps = n_comp, results_dir = sim_config$results_dir)
```

**Line 136: setup_esim call**
```r
# OLD:
setup_esim(w_var = chosen_var, ncomps = n_comp, results_dir)

# NEW:
setup_esim(ncomps = n_comp, results_dir = sim_config$results_dir)
```

**Line 164: full_data_fit call** (renamed from zfit)
```r
# OLD:
zfit(j, chosen_var, all_covs, all_data, sim_config$results_dir)

# NEW:
full_data_fit(j, sim_config$results_dir)
# Function loads X, z, d from comparison_XXX/ folder
# No more all_covs or all_data parameters!
# FOR NOW: Comment out this call - model fitting not ready yet
```

**Line 188: summary_oracle call**
```r
# OLD:
summary_oracle(x, truth, chosen_var, all_covs, sim_config$results_dir)

# NEW:
summary_oracle(x, sim_config$results_dir)
# Function loads X, z, d from comparison_XXX/ folder
# FOR NOW: Comment out this call - model fitting not ready yet
```

**Line 258: run_dt call**
```r
# OLD:
run_dt(comp_no = comp, k = 1, s = chosen_var, all_covs = all_covs, ...)

# NEW:
run_dt(comp_no = comp, k = 1, ...)
# Function loads X, d from comparison_XXX/ folder
# all_covs parameter removed - X already has correct columns
```

**Line 342: run_esim call**
```r
# OLD:
run_esim(comp, s = chosen_var, all_covs, ...)

# NEW:
run_esim(comp, ...)
# Function loads X, d from comparison_XXX/ folder
```

**Lines 393, 431, etc: summary_dt calls**
```r
# OLD:
summary_dt(x, n_folds = 1, cfg$chosen_var, ...)

# NEW:
summary_dt(x, n_folds = 1, ...)
# Function loads X, d from comparison_XXX/ folder
```

---

### 2.2 `sim_functions/sampling_and_setup.R` ⭐ START HERE - MAJOR CHANGES

**This is where the REAL data generation happens!** Each comparison gets unique X, z, and d via design-based sampling.

**Key functions to update:**

1. **`get_PS` function** - Looking into UPsystematic vs UPpoisson
```r
# CURRENT: Uses UPpoisson (line 20)
sample_idx <- which(UPpoisson(inclusion_probs) == 1)  # Random sample size

# ALTERNATIVE: UPsystematic (fixed-size sampling)
sample_idx <- UPsystematic(inclusion_probs)  # Guaranteed sample size

# Trade-offs under consideration:
# - UPpoisson: Can return 0 units (needs retry loop), but unbiased
# - UPsystematic: Fixed size, no retries, still design-unbiased
# For now: Keep current approach with retry loop in setup_comp()
```

2. **`setup_comp` function** - MAJOR REWRITE NEEDED
```r
# OLD signature:
setup_comp <- function(z_mean, z_var, ncomps, results_dir)

# NEW signature (same):
setup_comp <- function(ncomps, results_dir)

# What it now does:
# 1. Load acs_pop, predictor_names, response_var from sim_config.RDS
# 2. For each comparison (1:ncomps):
#    a. Draw stratified PPS sample via get_strat_PS() [FIXED SIZE]
#    b. Create survey design object with weights
#    c. Aggregate sample to PUMA level to create:
#       - X: predictor matrix (13 variables from predictor_names)
#       - z: direct estimates for response_var
#       - d: design-based variance
#    d. Check for zero sampled units in any PUMA (error if found)
#    e. Save X.RDS, z.RDS, d.RDS to comparison_XXX/
```

**New aggregation logic needed in setup_comp:**
```r
# After drawing sample (line ~78):
sample_df <- acs_pop[sample$idx, ]
sample_df$design_weight <- sample$weights

# Create predictor aggregations at PUMA level
predictor_names <- sim_config$predictor_names
X_agg <- sample_df %>%
  group_by(PUMA) %>%
  summarise(
    n = n(),  # Check for zero counts
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
    has_hicov = mean(HICOV == 1, na.rm = TRUE)
  )

# Check for PUMAs with zero samples
if (any(X_agg$n == 0)) {
  stop(sprintf("Comparison %d: %d PUMAs have zero sampled units",
               sim, sum(X_agg$n == 0)))
}

# Get z and d using survey design
sample_design <- svydesign(ids = ~1, weights = ~design_weight, data = sample_df)
formula_str <- paste0("~", response_var)
direst <- svyby(as.formula(formula_str), ~PUMA, sample_design, svymean, vartype = "var")

# Extract with PUMA alignment
puma_ids <- direst$PUMA
z_vals <- direst[[response_var]]
d_vals <- direst$var

# Ensure X_agg is in same order as z/d
X_agg <- X_agg %>% arrange(PUMA)
if (!all(X_agg$PUMA == puma_ids)) {
  stop("PUMA mismatch between X and z/d")
}

# Save with PUMA IDs
saveRDS(list(puma = puma_ids, values = z_vals), file.path(sim_folder, "z.RDS"))
saveRDS(list(puma = puma_ids, values = d_vals), file.path(sim_folder, "d.RDS"))
saveRDS(list(puma = puma_ids, X = X_agg %>% select(-PUMA, -n)),
        file.path(sim_folder, "X.RDS"))
```

3. **`esim_helper` function** - No changes needed (already loads d from folder)

4. **`setup_esim` function** - No changes needed

**Key statistical concepts:**
- PPS sampling: Probability proportional to size (using PWGTP weights)
- Stratified by PUMA: Each PUMA gets its own sample
- X, z, d all from SAME sample (no alignment issues!)
- Design-based variance via `survey` package
- Sampling method (Poisson vs systematic) still under evaluation

---

### 2.3 `sim_functions/full_data_fit.R` (renamed from zfit.R)

**Note:** This file fits models on the FULL observed data (X, z) for each comparison. Currently being rewritten for new spatial model.

**OLD signature:**
```r
zfit <- function(comp_no, chosen_var, all_covs, all_data, results_dir)
```

**NEW signature:**
```r
full_data_fit <- function(comp_no, results_dir)
```

**Inside function:**
```r
# Load comparison-specific data
comp_folder <- file.path(getwd(), results_dir, sprintf("comparison_%03d", comp_no))

# Load as lists with PUMA IDs
z_data <- readRDS(file.path(comp_folder, "z.RDS"))
d_data <- readRDS(file.path(comp_folder, "d.RDS"))
X_data <- readRDS(file.path(comp_folder, "X.RDS"))

# Extract values (already aligned by PUMA)
z <- z_data$values
d <- d_data$values
X_df <- X_data$X  # Data frame with predictor columns

# Create model matrix (adds intercept)
X <- model.matrix(~., X_df)

# TODO: Fit spatial basis model here (not implemented yet)
# For now, this function should be commented out in run_comparisons.R
```

**FOR NOW:** Comment out all calls to this function in run_comparisons.R - model code not ready yet

---

### 2.4 `sim_functions/summary_oracle.R` (Summary for Full Data Fits)

**Note:** This file calculates the Observed-Data Oracle (OD-Oracle) MSE, plus WAIC and DIC, for models fit on full data `z`. The "oracle" here means we know the true small area means and can calculate true MSE.

**OLD signature:**
```r
summary_oracle <- function(comp_no, truth, chosen_var, all_covs, results_dir)
```

**NEW signature:**
```r
summary_oracle <- function(comp_no, all_covs, results_dir)
```

**Inside function:**
```r
# ADD at start:
comp_folder <- file.path(getwd(), results_dir, sprintf("comparison_%03d", comp_no))
z <- readRDS(file.path(comp_folder, "z.RDS"))  # The "truth" for this comparison
d <- readRDS(file.path(comp_folder, "d.RDS"))

# Use z instead of truth
# Use d instead of chosen_var
```

---

### 2.5 `sim_functions/run_dt.R`

**OLD signature:**
```r
run_dt <- function(comp_no, k, s, all_covs, all_data, results_dir, eps = NULL, n_reps = NULL)
```

**NEW signature:**
```r
run_dt <- function(comp_no, k, results_dir, eps = NULL, n_reps = NULL)
```

**Inside function:**
```r
# ADD at start:
sim_folder <- file.path(results_dir, sprintf("comparison_%03d", comp_no))

# Load comparison-specific data
X_data <- readRDS(file.path(sim_folder, "X.RDS"))
d_data <- readRDS(file.path(sim_folder, "d.RDS"))

X_df <- X_data$X
d <- d_data$values

# No more all_covs or all_data parameters!
# X already has the correct predictors for this comparison
```

---

### 2.6 `sim_functions/summary_dt.R`

**OLD signature (approximate):**
```r
summary_dt <- function(comp_no, n_folds, chosen_var, loss_function, all_data, results_dir, ...)
```

**NEW signature:**
```r
summary_dt <- function(comp_no, n_folds, loss_function, results_dir, ...)
```

**Inside function:**
```r
# ADD at start:
sim_folder <- file.path(results_dir, sprintf("comparison_%03d", comp_no))

# Load comparison-specific data
X_data <- readRDS(file.path(sim_folder, "X.RDS"))
d_data <- readRDS(file.path(sim_folder, "d.RDS"))

# No more all_data parameter needed
```

---

### 2.7 `sim_functions/run_esim.R`

**OLD signature:**
```r
run_esim <- function(comp_no, s, all_covs, n_cores, all_data, results_dir, n_iters)
```

**NEW signature:**
```r
run_esim <- function(comp_no, n_cores, results_dir, n_iters)
```

**Inside function:**
```r
# ADD at start:
sim_folder <- file.path(results_dir, sprintf("comparison_%03d", comp_no))

# Load comparison-specific data
X_data <- readRDS(file.path(sim_folder, "X.RDS"))
d_data <- readRDS(file.path(sim_folder, "d.RDS"))

# No more all_covs, all_data, or s parameters
```

---

### 2.8 `sim_functions/summary_esim.R`

**Likely already loads data from folders, but check:**
```r
# Should load z and d from comparison folder if needed
sim_folder <- file.path(results_dir, sprintf("comparison_%03d", comp_no))
z <- readRDS(file.path(sim_folder, "z.RDS"))
d <- readRDS(file.path(sim_folder, "d.RDS"))
```

---

## PART 3: What Gets Saved Where

### Each comparison_XXX/ folder now contains:

```
comparison_001/
├── X.RDS              # NEW! Predictor matrix (PUMA-level) from this sample
│                      # Structure: list(puma = c(...), X = data.frame(...))
│                      # X has 13 columns: mean_age, pct_male, ..., has_hicov
├── z.RDS              # Direct estimates (PUMA-level) for response variable
│                      # Structure: list(puma = c(...), values = c(...))
├── d.RDS              # Design-based variances (PUMA-level)
│                      # Structure: list(puma = c(...), values = c(...))
├── dt_1fold/          # Data thinning results (1-fold)
├── dt_5fold/          # Data thinning results (5-fold)
├── emp_sim/           # Empirical simulation
│   ├── 001/
│   │   └── w.RDS      # Synthetic data generated using d from this comparison
│   ├── 002/
│   └── ...
└── fit_on_z/          # Full data fits (models fit on X and z)
```

**Key changes:**
- ✅ Added X.RDS (predictors vary per comparison!)
- ✅ All saved with PUMA IDs to ensure alignment
- ✅ X, z, d all from SAME sample (no mismatch possible)

### sim_config.RDS now contains:

```r
sim_config <- list(
  # NO LONGER CONTAINS:
  # - truth, chosen_var (never did)
  # - all_data (deleted!)
  # - all_covs (replaced)

  # NOW CONTAINS:
  acs_pop = acs_pop,                       # Individual-level population data
  samp_frac = 0.0225,                      # Sampling fraction (2.25% - tested optimal)
  response_var = "medi_cal_qualified",     # Which response to estimate
  X_approach = "population",               # How to compute X (population or estimated)
  predictor_names = predictor_names,       # Which variables to aggregate (13 vars)
  n_comp = n_comp,                         # Number of comparisons
  n_cores = n_cores,                       # Parallel processing
  results_dir = "_results"                 # Output directory
)
```

---

## PART 4: Key Principle

**STOP passing variances as parameters. START reading them from comparison folders.**

Every function that needs `d` or `z` should:
1. Build path: `sim_folder <- file.path(results_dir, sprintf("comparison_%03d", comp_no))`
2. Load data: `d <- readRDS(file.path(sim_folder, "d.RDS"))`
3. Load data: `z <- readRDS(file.path(sim_folder, "z.RDS"))`

---

## PART 5: Benefits of This Change

1. ✅ **More realistic** - Each comparison has different sampling variability
2. ✅ **Design-based** - Properly accounts for survey design
3. ✅ **Self-contained** - Each comparison folder has all its data
4. ✅ **Flexible** - Easy to add comparison-specific features later

---

## PART 6: Implementation Progress & Status

### ✅ COMPLETED PHASES (Architecture & Data Pipeline)

**Phase 0:** Data files ✅
- ✅ ca_pums_population.rds in place (1,761,769 individuals, 281 PUMAs)

**Phase 1:** `sampling_and_setup.R` implemented & tested ✅
- ✅ `create_X()` function with "population" and "estimated" approaches
- ✅ `setup_comp()` saves X, z, d as lists with PUMA IDs
- ✅ `esim_helper()` loads new list format
- ✅ Explicit PUMA sorting added (line 50-51) for guaranteed alignment
- ✅ Configurable via `sim_config$X_approach`
- ✅ Tested with 50 comparisons - all clean, median SE=3.4%
- ✅ Final params: samp_frac=0.0225, min_sample_size=30, response=medi_cal_qualified

**Phase 2:** Function signatures & data loading updated ✅
- ✅ run_comparisons.R: Updated data setup, removed all_data/all_covs, commented out model fitting
- ✅ run_dt.R: Loads X/z/d from comparison folder, removed old parameters
- ✅ summary_dt.R: Loads d from folder, updated all loss functions
- ✅ run_esim.R: Loads X/d from folder, removed old parameters
- ✅ summary_esim.R: Loads z from folder (new list format)
- ✅ full_data_fit.R: Renamed from zfit.R, loads X/z/d from folder
- ✅ summary_oracle.R: Loads z/d from folder (new list format)

**Phase 3:** Pipeline testing completed ✅
- ✅ test_setup_only.R: Verified setup_comp/setup_esim with n_comp=2
- ✅ test_pipeline_minimal.R: Verified DT + ESIM end-to-end
- ✅ **test_full_pipeline.R: Complete end-to-end test (Nov 7 2024)**
  - Configuration: 3 comparisons, 10 DT 1-fold configs, DT 5-fold, ESIM
  - Total runtime: ~37 seconds
  - Breakdown: Setup <1s, DT 1-fold 23.6s, DT 5-fold 2.4s, ESIM 9.2s
  - Results: 321 rows (297 DT + 24 ESIM across all configurations)
  - All methods completed successfully, no errors
  - Verified: Data generation, model fitting, summaries, result saving
- ✅ All data loading working correctly
- ✅ PUMA alignment maintained throughout
- ✅ Using placeholder Fay-Herriot model with minimal MCMC (`ndesired=10, nburn=10`) in:
  - `full_data_fit.R:42`
  - `run_dt.R:102-103`
  - `run_esim.R:64, 114`

**Architecture Migration Status:**
- All function signatures updated ✅
- All data loading working correctly ✅
- PUMA alignment guaranteed by explicit sorting ✅
- Tested end-to-end with DT and ESIM ✅

---

### 📋 TODO - MODEL IMPLEMENTATION (Phase 4)

**Current Blocker:** Spatial basis function model needs to be implemented before production runs.

**Phase 4: Spatial Model Implementation** 🔴 NOT YET STARTED
- ❌ Implement `models/spatial_basis_fit.R` (spatial basis function model)
- ❌ Update all three fitting functions to use spatial model:
  - `full_data_fit.R` (currently uses placeholder FH model)
  - `run_dt.R` (currently uses placeholder FH model)
  - `run_esim.R` (currently uses placeholder FH model)
- ❌ Increase MCMC iterations to production values in all three files:
  - Currently: `ndesired=10, nburn=10` (fast testing)
  - Production: `ndesired=1000, nburn=10000` (or appropriate for spatial model)
- ❌ Update `summary_oracle.R` if needed for spatial model output format
- ❌ Uncomment model fitting sections in `run_comparisons.R` (lines 119-193)
- ❌ Test full data fit with spatial model (n_comp=2)
- ❌ Run production comparisons (n_comp=50-70)

**What's Ready:**
- Data pipeline: Complete ✅
- All downstream functions (DT, ESIM): Complete ✅
- Sampling configuration: Finalized ✅
- Function signatures: All updated ✅

**What's Blocking Production:**
- Spatial model implementation: Required 🔴

**Key Test Files:**
- `test_setup_only.R`: Tests setup_comp/setup_esim with n_comp=2
- `test_pipeline_minimal.R`: Tests full pipeline (setup + DT + ESIM) with n_comp=2
- `test_two_configs.R`: Tests different sampling configurations with 50 comparisons

**Final Configuration:**
```r
sim_config <- list(
  acs_pop = acs_pop,
  samp_frac = 0.0225,                  # 2.25% sampling rate
  response_var = "medi_cal_qualified", # POVPIP < 138
  X_approach = "population",           # Fixed X across comparisons
  n_comp = 50,                         # Number of comparisons
  n_cores = 10,
  results_dir = "_results"
)
```

**Sampling Results (50 comparisons tested):**
- 0 problematic comparisons (no d=0, z=0, or z=1 cases)
- Median SE: 0.034 (3.4%)
- SE range: [0.009, 0.117]
- All 281 PUMAs sampled in every comparison
