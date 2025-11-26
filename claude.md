# Design-based Comparison Study For Data Thinning 

This is a repo for a research project I am conducting on data thinning for model selection. For the full background, see `_for_claude/background.tex`.

## Coding Context
I am working on this repo by myself. I would prioritize brevity / clarity and less on production-level software engineering concerns like error-catching.  Informative & concise comments would be appreciated (but avoid verbose comments).

Please be extra careful for making changes that involve nuanced statistical reasoning, evaluating results, etc. If you're unsure, just ask me (I don't mind). Ex: if you're unsure about what goes where in a model call please take a conservative approach. 

For tasks that are mostly software-enginnering based, I fully trust your judgement. 

### What needs to be done

Most of the code will be adapted from `dev/dt_choose_covs`. This repo differs in **two fundamental ways**:

1. **Data Generation:** Uses design-based sampling from California PUMS data
   - Each comparison = independent sample from `data/ca_pums_population.rds`
   - Predictors (X): Computed from full population (fixed across comparisons, configurable)
   - Response (z, d): Survey-weighted direct estimates from sample (vary across comparisons)
   - No pre-aggregated data - sampling happens on-the-fly

2. **Models:** Will use spatial basis function model instead of Fay-Herriot
   - New model code needs to be implemented: `spatial_basis_fit.R`
   - For now: Comment out all model fitting calls until model is ready

---

## Implementation Tasks

### ✅ Completed - Architecture & Data Pipeline

**Phase 1: Data Setup & Sampling** ✅
- ✅ Created `sim_functions/sampling_and_setup.R` with design-based sampling
- ✅ Implemented `create_X()` function with "population" and "estimated" approaches
- ✅ Made X configurable via `sim_config$X_approach`
- ✅ Save X.RDS, z.RDS, d.RDS with PUMA alignment (as lists with explicit sorting)
- ✅ Updated `esim_helper()` to load new list format
- ✅ Finalized sampling parameters: `samp_frac=0.0225`, `min_sample_size=30`
- ✅ Response variable: `medi_cal_qualified = ifelse(POVPIP < 138, 1, 0)`
- ✅ Tested with 50 comparisons - 0 problematic cases, median SE=3.4%

**Phase 2: Function Signatures & Data Loading** ✅
- ✅ Updated `run_comparisons.R` (data setup, function calls, commented out model fitting)
- ✅ Updated `run_dt.R` - loads X/z/d from comparison folder
- ✅ Updated `summary_dt.R` - loads d from comparison folder
- ✅ Updated `run_esim.R` - loads X/d from comparison folder
- ✅ Updated `summary_esim.R` - loads z from comparison folder
- ✅ Updated `full_data_fit.R` - loads X/z/d from comparison folder
- ✅ Updated `summary_oracle.R` - loads z/d from comparison folder

**Phase 3: Pipeline Testing** ✅
- ✅ Tested setup with `test_setup_only.R` (2 comparisons)
- ✅ Tested full pipeline with `test_pipeline_minimal.R` (DT + ESIM with n_comp=2)
- ✅ **Full pipeline test with `test_full_pipeline.R` (3 comparisons, Nov 7 2025):**
  - Configuration: 3 comparisons, 5 epsilon values (0.1, 0.3, 0.5, 0.7, 0.8), 2 repeat counts (1, 3)
  - Runtime: ~37 seconds total
    - Setup: <1 second
    - DT 1-fold: 23.6s (10 configs × 3 comparisons × 3 models)
    - DT 5-fold: 2.4s (3 comparisons × 5 folds × 3 models)
    - ESIM: 9.2s (3 comparisons × 100 iterations × 3 models)
    - Summaries: <2 seconds
  - Results: 321 rows (297 DT + 24 ESIM)
  - All methods completed successfully with no errors
  - Files created: X.RDS, z.RDS, d.RDS per comparison, all method results saved
- ✅ All data loading verified working correctly
- ✅ PUMA alignment maintained throughout pipeline
- ✅ Using placeholder Fay-Herriot model with minimal MCMC (`ndesired=10, nburn=10`) in all three fitting functions:
  - `full_data_fit.R:42`
  - `run_dt.R:102-103`
  - `run_esim.R:64, 114`

**Sampling Configuration (Updated Nov 14, 2025):**
- Sampling: Simple PPS using PWGTP (no informative sampling)
- **samp_frac: 0.01 (1%)** - Updated for production based on pilot tests
- min_sample_size: 30 per PUMA
- X_approach: "population" (fixed across comparisons)
- All 281 PUMAs sampled in each comparison
- **Response: WAGP (annual wages)** - Continuous variable with strong spatial signal

**Response Variable Selection (Updated Nov 2025):**
- **BREAKTHROUGH:** Continuous responses dominate binary responses for MSE differentiation
- **Winner:** **WAGP (annual wages)** - 27.7% OD-Oracle differentiation with weak priors 🏆
  - 2.5× better than best binary response (rents_home: 11%)
  - Clear U-shaped curve, interior optimum at nbasis=80
  - Source: `comprehensive_response_search.R` tested 11 responses
- **Runner-up:** **POVPIP (poverty ratio)** - 16.6% differentiation 🥈
  - 1.5× better than rents_home
  - Strong U-shaped curve, optimum at nbasis=60
- **Best binary:** rents_home - 11.1% differentiation with 2 covariates ⭐
- **Key insights:**
  1. **Weak priors essential** (a=b=c=d ≤ 0.001) to induce overfitting → U-shaped curves
  2. **Continuous > Binary:** Economic variables (wages, poverty) have natural spatial gradients
  3. **Fewer covariates** → More spatial signal (1-2 covariates optimal)
  4. **Sampling fraction:** 1% works well (0.5% also tested, more noise)
- **Test datasets created:** 5 responses in `exploration/testdata_*.RDA` (all 1% sampling)
  - povpip, rents_home, employed, owns_home, medi_cal_qualified
- **Notebook:** `exploration/explore_responses.Rmd` for interactive testing
- **Full documentation:** See `_for_claude/response_hacking.md` for complete exploration history

**Design-based Simulation Initial Analysis (Nov 2025):**
- **Completed:** Initial testing with `exploration/minimal_dt_test.Rmd`
- **Goal:** Test multiple response-covariate combinations across 3 direct estimates
- **Key findings:**
  1. **Priors:** Very diffuse priors (a=b=c=d ≤ 0.001) essential for non-flat curves
  2. **Covariates:** Dramatically affect curve shape, even for same response
  3. **Variability:** OD-Oracle curves vary across direct estimates (DT curves less so)
  4. **Epsilon:** DT nbasis selection not very sensitive to epsilon (tested 0.2 vs 0.4)
  5. **Response quality:** rents_home shows better spatial structure than owns_home
- **Advisor feedback:**
  - Favor rents_home with **fewer covariates** (2 or even 0)
  - More spatial variability → easier to see basis function effects
  - Consider nbasis 10-60 with finer grid (avoid high nbasis instability)
  - **Scale covariates** (good practice with diffuse priors)
- **Covariate analysis for rents_home** (`exploration/find_predictive_covariates_rents.R`):
  - **Strongest predictors:** mean_rooms (R²=0.81), pct_crowded (R²=0.79) ← TOO predictive
  - **Moderate predictors:** pct_married (R²=0.56), pct_citizen (R²=0.38)
  - **Weak predictors:** pct_hispanic (R²=0.19), mean_age (R²=0.15), pct_employed (R²=0.15)
  - **Trade-off:** Strong predictors → flat curves; weak/no predictors → strong spatial signal
- **Current recommendation for rents_home:**
  - **0 covariates** (intercept-only): Maximum spatial signal, clearest nbasis differentiation
  - **OR 1-2 weak covariates** (mean_age, pct_employed): Balance realism with signal
  - **Avoid:** Strong housing predictors (mean_rooms, pct_crowded) → leave little spatial variance

---

### ✅ Completed - Spatial Basis Function Model

**Phase 4: Spatial Model Implementation** ✅
- ✅ Implemented `models/spatial_basis_fh.R` - Gibbs sampler for spatial basis function Fay-Herriot model
- ✅ Implemented `models/spatial_basis_fh.stan` - Stan version for validation
- ✅ Fixed critical bug: Variable shadowing in Gibbs sampler (hyperparameters a,b were overwritten by temp vectors)
- ✅ Validated Gibbs sampler against Stan - estimates agree reasonably well
- ✅ Performance tested: Gibbs is ~3-3.5× faster than Stan (with nburn=18000 vs Stan warmup=1000)
- ✅ Convergence testing: 5-chain diagnostics confirm θ converges well (Rhat < 1.01) across nbasis=5-35

**Model Details:**
- **Specification:** θᵢ = Xᵢ'β + vᵢ + wᵢ
  - v = S·η, η ~ N(0, τ²Iᵣ) [spatial random effects via basis functions]
  - w ~ N(0, σ²Iₘ) [IID random effects]
  - S: m × r spatial basis matrix from Moran's I eigenvectors
- **Prior (default):** a=b=c=d=1
  - τ² ~ IG(1, 1) for spatial variance
  - σ² ~ IG(1, 1) for IID variance
  - Relatively non-informative prior for variance components
- **MCMC Settings:**
  - **NC test data (100 areas):**
    - Burn-in: 18,000 iterations (production)
    - Samples: 2,000 iterations
    - Runtime: ~3-5 seconds per fit (nbasis=5-35)
    - **Convergence:** Multi-chain diagnostics (5 chains) show θ Rhat < 1.01 across all nbasis values
    - **Note:** Variance components (τ², σ²) occasionally show Rhat > 1.1 at high nbasis, but small area estimates (θ) consistently converge well
  - **CA PUMA data (281 areas):**
    - **Recommended burn-in: 1,500-2,000 iterations** (validated with multi-chain diagnostics)
    - Samples: 2,000 iterations
    - **Convergence:** Excellent convergence even at nburn=1500 (all Rhat < 1.01, ESS 1800-10000)
    - Tested nbasis = 10, 20, 40, 60, 80, 100 (all converge well)
    - CA data converges much faster than NC - can use shorter burn-in for production
- **Basis functions:** Use `nbasis` parameter to control number of eigenvectors
  - **Note:** Maximum nbasis is limited by number of positive eigenvalues (>1e-10) from Moran's I decomposition
  - NC test data: 37 positive eigenvalues (100 areas) - tested nbasis = 5, 10, 20, 30, 35
  - CA PUMA data: ~110 positive eigenvalues (281 areas, 793 edges) - tested nbasis = 10-100
  - Requesting nbasis > available eigenvalues will use all available and trigger warning

**Files:**
- `models/spatial_basis_fh.R` - Main Gibbs sampler implementation
- `models/spatial_basis_fh.stan` - Stan implementation for validation
- `models/mcmc_helper.R` - Helper functions for MCMC
- **NC test data:**
  - `models/test_gibbs_convergence.R` - Multi-chain convergence diagnostics (5 chains)
  - `models/test_gibbs_vs_stan.R` - Timing and validation comparison with Stan
  - `models/diagnostics.Rmd` - Comprehensive diagnostic report with trace plots, ACF, spatial maps
- **CA PUMA test data:**
  - `models/generate_ca_test_sample.R` - Generate test sample from design-based comparison setup
  - `models/test_ca_gibbs_convergence.R` - Multi-chain convergence diagnostics (5 chains)
  - `models/test_ca_gibbs_vs_stan.R` - Timing and validation comparison with Stan
  - `models/diagnostics_ca.Rmd` - Comprehensive diagnostic report for CA data

---

### ✅ Completed - Adjacency Matrix for California PUMAs

**Adjacency matrix A has been created and saved**

Created `data/create_puma_adjacency.R` that automatically downloads PUMA shapefiles and generates the spatial adjacency matrix.

**Script:** `data/create_puma_adjacency.R`
- Tries multiple years (2022, 2021, 2020, 2019, 2018) to find matching vintage
- Found: **2022 (2020-based PUMAs)** perfectly aligns with population data
- Validates PUMA ID alignment before proceeding
- Creates binary adjacency matrix using rook contiguity

**Output:** `data/ca_puma_adjacency.RDA`
- 281×281 binary adjacency matrix
- 793 edges, average 5.64 neighbors per PUMA
- No island PUMAs (all have 2-15 neighbors)
- Perfectly aligned with `ca_pums_population.rds`

**Usage in models:**
```r
load("data/ca_puma_adjacency.RDA")  # Loads A (matrix) and puma_shape (sf object)
```

---

### ✅ Completed - nbasis Selection Analysis

**Single Dataset Analysis for Model Selection** ✅
- ✅ Created `one_data_analysis.Rmd` - Analysis to test data thinning for selecting optimal nbasis
- ✅ Compares three MSE metrics across nbasis values (10-100 by increments of 10):
  1. **OD-Oracle MSE:** Full-data fit evaluated against population truth (benchmark)
  2. **DT-Oracle MSE:** Thinned-data fit (eps=0.3) evaluated against truth (actual performance)
  3. **DT MSE:** Unbiased test set estimator (what you'd use in practice, no truth needed)
- ✅ Uses parallelization (10 cores) for data thinning fits
- ✅ Implements 3 repetitions per nbasis to reduce variance
- ✅ Includes both raw MSE plots and normalized [0,1] plots for comparison

**Key Components:**
- `create_dt_split()`: Single-fold data thinning (z_train, z_test)
- `calc_dt_mse()`: Unbiased MSE estimator from background.tex with bias correction
- De-scaling: DT models estimate eps*theta, so divide by eps before evaluating against truth
- Data alignment: Both test data and truth sorted by PUMA ID

**Purpose:** Validate that data thinning can identify the oracle-optimal nbasis using only train/test splits

---

### ✅ RESOLVED - Response Variable Selection

**Status:** COMPLETED - Found excellent configurations for data thinning testing

**Solution:** Comprehensive response search (Nov 2025) identified continuous responses as superior

**Final Configurations:**

| Response       | Type       | Covariates              | samp_frac  | MSE Diff  | Best nbasis | Quality       |
| -------------- | ---------- | ----------------------- | ---------- | --------- | ----------- | ------------- |
| **WAGP**       | Continuous | mean_age                | 0.005      | **27.7%** | 80          | ⭐⭐⭐ Excellent |
| **POVPIP**     | Continuous | mean_age                | 0.005-0.01 | **16.6%** | 60          | ⭐⭐ Excellent  |
| **rents_home** | Binary     | pct_asian, pct_bachelor | 0.01       | 11.1%     | 50          | ⭐ Good        |
| employed       | Binary     | mean_age                | 0.005      | 3.1%      | 70          | Moderate      |
| owns_home      | Binary     | 4 covariates            | 0.0075     | 2.1%      | 10          | Weak          |

**Key findings:**
- Continuous responses (WAGP, POVPIP) provide 2-3× better differentiation than binary
- Weak priors (a=b=c=d=0.001) essential for U-shaped curves
- All show clear interior optima (not at boundaries)
- **Recommendation:** Use WAGP or POVPIP for production runs

**Test datasets created:**
- `exploration/testdata_*.RDA` - 5 responses with 1% sampling, ready for testing
- `exploration/explore_responses.Rmd` - Interactive notebook for quick testing
- `exploration/TEST_DATASETS_README.md` - Documentation with recommendations

**Pilot Test Framework (Nov 2025):**
- ✅ Tested multiple response variables and model configurations (archived in `exploration/archive_pilots/`)
- ✅ Key finding: Fixed spatial effects work better than random effects for model selection
- **Status:** COMPLETED - Led to decision to use rents_home with fixed spatial basis

**Initial Analysis (Nov 16, 2025):**
- ✅ Ran 10 comparisons for both WAGP and rents_home with fixed spatial effects
- ✅ Analyzed oracle MSE consistency across methods (DT, WAIC, DIC, ESIM)
- **Key finding:** Random spatial effects produce inconsistent optimal nbasis across comparisons
- **Decision:** Proceed with rents_home + fixed spatial effects (NIKE swoosh pattern, interior minima)
- **Problem identified:** Optimal nbasis varies widely across comparisons (SD~26) - need more consistent signal
- **Results:** See `initial_conclusion/` for detailed analysis and plots

---

### 📁 Directory Cleanup (Nov 2025)

**Archived directories:**
- `archive/` - Old run_comparisons scripts for abandoned configurations
- `archive_results/` - Results from abandoned configurations (_results_wagp_random, etc.)
- `exploration/archive_pilots/` - Old pilot test directories
- `exploration/archive_scripts/` - Old exploration scripts and test datasets
- `exploration/archive_notebooks/` - Interactive Rmd files for exploratory testing
- `oracle_consistency_analysis/archive_response_search/` - Intermediate response screening results

**Removed files (Nov 2025 cleanup):**
- `models/nc_*.RDS/RDA` - NC test data (can regenerate if needed)
- `models/diagnostics.html` - Diagnostic report (can regenerate)
- `analyze_per_dataset_agreement.R` - Intermediate analysis script
- `analyze_per_dataset_proximity.R` - Intermediate analysis script

**Active directories:**
- `_results_pubcov/` - 10 comparisons with PUBCOV + equal_50 (for advisor report)
- `_results_rents/` - 10 comparisons with rents_home, 6 sampling configs (for advisor report)
- `initial_conclusion/` - Analysis from 10-comparison test (Nov 16, 2025)
- `oracle_consistency_analysis/` - Sampling design consistency tests (Nov 2025)

---

### ✅ Oracle Consistency Analysis (Nov 2025)

**Problem:** rents_home with fixed spatial effects shows interior minima (NIKE swoosh ✓) but **inconsistent optimal nbasis** across comparisons (SD~26, range 20-90)

**Goal:** Find response variable and sampling design that produces **consistent optimal nbasis** across comparisons to enable reliable method comparison

**Approach:** Expanded candidate response search, tested PUBCOV (Public Health Coverage) with multiple sampling schemes

**Status:** ✅ COMPLETED - Found winning configuration!

**Winner: PUBCOV + equal_50** 🏆

**Test configurations (10 comparisons each):**
1. **baseline:** samp_frac=0.01, min=30 (proportional allocation)
2. **equal_30:** equal allocation with n=30 per PUMA
3. **equal_50:** equal allocation with n=50 per PUMA
**Results:**

| Scheme    | MSE Variation | SD(nbasis) | Mean nbasis | Boundary % | **Result**        |
|-----------|---------------|------------|-------------|------------|-------------------|
| equal_50  | **19.5%**     | **13.2** ⭐ | 38          | **0%** ✓   | ✓✓✓ **PASS ALL**  |
| baseline  | 16.3%         | 14.0       | 38          | **0%** ✓   | ✓✓✓ **PASS ALL**  |
| equal_30  | 27.6% ⭐       | 13.5       | 25          | 30% ✗      | **FAIL**          |

**Key findings:**
- **PUBCOV (Public Health Coverage)** is first response variable to pass all filters
- **equal_50** provides best balance: 19.5% MSE differentiation, SD=13.2, no boundary selections
- equal_30 has highest MSE signal (27.6%) but 30% boundary selections (too noisy at n=30)
- Optimal nbasis range for equal_50: 20-50 (consistent across comparisons)

**Recommendation for production runs:**
- **Response variable:** PUBCOV (public health insurance coverage, binary)
- **Sampling scheme:** equal_50 (equal allocation, n=50 per PUMA)
- **Model config:** Intercept-only (no covariates), fixed spatial effects, weak priors (a=b=c=d=0.001)
- **nbasis range:** 10-100 by increments of 10

**Files:**
- `oracle_consistency_analysis/test_pubcov_detailed.R` - Main test script
- `oracle_consistency_analysis/test_pubcov_equal30.R` - Additional equal_30 test
- `oracle_consistency_analysis/analyze_pubcov_results.R` - Metrics computation
- `oracle_consistency_analysis/plot_pubcov_curves.R` - Visualization
- `oracle_consistency_analysis/response_search/pubcov_detailed/` - All results and plots

**Phase 2: Alternative Response Search** (Ongoing)
- Continue exploring alternative responses from other states (TX, NY, IL, LA, NJ)
- Test additional CA response variables with equal_50 sampling
- Goal: Find multiple response variables that pass filters for robustness testing

**Next Steps:**
- ~~Run full pipeline test: 10 comparisons with DT + ESIM methods to validate PUBCOV end-to-end~~
- ~~If successful, proceed to production: 50-70 comparisons for final method comparison~~

---

### ⚠️ CRITICAL FINDING: Flat Oracle MSE Curves (Nov 2025)

**Status:** ❌ **PROJECT PAUSED** - Response variables show insufficient steepness for meaningful method comparison

**What happened:**
1. Ran 10 comparisons for PUBCOV + equal_50 (passed v1 screening criteria)
2. Ran 10 comparisons for rents_home with 6 sampling configurations
3. **Discovered:** Oracle MSE curves are **flat** across all tested configurations

**Flatness Metrics (PUBCOV):**
- MSE penalty at ±10 nbasis: **3.3%** (target: >5%)
- MSE penalty at ±20 nbasis: **7.1%** (should be >10%)
- Models within 5% of optimal: **3.6 / 10** (too many)
- **Implication:** Methods can disagree by ±20-30 nbasis units yet have nearly identical MSE (~5% penalty)

**Method Performance (PUBCOV, 10 comparisons):**
- 1-fold DT (ε=0.5): 70% within ±20 of oracle, 5.1% MSE penalty
- 5-fold DT: 50% within ±20, 5.2% penalty
- DIC: 80% within ±20, 6.5% penalty
- ESIM: 50% within ±20, 5.3-5.6% penalty

**Why this is problematic:**
- Cannot meaningfully differentiate method quality when cost of being wrong is negligible
- Methods work "well enough" but impossible to identify which is better
- Defeats the purpose of method comparison study

**Results Available for Advisor Report:**
- `_results_pubcov/` - 10 comparisons, all methods
- `_results_rents/` - 10 comparisons, 6 sampling configurations
- `analyze_oracle_flatness.R` - Post-mortem analysis
- Plots show flat curves across all configurations

**Updated Screening Criteria:**
- See `_for_claude/screening_criteria_v2.md` for revised criteria
- **Added:** Steepness requirement (>5% penalty at ±10 units)
- **Removed:** Spatial correlation pre-screening (uninformative)

**Possible Explanations:**
1. Spatial basis function model is robust to nbasis mis-specification in 20-100 range
2. These response variables lack sufficient spatial signal
3. Need fundamentally different modeling approach

**Decision:**
- Report findings to advisors (flat curves = valuable negative result)
- Potentially pivot to studying robustness of spatial basis models (different research question)
- Or explore alternative model classes / response variables from other data sources

**Files:**
- `analyze_oracle_flatness.R` - Post-mortem flatness analysis
- `analyze_per_dataset_proximity.R` - REMOVED (housekeeping)
- `_results_pubcov/oracle_overlay_raw.png` - Visual proof of flatness
- `_results_rents/oracle_faceted_raw.png` - All 6 sampling configs (all flat)

---

### ✅ BREAKTHROUGH: CA Employed Response (Nov 2025)

**Status:** ✅ **PRODUCTION-READY CONFIGURATION IDENTIFIED**

After discovering flat curves with PUBCOV, conducted systematic testing of alternative response variables and found **CA employed** (employment-to-population ratio) produces excellent oracle properties.

**Critical Bug Found & Fixed (Nov 25, 2025):**
- **Issue:** Oracle screening scripts used wrong config keys (`population_data`, `adjacency_data`)
- Functions expected different keys (`population_file`, `adjacency_file`)
- Silent fallback to CA data masked the error → Oracle screening labeled "ny_phase2_employed" actually used **CA data**
- **Fix:** Removed all silent defaults, now requires explicit `population_file` and `adjacency_file` in all configs
- **Files changed:** `sim_functions/sampling_and_setup.R`, `sim_functions/full_data_fit.R`, `sim_functions/run_dt.R`, `sim_functions/run_esim.R`

**Systematic Audit Conducted:**
1. **Oracle screening pipeline:** ✓ Structurally valid (setup_comp → full_data_fit → summary_oracle)
2. **Plotting code:** ✓ Valid (correctly uses summary_oracle to read chains.RDS and theta_true.RDS)
3. **Bug fixes:** ✓ Removed CA fallbacks, standardized config keys
4. **Documentation:** Created `SYSTEMATIC_AUDIT_RESULTS.md`

**Sample Size Comparison (CA Employed, 20 comparisons each):**

| equal_n | Mean Opt | SD | Range | MSE Var | Penalty ±6 | Assessment |
|---------|----------|-----|--------|---------|------------|------------|
| 40 | 13.1 | 3.7 | [3,15] | 47.1% | 5.4% | ⚠️ Range touches boundary |
| 50 | 14.2 | 2.7 | [6,18] | 44.1% | 6.8% | ⚠️ Too close to boundary (6-3=3) |
| **75** | **16.4** | **2.7** | **[12,24]** | **31.5%** | **8.0%** | ✅ **SELECTED** |
| 100 | 18.9 | 8.3 | [15,51] | 24.7% | 8.7% | ❌ SD too high |

**Why equal_75 Selected:**
1. **Low SD (2.7)** - Consistent optimal across comparisons
2. **Safe range [12-24]** - Well away from lower boundary (nbasis=3), no boundary contamination
3. **Deep trough** - 8% average penalty at ±6, asymmetric (10.6% left, 5.3% right) reflects realistic cost structure
4. **Strong signal** - 31.5% total MSE variation, sufficient to differentiate methods
5. **Interior optimum** - Clear U-shaped curve with minimum at nbasis~15-16

**Final Production Configuration:**
```r
State: California (281 PUMAs)
Response: employed (employment-to-population ratio)
  - response_var = "employed"
  - response_filter = NULL  # Includes children ages 0-15 in denominator
  - Interpretation: Proportion of total population that is employed
Sampling: Equal allocation, n=75 per PUMA
Model: Fixed spatial basis (spatial_type="fixed")
Priors: c=d=0.001 (weak priors for U-shaped curves)
MCMC: nburn=1500, ndesired=2000
nbasis grid: [3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60]
Expected optimal: nbasis~15-16
```

**Test Script:**
- `test_ca_oracle_minimal.R` - Minimal oracle-only test (no excessive printing, easy to review)
- Validates oracle properties before running full pipeline

**Oracle Results Directories:**
- `_results_ca_employed_equal40_oracle/` - 20 comparisons (rejected: boundary issues)
- `_results_ca_employed_equal50_oracle/` - 20 comparisons (rejected: too close to boundary)
- `_results_ca_employed_equal75_oracle/` - 20 comparisons (**SELECTED**)
- `_results_ca_employed_equal100_oracle/` - 20 comparisons (rejected: high SD)

**Next Steps:**
1. Run full pilot (10-20 comparisons) with equal_75 + all methods (DT 1-fold, DT 5-fold, WAIC, DIC, ESIM)
2. Validate methods can identify optimal nbasis~15-16
3. If successful → scale to production (50-70 comparisons)

---

## Key Architecture Decisions

**No pre-aggregated data file:**
- ❌ DO NOT create `create_puma_data.R` or `data/puma_aggregated.rds`
- ✅ Each comparison samples individuals and aggregates on-the-fly

**Predictor matrix X - Fixed vs. Estimated:**
- **Current approach:** X computed from **true population** (Option 3)
  - X calculated once from full `acs_pop`, same for all comparisons
  - Aligns with SAE modeling assumption (X treated as known/fixed)
  - Cleaner for studying model selection methods
- **Alternative (for later):** X as direct estimates (Option 2)
  - Could compute X from survey sample using `svyby()` like we do for z
  - More conservative, both z and X vary across comparisons
  - Can make configurable via `sim_config$X_approach`
- **Response z:** Always survey-weighted direct estimates from sample (varies across comparisons)

**Function signatures changed:**
- OLD: Functions took `all_data`, `all_covs`, `chosen_var` as parameters
- NEW: Functions only take `comp_no` and `results_dir`, load X/z/d from comparison folders

**What gets saved per comparison:**
```
comparison_XXX/
├── X.RDS    # Predictors (13 vars) - currently from population, configurable
├── z.RDS    # Direct estimates from survey sample
├── d.RDS    # Design-based variances from survey sample
└── ...
```

See `_for_claude/architecture_change_analysis.md` for complete details.

---

### ✅ Full Pipeline Test Results (Nov 25, 2025)

**Status:** ✅ **COMPLETE** - Inline config approach validated, all methods tested

**Test Configuration:**
- Script: `test_ca_full_pipeline.R`
- Configuration: Inline config (saved as `model_config.RDS`)
- Response: CA employed (employment-to-population ratio)
- Sampling: Equal allocation, n=75 per PUMA
- Comparisons: 10
- Runtime: ~51 minutes (10 cores)

**Statistical Review Completed (Nov 25, 2025):**
- ✅ Reviewed all DT/ESIM/summary functions for misspecification
- ✅ **Bug fixed:** Predictive NLL variance formula in `summary_dt.R` (missing /eps² scaling)
- ✅ Verified: Plugin NLL, MSE bias correction, data thinning splits all correct
- ⚠️ **Note:** ESIM data fission formula (2z - w) increases noise - unclear if intentional, deferred

**Method Performance Results (10 comparisons):**

| Method | Mean nbasis | SD | Boundary % | Exact Match | Within ±3 | Ranking |
|--------|-------------|-----|------------|-------------|-----------|---------|
| **Oracle MSE** | 16.5 | 2.9 | 0% | - | - | Truth |
| **DIC** | 15.3 | 6.4 | 10% | **30%** | **60%** | 🥇 **Best** |
| **WAIC** | 15.3 | 6.2 | 10% | **20%** | **60%** | 🥈 2nd |
| DT eps=0.5, 3 reps | 12.9 | 5.3 | 10% | 20% | **70%** | 🥉 3rd |
| DT 5-fold | 14.4 | 8.1 | 10% | 10% | 30% | 5th |
| ESIM (50 iters) | 8.4 | 5.3 | **30%** | 30% | 40% | 4th |

**Key Findings:**

1. **WAIC/DIC excel:** 🏆
   - DIC: 30% exact agreement with oracle (best overall)
   - WAIC: 20% exact, 60% within ±3
   - Both correctly identify oracle mean (15.3 vs 16.5)
   - Only 10% boundary selections (very stable)

2. **DT shows systematic underestimation:**
   - All DT configs select fewer basis functions than oracle optimal
   - Best config (eps=0.5, 3 reps): 70% within ±3 but mean=12.9 vs oracle=16.5
   - Lower epsilon worse (eps=0.3: mean=9.0, 40% boundary)

3. **ESIM needs more iterations:**
   - 30% boundary selections (nbasis=3) with only 50 iters
   - When correct, exact match (30% rate)
   - Standard 100 iterations recommended for production

**Files:**
- `test_ca_full_pipeline.R` - Main test script with inline config
- `run_summary.R` - Summary script (reads `model_config.RDS`)
- `_results_ca_test_pipeline/` - Test results (10 comparisons)

**Workflow Established:**
```bash
./test_ca_full_pipeline.R              # Runs study, saves model_config.RDS
Rscript run_summary.R <results_dir>    # Analyzes results
```

**Next Steps:**
- Investigate DT systematic underestimation bias
- Scale to production (50-70 comparisons) if ready
- Consider testing more DT configurations or loss functions

---

## Statistical Context
This repo compares different model selection methods across **design-based synthetic datasets** that mimic real survey data.

**Data Generation Approach:**
- Population: Individual-level California PUMS data (`data/ca_pums_population.rds`)
- Sampling: Stratified PPS (probability proportional to size) by PUMA
- Per comparison: Draw independent sample → get z, d from sample; X from population
- Predictors (X): Computed once from full population (fixed across comparisons)
- Response (z, d): Survey-weighted direct estimates (vary across comparisons)
- Will run 50-70 comparisons (depends on computational time)

**What is new:** Using data thinning for model selection in small area estimation context. 

The methods we will compare are:

- Single-Fold Data Thinning (Primary method) with 1, 3, 5 repeats and evenly spaced out thinning/epsilon values (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), using two different loss functions (MSE vs. Negative Log-Likelihood)
- Multi-Fold (k=5) Data Thinning
- DIC / WAIC
- Empirical Simulation Study with two variants standard vs. Data Fission inspired approach 


### Evaluation

In this design-based simulation, we know the **true PUMA-level population means** from `ca_pums_population.rds`. This allows us to calculate several "Oracle" quantities:

- **Finite Population Oracle (FP-Oracle):** Score across datasets (expectation over sampling distribution). Estimated empirically across comparisons.
- **Observed-Data Oracle (OD-Oracle):** Computed knowing true population means. **Primary benchmark for this study.**
- **Thinned-Data Oracle (TD-Oracle):** Scores for specific data thinnings knowing true means.

**Primary metric: OD-Oracle MSE**
- Fit models on observed data (X, z) for a comparison
- Evaluate MSE against true population means
- Represents best achievable performance for that dataset
- Note: Some methods are likelihood-based; OD-Oracle for likelihood scores is less clear-cut 

Some ways to benchmark using OD-Oracle MSE:

1. Scatterplot of OD-Oracle MSE vs score from each method
2. Selected model: which model would be chosen based on each method
3. Ranking: did the ranking of the method agree with OD-Oracle based on Kendall's Tau or other similarity metrics?




