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
- ✅ **Full pipeline test with `test_full_pipeline.R` (3 comparisons, Nov 7 2024):**
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

**Sampling Configuration (finalized):**
- Sampling: Simple PPS using PWGTP (no informative sampling)
- samp_frac: 0.0225 (2.25%)
- min_sample_size: 30 per PUMA
- X_approach: "population" (fixed across comparisons)
- All 281 PUMAs sampled in each comparison

---

### ✅ Completed - Spatial Basis Function Model

**Phase 4: Spatial Model Implementation** ✅
- ✅ Implemented `models/spatial_basis_fh.R` - Gibbs sampler for spatial basis function Fay-Herriot model
- ✅ Implemented `models/spatial_basis_fh.stan` - Stan version for validation
- ✅ Fixed critical bug: Variable shadowing in Gibbs sampler (hyperparameters a,b were overwritten by temp vectors)
- ✅ Validated Gibbs sampler against Stan - estimates agree within ~10-30%
- ✅ Performance tested: Gibbs is ~2x faster than Stan even with 10x longer burn-in

**Model Details:**
- **Specification:** θᵢ = Xᵢ'β + vᵢ + wᵢ
  - v = S·η, η ~ N(0, τ²Iᵣ) [spatial random effects via basis functions]
  - w ~ N(0, σ²Iₘ) [IID random effects]
  - S: m × r spatial basis matrix from Moran's I eigenvectors
- **Prior (recommended):** a=c=2, b=d=0.25
  - τ² ~ IG(a, b) for spatial variance
  - σ² ~ IG(c, d) for IID variance
- **MCMC Settings (production):**
  - Burn-in: 10,000 iterations
  - Samples: 2,000 iterations
  - Runtime: ~2-3 seconds per fit (nbasis=5-30)
- **Basis functions:** Use `nbasis` parameter to control number of eigenvectors (tested: 5, 10, 20, 30)

**Files:**
- `models/spatial_basis_fh.R` - Main Gibbs sampler implementation
- `models/spatial_basis_fh.stan` - Stan implementation for validation
- `models/mcmc_helper.R` - Helper functions for MCMC

---

### 🔴 BLOCKER - Adjacency Matrix for California PUMAs

**Need to download PUMA shapefiles to create adjacency matrix A**

The spatial basis function model requires an adjacency matrix A for California PUMAs. This needs to be downloaded using the `tigris` package.

**Example code pattern** (adapted from `/Users/sho/dev/ssd_paper_code/sim_analysis_functions/load_county.R`):

```r
library(tigris)
library(spdep)

# Download PUMA shapefiles for California
puma_shape <- pumas(state = "CA", cb = FALSE, year = 2019)  # or appropriate year
rownames(puma_shape) <- puma_shape$GEOID10  # or GEOID depending on year
puma_shape <- puma_shape %>% arrange(GEOID10)

# Create adjacency matrix (rook contiguity)
A <- nb2mat(poly2nb(puma_shape), style = 'B', zero.policy = FALSE)

# Save for reuse
save(A, puma_shape, file = "data/ca_puma_adjacency.RDA")
```

**Required packages:**
- `tigris` - Download Census shapefiles
- `spdep` - Spatial dependence (poly2nb, nb2mat)

**Notes:**
- Use `cb = FALSE` for detailed boundaries (more accurate adjacency)
- `style = 'B'` gives binary adjacency (0/1)
- Make sure PUMA IDs align with those in `ca_pums_population.rds`
- This only needs to be done once and can be saved to `data/ca_puma_adjacency.RDA`

---

### 📋 TODO - Production Runs

**Phase 5: Integration & Production** 🔴 NOT YET STARTED
- ❌ **Download PUMA shapefiles and create adjacency matrix A** (blocker - see above)
- ❌ Create helper script (e.g., `data/create_puma_adjacency.R`) to download and save A
- ❌ Update all three fitting functions to use spatial basis function model:
  - `full_data_fit.R` (currently uses placeholder FH model)
  - `run_dt.R` (currently uses placeholder FH model)
  - `run_esim.R` (currently uses placeholder FH model)
- ❌ Update MCMC iterations to production values (ndesired=2000, nburn=10000)
- ❌ Uncomment model fitting sections in `run_comparisons.R`
- ❌ Test with spatial model on small number of comparisons
- ❌ Run full n_comp=50-70 comparisons

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




