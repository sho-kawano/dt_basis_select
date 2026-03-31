### Study Design Summary

**Section 6.3** in my paper compares data thinning against existing model selection methods on **proportional-to-population sampling designs**. See `draft_2026-02-04.tex` for the full paper. Previous documents may have used the term PPS but that's not the standard nomenclature as that stands for **Probability Proportional to Size**. We need to start using the term: **proportional allocation** (PA). 

**Designs needed:**

- PA at **0.75%** sampling rate
- PA at **1.25%** sampling rate
- PA at **1.75%** sampling rate

(Evenly spaced at 0.5% increments for clean experimental design)

**Number of samples:** S = 20 per design (for now -- will increase later)

**Methods to compare (5 total):**

1. Data thinning — MSE validation (ε = 0.5, 0.6, 0.7, R = 5)
2. Data thinning — Negative log-likelihood validation (ε = 0.5, 0.6, 0.7, R = 5)
3. DIC
4. WAIC
5. ESIM (with 100 simulation iterations -- may not have been done for these designs)

Note that the two Data Thinning methods here share ε and repeats so it's just a matter of post-processing / calculating the metrics differently.

**Epsilon choice:** Computing ε = 0.5, 0.6, 0.7 allows us to assess robustness across Section 6's recommended range. Main text can focus on one value (likely 0.6 or 0.7) with others relegated to appendix if results are consistent. 

**Candidate models:** p ∈ {3, 6, 9, ..., 60} spatial basis functions (20 models)

**Metrics needed per method × design:**

- **MAD**: Mean absolute deviation from per-sample oracle basis
- **Directional bias**: Signed deviation (negative = under-selection)
- **Oracle penalty**: MSE cost from suboptimal selection
This was analyzed in `analysis/10_config_comprehensive.Rmd`. Another document that would be great to look review again would be `Design-Based Sim Summary.md`.
---

### What Results Likely Exist

From prior simulation runs, you should have results for:

- **1.25% PA** ✓
- **1.5% PA**✓
- **2% PA**✓

And possibly additional PA designs (0.5%, 1.0%, 1.5%, 1.75%) from earlier exploratory work.

---

### What's Missing

- **0.75% PA design** — needs to be run (full simulations + DT + oracle)
- **ESIM for all PA designs** — emp_sim folders currently empty (computationally expensive)

---

### Status Check: 1.25% and 2% PA Designs

**Question 1: Do we have all 5 methods?**

For **1.25% PA** (`_results_prop1p25pct_comparison/`) and **2% PA** (`_results_prop2pct_comparison/`):

| Method | Status | Location |
|--------|--------|----------|
| **DT-MSE** | ✓ Ready | `results_multi_config/dt_all_10configs.RDS` |
| **DT-Likelihood** | ⚠️ Needs aggregation | Raw results exist in `dt_1fold/eps_*/` folders |
| **DIC** | ✓ Ready | `results_multi_config/oracle_all_10configs.RDS` |
| **WAIC** | ✓ Ready | `results_multi_config/oracle_all_10configs.RDS` |
| **ESIM** | ✗ Not run | `emp_sim/` folders are empty |

**Summary:** 3 methods ready, 1 needs post-processing, 1 missing entirely.

**DT details:**
- Available epsilon values: 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 (all in raw results)
- Available R values: 1, 3, 5 (but dt_all_10configs.RDS only has R=1 for MSE)
- Raw DT results contain model fits; `summary_dt.R` computes loss functions (MSE or plugin_NLL)
- Need to run aggregation with R=5 and plugin_NLL loss function

**Action items for 1.25% and 2% PA:**
1. Aggregate DT-Likelihood results (ε = 0.5, 0.6, 0.7, R = 5) using `summary_dt()` with `loss_function = "plugin_NLL"`
2. Decide whether to run ESIM (expensive) or proceed with 4 methods

---

---

### Question 2: Are results stored with the three metrics?

**Short answer:** No. The aggregated files contain **model-level** results (loss values for each model × method × sample). The three summary metrics are **computed on-the-fly** in analysis code.

**What's stored:**
- `oracle_all_10configs.RDS`: comp_no, model, nbasis, mse_true, WAIC, DIC, config
- `dt_all_10configs.RDS`: comp_no, method, model, nbasis, metric, loss_function, epsilon, n_reps_used, config

**Metric computation status:**

| Metric | Definition | Computed Where | Status |
|--------|-----------|----------------|--------|
| **MAD** | Mean absolute deviation from oracle nbasis | `analysis/10_config_comprehensive.Rmd` (lines 177-231) | ✓ Code exists |
| **Oracle penalty** | % MSE increase from suboptimal selection | `analysis/10_config_comprehensive.Rmd` (lines 236-296) | ✓ Code exists |
| **Directional bias** | Signed deviation (negative = under-selection) | Not yet computed | ⚠️ Needs to be added |

**How they're computed:**
1. Find oracle-selected model per sample: `min(mse_true)` → `oracle_nbasis`
2. Find method-selected model per sample: `min(DIC/WAIC/metric)` → `method_nbasis`
3. Compute deviations per sample, then aggregate across samples per config

**Directional bias** would be identical to MAD computation but without `abs()`:
```r
mutate(directional_bias = method_nbasis - oracle_nbasis)  # instead of abs(...)
```

**Action item:** Add directional bias computation to analysis workflow (straightforward extension of existing MAD code).

---

### Question 3: Architecture improvements

**Question:** Is there a way we can improve the architecture without breaking the other dependencies (specifically the `sectionX.Rmd` in `analysis/`)?

**Current architecture:**
```
Raw results (per comparison)
  → Aggregation scripts (aggregate_*.R)
  → Aggregated files (results_multi_config/*.RDS)
  → Analysis Rmds (analysis/*.Rmd)
  → Metrics computed on-the-fly in each Rmd
```

**Issues identified:**
1. `dt_all_10configs.RDS` only has **MSE loss**, not plugin_NLL
2. Metrics computation code is **duplicated** across multiple Rmds
3. No aggregated results specifically for PA designs needed for Section 6.3/7
4. `dt_all_10configs.RDS` has R=1 and R=5, but section7 analysis may need to filter specifically

**Proposed approach: New standalone aggregation system**

Philosophy: Create a clean, purpose-built aggregation optimized for methods comparison. Keep old system (`aggregate_all_10configs.R` + existing files) for backwards compatibility. Migrate notebooks to new system later during repo polishing.

---

#### New Aggregation System Design

**File:** `aggregate_methods_comparison.R`

**What it aggregates:**
- **PA designs:** 0.75%, 1.25%, 2% (the 3 designs for Section 7)
- **DT results:** ε ∈ {0.5, 0.6, 0.7}, R = 5, **both MSE and plugin_NLL**
- **Oracle/DIC/WAIC:** From full data fits
- **20 samples per design** (expandable to 50+ later)

**Output structure:**
```
results_multi_config/
  methods_comparison/
    dt_results.RDS           # All DT results (both loss functions)
    oracle_results.RDS       # Oracle MSE, DIC, WAIC
    metrics_summary.RDS      # Pre-computed MAD, bias, penalty
```

**Key features:**
1. **Both loss functions by default** - No need to re-aggregate later
2. **Focused scope** - Only configs needed for Section 7, faster to run
3. **Pre-computed metrics** - MAD, directional bias, oracle penalty included
4. **Clean organization** - Separate folder for new system
5. **Self-contained** - Includes helper functions for metric computation
6. **Extensible** - Easy to add 0.75% PA results once available

---

#### Supporting Infrastructure

**File:** `sim_functions/compute_metrics.R`

Shared function for computing evaluation metrics (used by aggregation script and can be used by analysis Rmds):

```r
compute_metrics <- function(oracle_data, method_selections, method_name) {
  # oracle_data: comp_no, model, nbasis, mse_true
  # method_selections: comp_no, selected_nbasis
  # Returns: data frame with MAD, directional_bias, oracle_penalty
}
```

---

#### Migration Plan

**Phase 1 (Now):**
Build new system for Section 7, old system untouched
- Create `aggregate_methods_comparison.R` with new structure
- Create `compute_metrics.R` helper function
- Create Section 7 analysis Rmd using new aggregation
- Old system remains intact (Sections 2-6 continue working)

**Phase 2 (Post-Draft Completion):**
Migrate all notebooks during repo polishing
- Update all analysis notebooks to use new system
- Create equivalent aggregations for other sections
- Deprecate old aggregation files
- Clean unified aggregation system
- Well-documented for reproducibility

---

---

### Current Status

**Completed (Phase 1):**
- ✅ Created `sim_functions/compute_metrics.R` with helper functions
- ✅ Created `aggregate_methods_comparison.R` for new aggregation system
- ✅ Tested with existing 1.25% and 2% PA data
- ✅ New metrics include: mean_selected_nbasis, mean_oracle_nbasis, MAD, directional_bias, oracle_penalty

**Output available:**
- `results_multi_config/methods_comparison/dt_results.RDS` - DT with both MSE and plugin_NLL
- `results_multi_config/methods_comparison/oracle_results.RDS` - Oracle MSE, DIC, WAIC
- `results_multi_config/methods_comparison/metrics_summary.RDS` - Pre-computed evaluation metrics

**Remaining work:**
1. Run 0.75% PA design simulations (lowest sample size needed)
2. Run ESIM for PA designs (required for methods comparison)
3. Update aggregation to include 0.75% PA once available
4. Create Section 6.3 methods comparison notebook using new infrastructure

This approach gives you a clean foundation to build on without breaking existing work.

---

## Task 5 & 6 Scoping

### Task 5: Run 0.75% PA Design Simulations

**What:** Generate full simulation pipeline for prop_0.75pct design (lowest sample size PA design).

**Pipeline steps:**
1. **Setup comparisons** (`setup_comp`)
   - Generate 20 samples using proportional allocation at 0.75% sampling rate
   - Target: ~47 mean n per PUMA (range: ~25-82 based on population distribution)
   - Compute HT direct estimates (z) and design variances (d) per comparison

2. **Full data fits** (`full_data_fit`)
   - Fit all 20 models (p = 3, 6, 9, ..., 60) on observed z
   - Compute oracle MSE, DIC, WAIC per model
   - ~10-15 min per comparison × 20 = **~3-5 hours**

3. **Data thinning** (`run_dt`)
   - Run 1-fold DT for ε ∈ {0.5, 0.6, 0.7} with R=5 repetitions
   - 3 epsilon × 5 reps × 20 models per comparison
   - ~15-20 min per comparison × 20 = **~5-7 hours**

**Configuration:**
```r
samp_frac = 0.0075       # 0.75% sampling rate
equal_allocation = FALSE  # Proportional allocation
min_sample_size = 10      # Minimum n per PUMA
n_comp = 20               # 20 comparisons
results_dir = "_results_prop0.75pct_comparison"
```

**Total estimated time:** ~8-12 hours (parallelized with 11 cores)

**Output:** `_results_prop0.75pct_comparison/` with structure matching existing PA designs

---

### Task 6: Run ESIM for PA Designs

**What:** Empirical simulation - fit all candidate models on synthetic datasets to evaluate selection performance.

**Why ESIM is needed:**
- Benchmark method that doesn't use information criteria (DIC/WAIC) or data thinning
- Tests "if we could resample from the true distribution, which model would we select?"
- Gold standard for comparing selection methods (used in existing equal_75 analysis)

**How ESIM works:**
1. **Setup** (`setup_esim`): Generate synthetic direct estimates for each iteration
   - For each comparison: w_i ~ N(z_i, d_i) where z = observed estimate, d = design variance
   - Create 100 synthetic datasets per comparison
   - Fast: ~1 second per comparison

2. **Fitting** (`run_esim`): Fit all 20 models on each synthetic dataset
   - 100 iterations × 20 models × 20 comparisons = **40,000 model fits**
   - Each fit: MCMC with 1500 burn-in + 2000 samples
   - This is the computationally expensive step

**Designs to run:**
1. prop_0.75pct (20 comparisons × 100 iterations)
2. prop_1p25pct (20 comparisons × 100 iterations)
3. prop_1p75pct (20 comparisons × 100 iterations)

**Computational cost per design:**
- Setup: ~20 seconds (trivial)
- Fitting: 100 iterations × 20 models × ~1-2 min per model = **~33-67 hours per design**
  - Can parallelize across comparisons (20 parallel jobs, each doing 100×20 fits sequentially)
  - With 11 cores + efficient scheduling: **~6-12 hours wall time per design**

**Total for all 3 designs:** ~18-36 hours wall time (parallelized)

**Can we reduce iterations?**
- Original analysis used 100 iterations
- Could test with 50 iterations first (~half the time)
- User wants ESIM included (not optional), so need to decide on iteration count

**Configuration:**
```r
n_iter_esim = 100  # Number of synthetic datasets per comparison
n_cores = 11       # Parallel processing
# For each design:
setup_esim(ncomps = 20, results_dir = "...", niters = 100)
run_esim(comp_no, n_cores = 1, results_dir = "...",
         model_config, n_iters = 100)
```

**Output:** `emp_sim/XXX/results.RDS` for each iteration, containing model fits

---

### Combined Timeline Estimate

**Sequential approach:**
- Task 5 (0.75% PA sims): 8-12 hours
- Task 6 (ESIM all designs): 18-36 hours
- **Total: 26-48 hours**

**Can overlap:**
- Run Task 6 for prop_1p25pct and prop_2pct while Task 5 runs (data already exists)
- Then run Task 6 for prop_0.75pct after Task 5 completes
- **Overlapped total: ~18-36 hours** (dominated by ESIM)

