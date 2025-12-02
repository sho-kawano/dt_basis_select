# Reproduction Guide - 10-Config Comprehensive Analysis

## Overview

This guide documents all files needed to reproduce the 10-config sampling design comparison analysis that compares Data Thinning (DT) against DIC/WAIC for spatial basis function model selection.

---

## Required Files

### 1. **Core Simulation Functions** (in `sim_functions/`)

```
sim_functions/
├── sampling_and_setup.R         # Sampling from CA PUMS + setup
├── full_data_fit.R              # Oracle fits (DIC/WAIC computation)
├── run_dt.R                     # Data thinning implementation
├── summary_dt.R                 # DT results summarization
├── summary_oracle.R             # Oracle/DIC/WAIC summarization
└── spatial_basis_fit.R          # Spatial basis function Gibbs sampler
```

### 2. **Data Files** (in `data/`)

```
data/
├── ca_pums_population.rds       # CA PUMS individual-level data
└── ca_puma_adjacency.RDA        # CA PUMA adjacency matrix (281 PUMAs, 793 edges)
```

### 3. **Main Execution Script**

```
test_sampling_designs_10configs.R  # Runs all 10 sampling designs
```

This script:
- Defines base configuration (response var, priors, MCMC settings)
- Loops through 10 configs
- For each config: sampling → full-data fits → DT (1-fold, 3-fold, 5-fold)
- Saves results to `_results_*_comparison/` directories

### 4. **Summary/Aggregation Scripts**

```
aggregate_additional_configs.R   # Aggregates 4 additional configs (skip ESIM)
analyze_all_dt_configs.R         # Computes MAD, top 3, DT vs DIC gap
```

### 5. **Analysis Output Files**

```
results_multi_config/
├── dt_all_10configs.RDS         # All DT results (10 configs)
└── oracle_all_10configs.RDS     # All oracle/DIC/WAIC results (10 configs)
```

---

## Reproduction Steps

### Step 1: Run All 10 Sampling Designs

```bash
./test_sampling_designs_10configs.R
```

**Runtime:** ~10-15 hours (depends on cores)

**Output:** 10 results directories:
```
_results_equal40_comparison/
_results_equal50_comparison/
_results_equal75_comparison/
_results_equal100_comparison/
_results_prop0.5pct_comparison/
_results_prop1pct_comparison/
_results_prop1p25pct_comparison/
_results_prop1p5pct_comparison/
_results_prop1p75pct_comparison/
_results_prop2pct_comparison/
```

Each directory contains:
```
comparison_XXX/
├── X.RDS                        # Predictors (from population)
├── z.RDS                        # Direct estimates (from sample)
├── d.RDS                        # Design variances (from sample)
├── fit_on_z/                    # Full-data fit (oracle)
│   └── chains.RDS
├── dt_1fold/                    # DT 1-fold results
│   ├── eps_0.30/
│   ├── eps_0.50/
│   └── eps_0.70/
├── dt_3fold/                    # DT 3-fold results
└── dt_5fold/                    # DT 5-fold results
```

### Step 2: Aggregate Results

For the 6 original configs (already done):
```
results_multi_config/dt_all.RDS
results_multi_config/oracle_all.RDS
```

For all 10 configs:
```bash
./aggregate_additional_configs.R
```

**Output:**
```
results_multi_config/dt_all_10configs.RDS
results_multi_config/oracle_all_10configs.RDS
```

### Step 3: Comprehensive Analysis

```bash
./analyze_all_dt_configs.R > analysis_output.txt
```

**Output:**
- Top 3 methods per config
- Best DT vs DIC gap per config
- Winner summary table

---

## Configuration Details

### The 10 Sampling Designs

| Config | Type | Parameter | Mean n/PUMA | Interpretation |
|--------|------|-----------|-------------|----------------|
| **equal_40** | Equal | n=40 | 40 | Low equal sample |
| **equal_50** | Equal | n=50 | 50 | Medium equal sample |
| **equal_75** | Equal | n=75 | 75 | Baseline |
| **equal_100** | Equal | n=100 | 100 | High equal sample |
| **prop_0.5pct** | PPS | 0.5% | ~85 | Very low % |
| **prop_1pct** | PPS | 1.0% | ~170 | Low % |
| **prop_1p25pct** | PPS | 1.25% | ~212 | Medium-low % |
| **prop_1p5pct** | PPS | 1.5% | ~255 | Medium % |
| **prop_1p75pct** | PPS | 1.75% | ~297 | Medium-high % |
| **prop_2pct** | PPS | 2.0% | ~340 | High % |

### Model Configuration

**Response:** Employment-to-population ratio (`employed`)
- Binary indicator (1 if employed, 0 otherwise)
- Includes all ages (children 0-15 in denominator)

**Model:** Spatial basis function model
- Intercept-only (no covariates)
- Fixed spatial basis (`spatial_type = "fixed"`)
- nbasis grid: [3, 6, 9, ..., 60] (20 values)

**Priors:** Weak priors
- c = 0.001
- d = 0.001

**MCMC:**
- nburn = 1500
- ndesired = 2000
- nthin = 1

**DT Configurations:**
- 1-fold: ε ∈ {0.3, 0.5, 0.7}, n_reps = 5
- 3-fold: 1 rep
- 5-fold: 1 rep

---

## Key Results (From Current Run)

### DT vs DIC Win Rate
- **DT Wins:** 7/10 configs (70%)
- **DIC Wins:** 3/10 configs (30%)

### Largest DT Advantages
1. **prop_1p75pct:** DT wins by 4.65 MAD
2. **equal_100:** DT wins by 3.90 MAD
3. **prop_2pct:** DT wins by 1.65 MAD

### Best Overall Performance
- **prop_1p25pct:** DT ε=0.7 n=5 achieves MAD=4.35 (best across all configs)

---

## File Dependencies

```
test_sampling_designs_10configs.R
  ├─ sim_functions/sampling_and_setup.R
  │    └─ data/ca_pums_population.rds
  ├─ sim_functions/full_data_fit.R
  │    ├─ sim_functions/spatial_basis_fit.R
  │    └─ data/ca_puma_adjacency.RDA
  └─ sim_functions/run_dt.R
       └─ sim_functions/spatial_basis_fit.R

aggregate_additional_configs.R
  ├─ sim_functions/summary_dt.R
  ├─ sim_functions/summary_oracle.R
  ├─ results_multi_config/dt_all.RDS (6 configs)
  └─ results_multi_config/oracle_all.RDS (6 configs)

analyze_all_dt_configs.R
  ├─ results_multi_config/dt_all_10configs.RDS
  └─ results_multi_config/oracle_all_10configs.RDS
```

---

## Notes

### Why No ESIM?
None of the 10 configs have ESIM results. The aggregation script (`aggregate_additional_configs.R`) explicitly skips ESIM to avoid errors during summary computation.

### Why Only MSE Loss?
The 4 additional configs (prop_1p25pct, prop_1p5pct, prop_1p75pct, equal_100) only have MSE loss computed to save time. The comprehensive analysis focuses on MSE-based DT.

### Missing ε=0.3 in Some Configs
The 4 additional configs only have ε ∈ {0.5, 0.7} to reduce computational burden. This doesn't affect main conclusions since ε=0.3 rarely wins.

---

## Troubleshooting

### Issue: "n_comp not found" error
**Cause:** Model config missing `n_comp` parameter
**Fix:** Ensure `model_config.RDS` has `n_comp` field

### Issue: ESIM errors during summary
**Cause:** No ESIM results exist
**Fix:** Use `aggregate_additional_configs.R` which skips ESIM

### Issue: "epsilon == 0.6" finds 0 rows
**Cause:** Floating point precision
**Fix:** Use `abs(epsilon - 0.6) < 0.001` for comparisons

---

## Citation

If using this analysis framework:
- See `CLAUDE.md` for project background
- See `COMPREHENSIVE_SUMMARY_10_CONFIGS.md` for detailed results
