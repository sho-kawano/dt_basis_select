# Spatial Basis Function Selection via Data Thinning

A design-based simulation study comparing methods for selecting the number of
spatial basis functions in Fay-Herriot small area estimation models.

## Overview

We use California PUMS 2019–2023 (~1.76M person records, 281 PUMAs) as a
synthetic finite population, with employment-to-population rate as the target
parameter. We draw S=50 independent samples per design using stratified Poisson
PPS sampling (`sim_functions/sampling_and_setup.R`) and fit candidate
Fay-Herriot models with p ∈ {3, 6, …, 60} Moran's I spatial basis functions
via Gibbs sampling (`models/spatial_basis_fh.R`).

Two sampling designs are considered: equal allocation (target n_i per area: 30,
50, 75, 100) for tuning thinning parameters, and proportional-to-population
allocation (0.75%, 1.25%, 1.75% sampling rates) for method comparison. Horvitz-Thompson
direct estimates and design-based variance estimates are computed via
`survey::svymean()`.

Methods compared:
- **Data Thinning (DT):** out-of-sample validation via Gaussian data thinning (ε=0.7, R=5)
- **DIC / WAIC:** likelihood-based information criteria
- **ESIM:** empirical simulation / parametric bootstrap (L=100)

## Structure

| Directory | Contents |
|-----------|----------|
| `sim_functions/` | Core pipeline: sampling, model fitting, DT, ESIM, metrics |
| `models/` | Spatial basis Fay-Herriot Gibbs sampler |
| `p1_analysis/` | Paper analysis notebooks and figures (Sections 2–6) |
| `p2_analysis/` | Exploratory analysis (IS oracle) |
| `results_summary/` | Aggregated results (S=50) |
| `run_equal_allocation.R` | Driver: equal-allocation designs (Section 6.2) |
| `run_methodcomp.R` | Driver: proportional-allocation method comparison (Section 6.3) |

## Data

Population file and adjacency matrix are in `data/` (not tracked in git due to size).
