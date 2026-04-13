# Spatial Basis Function Selection via Data Thinning

A design-based simulation study comparing methods for selecting the number of
spatial basis functions in Fay-Herriot small area estimation models.

**Paper:** Kawano, Parker, Li (2026). *On Data Thinning for Model Validation in Small Area Estimation.* [arXiv:2604.04141](https://arxiv.org/abs/2604.04141)

## Overview

We use California PUMS 2019–2023 (~1.76M person records, 281 PUMAs) as a
synthetic finite population, with employment-to-population rate as the target
parameter. We draw S=50 independent samples per design using stratified Poisson
PPS sampling (`sim_functions/sampling_and_setup.R`) and fit candidate
Fay-Herriot models with p ∈ {3, 6, …, 60} Moran's I spatial basis functions
via Gibbs sampling (`models/spatial_basis_fh.R`).

Two sampling designs are considered. Equal allocation designs (target n_i per
area: 30, 40, 50, 75, 100, 125) are used to study how the thinning fraction ε and
number of repeats R affect model selection performance. Proportional-to-population
allocation designs (0.75%, 1.25%, 1.75% sampling rates) are used to benchmark
DT against existing methods under realistic sample size variation.
Horvitz-Thompson direct estimates and design-based variance estimates are
computed via `survey::svymean()`.

Benchmarks in method comparison:
- **DIC / WAIC:** likelihood-based information criteria
- **ESIM:** empirical simulation method adding noise to direct estimates treated as truth (L=100)

## Structure

| Directory / File | Contents |
|------------------|----------|
| `sim_functions/` | Core pipeline: sampling, model fitting, DT, ESIM, metrics |
| `models/` | Spatial basis Fay-Herriot Gibbs sampler |
| `core_analysis/` | Paper analysis notebooks and figures (Sections 2–6) |
| `results_summary/` | Aggregated results (S=50) |
| `run_equal_allocation.R` | Driver: equal-allocation designs (Section 6.2) |
| `aggregate_equal_allocation.R` | Aggregates equal-allocation results → `results_summary/` |
| `run_methodcomp.R` | Driver: proportional-allocation method comparison (Section 6.3) |
| `aggregate_methodcomp.R` | Aggregates method comparison results → `results_summary/` |

## Data

Population file and adjacency matrix are in `data/` (not tracked in git due to size).
