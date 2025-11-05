# Model Selection via Data Thinning

Comparison of model selection methods for small area estimation using Fay-Herriot models with Illinois county data.

## Quick Start

### 1. Setup

```bash
# Install R packages
Rscript -e 'install.packages(c("tidyverse", "parallel", "doParallel", "Matrix", "SUMMER", "tidycensus", "mvtnorm", "LaplacesDemon"))'

# Configure Census API key
cp .Renviron.example .Renviron
# Edit .Renviron and add: CENSUS_API_KEY=your_key_here
# Get a key at: https://api.census.gov/data/key_signup.html
```

### 2. Run

```bash
# Full pipeline (50 comparisons, ~2-3 hours)
Rscript run_comparisons.R

# Or just recompute summaries from existing results (~5 min)
Rscript run_summaries_only.R
```

### 3. Analyze

```r
# Load results
results <- readRDS("results.RDS")

# Open analysis notebook
# analysis.Rmd
```

## What This Does

Compares model selection methods to choose the best Fay-Herriot model (1, 4, or 7 covariates):

**Methods compared:**
- **Data Thinning (DT)** - 1-fold and 5-fold variants with MSE, plugin_NLL, predictive_NLL losses
- **Empirical Simulation (ESIM)** - Standard and data fission validation
- **Information Criteria** - DIC, WAIC (full-data methods)

**Key finding:** DT 1-fold with ε=0.05 achieves 100% accuracy (50/50 comparisons) and perfect ranking correlation (τ=1.0) with the true model ordering.

## Saving Results

To save results from a run, set `SAVE_AS` at the top of the script:

```r
# In run_comparisons.R or run_summaries_only.R (line 17)
SAVE_AS <- "baseline_eps0.01-0.80"  # or NULL to skip
```

This copies results to `saved_results/<name>/`:
- `results.RDS`, `results_list.RDS`, `sim_config.RDS`
- `info.txt` (timestamp and basic info)

Load later:
```r
results <- readRDS("saved_results/baseline_eps0.01-0.80/results.RDS")
config <- readRDS("saved_results/baseline_eps0.01-0.80/sim_config.RDS")
```

## File Structure

```
├── run_comparisons.R          # Main pipeline (creates config + runs all methods)
├── run_summaries_only.R       # Recompute summaries only
├── analysis.Rmd               # Analysis notebook
├── sim_functions/             # Core implementations
│   ├── run_dt.R               # Data thinning
│   ├── run_esim.R             # Empirical simulation
│   ├── summary_dt.R           # DT metrics
│   ├── summary_esim.R         # ESIM metrics
│   └── zfit.R                 # Oracle fits
├── samplers/                  # MCMC samplers (FH, BYM, CAR, etc.)
├── data/                      # Illinois county data
├── _results/                  # Per-comparison results
└── saved_results/             # Archived runs
```

## Configuration

Edit parameters in `run_comparisons.R` (lines 82-88):

```r
n_comp <- 50                   # Number of comparisons
n_cores <- 11                  # Parallel cores
chosen_var <- d / y^2          # Noise variance (line 88)
```

For summaries, edit `run_summaries_only.R` (lines 25-27):

```r
epsilon_values <- c(0.01, 0.02, ..., 0.80)  # DT thinning parameters
repeat_counts <- c(1, 3, 5)                 # DT repetitions
loss_functions <- c("MSE", "plugin_NLL", "predictive_NLL")
```

## Troubleshooting

**"CENSUS_API_KEY not set"**
- Create `.Renviron` from `.Renviron.example`
- Restart R after editing

**Out of memory**
- Reduce `n_cores` in config
- Run fewer comparisons for testing

**Port conflicts**
- Don't create nested parallel clusters
- Pipeline handles parallelization automatically

## System Requirements

- R 4.0+
- 16GB+ RAM (for 11 cores)
- ~5GB disk space
- ~2-3 hours runtime (50 comparisons)

## Citation

For methodology details, see associated paper on data thinning for model selection.
