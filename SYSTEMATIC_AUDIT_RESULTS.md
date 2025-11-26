# Systematic Audit Results

## 1. Oracle Screening Code - VALID

**Files:**
- `oracle_consistency_analysis/phase2_employed_fine.R` (main)
- `sim_functions/sampling_and_setup.R`
- `sim_functions/full_data_fit.R`
- `sim_functions/summary_oracle.R`
- `models/spatial_basis_fh.R`

**Data flow:**
1. `setup_comp()` → Creates theta_true.RDS, X.RDS, z.RDS, d.RDS
2. `full_data_fit()` → Creates fit_on_z/chains.RDS (one per nbasis)
3. `summary_oracle()` → Reads chains.RDS + theta_true.RDS → Returns MSE

**✓ Structurally valid**

## 2. Plotting Code - VALID

**Files:**
- `oracle_consistency_analysis/plot_phase2_results.R`
- Uses `summary_oracle()` to extract results
- Creates 3 plots: overlay, faceted, distribution

**✓ Correctly reads chains.RDS and theta_true.RDS**

## 3. Components Needing Attention

**HIGH PRIORITY:**
- Config key mismatch: phase2_employed_fine.R uses `population_data/adjacency_data`, but functions expect `population_file/adjacency_file`

**MEDIUM PRIORITY:**
- Oracle screening mislabeled as "ny_phase2" but uses CA data (281 PUMAs)
- Should rename directory for clarity

**LOW PRIORITY:**
- Document response_filter=NULL behavior (includes children)
- Add MCMC convergence checks (optional)

**✓ FIXED:**
- Removed silent CA data fallbacks from all functions

## 4. Minimal CA Oracle Test

**File:** `test_ca_oracle_minimal.R`

**What it does:**
1. Setup 20 comparisons with CA data
2. Run oracle fits for nbasis=[5,10,15,...,100]
3. Analyze and plot results
4. Minimal output - easy to review

**Config:**
- CA data (281 PUMAs)
- response_filter=NULL (employment-to-population ratio)
- equal_n=50
- No unnecessary printing

**To run:**
```bash
Rscript test_ca_oracle_minimal.R
```

**Expected:** Should match oracle screening (optimal nbasis~15, U-shaped curve)
