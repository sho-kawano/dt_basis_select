# Data Thinning for Model Validation in Small Area Estimation

Research project applying data thinning to model selection for Fay-Herriot models. See `_for_claude/draft_2026-03-26.tex` for the paper and `_for_claude/PIPELINE_REFERENCE.md` for pipeline details.

**Status:** Full draft complete. All empirical work done (ESIM backfill, aggregation, notebook, paper text updated).

**Paper narrative (Section 6.3):** DT is the most *stable* method — consistent conservative bias across all designs, no sign-flip. WAIC wins at 0.75%, DIC wins at 1.25%, ESIM wins at 1.75%. DT stays consistently negative bias — never over-selects. Do NOT frame as "DT wins overall." Do NOT crown an overall winner.

**START HERE for new sessions:** Read paper Sections 2.3 and 6 in `_for_claude/draft_2026-03-26.tex`.

## Coding Context

Solo project. Prioritize brevity/clarity over production-style error handling. Be extra careful with statistical reasoning -- ask if unsure.

**Data analysis:** Be skeptical, avoid overstatement, verify quantitative claims, prefer conservative interpretations.

**For parallel processing:** ALWAYS use FORK clusters (`type = "FORK"`).

**Floating point comparisons:** Always use tolerance with epsilon values:
```r
filter(abs(epsilon - 0.6) < 0.001)  # NOT: filter(epsilon == 0.6)
```

**Testing analysis notebooks:** Use `knitr::purl()` to extract code, not direct render:
```bash
Rscript -e "knitr::purl('notebook.Rmd', output='/tmp/test.R', quiet=TRUE); source('/tmp/test.R')"
```

## CRITICAL RULE: NEVER DELETE RESULTS DIRECTORIES

**NEVER use `unlink()`, `rm -rf`, or any method to delete `_results_*` directories.** They take hours to generate. Always create new, uniquely named directories. Ask if uncertain.

---

## Paper Structure (`_for_claude/draft_2026-02-22.tex`)

| Section | Title | Status |
|---------|-------|--------|
| 1 | Introduction | Written |
| 2 | Background (FH model, Gaussian DT, Motivating Example) | Written |
| 3 | MSE-Based Validation with Data Thinning | Written (theory + figures) |
| 4 | Repeated Thinning | Written |
| 5 | Likelihood-Based Validation | Written |
| 6 | Empirical Analysis and Model Comparison | Written |
| — | Discussion and Future Work | Written (last section, no number) |

**No Section 7.** What was previously called Section 7 is now Section 6.3.

### Paper framing
SAE -> DT. Solving a practical need in SAE, with generalizable insights about DT as secondary contribution.

### Section 6 structure
- **6.1** Simulation Framework (shared setup for all empirical sections)
- **6.2** Effect of epsilon and Repeats R — equal allocation designs, S=50, figures/plot4_basis.png
- **6.3** Comparison with Existing Methods — PA designs, DT vs DIC/WAIC/ESIM, figures/plot5_comparison.png

### Key theoretical results (Section 3)
- Thinning gap is monotonically decreasing in epsilon (Corollary on monotonicity)
- Thinning gap is model-dependent: complex models have larger gaps
- Fundamental tension: gap favors high epsilon, variance favors low epsilon
- No single epsilon is optimal across candidate models

### Ground truth / evaluation metrics (finalized)
- **Oracle basis:** per-sample argmin_p sum_i (theta_hat_i^(p) - theta_i)^2
- **MAE:** mean absolute error from oracle basis count (primary metric)
- **Mean Bias:** mean signed deviation (under-selection = negative)
- Candidate grid: p ∈ {3, 6, 9, ..., 60} (20 models)

### Recommendation (from Section 6.2)
epsilon ≈ 0.6–0.7 with R ≥ 5 repeats. Paper uses epsilon=0.7 for method comparison.

### Section 6.3 results (S=50, seeds 1–50, new trio)

| | 0.75% | 1.25% | 1.75% | Overall |
|---|---|---|---|---|
| DIC | 6.18 | **5.70** | 10.7 | 7.52 |
| DT-MSE | 7.74 | 5.82 | 9.54 | 7.70 |
| DT-NLL | 7.68 | 5.88 | 10.1 | 7.90 |
| ESIM | 9.60 | 7.44 | **7.38** | 8.14 |
| WAIC | **5.94** | 6.78 | 13.0 | 8.58 |

### Figure scripts
- `analysis/section2_spatial_illustration.Rmd` -- Section 2 figure
- `analysis/section3_figures.Rmd` -- Sections 3/4/5 gap/variance/tradeoff figures
- `analysis/section4_analysis.Rmd` -- Section 4 variance ratio analysis
- `analysis/section6_epsilon_analysis.Rmd` -- Section 6.2 epsilon/R analysis (uses equal alloc)
- `analysis/section6_methods_comparison.Rmd` -- Section 6.3 method comparison (uses PA alloc)
- `analysis/oracle_margin_analysis.R` -- Figure 3 (plot6_oracle_margin.png): oracle signal vs accuracy at 1.75% PA
- `analysis/supplement_thinning_maps.Rmd` -- Supplement figure (thinning maps)

---

## Empirical Setup

- **Population:** California PUMS, ~1.76M person records across 281 PUMAs (`data/ca_pums_population.rds`)
- **Response:** Employment-to-population ratio (`employed`), binary indicator
- **Adjacency:** `data/ca_puma_adjacency.RDA` (793 edges, rook contiguity)
- **Model:** Intercept-only Fay-Herriot with Moran's I spatial basis functions (fixed basis, `spatial_type="fixed"`, `X_covariates = NULL`)
- **Candidate models:** p ∈ {3, 6, 9, ..., 60} (20 models)
- **Priors:** c=d=0.001 (weak)
- **MCMC:** nburn=1500, ndesired=2000
- **X:** Intercept + spatial basis functions (no auxiliary covariates)
- **z, d:** Horvitz-Thompson direct estimates and design-based sampling variances via `survey::svymean()` (vary per sample)

### Sampling mechanism (`sim_functions/sampling_and_setup.R`)

Both equal and proportional allocation use the **same within-PUMA mechanism**: stratified Poisson PPS sampling with inclusion probabilities proportional to person weights (PWGTP), implemented via `sampling::UPpoisson()`. Design weights = 1/inclusion_prob.

- **Equal allocation:** Target sample size fixed at n_i = n for all PUMAs. Realized sizes vary around target (e.g., ±25% for n=40) due to Poisson mechanism.
- **Proportional (PPS) allocation:** Target n_i = max(floor(N_i × samp_frac), min_sample_size), so larger PUMAs get more samples.

### Survey designs used in paper

**Sections 2–6.2 (equal allocation, S=50):**
equal_30, equal_40, equal_50, equal_75, equal_100, equal_125

**Section 6.3 (PA allocation, S=50, new trio):**
prop_0.75pct (~47 avg n), prop_1p25pct (~78 avg n), prop_1p75pct (~110 avg n)

### Methods compared
- **Data Thinning (DT):** Repeated single-set thinning, eps=0.7 R=5
- **DIC / WAIC:** Information criteria benchmarks
- **ESIM:** Empirical simulation (L=100 iterations)
- **OD-Oracle:** Per-sample oracle using true population means (ground truth)

---

## Results & Data

### Paper final data (`results_summary/`)
- `equal_allocation_results.RDS` — S=50, 6 equal designs, DT + oracle + DIC/WAIC
- `thinned_oracle_s50.RDS` — S=50, thinned oracle MSE for gap analysis
- `methodcomp_results.RDS` — S=50, new trio (0.75%/1.25%/1.75%), DT + ESIM + oracle + DIC/WAIC
- `is_oracle_equal_alloc.RDS` — IS-oracle results for equal allocation (exploratory, Paper 2)

### Per-design results directories
```
# Equal allocation (S=50) — Sections 2-6.2
_results_equal30/    _results_equal40/    _results_equal50/
_results_equal75/    _results_equal100/   _results_equal125/

# PA designs — Section 6.3 trio (seeds 1-50)
_results_prop0.75pct/    _results_prop1p25pct/    _results_prop1p75pct/

# PA designs — kept for reference (not in trio)
_results_prop1p5pct/    _results_prop2p25pct/
```

NEVER DELETE any of these.

---

## Key Findings (Section 6.3, S=50, new trio, eps=0.7 R=5)

- **Per-design winners:** WAIC at 0.75%, DIC at 1.25%, ESIM at 1.75%. No single method dominates.
- **Stability:** DT-NLL stays consistently negative bias — never over-selects. DIC/WAIC bias shifts across designs.
- **ESIM:** Bad at low n (9.60 at 0.75%), best at high n (7.38 at 1.75%) — opposite arc from info criteria.
- **Oracle margin analysis:** At 1.75%, filtering to decisive oracle samples shows ESIM and DT improve dramatically while DIC/WAIC are flat — out-of-sample methods track the oracle signal, info criteria don't.

## Architecture

**Pipeline flow:** `setup_comp()` -> `full_data_fit()` -> `run_dt()` / `run_esim()` -> `summary_*()` functions

**Run + aggregate scripts:**
- `run_equal_allocation.R` + `aggregate_equal_allocation.R` → Sections 2–6.2
- `run_methodcomp.R` + `aggregate_methodcomp.R` → Section 6.3

**Per comparison folder:**
```
comparison_XXX/
  X.RDS, z.RDS, d.RDS     # Data
  fit_on_z/chains.RDS      # Full-data model fits
  dt_1fold/eps_*/           # DT results by epsilon
  emp_sim/                  # ESIM iterations
```

**Core functions:**
- `sim_functions/sampling_and_setup.R` -- Sample generation, setup, setup_esim
- `sim_functions/full_data_fit.R` -- Fit models on full data
- `sim_functions/run_dt.R` -- Data thinning
- `sim_functions/run_esim.R` -- Empirical simulation
- `sim_functions/summary_oracle.R` -- Oracle MSE, DIC, WAIC computation
- `sim_functions/summary_dt.R` -- DT loss computation
- `sim_functions/summary_esim.R` -- ESIM loss computation
- `models/spatial_basis_fh.R` -- Spatial basis FH model (Gibbs sampler)

**Analysis notebooks:**
- `analysis/section2_spatial_illustration.Rmd` -- Section 2 figure
- `analysis/section3_figures.Rmd` -- Section 3/4/5 gap/variance/tradeoff figures
- `analysis/section4_analysis.Rmd` -- Section 4 variance ratio
- `analysis/section6_epsilon_analysis.Rmd` -- Section 6.2 epsilon/R analysis (equal allocation)
- `analysis/section6_methods_comparison.Rmd` -- Section 6.3 method comparison (PA allocation)
- `analysis/oracle_margin_analysis.R` -- Figure plot6 oracle margin (1.75% PA)

See `_for_claude/PIPELINE_REFERENCE.md` for full pipeline documentation.
