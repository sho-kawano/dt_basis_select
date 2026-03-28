# Response Variable Screening Criteria (Version 2)

**Updated:** November 2025
**Status:** Post-mortem revision after PUBCOV/rents_home flatness discovery

## Motivation

To evaluate data thinning for model selection in small area estimation, we need response variables that:
1. Have clear optimal nbasis values (interior minima)
2. Show **consistent optima** across independent samples
3. Have **steep enough MSE surfaces** that being wrong matters
4. Provide meaningful differentiation between methods

Version 1 criteria (MSE differentiation + consistency) were **insufficient** - they allowed flat curves where methods could disagree substantially without meaningful MSE penalty.

---

## Screening Pipeline (Updated)

### Phase 1: Single-Dataset Analysis (Per Response)

For each candidate response variable, generate **one test dataset** and fit models across nbasis range (10-100 by 10).

**Required Filters:**

#### Filter 1: Interior Minimum ✓
- **Metric:** Best nbasis location
- **Criterion:** Must be in range [20, 90]
- **Rationale:** Avoid boundary selections (indicates model mis-specification)

#### Filter 2: MSE Differentiation ✓
- **Metric:** `(max_MSE - min_MSE) / min_MSE × 100`
- **Criterion:** 15-30%
- **Rationale:** Overall signal must exist for methods to detect

#### Filter 3: **Curve Steepness** ⭐ NEW
- **Metric (Option A):** Average MSE penalty at ±10 nbasis units from optimal
- **Criterion:** **>5%** (strict) or **>3%** (moderate)
- **Metric (Option B):** Number of nbasis values within 5% of optimal MSE
- **Criterion:** **<3 models** (strict) or **<4** (moderate)
- **Rationale:** Ensures meaningful cost of selecting wrong nbasis
- **Why this matters:** Flat curves allow methods to disagree by ±20-30 units with only ~5% MSE penalty, making method comparison uninformative

**Output:** List of candidate response variables passing all filters

---

### Phase 2: Multi-Dataset Consistency (Per Candidate)

For each response that passed Phase 1, generate **10 independent datasets** to check consistency.

**Required Filters:**

#### Filter 4: Optimal nbasis Consistency ✓
- **Metric:** SD of optimal nbasis across 10 datasets
- **Criterion:** SD < 15
- **Rationale:** Oracle optimal should be stable across sampling variation

#### Filter 5: No Boundary Selections ✓
- **Metric:** % of datasets with boundary optimal (nbasis=10 or 100)
- **Criterion:** 0%
- **Rationale:** Boundary selections indicate dataset-specific pathologies

**Output:** Final response variable(s) for production runs

---

## Deprecated Criteria

### ❌ Removed: Step 0 - Spatial Correlation Screening

**Previous approach:**
- Compute spatial correlation of response at population level
- Only test responses with "sufficient" spatial structure

**Why removed:**
- Did not predict flat/steep curve behavior
- Spatial correlation doesn't guarantee the model will benefit from nbasis variation
- Empirical testing (Phase 1) is more direct and informative

---

## Lessons from PUBCOV & rents_home

### What Went Wrong:

Both PUBCOV and rents_home passed Version 1 criteria but failed in practice:

| Response   | Config    | MSE Diff | SD(nbasis) | Boundary | **Penalty @±10** | **Outcome** |
|------------|-----------|----------|------------|----------|------------------|-------------|
| PUBCOV     | equal_50  | 19.5%    | 13.2 ✓     | 0% ✓     | **3.3%** ✗       | Too flat    |
| rents_home | All 6     | ~20%     | ~15 ✓      | 0% ✓     | **~2-3%** ✗      | Too flat    |

### Key Insight:

**High MSE differentiation ≠ Steep curves**

- Differentiation measures **global range** (best vs. worst)
- Steepness measures **local penalty** (cost of being wrong by 10-20 units)
- Can have 20% global range but only 3% penalty at ±10 units (flat)

### What This Means:

For spatial basis function Fay-Herriot models with these responses:
- Oracle MSE surfaces are **relatively flat** in nbasis range 20-60
- Methods (DT, WAIC, ESIM) can disagree substantially (±20-30 nbasis)
- Yet achieve similar MSE performance (~5% penalty)
- **Difficult to meaningfully compare method quality**

---

## Recommended Strict Criteria (Summary)

| Filter | Metric | Criterion | Phase |
|--------|--------|-----------|-------|
| 1. Interior minimum | Best nbasis | ∈ [20, 90] | 1 |
| 2. MSE differentiation | (max-min)/min × 100 | 15-30% | 1 |
| **3. Steepness** ⭐ | **Penalty @±10 units** | **>5%** | **1** |
| 4. Consistency | SD(optimal nbasis) | <15 | 2 |
| 5. No boundaries | % boundary optimal | 0% | 2 |

**Note:** Can relax Filter 3 to >3% penalty if needed, but expect weaker method differentiation.

---

## Implementation Notes

**Steepness calculation:**
```r
# For each dataset:
oracle_results %>%
  group_by(comp_no) %>%
  mutate(
    min_mse = min(mse_true),
    optimal_nbasis = nbasis[which.min(mse_true)],
    dist_from_optimal = abs(nbasis - optimal_nbasis),
    pct_from_min = (mse_true / min_mse - 1) * 100
  ) %>%
  filter(dist_from_optimal == 10) %>%
  summarise(penalty_at_10 = mean(pct_from_min))
```

**Alternative steepness (width method):**
```r
# Count models within 5% of optimal
oracle_results %>%
  group_by(comp_no) %>%
  mutate(
    min_mse = min(mse_true),
    pct_from_min = (mse_true / min_mse - 1) * 100
  ) %>%
  summarise(n_within_5pct = sum(pct_from_min <= 5))
```

---

## Future Directions

If no response variables pass strict criteria, consider:

1. **Different model class:** Try non-spatial models or different spatial structures
2. **Different outcomes:** Test responses from other states or data sources
3. **Different nbasis range:** Test finer grid (nbasis=5-50 by 5) or different range
4. **Accept moderate steepness:** Use >3% threshold and report limitations
5. **Alternative research question:** Study robustness of spatial basis models to nbasis mis-specification (could be interesting finding!)

---

## References

- Original screening criteria: `_for_claude/response_hacking.md`
- Post-mortem analysis: `analyze_oracle_flatness.R`
- PUBCOV results: `_results_pubcov/`
- rents_home results: `_results_rents/`
