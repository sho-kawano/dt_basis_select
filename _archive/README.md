# Documentation for Claude

This folder contains background information and analysis documentation for the project.

## Active Documents

### 📄 background.tex
- Mathematical background on data thinning methodology
- Oracle definitions (FP-Oracle, OD-Oracle, TD-Oracle)
- Model specification and evaluation metrics
- **Status:** Current reference

### 📄 response_hacking.md ⭐
- **MOST COMPREHENSIVE** - Complete history of response variable exploration
- Documents discovery of weak priors + continuous responses
- Final recommendations: WAGP (27.7%), POVPIP (16.6%), rents_home (11.1%)
- Includes all experiments from Jan 2025
- **Status:** **PRIMARY REFERENCE** for response variable selection

### 📄 Data Thinning Paper.pdf
- Original data thinning methodology paper
- Theoretical foundation for the project
- **Status:** Reference material

## Archive / Historical

### 📄 response_variable_analysis.md
- Initial analysis from early exploration (Jan 12, 2025)
- Analyzed population heterogeneity (PUMA SD)
- **Status:** SUPERSEDED by response_hacking.md
- Kept for historical reference only

### 📄 architecture_change_analysis.md
- Documents transition from pre-aggregated to on-the-fly sampling
- Explains predictor matrix X decisions (population vs. estimated)
- **Status:** Historical record of architecture decisions

---

## Quick Reference

**For response variable selection:** Read `response_hacking.md` (comprehensive)

**For methodology background:** Read `background.tex`

**For data generation approach:** Read `architecture_change_analysis.md`

---

## Key Findings Summary

From `response_hacking.md`:

1. **Continuous responses dominate binary** (2-3× better differentiation)
2. **Weak priors essential** (a=b=c=d ≤ 0.001) to induce overfitting
3. **Fewer covariates = more spatial signal** (1-2 optimal)
4. **Winner:** WAGP (wages) - 27.7% OD-Oracle MSE differentiation
5. **Runner-up:** POVPIP (poverty) - 16.6% differentiation
6. **Best binary:** rents_home - 11.1% with 2 covariates

All U-shaped curves with clear interior optima.
