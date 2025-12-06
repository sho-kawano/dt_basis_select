#!/usr/bin/env Rscript
# ==============================================================================
# Six-Dimensional Analysis Framework
# ==============================================================================
# Unifying lenses across equal allocation and proportional designs

library(tidyverse)

cat("\n=== SIX-DIMENSIONAL ANALYSIS FRAMEWORK ===\n\n")

# Load data
oracle_all <- readRDS("results_multi_config/oracle_all.RDS")
dt_all <- readRDS("results_multi_config/dt_all.RDS")

configs <- c("equal_40", "equal_50", "prop_0.5pct", "prop_1pct", "prop_2pct")

config_to_dir <- c(
  "equal_40" = "_results_equal40_comparison",
  "equal_50" = "_results_equal50_comparison",
  "prop_0.5pct" = "_results_prop0.5pct_comparison",
  "prop_1pct" = "_results_prop1pct_comparison",
  "prop_2pct" = "_results_prop2pct_comparison"
)

# ==============================================================================
# DIMENSION 1: Signal-to-Noise (CV)
# ==============================================================================

cv_stats <- map_df(configs, function(cfg) {
  results_dir <- config_to_dir[cfg]
  if (!dir.exists(results_dir)) return(NULL)

  comp_dirs <- list.files(results_dir, pattern = "^comparison_", full.names = TRUE)

  comp_data <- map_df(comp_dirs, function(comp_dir) {
    z_list <- readRDS(file.path(comp_dir, "z.RDS"))
    d_list <- readRDS(file.path(comp_dir, "d.RDS"))
    z <- z_list$values
    d <- d_list$values

    tibble(
      config = cfg,
      cv = sqrt(d) / abs(z)
    )
  })

  comp_data %>%
    group_by(config) %>%
    summarise(
      mean_cv = mean(cv, na.rm = TRUE),
      median_cv = median(cv, na.rm = TRUE),
      .groups = "drop"
    )
})

# ==============================================================================
# DIMENSION 2: Oracle Signal (MSE variation)
# ==============================================================================

oracle_signal <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  summarise(
    mse_range = max(mse_true) - min(mse_true),
    mse_min = min(mse_true),
    .groups = "drop"
  ) %>%
  group_by(config) %>%
  summarise(
    mse_variation_pct = mean(mse_range / mse_min) * 100,
    .groups = "drop"
  )

# ==============================================================================
# DIMENSION 3: Variance Heterogeneity (d_ratio)
# ==============================================================================

heterogeneity <- map_df(configs, function(cfg) {
  results_dir <- config_to_dir[cfg]
  if (!dir.exists(results_dir)) return(NULL)

  comp_dirs <- list.files(results_dir, pattern = "^comparison_", full.names = TRUE)

  # Compute d_ratio per comparison, then average
  d_ratios <- map_dbl(comp_dirs, function(comp_dir) {
    d_list <- readRDS(file.path(comp_dir, "d.RDS"))
    d <- d_list$values
    max(d) / min(d)
  })

  tibble(
    config = cfg,
    mean_d_ratio = mean(d_ratios),
    median_d_ratio = median(d_ratios)
  )
})

# ==============================================================================
# DIMENSION 4: Oracle Stability (oracle SD)
# ==============================================================================

oracle_stability <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  group_by(config) %>%
  summarise(
    oracle_sd = sd(nbasis),
    .groups = "drop"
  )

# ==============================================================================
# DIMENSION 5: Model Complexity Support
# ==============================================================================

complexity_support <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  group_by(config) %>%
  summarise(
    mean_optimal_nbasis = mean(nbasis),
    median_optimal_nbasis = median(nbasis),
    .groups = "drop"
  )

# ==============================================================================
# DIMENSION 6: Selection Task Difficulty (composite)
# ==============================================================================

# Compute as: oracle_sd / (mse_variation_pct / 100)
# High value = unstable oracle despite having signal
task_difficulty <- oracle_stability %>%
  left_join(oracle_signal, by = "config") %>%
  mutate(
    difficulty_score = oracle_sd / (mse_variation_pct / 100)
  ) %>%
  select(config, difficulty_score)

# ==============================================================================
# COMBINE ALL DIMENSIONS
# ==============================================================================

all_dims <- cv_stats %>%
  left_join(oracle_signal, by = "config") %>%
  left_join(heterogeneity, by = "config") %>%
  left_join(oracle_stability, by = "config") %>%
  left_join(complexity_support, by = "config") %>%
  left_join(task_difficulty, by = "config")

# Add effective sample size
eff_n <- map_df(configs, function(cfg) {
  results_dir <- config_to_dir[cfg]
  if (!dir.exists(results_dir)) return(NULL)

  comp_dirs <- list.files(results_dir, pattern = "^comparison_", full.names = TRUE)

  eff_ns <- map_dbl(comp_dirs, function(comp_dir) {
    d_list <- readRDS(file.path(comp_dir, "d.RDS"))
    d <- d_list$values
    mean(1/d)
  })

  tibble(
    config = cfg,
    mean_eff_n = mean(eff_ns)
  )
})

all_dims <- all_dims %>%
  left_join(eff_n, by = "config")

# ==============================================================================
# METHOD PERFORMANCE
# ==============================================================================

# Get oracle selections
oracle_selections <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, oracle_nbasis = nbasis)

# DIC
dic_perf <- oracle_all %>%
  filter(!is.na(DIC)) %>%
  group_by(config, comp_no) %>%
  slice_min(DIC, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, nbasis) %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(dev = nbasis - oracle_nbasis) %>%
  group_by(config) %>%
  summarise(dic_mad = mean(abs(dev)), .groups = "drop")

# WAIC
waic_perf <- oracle_all %>%
  filter(!is.na(WAIC)) %>%
  group_by(config, comp_no) %>%
  slice_min(WAIC, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, nbasis) %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(dev = nbasis - oracle_nbasis) %>%
  group_by(config) %>%
  summarise(waic_mad = mean(abs(dev)), .groups = "drop")

# Best DT 1-fold
dt_best <- dt_all %>%
  filter(method_name == "dt_1fold") %>%
  group_by(config, comp_no, epsilon, n_reps_used) %>%
  slice_min(metric, n = 1, with_ties = FALSE) %>%
  left_join(oracle_selections, by = c("config", "comp_no")) %>%
  mutate(dev = nbasis - oracle_nbasis) %>%
  group_by(config, epsilon, n_reps_used) %>%
  summarise(mad = mean(abs(dev)), .groups = "drop") %>%
  group_by(config) %>%
  slice_min(mad, n = 1, with_ties = FALSE) %>%
  select(config, best_dt_mad = mad, best_dt_config = epsilon)

# Add performance
all_dims <- all_dims %>%
  left_join(dic_perf, by = "config") %>%
  left_join(waic_perf, by = "config") %>%
  left_join(dt_best, by = "config")

# ==============================================================================
# DISPLAY RESULTS
# ==============================================================================

cat("=== SIX DIMENSIONS (sorted by difficulty score) ===\n\n")
print(all_dims %>%
        arrange(difficulty_score) %>%
        select(config, mean_cv, mse_variation_pct, mean_d_ratio,
               oracle_sd, mean_optimal_nbasis, difficulty_score) %>%
        knitr::kable(digits = 2))

cat("\n\n=== WITH METHOD PERFORMANCE ===\n\n")
print(all_dims %>%
        arrange(difficulty_score) %>%
        select(config, difficulty_score, dic_mad, waic_mad, best_dt_mad) %>%
        knitr::kable(digits = 2))

cat("\n\n=== CORRELATIONS WITH DIFFICULTY SCORE ===\n\n")
cors <- tibble(
  dimension = c("CV", "MSE_variation", "d_ratio", "oracle_SD", "optimal_nbasis", "eff_n"),
  correlation_with_difficulty = c(
    cor(all_dims$mean_cv, all_dims$difficulty_score),
    cor(all_dims$mse_variation_pct, all_dims$difficulty_score),
    cor(all_dims$mean_d_ratio, all_dims$difficulty_score),
    cor(all_dims$oracle_sd, all_dims$difficulty_score),
    cor(all_dims$mean_optimal_nbasis, all_dims$difficulty_score),
    cor(all_dims$mean_eff_n, all_dims$difficulty_score)
  )
)
print(cors %>% knitr::kable(digits = 3))

cat("\n\n=== CORRELATIONS: DIMENSIONS vs. METHOD PERFORMANCE ===\n\n")

# Correlations with DIC MAD
cat("Correlations with DIC MAD:\n")
dic_cors <- tibble(
  dimension = c("CV", "MSE_var", "d_ratio", "oracle_SD", "difficulty"),
  correlation = c(
    cor(all_dims$mean_cv, all_dims$dic_mad),
    cor(all_dims$mse_variation_pct, all_dims$dic_mad),
    cor(all_dims$mean_d_ratio, all_dims$dic_mad),
    cor(all_dims$oracle_sd, all_dims$dic_mad),
    cor(all_dims$difficulty_score, all_dims$dic_mad)
  )
)
print(dic_cors %>% knitr::kable(digits = 3))

cat("\n\nCorrelations with WAIC MAD:\n")
waic_cors <- tibble(
  dimension = c("CV", "MSE_var", "d_ratio", "oracle_SD", "difficulty"),
  correlation = c(
    cor(all_dims$mean_cv, all_dims$waic_mad),
    cor(all_dims$mse_variation_pct, all_dims$waic_mad),
    cor(all_dims$mean_d_ratio, all_dims$waic_mad),
    cor(all_dims$oracle_sd, all_dims$waic_mad),
    cor(all_dims$difficulty_score, all_dims$waic_mad)
  )
)
print(waic_cors %>% knitr::kable(digits = 3))

cat("\n\nCorrelations with best DT MAD:\n")
dt_cors <- tibble(
  dimension = c("CV", "MSE_var", "d_ratio", "oracle_SD", "difficulty"),
  correlation = c(
    cor(all_dims$mean_cv, all_dims$best_dt_mad),
    cor(all_dims$mse_variation_pct, all_dims$best_dt_mad),
    cor(all_dims$mean_d_ratio, all_dims$best_dt_mad),
    cor(all_dims$oracle_sd, all_dims$best_dt_mad),
    cor(all_dims$difficulty_score, all_dims$best_dt_mad)
  )
)
print(dt_cors %>% knitr::kable(digits = 3))

cat("\n\n=== INTERPRETATION ===\n\n")
cat("Difficulty Score = oracle_SD / (MSE_variation / 100)\n")
cat("  - High score = unstable oracle despite having signal\n")
cat("  - Measures: 'How hard is selection given the data?'\n\n")

# Regime classification
all_dims_regime <- all_dims %>%
  mutate(
    regime = case_when(
      difficulty_score < 2.0 ~ "Easy",
      difficulty_score < 3.5 ~ "Medium",
      TRUE ~ "Hard"
    )
  )

cat("Selection Difficulty Regimes:\n\n")
print(all_dims_regime %>%
        select(config, difficulty_score, regime, dic_mad, waic_mad, best_dt_mad) %>%
        arrange(difficulty_score) %>%
        knitr::kable(digits = 2))

cat("\n")
