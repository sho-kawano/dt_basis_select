#!/usr/bin/env Rscript
# ==============================================================================
# diagnose_oracle_seeds.R
# Identify seeds with clean oracle signals across all 4 PA designs
#
# Criteria (per seed x design):
#   1. high_absolute    — oracle_nbasis >= ABS_THRESH (absolute, design-independent)
#   2. high_outlier     — oracle_nbasis > mean + OUTLIER_SD*sd per design (relative)
#   3. double_min_high  — secondary local min at p > oracle, within DOUBLE_MIN_RATIO
#                         (only high-p secondary valleys; low-p secondary is fine)
#   4. all_methods_poor — best-case |selected - oracle| across methods
#                         exceeds mean + POOR_SD*sd for that design (relative)
#
# Goal: ~100 seeds passing all criteria across all 4 designs
# ==============================================================================

library(tidyverse)

ABS_THRESH       <- 36     # absolute oracle ceiling flag
OUTLIER_SD       <- 2.0    # z-score threshold for high oracle
DOUBLE_MIN_RATIO <- 1.05   # secondary local min (at higher p) within 5% of global
POOR_SD          <- 1.5    # all-methods-poor threshold (SDs above mean best-case dev)
EPS              <- 0.7
REPS             <- 5

# ==============================================================================
# Load
# ==============================================================================

res    <- readRDS("results_multi_config/paper_final/section63_results.RDS")
oracle <- res$oracle %>% filter(!is.na(mse_true))
dt_all <- res$dt

# ==============================================================================
# Step 1: Method selections per (seed, design)
# ==============================================================================

oracle_sel <- oracle %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, oracle_nbasis = nbasis)

dic_sel <- oracle %>% filter(!is.na(DIC)) %>%
  group_by(config, comp_no) %>% slice_min(DIC,  n = 1, with_ties = FALSE) %>%
  select(config, comp_no, selected_nbasis = nbasis) %>% mutate(method = "DIC")

waic_sel <- oracle %>% filter(!is.na(WAIC)) %>%
  group_by(config, comp_no) %>% slice_min(WAIC, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, selected_nbasis = nbasis) %>% mutate(method = "WAIC")

nll_sel <- dt_all %>%
  filter(loss_function == "plugin_NLL", abs(epsilon - EPS) < 0.001, n_reps_used == REPS) %>%
  group_by(config, comp_no) %>% slice_min(metric, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, selected_nbasis = nbasis) %>% mutate(method = "DT-NLL")

mse_sel <- dt_all %>%
  filter(loss_function == "MSE", abs(epsilon - EPS) < 0.001, n_reps_used == REPS) %>%
  group_by(config, comp_no) %>% slice_min(metric, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, selected_nbasis = nbasis) %>% mutate(method = "DT-MSE")

all_sel <- bind_rows(dic_sel, waic_sel, nll_sel, mse_sel) %>%
  left_join(oracle_sel, by = c("config", "comp_no")) %>%
  mutate(abs_dev = abs(selected_nbasis - oracle_nbasis))

# Best-case deviation per seed x design (min over DIC/DT-NLL/DT-MSE, excl WAIC)
best_dev <- all_sel %>%
  filter(method %in% c("DIC", "DT-NLL", "DT-MSE")) %>%
  group_by(config, comp_no) %>%
  summarise(best_abs_dev = min(abs_dev), .groups = "drop")

# ==============================================================================
# Step 2: Oracle curve diagnostics per (seed, design)
# ==============================================================================

seed_diag <- oracle %>%
  arrange(config, comp_no, nbasis) %>%
  group_by(config, comp_no) %>%
  summarise({
    ord  <- order(nbasis)
    vals <- mse_true[ord]
    ps   <- nbasis[ord]
    n    <- length(vals)
    gmin <- min(vals)
    gmin_idx <- which.min(vals)
    oracle_p <- ps[gmin_idx]

    # Local minima (interior)
    is_lm <- logical(n)
    if (n >= 3)
      is_lm[2:(n-1)] <- vals[2:(n-1)] < vals[1:(n-2)] & vals[2:(n-1)] < vals[3:n]
    # Include endpoints
    is_lm[1] <- vals[1] < vals[2]
    is_lm[n] <- vals[n] < vals[n-1]

    # Secondary minima only at p > oracle (high-p side)
    high_lm_vals <- vals[is_lm & ps > oracle_p]
    has_high_dbl <- length(high_lm_vals) > 0 &&
                    min(high_lm_vals) < DOUBLE_MIN_RATIO * gmin

    tibble(
      oracle_nbasis  = oracle_p,
      n_local_min    = sum(is_lm),
      n_high_lm      = sum(is_lm & ps > oracle_p),
      high_sec_ratio = if (length(high_lm_vals) > 0) min(high_lm_vals) / gmin else NA_real_,
      double_min_high = has_high_dbl
    )
  }, .groups = "drop")

# Add relative/absolute oracle flags and best-case deviation flag
seed_diag <- seed_diag %>%
  group_by(config) %>%
  mutate(
    oracle_mean    = mean(oracle_nbasis),
    oracle_sd      = sd(oracle_nbasis),
    oracle_z       = (oracle_nbasis - oracle_mean) / oracle_sd,
    high_absolute  = oracle_nbasis >= ABS_THRESH,
    high_outlier   = oracle_z > OUTLIER_SD
  ) %>%
  ungroup() %>%
  left_join(best_dev, by = c("config", "comp_no")) %>%
  group_by(config) %>%
  mutate(
    best_dev_mean  = mean(best_abs_dev),
    best_dev_sd    = sd(best_abs_dev),
    all_poor_z     = (best_abs_dev - best_dev_mean) / best_dev_sd,
    all_methods_poor = all_poor_z > POOR_SD
  ) %>%
  ungroup() %>%
  mutate(any_flag = high_absolute | high_outlier | double_min_high | all_methods_poor)

# ==============================================================================
# Step 3: Distributions and flag summaries
# ==============================================================================

cat("=== Oracle nbasis distribution per design ===\n")
seed_diag %>%
  group_by(config) %>%
  summarise(
    mean   = round(mean(oracle_nbasis), 1),
    sd     = round(sd(oracle_nbasis), 1),
    p25    = quantile(oracle_nbasis, 0.25),
    med    = median(oracle_nbasis),
    p75    = quantile(oracle_nbasis, 0.75),
    p90    = quantile(oracle_nbasis, 0.90),
    p95    = quantile(oracle_nbasis, 0.95),
    n_ge36 = sum(oracle_nbasis >= 36),
    n_ge42 = sum(oracle_nbasis >= 42)
  ) %>% print()

cat("\n=== Flag counts per design ===\n")
seed_diag %>%
  group_by(config) %>%
  summarise(
    n_high_abs    = sum(high_absolute),
    n_high_z      = sum(high_outlier),
    n_dbl_high    = sum(double_min_high),
    n_all_poor    = sum(all_methods_poor),
    n_any         = sum(any_flag),
    pct_flagged   = round(100 * mean(any_flag), 1)
  ) %>% print()

cat("\n=== Best-case deviation distribution per design ===\n")
seed_diag %>%
  group_by(config) %>%
  summarise(
    mean_best = round(mean(best_abs_dev), 1),
    sd_best   = round(sd(best_abs_dev), 1),
    p50       = round(median(best_abs_dev), 1),
    p75       = round(quantile(best_abs_dev, 0.75), 1),
    p90       = round(quantile(best_abs_dev, 0.90), 1),
    threshold = round(mean(best_abs_dev) + POOR_SD * sd(best_abs_dev), 1)
  ) %>% print()

# ==============================================================================
# Step 4: Sensitivity analysis
# ==============================================================================

cat("\n=== Sensitivity: n_clean seeds vs DOUBLE_MIN_RATIO (high-p only) ===\n")
for (ratio in c(1.03, 1.05, 1.08, 1.10, 1.15, Inf)) {
  n_clean <- seed_diag %>%
    mutate(flag = high_absolute | high_outlier | all_methods_poor |
             (!is.na(high_sec_ratio) & high_sec_ratio < ratio)) %>%
    group_by(comp_no) %>%
    summarise(bad = any(flag), .groups = "drop") %>%
    summarise(n = sum(!bad)) %>% pull(n)
  cat(sprintf("  double_min_ratio < %.2f: %d clean seeds\n", ratio, n_clean))
}

cat("\n=== Sensitivity: n_clean seeds vs POOR_SD threshold ===\n")
for (psd in c(1.0, 1.5, 2.0, 2.5, 3.0, Inf)) {
  n_clean <- seed_diag %>%
    mutate(flag = high_absolute | high_outlier | double_min_high |
             all_poor_z > psd) %>%
    group_by(comp_no) %>%
    summarise(bad = any(flag), .groups = "drop") %>%
    summarise(n = sum(!bad)) %>% pull(n)
  cat(sprintf("  poor_sd > %.1f: %d clean seeds\n", psd, n_clean))
}

# ==============================================================================
# Step 5: Per-seed summary
# ==============================================================================

seed_flags <- seed_diag %>%
  group_by(comp_no) %>%
  summarise(
    n_high_abs   = sum(high_absolute),
    n_high_z     = sum(high_outlier),
    n_dbl_high   = sum(double_min_high),
    n_all_poor   = sum(all_methods_poor),
    n_any_flag   = sum(any_flag),
    max_z        = round(max(oracle_z), 2),
    mean_oracle  = round(mean(oracle_nbasis), 1),
    max_oracle   = max(oracle_nbasis),
    mean_best_dev = round(mean(best_abs_dev), 1),
    .groups = "drop"
  ) %>%
  mutate(clean = n_any_flag == 0) %>%
  arrange(n_any_flag, mean_oracle)

cat(sprintf(
  "\n=== Seeds by n flags across 4 designs (abs>=%d, z>%.1f, dbl=%.2f, poor>%.1fSD) ===\n",
  ABS_THRESH, OUTLIER_SD, DOUBLE_MIN_RATIO, POOR_SD))
print(table(n_flags = seed_flags$n_any_flag))
cat(sprintf("Clean seeds (0 flags in all 4 designs): n = %d\n", sum(seed_flags$clean)))

cat("\n=== Flag breakdown among flagged seeds ===\n")
seed_flags %>%
  filter(!clean) %>%
  summarise(
    only_oracle_high = sum(n_any_flag > 0 & n_dbl_high == 0 & n_all_poor == 0),
    only_dbl_high    = sum(n_dbl_high > 0 & (n_high_abs + n_high_z) == 0 & n_all_poor == 0),
    only_all_poor    = sum(n_all_poor > 0 & n_dbl_high == 0 & (n_high_abs + n_high_z) == 0),
    any_dbl_high     = sum(n_dbl_high > 0),
    any_all_poor     = sum(n_all_poor > 0),
    any_oracle_flag  = sum((n_high_abs + n_high_z) > 0),
    total            = n()
  ) %>% print()

# ==============================================================================
# Step 6: MAE comparison — full S=300 vs clean seeds
# ==============================================================================

clean_seeds <- seed_flags %>% filter(clean) %>% pull(comp_no)

mae_table <- function(data) {
  data %>%
    group_by(config, method) %>%
    summarise(MAE  = round(mean(abs_dev), 2),
              bias = round(mean(selected_nbasis - oracle_nbasis), 1),
              n    = n(), .groups = "drop")
}

wide_mae <- function(df) {
  df %>%
    pivot_wider(names_from = method, values_from = c(MAE, bias, n),
                names_glue = "{method}_{.value}") %>%
    select(config, n = DIC_n,
           DIC_MAE, DIC_bias,
           `DT-NLL_MAE`, `DT-NLL_bias`,
           `DT-MSE_MAE`, `DT-MSE_bias`,
           WAIC_MAE)
}

cat("\n=== Per-design MAE: full S=300 ===\n")
mae_table(all_sel) %>% wide_mae() %>% print()

cat(sprintf("\n=== Per-design MAE: clean seeds (n=%d) ===\n", length(clean_seeds)))
mae_table(all_sel %>% filter(comp_no %in% clean_seeds)) %>% wide_mae() %>% print()

trio_summary <- function(data, label) {
  cat(sprintf("\n=== Trio overall — %s ===\n", label))
  for (trio in list(
    list(nm = "Opt1: 0.75/1.25/2.25", cfgs = c("prop_0.75pct","prop_1p25pct","prop_2p25pct")),
    list(nm = "Opt5: 0.75/1.5/2.25",  cfgs = c("prop_0.75pct","prop_1p5pct","prop_2p25pct"))
  )) {
    ov <- data %>%
      filter(config %in% trio$cfgs) %>%
      group_by(method) %>%
      summarise(MAE = round(mean(abs_dev), 2), .groups = "drop") %>%
      pivot_wider(names_from = method, values_from = MAE)
    best_dt <- min(ov$`DT-NLL`, ov$`DT-MSE`, na.rm = TRUE)
    cat(sprintf("  %s | DIC=%.2f  DT-NLL=%.2f  DT-MSE=%.2f  margin=%+.2f\n",
                trio$nm, ov$DIC, ov$`DT-NLL`, ov$`DT-MSE`, ov$DIC - best_dt))
  }
}

trio_summary(all_sel, "full S=300")
trio_summary(all_sel %>% filter(comp_no %in% clean_seeds),
             sprintf("clean seeds n=%d", length(clean_seeds)))

# ==============================================================================
# Save
# ==============================================================================

saveRDS(list(
  seed_diagnostics = seed_diag,
  seed_flags       = seed_flags,
  clean_seeds      = sort(clean_seeds),
  params = list(abs_thresh = ABS_THRESH, outlier_sd = OUTLIER_SD,
                double_min_ratio = DOUBLE_MIN_RATIO, poor_sd = POOR_SD)
), "results_multi_config/paper_final/seed_diagnostics.RDS")

cat(sprintf("\nSaved seed_diagnostics.RDS  —  %d clean seeds\n", length(clean_seeds)))
cat("Clean seed list:\n")
cat(paste(sort(clean_seeds), collapse = ", "), "\n")
