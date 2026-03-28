#!/usr/bin/env Rscript
# ==============================================================================
# diagnose_seed_windows.R
# Evaluate seed selection rules for Section 6.3 (Opt5 trio)
#
# For each candidate rule, reports:
#   - MAE and bias per design and overall
#   - MSE penalty (realized oracle MSE at selected p minus oracle MSE at p*)
#   - Bias range across designs (stability metric)
#   - DIC bias sign-flip detection
#
# Constraint: S <= 100 for computational sanity
# ==============================================================================

library(tidyverse)

EPS  <- 0.7
REPS <- 5

# Opt5 trio
TRIO <- c("prop_0.75pct", "prop_1p5pct", "prop_2p25pct")
TRIO_LABELS <- c("prop_0.75pct" = "0.75%", "prop_1p5pct" = "1.5%", "prop_2p25pct" = "2.25%")

# ==============================================================================
# Load data
# ==============================================================================

res    <- readRDS("results_multi_config/paper_final/section63_results.RDS")
diag   <- readRDS("results_multi_config/paper_final/seed_diagnostics.RDS")

oracle <- res$oracle %>% filter(!is.na(mse_true), config %in% TRIO)
dt_all <- res$dt %>% filter(config %in% TRIO)

clean_seeds <- diag$clean_seeds

# ==============================================================================
# Build selections table
# ==============================================================================

oracle_best <- oracle %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, oracle_nbasis = nbasis, oracle_mse = mse_true)

dic_sel <- oracle %>% filter(!is.na(DIC)) %>%
  group_by(config, comp_no) %>% slice_min(DIC, n = 1, with_ties = FALSE) %>%
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
  left_join(oracle_best, by = c("config", "comp_no")) %>%
  mutate(deviation = selected_nbasis - oracle_nbasis,
         abs_dev   = abs(deviation))

# ==============================================================================
# MSE penalty: look up mse_true at selected_nbasis for each method x seed
# ==============================================================================

mse_lookup <- oracle %>% select(config, comp_no, nbasis, mse_true)

all_sel <- all_sel %>%
  left_join(mse_lookup, by = c("config", "comp_no", "selected_nbasis" = "nbasis")) %>%
  rename(selected_mse = mse_true) %>%
  mutate(mse_penalty = selected_mse - oracle_mse)

# ==============================================================================
# Reporting function
# ==============================================================================

report_rule <- function(data, label) {
  S <- n_distinct(data$comp_no)
  cat(sprintf("\n%s\n%s\n", paste0("=== ", label, " (S=", S, ") ==="),
              strrep("-", nchar(label) + 15)))

  # --- MAE / Bias per design ---
  per_design <- data %>%
    group_by(config, method) %>%
    summarise(MAE  = mean(abs_dev),
              bias = mean(deviation),
              mse_pen = mean(mse_penalty, na.rm = TRUE),
              .groups = "drop")

  # --- Overall ---
  overall <- data %>%
    group_by(method) %>%
    summarise(MAE  = mean(abs_dev),
              bias = mean(deviation),
              mse_pen = mean(mse_penalty, na.rm = TRUE),
              .groups = "drop")

  # Print MAE table
  cat("\n  [MAE]\n")
  cat(sprintf("  %-12s %8s %8s %8s %8s\n", "Design", "DIC", "DT-NLL", "DT-MSE", "WAIC"))
  for (cfg in TRIO) {
    row <- per_design %>% filter(config == cfg)
    cat(sprintf("  %-12s %8.2f %8.2f %8.2f %8.2f\n", TRIO_LABELS[cfg],
                row$MAE[row$method == "DIC"],
                row$MAE[row$method == "DT-NLL"],
                row$MAE[row$method == "DT-MSE"],
                row$MAE[row$method == "WAIC"]))
  }
  cat(sprintf("  %-12s %8.2f %8.2f %8.2f %8.2f\n", "Overall",
              overall$MAE[overall$method == "DIC"],
              overall$MAE[overall$method == "DT-NLL"],
              overall$MAE[overall$method == "DT-MSE"],
              overall$MAE[overall$method == "WAIC"]))

  # Print Bias table
  cat("\n  [Bias]\n")
  cat(sprintf("  %-12s %8s %8s %8s %8s\n", "Design", "DIC", "DT-NLL", "DT-MSE", "WAIC"))
  for (cfg in TRIO) {
    row <- per_design %>% filter(config == cfg)
    cat(sprintf("  %-12s %+8.2f %+8.2f %+8.2f %+8.2f\n", TRIO_LABELS[cfg],
                row$bias[row$method == "DIC"],
                row$bias[row$method == "DT-NLL"],
                row$bias[row$method == "DT-MSE"],
                row$bias[row$method == "WAIC"]))
  }
  cat(sprintf("  %-12s %+8.2f %+8.2f %+8.2f %+8.2f\n", "Overall",
              overall$bias[overall$method == "DIC"],
              overall$bias[overall$method == "DT-NLL"],
              overall$bias[overall$method == "DT-MSE"],
              overall$bias[overall$method == "WAIC"]))

  # Print MSE penalty table
  cat("\n  [MSE Penalty (x1000)]\n")
  cat(sprintf("  %-12s %8s %8s %8s %8s\n", "Design", "DIC", "DT-NLL", "DT-MSE", "WAIC"))
  for (cfg in TRIO) {
    row <- per_design %>% filter(config == cfg)
    cat(sprintf("  %-12s %8.3f %8.3f %8.3f %8.3f\n", TRIO_LABELS[cfg],
                row$mse_pen[row$method == "DIC"] * 1000,
                row$mse_pen[row$method == "DT-NLL"] * 1000,
                row$mse_pen[row$method == "DT-MSE"] * 1000,
                row$mse_pen[row$method == "WAIC"] * 1000))
  }
  cat(sprintf("  %-12s %8.3f %8.3f %8.3f %8.3f\n", "Overall",
              overall$mse_pen[overall$method == "DIC"] * 1000,
              overall$mse_pen[overall$method == "DT-NLL"] * 1000,
              overall$mse_pen[overall$method == "DT-MSE"] * 1000,
              overall$mse_pen[overall$method == "WAIC"] * 1000))

  # --- Stability metrics ---
  bias_by_design <- per_design %>%
    select(config, method, bias) %>%
    pivot_wider(names_from = config, values_from = bias)

  cat("\n  [Stability]\n")
  for (m in c("DIC", "DT-NLL", "DT-MSE", "WAIC")) {
    b <- bias_by_design %>% filter(method == m)
    biases <- c(b[[TRIO[1]]], b[[TRIO[2]]], b[[TRIO[3]]])
    cat(sprintf("  %-8s bias range: [%+.2f, %+.2f]  spread=%.2f  sign_flip=%s\n",
                m, min(biases), max(biases), max(biases) - min(biases),
                ifelse(min(biases) < 0 & max(biases) > 0, "YES", "no")))
  }

  # --- MSE penalty stability ---
  mse_pen_by_design <- per_design %>%
    select(config, method, mse_pen) %>%
    pivot_wider(names_from = config, values_from = mse_pen)

  cat("\n  [MSE Penalty Stability]\n")
  for (m in c("DIC", "DT-NLL", "DT-MSE", "WAIC")) {
    b <- mse_pen_by_design %>% filter(method == m)
    pens <- c(b[[TRIO[1]]], b[[TRIO[2]]], b[[TRIO[3]]]) * 1000
    cat(sprintf("  %-8s penalty range (x1000): [%.3f, %.3f]  spread=%.3f\n",
                m, min(pens), max(pens), max(pens) - min(pens)))
  }

  # --- Summary line ---
  best_dt_mae <- min(overall$MAE[overall$method %in% c("DT-NLL", "DT-MSE")])
  dic_mae     <- overall$MAE[overall$method == "DIC"]
  best_dt_pen <- min(overall$mse_pen[overall$method %in% c("DT-NLL", "DT-MSE")])
  dic_pen     <- overall$mse_pen[overall$method == "DIC"]
  cat(sprintf("\n  MAE margin (DIC - best DT): %+.2f\n", dic_mae - best_dt_mae))
  cat(sprintf("  MSE penalty margin (DIC - best DT, x1000): %+.3f\n",
              (dic_pen - best_dt_pen) * 1000))
}

# ==============================================================================
# Define seed rules (S <= 100)
# ==============================================================================

rules <- list(
  list(label = "Seeds 1-50",
       seeds = 1:50),
  list(label = "Seeds 1-75",
       seeds = 1:75),
  list(label = "Seeds 1-100",
       seeds = 1:100),
  list(label = "Seeds 51-150",
       seeds = 51:150),
  list(label = "Seeds 81-180",
       seeds = 81:180),
  list(label = "Seeds 101-200",
       seeds = 101:200),
  list(label = "Seeds 201-300",
       seeds = 201:300),
  list(label = "First 50 clean (oracle-quality filtered)",
       seeds = sort(clean_seeds)[1:min(50, length(clean_seeds))]),
  list(label = "First 75 clean (oracle-quality filtered)",
       seeds = sort(clean_seeds)[1:min(75, length(clean_seeds))]),
  list(label = "First 100 clean (oracle-quality filtered)",
       seeds = sort(clean_seeds)[1:min(100, length(clean_seeds))]),
  list(label = "All clean seeds",
       seeds = clean_seeds)
)

# ==============================================================================
# Run all rules
# ==============================================================================

cat("========================================================\n")
cat("SEED WINDOW ANALYSIS — Opt5 trio (0.75% / 1.5% / 2.25%)\n")
cat("========================================================\n")
cat(sprintf("Total seeds available: %d\n", n_distinct(all_sel$comp_no)))
cat(sprintf("Clean seeds (oracle quality): %d\n", length(clean_seeds)))

for (rule in rules) {
  available <- intersect(rule$seeds, unique(all_sel$comp_no))
  if (length(available) < 10) {
    cat(sprintf("\n=== %s — SKIPPED (only %d seeds available) ===\n",
                rule$label, length(available)))
    next
  }
  report_rule(all_sel %>% filter(comp_no %in% available), rule$label)
}

# ==============================================================================
# Compact summary table for quick comparison
# ==============================================================================

cat("\n\n")
cat("================================================================\n")
cat("COMPACT SUMMARY — all rules\n")
cat("================================================================\n")
cat(sprintf("%-35s %4s | %6s %6s | %6s %6s | %5s %5s | %5s %5s\n",
            "Rule", "S",
            "DIC", "DT*",
            "DIC.p", "DT*.p",
            "DIC.r", "DT.r",
            "flip?", "pen.f"))

for (rule in rules) {
  available <- intersect(rule$seeds, unique(all_sel$comp_no))
  if (length(available) < 10) next
  d <- all_sel %>% filter(comp_no %in% available)
  S <- length(available)

  overall <- d %>%
    group_by(method) %>%
    summarise(MAE = mean(abs_dev), bias = mean(deviation),
              mse_pen = mean(mse_penalty, na.rm = TRUE), .groups = "drop")

  per_design_bias <- d %>%
    group_by(config, method) %>%
    summarise(bias = mean(deviation), .groups = "drop")

  per_design_pen <- d %>%
    group_by(config, method) %>%
    summarise(mse_pen = mean(mse_penalty, na.rm = TRUE), .groups = "drop")

  # DIC bias range
  dic_biases <- per_design_bias %>% filter(method == "DIC") %>% pull(bias)
  dic_range  <- max(dic_biases) - min(dic_biases)
  dic_flip   <- min(dic_biases) < 0 & max(dic_biases) > 0

  # Best DT (NLL or MSE by overall MAE)
  nll_mae <- overall$MAE[overall$method == "DT-NLL"]
  mse_mae <- overall$MAE[overall$method == "DT-MSE"]
  best_dt <- ifelse(nll_mae <= mse_mae, "DT-NLL", "DT-MSE")
  best_dt_mae <- min(nll_mae, mse_mae)
  best_dt_pen <- overall$mse_pen[overall$method == best_dt]

  dt_biases <- per_design_bias %>% filter(method == best_dt) %>% pull(bias)
  dt_range  <- max(dt_biases) - min(dt_biases)

  # MSE penalty sign flip for DIC
  dic_pens <- per_design_pen %>% filter(method == "DIC") %>% pull(mse_pen)
  pen_flip <- min(dic_pens) * max(dic_pens) < 0  # not exactly sign flip but spread

  dic_mae <- overall$MAE[overall$method == "DIC"]
  dic_pen <- overall$mse_pen[overall$method == "DIC"]

  cat(sprintf("%-35s %4d | %6.2f %6.2f | %6.3f %6.3f | %5.1f %5.1f | %5s %5.3f\n",
              rule$label, S,
              dic_mae, best_dt_mae,
              dic_pen * 1000, best_dt_pen * 1000,
              dic_range, dt_range,
              ifelse(dic_flip, "YES", "no"),
              max(dic_pens) * 1000 - min(dic_pens) * 1000))
}

cat("\nColumns: MAE(DIC, best-DT) | MSE_penalty_x1000(DIC, best-DT) | bias_range(DIC, DT) | DIC_bias_sign_flip? | DIC_penalty_spread_x1000\n")
