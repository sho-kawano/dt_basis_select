#!/usr/bin/env Rscript
# Quick check: what do seeds 1-75 + 201-275 look like combined?

library(tidyverse)

EPS  <- 0.7; REPS <- 5
TRIO <- c("prop_0.75pct", "prop_1p5pct", "prop_2p25pct")
LABELS <- c("prop_0.75pct" = "0.75%", "prop_1p5pct" = "1.5%", "prop_2p25pct" = "2.25%")

res    <- readRDS("results_multi_config/paper_final/section63_results.RDS")
oracle <- res$oracle %>% filter(!is.na(mse_true), config %in% TRIO)
dt_all <- res$dt %>% filter(config %in% TRIO)

oracle_best <- oracle %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, oracle_nbasis = nbasis, oracle_mse = mse_true)

make_sel <- function(data, col, method_name) {
  data %>% group_by(config, comp_no) %>%
    slice_min({{col}}, n = 1, with_ties = FALSE) %>%
    select(config, comp_no, selected_nbasis = nbasis) %>%
    mutate(method = method_name)
}

all_sel <- bind_rows(
  make_sel(oracle %>% filter(!is.na(DIC)), DIC, "DIC"),
  make_sel(oracle %>% filter(!is.na(WAIC)), WAIC, "WAIC"),
  make_sel(dt_all %>% filter(loss_function == "plugin_NLL", abs(epsilon-EPS)<0.001, n_reps_used==REPS), metric, "DT-NLL"),
  make_sel(dt_all %>% filter(loss_function == "MSE", abs(epsilon-EPS)<0.001, n_reps_used==REPS), metric, "DT-MSE")
) %>%
  left_join(oracle_best, by = c("config", "comp_no")) %>%
  mutate(deviation = selected_nbasis - oracle_nbasis, abs_dev = abs(deviation))

methods <- c("DIC", "DT-NLL", "DT-MSE", "WAIC")

report <- function(seeds, label) {
  d <- all_sel %>% filter(comp_no %in% seeds)
  S <- n_distinct(d$comp_no)
  cat(sprintf("\n=== %s (S=%d) ===\n", label, S))

  per_design <- d %>% group_by(config, method) %>%
    summarise(MAE = mean(abs_dev), bias = mean(deviation), .groups = "drop")
  overall <- d %>% group_by(method) %>%
    summarise(MAE = mean(abs_dev), bias = mean(deviation), .groups = "drop")

  cat("\n  MAE:\n")
  cat(sprintf("  %-8s %8s %8s %8s %8s\n", "", "DIC", "DT-NLL", "DT-MSE", "WAIC"))
  for (cfg in TRIO) {
    vals <- sapply(methods, function(m) per_design$MAE[per_design$config == cfg & per_design$method == m])
    cat(sprintf("  %-8s %8.2f %8.2f %8.2f %8.2f\n", LABELS[cfg], vals[1], vals[2], vals[3], vals[4]))
  }
  ov <- sapply(methods, function(m) overall$MAE[overall$method == m])
  cat(sprintf("  %-8s %8.2f %8.2f %8.2f %8.2f\n", "Overall", ov[1], ov[2], ov[3], ov[4]))

  cat("\n  Bias:\n")
  cat(sprintf("  %-8s %8s %8s %8s %8s\n", "", "DIC", "DT-NLL", "DT-MSE", "WAIC"))
  for (cfg in TRIO) {
    vals <- sapply(methods, function(m) per_design$bias[per_design$config == cfg & per_design$method == m])
    cat(sprintf("  %-8s %+8.2f %+8.2f %+8.2f %+8.2f\n", LABELS[cfg], vals[1], vals[2], vals[3], vals[4]))
  }
  ov_b <- sapply(methods, function(m) overall$bias[overall$method == m])
  cat(sprintf("  %-8s %+8.2f %+8.2f %+8.2f %+8.2f\n", "Overall", ov_b[1], ov_b[2], ov_b[3], ov_b[4]))

  # Stability
  cat("\n  Stability:\n")
  for (m in methods) {
    biases <- per_design %>% filter(method == m) %>% pull(bias)
    cat(sprintf("  %-8s range: [%+.2f, %+.2f]  spread=%.2f  flip=%s\n",
                m, min(biases), max(biases), max(biases) - min(biases),
                ifelse(min(biases) < 0 & max(biases) > 0, "YES", "no")))
  }
  cat(sprintf("\n  MAE margin (DIC - DT-NLL): %+.2f\n", ov[1] - ov[2]))
}

report(1:75, "Primary: seeds 1-75")
report(201:275, "Replication: seeds 201-275")
report(c(1:75, 201:275), "Combined: seeds 1-75 + 201-275")
