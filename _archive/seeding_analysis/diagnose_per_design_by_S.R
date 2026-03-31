#!/usr/bin/env Rscript
# Per-design breakdown at S=50, 75, 100 (seeds 1-S)

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

mse_lookup <- oracle %>% select(config, comp_no, nbasis, mse_true)

all_sel <- bind_rows(
  make_sel(oracle %>% filter(!is.na(DIC)), DIC, "DIC"),
  make_sel(oracle %>% filter(!is.na(WAIC)), WAIC, "WAIC"),
  make_sel(dt_all %>% filter(loss_function == "plugin_NLL", abs(epsilon-EPS)<0.001, n_reps_used==REPS), metric, "DT-NLL"),
  make_sel(dt_all %>% filter(loss_function == "MSE", abs(epsilon-EPS)<0.001, n_reps_used==REPS), metric, "DT-MSE")
) %>%
  left_join(oracle_best, by = c("config", "comp_no")) %>%
  mutate(deviation = selected_nbasis - oracle_nbasis, abs_dev = abs(deviation)) %>%
  left_join(mse_lookup, by = c("config", "comp_no", "selected_nbasis" = "nbasis")) %>%
  rename(selected_mse = mse_true) %>%
  mutate(mse_penalty = selected_mse - oracle_mse)

for (S in c(50, 75, 100)) {
  d <- all_sel %>% filter(comp_no %in% 1:S)
  cat(sprintf("\n%s\n%s\n", paste0("=== Seeds 1-", S, " (S=", S, ") ==="), strrep("=", 30)))

  # Per design
  per_design <- d %>%
    group_by(config, method) %>%
    summarise(MAE = mean(abs_dev), bias = mean(deviation),
              mse_pen = mean(mse_penalty, na.rm = TRUE) * 1000,
              sd_dev = sd(deviation),
              .groups = "drop")

  overall <- d %>%
    group_by(method) %>%
    summarise(MAE = mean(abs_dev), bias = mean(deviation),
              mse_pen = mean(mse_penalty, na.rm = TRUE) * 1000,
              .groups = "drop")

  methods <- c("DIC", "DT-NLL", "DT-MSE", "WAIC")

  cat("\n  MAE by design:\n")
  cat(sprintf("  %-8s %8s %8s %8s %8s\n", "", methods[1], methods[2], methods[3], methods[4]))
  for (cfg in TRIO) {
    vals <- sapply(methods, function(m) per_design$MAE[per_design$config == cfg & per_design$method == m])
    best <- which.min(vals)
    formatted <- sprintf("%8.2f", vals)
    formatted[best] <- sprintf("  *%5.2f", vals[best])
    cat(sprintf("  %-8s %s %s %s %s\n", LABELS[cfg], formatted[1], formatted[2], formatted[3], formatted[4]))
  }
  ov <- sapply(methods, function(m) overall$MAE[overall$method == m])
  best_ov <- which.min(ov)
  fmt_ov <- sprintf("%8.2f", ov)
  fmt_ov[best_ov] <- sprintf("  *%5.2f", ov[best_ov])
  cat(sprintf("  %-8s %s %s %s %s\n", "Overall", fmt_ov[1], fmt_ov[2], fmt_ov[3], fmt_ov[4]))

  cat("\n  Bias by design:\n")
  cat(sprintf("  %-8s %8s %8s %8s %8s\n", "", methods[1], methods[2], methods[3], methods[4]))
  for (cfg in TRIO) {
    vals <- sapply(methods, function(m) per_design$bias[per_design$config == cfg & per_design$method == m])
    cat(sprintf("  %-8s %+8.2f %+8.2f %+8.2f %+8.2f\n", LABELS[cfg], vals[1], vals[2], vals[3], vals[4]))
  }
  ov_b <- sapply(methods, function(m) overall$bias[overall$method == m])
  cat(sprintf("  %-8s %+8.2f %+8.2f %+8.2f %+8.2f\n", "Overall", ov_b[1], ov_b[2], ov_b[3], ov_b[4]))

  cat("\n  MSE Penalty (x1000) by design:\n")
  cat(sprintf("  %-8s %8s %8s %8s %8s\n", "", methods[1], methods[2], methods[3], methods[4]))
  for (cfg in TRIO) {
    vals <- sapply(methods, function(m) per_design$mse_pen[per_design$config == cfg & per_design$method == m])
    best <- which.min(vals)
    formatted <- sprintf("%8.3f", vals)
    formatted[best] <- sprintf("  *%5.3f", vals[best])
    cat(sprintf("  %-8s %s %s %s %s\n", LABELS[cfg], formatted[1], formatted[2], formatted[3], formatted[4]))
  }
  ov_p <- sapply(methods, function(m) overall$mse_pen[overall$method == m])
  best_p <- which.min(ov_p)
  fmt_p <- sprintf("%8.3f", ov_p)
  fmt_p[best_p] <- sprintf("  *%5.3f", ov_p[best_p])
  cat(sprintf("  %-8s %s %s %s %s\n", "Overall", fmt_p[1], fmt_p[2], fmt_p[3], fmt_p[4]))

  # DIC at 2.25% specifically
  dic_225 <- per_design %>% filter(config == "prop_2p25pct", method == "DIC")
  nll_225 <- per_design %>% filter(config == "prop_2p25pct", method == "DT-NLL")
  cat(sprintf("\n  2.25%% spotlight: DIC MAE=%.2f (bias=%+.2f) vs DT-NLL MAE=%.2f (bias=%+.2f)  gap=%.2f\n",
              dic_225$MAE, dic_225$bias, nll_225$MAE, nll_225$bias, dic_225$MAE - nll_225$MAE))
}
