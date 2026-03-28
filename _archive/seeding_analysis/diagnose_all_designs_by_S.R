#!/usr/bin/env Rscript
# All 3 designs × 3 methods × 3 S values — MAE, bias, MSE penalty

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
  make_sel(dt_all %>% filter(loss_function == "plugin_NLL", abs(epsilon-EPS)<0.001, n_reps_used==REPS), metric, "DT-NLL"),
  make_sel(dt_all %>% filter(loss_function == "MSE", abs(epsilon-EPS)<0.001, n_reps_used==REPS), metric, "DT-MSE")
) %>%
  left_join(oracle_best, by = c("config", "comp_no")) %>%
  mutate(deviation = selected_nbasis - oracle_nbasis, abs_dev = abs(deviation)) %>%
  left_join(mse_lookup, by = c("config", "comp_no", "selected_nbasis" = "nbasis")) %>%
  rename(selected_mse = mse_true) %>%
  mutate(mse_penalty = selected_mse - oracle_mse)

methods <- c("DIC", "DT-NLL", "DT-MSE")

for (metric_name in c("MAE", "Bias", "MSE Penalty (x1000)")) {
  cat(sprintf("\n=== %s ===\n", metric_name))
  cat(sprintf("%-8s %-8s | %8s %8s %8s | %8s %8s %8s | %8s %8s %8s\n",
              "", "", "S=50", "", "", "S=75", "", "", "S=100", "", ""))
  cat(sprintf("%-8s %-8s | %8s %8s %8s | %8s %8s %8s | %8s %8s %8s\n",
              "Design", "", "DIC", "DT-NLL", "DT-MSE", "DIC", "DT-NLL", "DT-MSE", "DIC", "DT-NLL", "DT-MSE"))
  cat(strrep("-", 105), "\n")

  for (cfg in c(TRIO, "Overall")) {
    vals <- matrix(NA, 3, 3)  # S x method
    for (si in 1:3) {
      S <- c(50, 75, 100)[si]
      d <- all_sel %>% filter(comp_no %in% 1:S)
      if (cfg == "Overall") {
        summ <- d %>% group_by(method) %>%
          summarise(MAE = mean(abs_dev), bias = mean(deviation),
                    mse_pen = mean(mse_penalty, na.rm = TRUE) * 1000, .groups = "drop")
      } else {
        summ <- d %>% filter(config == cfg) %>% group_by(method) %>%
          summarise(MAE = mean(abs_dev), bias = mean(deviation),
                    mse_pen = mean(mse_penalty, na.rm = TRUE) * 1000, .groups = "drop")
      }
      for (mi in 1:3) {
        row <- summ %>% filter(method == methods[mi])
        vals[si, mi] <- switch(metric_name,
          "MAE" = row$MAE,
          "Bias" = row$bias,
          "MSE Penalty (x1000)" = row$mse_pen)
      }
    }

    label <- if (cfg == "Overall") "Overall" else LABELS[cfg]
    if (metric_name == "Bias") {
      cat(sprintf("%-17s | %+8.2f %+8.2f %+8.2f | %+8.2f %+8.2f %+8.2f | %+8.2f %+8.2f %+8.2f\n",
                  label, vals[1,1], vals[1,2], vals[1,3], vals[2,1], vals[2,2], vals[2,3], vals[3,1], vals[3,2], vals[3,3]))
    } else if (metric_name == "MSE Penalty (x1000)") {
      cat(sprintf("%-17s | %8.3f %8.3f %8.3f | %8.3f %8.3f %8.3f | %8.3f %8.3f %8.3f\n",
                  label, vals[1,1], vals[1,2], vals[1,3], vals[2,1], vals[2,2], vals[2,3], vals[3,1], vals[3,2], vals[3,3]))
    } else {
      cat(sprintf("%-17s | %8.2f %8.2f %8.2f | %8.2f %8.2f %8.2f | %8.2f %8.2f %8.2f\n",
                  label, vals[1,1], vals[1,2], vals[1,3], vals[2,1], vals[2,2], vals[2,3], vals[3,1], vals[3,2], vals[3,3]))
    }
  }
}
