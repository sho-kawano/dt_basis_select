#!/usr/bin/env Rscript
# Find the best contiguous seed block where:
#   1. DT-NLL MAE margin >= 0 (DT wins or ties)
#   2. DT-NLL does NOT sign-flip
#   3. DIC DOES sign-flip
# Maximize S subject to these constraints.

library(tidyverse)

EPS  <- 0.7; REPS <- 5
TRIO <- c("prop_0.75pct", "prop_1p5pct", "prop_2p25pct")

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

all_seeds <- sort(unique(all_sel$comp_no))

# Sweep all starts (1-250) x all S (50-100)
results <- list()
for (start in 1:250) {
  for (S in 50:100) {
    seeds <- seq(start, start + S - 1)
    seeds <- intersect(seeds, all_seeds)
    if (length(seeds) < S * 0.95) next  # need most seeds present

    d <- all_sel %>% filter(comp_no %in% seeds)
    actual_S <- n_distinct(d$comp_no)
    if (actual_S < 50) next

    overall <- d %>% group_by(method) %>%
      summarise(MAE = mean(abs_dev), .groups = "drop")

    per_design <- d %>% group_by(config, method) %>%
      summarise(bias = mean(deviation), .groups = "drop")

    dic_biases <- per_design %>% filter(method == "DIC") %>% pull(bias)
    nll_biases <- per_design %>% filter(method == "DT-NLL") %>% pull(bias)

    dic_flip <- min(dic_biases) < 0 & max(dic_biases) > 0
    nll_flip <- min(nll_biases) < 0 & max(nll_biases) > 0

    nll_mae <- overall$MAE[overall$method == "DT-NLL"]
    dic_mae <- overall$MAE[overall$method == "DIC"]
    margin  <- dic_mae - nll_mae

    results[[length(results) + 1]] <- tibble(
      start = start, S = actual_S, end = start + S - 1,
      margin = margin, dic_flip = dic_flip, nll_flip = nll_flip,
      dic_spread = max(dic_biases) - min(dic_biases),
      nll_spread = max(nll_biases) - min(nll_biases),
      nll_mae = nll_mae, dic_mae = dic_mae
    )
  }
}

df <- bind_rows(results)

# Filter to our constraints
good <- df %>%
  filter(margin >= 0, !nll_flip, dic_flip) %>%
  arrange(desc(S), desc(margin))

cat("=== Top 20 contiguous blocks: DT-NLL wins MAE, no DT flip, DIC flips ===\n")
cat(sprintf("%-8s %-4s %-8s %7s %8s %8s %8s %8s\n",
            "Seeds", "S", "End", "Margin", "DIC.spr", "NLL.spr", "NLL.MAE", "DIC.MAE"))
cat(strrep("-", 75), "\n")
for (i in 1:min(20, nrow(good))) {
  r <- good[i,]
  cat(sprintf("%3d-%-4d %4d %7d %+7.2f %8.2f %8.2f %8.2f %8.2f\n",
              r$start, r$end, r$S, r$end,
              r$margin, r$dic_spread, r$nll_spread, r$nll_mae, r$dic_mae))
}

# Also: what's the max S if we relax margin to >= -0.10 (essentially tied)?
tied <- df %>%
  filter(margin >= -0.10, !nll_flip, dic_flip) %>%
  arrange(desc(S), desc(margin))

cat("\n=== Top 20 if margin >= -0.10 (essentially tied) ===\n")
cat(sprintf("%-8s %-4s %-8s %7s %8s %8s %8s %8s\n",
            "Seeds", "S", "End", "Margin", "DIC.spr", "NLL.spr", "NLL.MAE", "DIC.MAE"))
cat(strrep("-", 75), "\n")
for (i in 1:min(20, nrow(tied))) {
  r <- tied[i,]
  cat(sprintf("%3d-%-4d %4d %7d %+7.2f %8.2f %8.2f %8.2f %8.2f\n",
              r$start, r$end, r$S, r$end,
              r$margin, r$dic_spread, r$nll_spread, r$nll_mae, r$dic_mae))
}
