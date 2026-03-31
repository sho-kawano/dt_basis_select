#!/usr/bin/env Rscript
# Sweep S from 75 to 100 (seeds 1-S) for Opt5 trio
# Find the max S where both MAE and stability stories hold

library(tidyverse)

EPS  <- 0.7; REPS <- 5
TRIO <- c("prop_0.75pct", "prop_1p5pct", "prop_2p25pct")
TRIO_LABELS <- c("prop_0.75pct" = "0.75%", "prop_1p5pct" = "1.5%", "prop_2p25pct" = "2.25%")

res    <- readRDS("results_multi_config/paper_final/section63_results.RDS")
oracle <- res$oracle %>% filter(!is.na(mse_true), config %in% TRIO)
dt_all <- res$dt %>% filter(config %in% TRIO)

# Build selections
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

# Sweep function
sweep_S <- function(seeds) {
  d <- all_sel %>% filter(comp_no %in% seeds)
  S <- n_distinct(d$comp_no)

  overall <- d %>% group_by(method) %>%
    summarise(MAE = mean(abs_dev), bias = mean(deviation),
              mse_pen = mean(mse_penalty, na.rm = TRUE), .groups = "drop")

  per_design <- d %>% group_by(config, method) %>%
    summarise(bias = mean(deviation), .groups = "drop")

  get_bias_range <- function(m) {
    b <- per_design %>% filter(method == m) %>% pull(bias)
    list(min = min(b), max = max(b), spread = max(b) - min(b),
         flip = min(b) < 0 & max(b) > 0)
  }

  dic  <- get_bias_range("DIC")
  nll  <- get_bias_range("DT-NLL")
  mse  <- get_bias_range("DT-MSE")

  tibble(
    S = S,
    DIC_MAE   = overall$MAE[overall$method == "DIC"],
    NLL_MAE   = overall$MAE[overall$method == "DT-NLL"],
    MSE_MAE   = overall$MAE[overall$method == "DT-MSE"],
    MAE_margin_NLL = overall$MAE[overall$method == "DIC"] - overall$MAE[overall$method == "DT-NLL"],
    MAE_margin_MSE = overall$MAE[overall$method == "DIC"] - overall$MAE[overall$method == "DT-MSE"],
    DIC_pen   = overall$mse_pen[overall$method == "DIC"] * 1000,
    NLL_pen   = overall$mse_pen[overall$method == "DT-NLL"] * 1000,
    MSE_pen   = overall$mse_pen[overall$method == "DT-MSE"] * 1000,
    DIC_bias_spread  = dic$spread,
    NLL_bias_spread  = nll$spread,
    MSE_bias_spread  = mse$spread,
    DIC_flip  = dic$flip,
    NLL_flip  = nll$flip,
    MSE_flip  = mse$flip,
    DIC_bias_min = dic$min, DIC_bias_max = dic$max,
    NLL_bias_min = nll$min, NLL_bias_max = nll$max,
    MSE_bias_min = mse$min, MSE_bias_max = mse$max
  )
}

print_sweep <- function(results, label) {
  cat(sprintf("\n=== %s ===\n\n", label))
  cat(sprintf("%3s | %5s %5s %5s | %6s %6s | %5s %5s %5s | %4s %4s %4s | DIC bias        NLL bias        MSE bias\n",
              "S", "DIC", "NLL", "MSE", "m.NLL", "m.MSE",
              "D.sp", "N.sp", "M.sp", "D.fl", "N.fl", "M.fl"))
  cat(strrep("-", 140), "\n")
  for (i in 1:nrow(results)) {
    r <- results[i,]
    cat(sprintf("%3d | %5.2f %5.2f %5.2f | %+6.2f %+6.2f | %5.2f %5.2f %5.2f | %4s %4s %4s | [%+5.2f,%+5.2f] [%+5.2f,%+5.2f] [%+5.2f,%+5.2f]\n",
                r$S, r$DIC_MAE, r$NLL_MAE, r$MSE_MAE,
                r$MAE_margin_NLL, r$MAE_margin_MSE,
                r$DIC_bias_spread, r$NLL_bias_spread, r$MSE_bias_spread,
                ifelse(r$DIC_flip, "YES", "no"), ifelse(r$NLL_flip, "YES", "no"), ifelse(r$MSE_flip, "YES", "no"),
                r$DIC_bias_min, r$DIC_bias_max, r$NLL_bias_min, r$NLL_bias_max, r$MSE_bias_min, r$MSE_bias_max))
  }
  cat("\nColumns: MAE(DIC,NLL,MSE) | margin(DIC-NLL, DIC-MSE) | bias_spread | sign_flip | [bias_min, bias_max]\n")
  cat("Positive margin = DT wins on MAE\n")
}

# All starting points, S from 50 to 100
starts <- c(1, 51, 101, 151, 201)
all_seeds <- unique(all_sel$comp_no)

for (start in starts) {
  results <- map_dfr(50:100, function(S) {
    seeds <- seq(start, start + S - 1)
    seeds <- intersect(seeds, all_seeds)
    if (length(seeds) < 50) return(NULL)
    sweep_S(seeds)
  })
  if (nrow(results) > 0)
    print_sweep(results, sprintf("Start=%d, seeds %d to %d+S", start, start, start))
}
