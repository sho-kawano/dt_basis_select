#!/usr/bin/env Rscript
# ==============================================================================
# Comprehensive DIC vs DT diagnostic across all available S=50 designs
# Per-comparison AND per-area analysis
# ==============================================================================

library(tidyverse)

EPS <- 0.7; REPS <- 5

# ==============================================================================
# HELPERS
# ==============================================================================

# Compute per-comparison selections + MSE excess for all methods
per_comp_metrics <- function(dt_data, oracle_data, config_label) {
  oracle_sel <- oracle_data %>%
    filter(!is.na(mse_true)) %>%
    group_by(config, comp_no) %>%
    slice_min(mse_true, n=1, with_ties=FALSE) %>%
    select(config, comp_no, oracle_nbasis=nbasis, mse_oracle=mse_true)

  mse_lkp <- oracle_data %>%
    filter(!is.na(mse_true)) %>%
    select(config, comp_no, nbasis, mse_true)

  dt_filt <- dt_data %>%
    filter(abs(epsilon - EPS) < 0.001, n_reps_used == REPS)

  dt_mse <- dt_filt %>% filter(loss_function == "MSE") %>%
    group_by(config, comp_no) %>%
    slice_min(metric, n=1, with_ties=FALSE) %>%
    select(config, comp_no, dt_mse_nbasis=nbasis)

  dt_nll <- dt_filt %>% filter(loss_function == "plugin_NLL") %>%
    group_by(config, comp_no) %>%
    slice_min(metric, n=1, with_ties=FALSE) %>%
    select(config, comp_no, dt_nll_nbasis=nbasis)

  dic_sel <- oracle_data %>%
    filter(!is.na(DIC)) %>%
    group_by(config, comp_no) %>%
    slice_min(DIC, n=1, with_ties=FALSE) %>%
    select(config, comp_no, dic_nbasis=nbasis)

  waic_sel <- oracle_data %>%
    filter(!is.na(WAIC)) %>%
    group_by(config, comp_no) %>%
    slice_min(WAIC, n=1, with_ties=FALSE) %>%
    select(config, comp_no, waic_nbasis=nbasis)

  oracle_sel %>%
    left_join(dic_sel,  by=c("config","comp_no")) %>%
    left_join(dt_mse,   by=c("config","comp_no")) %>%
    left_join(dt_nll,   by=c("config","comp_no")) %>%
    left_join(waic_sel, by=c("config","comp_no")) %>%
    left_join(mse_lkp %>% rename(dic_mse=mse_true),
              by=c("config","comp_no","dic_nbasis"="nbasis")) %>%
    left_join(mse_lkp %>% rename(dt_mse_mse=mse_true),
              by=c("config","comp_no","dt_mse_nbasis"="nbasis")) %>%
    left_join(mse_lkp %>% rename(dt_nll_mse=mse_true),
              by=c("config","comp_no","dt_nll_nbasis"="nbasis")) %>%
    left_join(mse_lkp %>% rename(waic_mse=mse_true),
              by=c("config","comp_no","waic_nbasis"="nbasis")) %>%
    mutate(
      dic_excess     = (dic_mse     - mse_oracle) / mse_oracle * 100,
      dt_mse_excess  = (dt_mse_mse  - mse_oracle) / mse_oracle * 100,
      dt_nll_excess  = (dt_nll_mse  - mse_oracle) / mse_oracle * 100,
      waic_excess    = (waic_mse    - mse_oracle) / mse_oracle * 100,
      dt_nll_wins    = dt_nll_excess < dic_excess,
      dic_wins       = dic_excess < dt_nll_excess,
      dt_nll_fail    = dt_nll_excess > 10,
      dic_fail       = dic_excess > 10,
      both_fail      = dt_nll_fail & dic_fail,
      dataset        = config_label
    )
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading data...\n")
pa_old  <- readRDS("results_multi_config/paper_final/pa_method_comparison.RDS")
pa_new  <- readRDS("results_multi_config/paper_final/pa_method_comparison_final.RDS")
eq      <- readRDS("results_multi_config/paper_final/equal_allocation_results.RDS")

# Combine PA datasets (old: 0.75/1.25/1.75; new: 0.75/1.25/1.5)
dt_all <- bind_rows(
  pa_old$dt  %>% mutate(dataset="PA-old"),
  pa_new$dt  %>% mutate(dataset="PA-new"),
  eq$dt_1fold %>% mutate(dataset="Equal")
)
ora_all <- bind_rows(
  pa_old$oracle %>% mutate(dataset="PA-old"),
  pa_new$oracle %>% mutate(dataset="PA-new"),
  eq$oracle     %>% mutate(dataset="Equal")
)

# ==============================================================================
# PER-COMPARISON METRICS
# ==============================================================================

cat("Computing per-comparison metrics...\n")

pc_pa_old <- per_comp_metrics(pa_old$dt,    pa_old$oracle,    "PA-old")
pc_pa_new <- per_comp_metrics(pa_new$dt,    pa_new$oracle,    "PA-new")
pc_eq     <- per_comp_metrics(eq$dt_1fold,  eq$oracle,        "Equal")

all_pc <- bind_rows(pc_pa_old, pc_pa_new, pc_eq)

# ==============================================================================
# SUMMARY TABLE BY DESIGN
# ==============================================================================

cat("\n=== Summary by design: MSE excess% (DIC | DT-NLL) ===\n")
all_pc %>%
  group_by(dataset, config) %>%
  summarise(
    n          = n(),
    dic_excess = mean(dic_excess, na.rm=TRUE),
    dtnll_excess = mean(dt_nll_excess, na.rm=TRUE),
    dic_mae    = mean(abs(dic_nbasis - oracle_nbasis), na.rm=TRUE),
    dtnll_mae  = mean(abs(dt_nll_nbasis - oracle_nbasis), na.rm=TRUE),
    pct_dtnll_wins = mean(dt_nll_wins, na.rm=TRUE)*100,
    pct_dic_wins   = mean(dic_wins, na.rm=TRUE)*100,
    pct_dtnll_fail = mean(dt_nll_fail, na.rm=TRUE)*100,
    pct_dic_fail   = mean(dic_fail, na.rm=TRUE)*100,
    .groups="drop"
  ) %>%
  mutate(
    winner = ifelse(dic_excess < dtnll_excess, "DIC", "DT-NLL"),
    cell   = sprintf("%.1f%% | %.1f%%  [DT wins %.0f%%]",
                     dic_excess, dtnll_excess, pct_dtnll_wins)
  ) %>%
  select(dataset, config, n, cell, pct_dtnll_fail, pct_dic_fail) %>%
  arrange(dataset, config) %>%
  print(n=50)

# ==============================================================================
# WHERE DT FAILS (excess > 10%) — what's the oracle nbasis?
# ==============================================================================

cat("\n=== DT-NLL failures (excess > 10%): oracle nbasis distribution ===\n")
all_pc %>%
  filter(dt_nll_fail) %>%
  group_by(dataset, config) %>%
  summarise(
    n_fail       = n(),
    oracle_mean  = mean(oracle_nbasis),
    oracle_max   = max(oracle_nbasis),
    dt_nll_sel   = mean(dt_nll_nbasis),
    dic_sel      = mean(dic_nbasis),
    dic_excess_when_dt_fails = mean(dic_excess),
    .groups="drop"
  ) %>%
  arrange(dataset, config) %>%
  print(n=50)

# ==============================================================================
# WHERE DIC FAILS (excess > 10%) — what's the oracle nbasis?
# ==============================================================================

cat("\n=== DIC failures (excess > 10%): oracle nbasis distribution ===\n")
all_pc %>%
  filter(dic_fail) %>%
  group_by(dataset, config) %>%
  summarise(
    n_fail       = n(),
    oracle_mean  = mean(oracle_nbasis),
    oracle_max   = max(oracle_nbasis),
    dic_sel      = mean(dic_nbasis),
    dt_nll_sel   = mean(dt_nll_nbasis),
    dt_excess_when_dic_fails = mean(dt_nll_excess),
    .groups="drop"
  ) %>%
  arrange(dataset, config) %>%
  print(n=50)

# ==============================================================================
# WHERE DT WINS BIG (DT excess < DIC excess by > 5pp)
# ==============================================================================

cat("\n=== DT-NLL wins by >5pp over DIC ===\n")
all_pc %>%
  mutate(dt_advantage = dic_excess - dt_nll_excess) %>%
  filter(dt_advantage > 5) %>%
  group_by(dataset, config) %>%
  summarise(
    n            = n(),
    avg_advantage = mean(dt_advantage),
    oracle_mean  = mean(oracle_nbasis),
    dic_sel      = mean(dic_nbasis),
    dt_nll_sel   = mean(dt_nll_nbasis),
    .groups="drop"
  ) %>%
  arrange(desc(avg_advantage)) %>%
  print(n=50)

# ==============================================================================
# WHERE DIC WINS BIG (DIC excess < DT excess by > 5pp)
# ==============================================================================

cat("\n=== DIC wins by >5pp over DT-NLL ===\n")
all_pc %>%
  mutate(dic_advantage = dt_nll_excess - dic_excess) %>%
  filter(dic_advantage > 5) %>%
  group_by(dataset, config) %>%
  summarise(
    n             = n(),
    avg_advantage = mean(dic_advantage),
    oracle_mean   = mean(oracle_nbasis),
    dic_sel       = mean(dic_nbasis),
    dt_nll_sel    = mean(dt_nll_nbasis),
    .groups="drop"
  ) %>%
  arrange(desc(avg_advantage)) %>%
  print(n=50)

# ==============================================================================
# PER-AREA ANALYSIS — targeted at worst DT failure comps
# ==============================================================================

cat("\n=== Worst DT-NLL failures (top 10 by excess) ===\n")
worst_dt <- all_pc %>%
  arrange(desc(dt_nll_excess)) %>%
  select(dataset, config, comp_no, oracle_nbasis, dt_nll_nbasis, dic_nbasis,
         dt_nll_excess, dic_excess) %>%
  head(10)
print(worst_dt)

cat("\n=== Worst DIC failures (top 10 by excess) ===\n")
worst_dic <- all_pc %>%
  arrange(desc(dic_excess)) %>%
  select(dataset, config, comp_no, oracle_nbasis, dic_nbasis, dt_nll_nbasis,
         dic_excess, dt_nll_excess) %>%
  head(10)
print(worst_dic)

# ==============================================================================
# PER-AREA: load chains for a few worst DT failures and check PUMA-level MSE
# ==============================================================================

cat("\n=== Per-area MSE: worst DT failures ===\n")

# Map config to results directory
dir_map <- c(
  prop_0.75pct = "_results_prop0.75pct",
  prop_1p25pct = "_results_prop1p25pct",
  prop_1p75pct = "_results_prop1p75pct_comparison",
  prop_1p5pct  = "_results_prop1p5pct_final",
  equal_30     = "_results_equal30",
  equal_40     = "_results_equal40",
  equal_50     = "_results_equal50",
  equal_75     = "_results_equal75",
  equal_100    = "_results_equal100",
  equal_125    = "_results_equal125"
)

examine_comp <- function(config, comp_no, dt_sel, dic_sel, oracle_sel) {
  res_dir <- dir_map[config]
  if (is.na(res_dir)) return(NULL)

  comp_folder <- file.path(res_dir, sprintf("comparison_%03d", comp_no))
  chains_file <- file.path(comp_folder, "fit_on_z", "chains.RDS")
  if (!file.exists(chains_file)) return(NULL)

  chains    <- readRDS(chains_file)
  theta_true <- readRDS(file.path(res_dir, "theta_true.RDS"))

  nbasis_vals <- sapply(chains, function(x) x$nbasis)
  truth <- theta_true$values

  get_puma_mse <- function(p) {
    idx <- which(nbasis_vals == p)
    if (length(idx) == 0) return(NULL)
    est <- colMeans(chains[[idx]]$chain$theta)
    (est - truth)^2
  }

  oracle_puma <- get_puma_mse(oracle_sel)
  dt_puma     <- get_puma_mse(dt_sel)
  dic_puma    <- get_puma_mse(dic_sel)

  if (is.null(oracle_puma) || is.null(dt_puma) || is.null(dic_puma)) return(NULL)

  tibble(
    config=config, comp_no=comp_no,
    puma=seq_along(truth),
    oracle_nbasis=oracle_sel, dt_nbasis=dt_sel, dic_nbasis=dic_sel,
    dt_excess_puma  = (dt_puma  - oracle_puma) / mean(oracle_puma) * 100,
    dic_excess_puma = (dic_puma - oracle_puma) / mean(oracle_puma) * 100
  )
}

# Run on worst 5 DT failures
puma_results <- map(1:min(5, nrow(worst_dt)), function(i) {
  row <- worst_dt[i,]
  examine_comp(row$config, row$comp_no, row$dt_nll_nbasis,
               row$dic_nbasis, row$oracle_nbasis)
}) %>% bind_rows()

if (nrow(puma_results) > 0) {
  cat("\nPer-PUMA excess summary (DT-NLL vs DIC) for worst DT failures:\n")
  puma_results %>%
    group_by(config, comp_no, oracle_nbasis, dt_nbasis, dic_nbasis) %>%
    summarise(
      dt_puma_excess_mean  = mean(dt_excess_puma),
      dt_puma_excess_max   = max(dt_excess_puma),
      dic_puma_excess_mean = mean(dic_excess_puma),
      dic_puma_excess_max  = max(dic_excess_puma),
      n_pumas_dt_worse     = sum(dt_excess_puma > dic_excess_puma),
      n_pumas              = n(),
      .groups="drop"
    ) %>%
    print()
}
