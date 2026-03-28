#!/usr/bin/env Rscript
# ==============================================================================
# Diagnostic: Compare S=20 vs S=50 results with consistent metric computation
# ==============================================================================
# Checks whether the DT win at S=20 is:
#   (a) A seed artifact (S=50 comps 1-20 should match S=20)
#   (b) A metric computation difference
#   (c) Genuine sampling variability
# ==============================================================================

library(tidyverse)

# Consistent metric computation (same as all other scripts)
compute_summary <- function(dt_data, oracle_data, label,
                            eps_val = 0.7, n_reps_val = 5,
                            comp_filter = NULL) {

  if (!is.null(comp_filter))
    oracle_data <- oracle_data %>% filter(comp_no %in% comp_filter)

  oracle_sel <- oracle_data %>%
    filter(!is.na(mse_true)) %>%
    group_by(config, comp_no) %>%
    slice_min(mse_true, n = 1, with_ties = FALSE) %>%
    select(config, comp_no, oracle_nbasis = nbasis, mse_oracle = mse_true)

  mse_lookup <- oracle_data %>%
    select(config, comp_no, nbasis, mse_true) %>%
    filter(!is.na(mse_true))

  # DT selections
  dt_filtered <- dt_data %>%
    filter(abs(epsilon - eps_val) < 0.001, n_reps_used == n_reps_val)

  if (!is.null(comp_filter))
    dt_filtered <- dt_filtered %>% filter(comp_no %in% comp_filter)

  dt_mse_sel <- dt_filtered %>%
    filter(loss_function == "MSE") %>%
    group_by(config, comp_no) %>%
    slice_min(metric, n = 1, with_ties = FALSE) %>%
    select(config, comp_no, selected_nbasis = nbasis) %>%
    mutate(method = "DT-MSE")

  dt_nll_sel <- dt_filtered %>%
    filter(loss_function == "plugin_NLL") %>%
    group_by(config, comp_no) %>%
    slice_min(metric, n = 1, with_ties = FALSE) %>%
    select(config, comp_no, selected_nbasis = nbasis) %>%
    mutate(method = "DT-NLL")

  # DIC/WAIC selections
  ic_sel <- oracle_data %>%
    { if (!is.null(comp_filter)) filter(., comp_no %in% comp_filter) else . } %>%
    pivot_longer(cols = c(DIC, WAIC), names_to = "method", values_to = "ic") %>%
    filter(!is.na(ic)) %>%
    group_by(config, comp_no, method) %>%
    slice_min(ic, n = 1, with_ties = FALSE) %>%
    select(config, comp_no, method, selected_nbasis = nbasis)

  all_sel <- bind_rows(dt_mse_sel, dt_nll_sel, ic_sel) %>%
    left_join(oracle_sel, by = c("config", "comp_no")) %>%
    left_join(mse_lookup, by = c("config", "comp_no", "selected_nbasis" = "nbasis")) %>%
    mutate(
      rel_excess = (mse_true - mse_oracle) / mse_oracle * 100,
      deviation  = selected_nbasis - oracle_nbasis
    )

  result <- all_sel %>%
    group_by(config, method) %>%
    summarise(
      excess = mean(rel_excess, na.rm = TRUE),
      mae    = mean(abs(deviation), na.rm = TRUE),
      n      = n(),
      .groups = "drop"
    ) %>%
    mutate(label = label)

  result
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading data...\n")

# S=50 PA (old seeds 1-50, designs: 0.75%, 1.25%, 1.75%)
pa_s50 <- readRDS("results_multi_config/paper_final/pa_method_comparison.RDS")
dt_pa   <- pa_s50$dt
ora_pa  <- pa_s50$oracle

# S=20 pre-computed metrics (PA designs, comps 1-20)
s20_metrics <- readRDS("results_multi_config/methods_comparison/metrics_summary.RDS")

# Equal allocation S=50
eq_s50   <- readRDS("results_multi_config/paper_final/equal_allocation_results.RDS")
dt_eq    <- eq_s50$dt_1fold
ora_eq   <- eq_s50$oracle

cat(sprintf("PA S=50: %d DT rows, %d oracle rows\n", nrow(dt_pa), nrow(ora_pa)))
cat(sprintf("Equal S=50: %d DT rows, %d oracle rows\n", nrow(dt_eq), nrow(ora_eq)))
cat(sprintf("S=20 pre-computed: %d rows\n\n", nrow(s20_metrics)))

# ==============================================================================
# CHECK: What seeds/comps are in S=20 legacy data?
# ==============================================================================

cat("=== Comp numbers in S=20 legacy (oracle data) ===\n")
s20_oracle <- readRDS("results_multi_config/methods_comparison/oracle_results.RDS")
cat(sprintf("  Configs: %s\n", paste(unique(s20_oracle$config), collapse = ", ")))
cat(sprintf("  Comp nos: %s\n\n", paste(sort(unique(s20_oracle$comp_no)), collapse = ", ")))

# ==============================================================================
# RECOMPUTE: PA S=50 — full 50 comps
# ==============================================================================

cat("=== PA S=50, all 50 comparisons ===\n")
res_pa_s50 <- compute_summary(dt_pa, ora_pa, "PA S=50 (all)")
res_pa_s50 %>%
  arrange(config, method) %>%
  mutate(cell = sprintf("%.1f%% (%.1f)", excess, mae)) %>%
  select(config, method, cell) %>%
  pivot_wider(names_from = method, values_from = cell) %>%
  print()

# ==============================================================================
# RECOMPUTE: PA S=50 — restrict to comps 1-20 (same as S=20)
# ==============================================================================

cat("\n=== PA S=50, comps 1-20 only (matches S=20 seed set) ===\n")
res_pa_s20slice <- compute_summary(dt_pa, ora_pa, "PA S=50 slice (1-20)",
                                    comp_filter = 1:20)
res_pa_s20slice %>%
  arrange(config, method) %>%
  mutate(cell = sprintf("%.1f%% (%.1f)", excess, mae)) %>%
  select(config, method, cell) %>%
  pivot_wider(names_from = method, values_from = cell) %>%
  print()

# ==============================================================================
# S=20 PRE-COMPUTED (from methods_comparison/metrics_summary.RDS)
# ==============================================================================

cat("\n=== S=20 pre-computed metrics (MAD + oracle_penalty) ===\n")
s20_metrics %>%
  filter(grepl("ε=0.7", method) | method %in% c("DIC", "WAIC")) %>%
  filter(grepl("R=5", method) | method %in% c("DIC", "WAIC")) %>%
  arrange(config, method) %>%
  mutate(cell = sprintf("%.1f%% (%.1f)", oracle_penalty, MAD)) %>%
  select(config, method, cell) %>%
  pivot_wider(names_from = method, values_from = cell) %>%
  print()

# ==============================================================================
# EQUAL ALLOCATION: Full S=50 by design
# ==============================================================================

cat("\n=== Equal allocation S=50 by design ===\n")
res_eq <- compute_summary(dt_eq, ora_eq, "Equal S=50")
res_eq %>%
  mutate(n_label = as.integer(gsub("equal_", "", config))) %>%
  arrange(n_label, method) %>%
  mutate(cell = sprintf("%.1f%% (%.1f)", excess, mae)) %>%
  select(n_label, method, cell) %>%
  pivot_wider(names_from = method, values_from = cell) %>%
  arrange(n_label) %>%
  print()

# ==============================================================================
# ORACLE STATS by dataset
# ==============================================================================

cat("\n=== Oracle nbasis stats (PA S=50) ===\n")
ora_pa %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  group_by(config) %>%
  summarise(n = n(), mean = mean(nbasis), sd = sd(nbasis),
            min = min(nbasis), max = max(nbasis), .groups = "drop") %>%
  print()

cat("\n=== Oracle nbasis stats (PA S=50, comps 1-20) ===\n")
ora_pa %>%
  filter(!is.na(mse_true), comp_no <= 20) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  group_by(config) %>%
  summarise(n = n(), mean = mean(nbasis), sd = sd(nbasis),
            min = min(nbasis), max = max(nbasis), .groups = "drop") %>%
  print()

cat("\n=== DT selected nbasis stats (PA S=50, eps=0.7, MSE) ===\n")
dt_pa %>%
  filter(abs(epsilon - 0.7) < 0.001, n_reps_used == 5,
         loss_function == "MSE") %>%
  group_by(config, comp_no) %>%
  slice_min(metric, n = 1, with_ties = FALSE) %>%
  group_by(config) %>%
  summarise(n = n(), mean = mean(nbasis), sd = sd(nbasis), .groups = "drop") %>%
  print()
