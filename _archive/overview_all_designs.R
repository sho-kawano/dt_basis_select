#!/usr/bin/env Rscript
# Overview: DT-NLL eps=0.7 R=5, DIC, WAIC across all 13 designs (S=50)
library(tidyverse)

EPS <- 0.7; REPS <- 5

# ==============================================================================
# LOAD
# ==============================================================================
eq  <- readRDS("results_multi_config/paper_final/equal_allocation_results.RDS")
pa  <- readRDS("results_multi_config/paper_final/pa_method_comparison.RDS")
paf <- readRDS("results_multi_config/paper_final/pa_method_comparison_final.RDS")
ph  <- readRDS("results_multi_config/paper_final/phase1_results.RDS")

# Combine oracle (paf contributes only prop_1p5pct, not the overlap designs)
oracle <- bind_rows(
  eq$oracle,
  pa$oracle,
  paf$oracle %>% filter(config == "prop_1p5pct"),
  ph$oracle
)

# Combine DT (same dedup rule)
dt <- bind_rows(
  eq$dt_1fold,
  pa$dt,
  paf$dt %>% filter(config == "prop_1p5pct"),
  ph$dt
)

cat(sprintf("Designs in oracle: %s\n", paste(sort(unique(oracle$config)), collapse=", ")))

# ==============================================================================
# SELECTIONS: argmin per method per (config, comp_no)
# ==============================================================================

oracle_sel <- oracle %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n=1, with_ties=FALSE) %>%
  select(config, comp_no, oracle_nbasis=nbasis)

dic_sel <- oracle %>%
  filter(!is.na(DIC)) %>%
  group_by(config, comp_no) %>%
  slice_min(DIC, n=1, with_ties=FALSE) %>%
  select(config, comp_no, selected_nbasis=nbasis) %>%
  mutate(method="DIC")

waic_sel <- oracle %>%
  filter(!is.na(WAIC)) %>%
  group_by(config, comp_no) %>%
  slice_min(WAIC, n=1, with_ties=FALSE) %>%
  select(config, comp_no, selected_nbasis=nbasis) %>%
  mutate(method="WAIC")

dtnll_sel <- dt %>%
  filter(loss_function == "plugin_NLL",
         abs(epsilon - EPS) < 0.001,
         n_reps_used == REPS) %>%
  group_by(config, comp_no) %>%
  slice_min(metric, n=1, with_ties=FALSE) %>%
  select(config, comp_no, selected_nbasis=nbasis) %>%
  mutate(method="DT-NLL")

dtmse_sel <- dt %>%
  filter(loss_function == "MSE",
         abs(epsilon - EPS) < 0.001,
         n_reps_used == REPS) %>%
  group_by(config, comp_no) %>%
  slice_min(metric, n=1, with_ties=FALSE) %>%
  select(config, comp_no, selected_nbasis=nbasis) %>%
  mutate(method="DT-MSE")

all_sel <- bind_rows(dic_sel, waic_sel, dtnll_sel, dtmse_sel) %>%
  left_join(oracle_sel, by=c("config","comp_no")) %>%
  mutate(deviation = selected_nbasis - oracle_nbasis)

# ==============================================================================
# SUMMARY TABLE
# ==============================================================================

design_order <- c(
  "equal_30","equal_40","equal_50","equal_75","equal_100","equal_125","equal_150",
  "prop_0.75pct","prop_1p25pct","prop_1p5pct","prop_1p75pct","prop_2pct","prop_2p25pct"
)
avg_n_map <- c(
  equal_30=30, equal_40=40, equal_50=50, equal_75=75, equal_100=100,
  equal_125=125, equal_150=150,
  prop_0.75pct=47, prop_1p25pct=78, prop_1p5pct=94, prop_1p75pct=110,
  prop_2pct=125, prop_2p25pct=141
)

summary_tbl <- all_sel %>%
  group_by(config, method) %>%
  summarise(
    MAE  = mean(abs(deviation)),
    bias = mean(deviation),
    n    = n(),
    .groups="drop"
  ) %>%
  mutate(
    avg_n  = avg_n_map[config],
    config = factor(config, levels=design_order)
  ) %>%
  arrange(config) %>%
  select(config, avg_n, method, MAE, bias, n) %>%
  pivot_wider(names_from=method, values_from=c(MAE, bias),
              names_glue="{method}_{.value}") %>%
  select(config, avg_n,
         `DIC_MAE`, `DIC_bias`,
         `DT-NLL_MAE`, `DT-NLL_bias`,
         `DT-MSE_MAE`, `DT-MSE_bias`,
         `WAIC_MAE`, `WAIC_bias`)

cat("\n=== DT-NLL (eps=0.7, R=5) vs DIC vs WAIC — all designs, S=50 ===\n")
cat("MAE = mean |selected - oracle|;  bias = mean (selected - oracle)\n\n")
print(summary_tbl %>% mutate(across(where(is.numeric), ~round(., 1))), n=20)

# ==============================================================================
# ORACLE STATS
# ==============================================================================

cat("\n=== Oracle basis stats ===\n")
oracle_sel %>%
  mutate(avg_n = avg_n_map[config],
         config = factor(config, levels=design_order)) %>%
  group_by(config, avg_n) %>%
  summarise(n=n(), mean=mean(oracle_nbasis), sd=sd(oracle_nbasis),
            min=min(oracle_nbasis), max=max(oracle_nbasis), .groups="drop") %>%
  arrange(config) %>%
  print(n=20)

# ==============================================================================
# WINNER TABLE (which method has lower MAE per design?)
# ==============================================================================

cat("\n=== Winner (DIC vs DT-NLL, by MAE) ===\n")
summary_tbl %>%
  mutate(winner = case_when(
           `DIC_MAE` == pmin(`DIC_MAE`, `DT-NLL_MAE`, `DT-MSE_MAE`) ~ "DIC",
           `DT-NLL_MAE` <= `DT-MSE_MAE` ~ "DT-NLL",
           TRUE ~ "DT-MSE")) %>%
  select(config, avg_n, DIC=`DIC_MAE`, `DT-NLL`=`DT-NLL_MAE`, `DT-MSE`=`DT-MSE_MAE`, winner) %>%
  mutate(across(where(is.numeric), ~round(., 1))) %>%
  print(n=20)

# ==============================================================================
# BIAS CONSISTENCY: how stable is each method's bias across designs?
# ==============================================================================

# ==============================================================================
# TRIO COMPARISON: which option makes DT "best overall"?
# ==============================================================================

cat("\n=== Trio-level averages (MAE) ===\n")
trios <- list(
  # Original candidates
  "Opt1: PA 0.75/1.25/2.25"  = c("prop_0.75pct","prop_1p25pct","prop_2p25pct"),
  "Opt2: Equal 30/75/150"    = c("equal_30","equal_75","equal_150"),
  "Opt3: PA 1.25/1.75/2.25"  = c("prop_1p25pct","prop_1p75pct","prop_2p25pct"),
  "Opt4: Equal 30/100/150"   = c("equal_30","equal_100","equal_150"),
  # New options
  "Opt5: PA 0.75/1.5/2.25"   = c("prop_0.75pct","prop_1p5pct","prop_2p25pct"),
  "Opt6: PA 1.5/1.75/2.25"   = c("prop_1p5pct","prop_1p75pct","prop_2p25pct"),
  "Opt7: Equal 50/100/150"   = c("equal_50","equal_100","equal_150"),
  "Opt8: Equal 40/100/150"   = c("equal_40","equal_100","equal_150"),
  # Original paper trio (S=20 baseline)
  "ORIG: PA 0.75/1.25/1.75"  = c("prop_0.75pct","prop_1p25pct","prop_1p75pct")
)

# Build full per-design summary for each trio
trio_results <- list()
for (trio_name in names(trios)) {
  trio_designs <- trios[[trio_name]]

  per_design <- all_sel %>%
    filter(config %in% trio_designs) %>%
    group_by(config, method) %>%
    summarise(MAE=mean(abs(deviation)), bias=mean(deviation), .groups="drop") %>%
    mutate(avg_n = avg_n_map[config]) %>%
    arrange(avg_n, method)

  overall <- all_sel %>%
    filter(config %in% trio_designs) %>%
    group_by(method) %>%
    summarise(MAE=round(mean(abs(deviation)),2), .groups="drop") %>%
    pivot_wider(names_from=method, values_from=MAE)

  best_dt_mae <- min(overall$`DT-NLL`, overall$`DT-MSE`, na.rm=TRUE)
  margin <- overall$DIC - best_dt_mae
  dt_winner <- ifelse(overall$`DT-NLL` <= overall$`DT-MSE`, "DT-NLL", "DT-MSE")

  trio_results[[trio_name]] <- list(
    per_design=per_design, overall=overall,
    margin=margin, dt_winner=dt_winner
  )
}

# Summary ranking table
cat("\n=== TRIO RANKING by DT margin over DIC ===\n")
cat("(positive margin = DT wins overall; larger = more robust)\n\n")
ranking <- data.frame(
  trio     = names(trios),
  DIC      = sapply(trio_results, function(x) x$overall$DIC),
  DT_NLL   = sapply(trio_results, function(x) x$overall$`DT-NLL`),
  DT_MSE   = sapply(trio_results, function(x) x$overall$`DT-MSE`),
  WAIC     = sapply(trio_results, function(x) x$overall$WAIC),
  margin   = sapply(trio_results, function(x) round(x$margin, 2)),
  dt_winner= sapply(trio_results, function(x) x$dt_winner)
) %>% arrange(desc(margin))
print(ranking, row.names=FALSE)

# Per-design detail for each trio
cat("\n\n=== PER-DESIGN DETAIL (MAE | bias) ===\n")
for (trio_name in names(trios)) {
  cat(sprintf("\n--- %s ---\n", trio_name))
  trio_results[[trio_name]]$per_design %>%
    select(config, avg_n, method, MAE, bias) %>%
    mutate(across(c(MAE, bias), ~round(., 1))) %>%
    pivot_wider(names_from=method, values_from=c(MAE, bias),
                names_glue="{method}_{.value}") %>%
    select(config, avg_n,
           DIC_MAE, DIC_bias,
           `DT-NLL_MAE`, `DT-NLL_bias`,
           `DT-MSE_MAE`, `DT-MSE_bias`) %>%
    print()
}

# ==============================================================================
# GRID TRUNCATION ANALYSIS — Opt1 (PA 0.75/1.25/2.25)
# ==============================================================================

OPT1 <- c("prop_0.75pct", "prop_1p25pct", "prop_2p25pct")

trio_stats <- function(oracle, dt, designs, lo, hi) {
  o <- oracle %>% filter(config %in% designs, nbasis >= lo, nbasis <= hi)
  d <- dt     %>% filter(config %in% designs, nbasis >= lo, nbasis <= hi)

  osel <- o %>% filter(!is.na(mse_true)) %>%
    group_by(config, comp_no) %>%
    slice_min(mse_true, n=1, with_ties=FALSE) %>%
    select(config, comp_no, oracle_nbasis=nbasis)

  dic_s  <- o %>% filter(!is.na(DIC))  %>% group_by(config, comp_no) %>%
    slice_min(DIC,  n=1, with_ties=FALSE) %>% select(config, comp_no, selected_nbasis=nbasis) %>% mutate(method="DIC")
  waic_s <- o %>% filter(!is.na(WAIC)) %>% group_by(config, comp_no) %>%
    slice_min(WAIC, n=1, with_ties=FALSE) %>% select(config, comp_no, selected_nbasis=nbasis) %>% mutate(method="WAIC")
  nll_s  <- d %>% filter(loss_function=="plugin_NLL", abs(epsilon-0.7)<0.001, n_reps_used==5) %>%
    group_by(config, comp_no) %>% slice_min(metric, n=1, with_ties=FALSE) %>%
    select(config, comp_no, selected_nbasis=nbasis) %>% mutate(method="DT-NLL")
  mse_s  <- d %>% filter(loss_function=="MSE", abs(epsilon-0.7)<0.001, n_reps_used==5) %>%
    group_by(config, comp_no) %>% slice_min(metric, n=1, with_ties=FALSE) %>%
    select(config, comp_no, selected_nbasis=nbasis) %>% mutate(method="DT-MSE")

  sel <- bind_rows(dic_s, waic_s, nll_s, mse_s) %>%
    left_join(osel, by=c("config","comp_no")) %>%
    mutate(deviation = selected_nbasis - oracle_nbasis)

  # Oracle stats (to see how truncation shifts oracle distribution)
  ostats <- osel %>% group_by(config) %>%
    summarise(n=n(), mean=round(mean(oracle_nbasis),1), sd=round(sd(oracle_nbasis),1),
              max=max(oracle_nbasis), .groups="drop")

  # Per-design MAE + bias
  per_d <- sel %>% group_by(config, method) %>%
    summarise(MAE=round(mean(abs(deviation)),2), bias=round(mean(deviation),1),
              n=n(), .groups="drop") %>%
    mutate(avg_n=avg_n_map[config])

  # Overall
  overall <- sel %>% group_by(method) %>%
    summarise(MAE=round(mean(abs(deviation)),2), .groups="drop") %>%
    pivot_wider(names_from=method, values_from=MAE)

  list(ostats=ostats, per_design=per_d, overall=overall)
}

grids <- list(
  "(3,60) full"  = c(3, 60),
  "(6,60)"       = c(6, 60),
  "(3,48)"       = c(3, 48),
  "(6,48)"       = c(6, 48),
  "(3,45)"       = c(3, 45),
  "(6,45)"       = c(6, 45)
)

cat("\n")
cat("================================================================================\n")
cat("  GRID TRUNCATION — Opt1: PA 0.75 / 1.25 / 2.25\n")
cat("================================================================================\n")

for (gname in names(grids)) {
  lo <- grids[[gname]][1]; hi <- grids[[gname]][2]
  res <- trio_stats(oracle, dt, OPT1, lo, hi)

  cat(sprintf("\n--- Grid [%d, %d] (%d models) ---\n",
              lo, hi, length(seq(lo, hi, by=3))))

  cat("Oracle stats:\n")
  print(res$ostats, row.names=FALSE)

  cat("Overall MAE:\n")
  print(res$overall)

  # Compute DT margin
  best_dt <- min(res$overall$`DT-NLL`, res$overall$`DT-MSE`, na.rm=TRUE)
  dt_w    <- ifelse(res$overall$`DT-NLL` <= res$overall$`DT-MSE`, "DT-NLL", "DT-MSE")
  cat(sprintf("DT margin over DIC: %+.2f  (%s wins)\n", res$overall$DIC - best_dt, dt_w))

  cat("Per-design (MAE | bias):\n")
  res$per_design %>%
    select(config, avg_n, method, MAE, bias) %>%
    pivot_wider(names_from=method, values_from=c(MAE, bias), names_glue="{method}_{.value}") %>%
    select(config, avg_n, DIC_MAE, DIC_bias, `DT-NLL_MAE`, `DT-NLL_bias`,
           `DT-MSE_MAE`, `DT-MSE_bias`) %>%
    print()
}

cat("\n=== Bias by method across all designs ===\n")
cat("(Consistent = same sign and similar magnitude across designs)\n\n")

all_sel %>%
  mutate(avg_n = avg_n_map[config],
         config = factor(config, levels=design_order)) %>%
  group_by(method) %>%
  summarise(
    mean_bias   = mean(deviation),
    sd_bias     = sd(deviation),
    min_bias    = min(deviation),
    max_bias    = max(deviation),
    pct_neg     = mean(deviation < 0) * 100,   # % of comparisons under-selecting
    .groups="drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 2))) %>%
  print()

cat("\n=== Bias per design per method (rows = designs, cols = methods) ===\n")
all_sel %>%
  mutate(avg_n = avg_n_map[config],
         config = factor(config, levels=design_order)) %>%
  group_by(config, avg_n, method) %>%
  summarise(bias = round(mean(deviation), 1), .groups="drop") %>%
  pivot_wider(names_from=method, values_from=bias) %>%
  arrange(config) %>%
  select(config, avg_n, DIC, `DT-NLL`, `DT-MSE`, WAIC) %>%
  print(n=20)
