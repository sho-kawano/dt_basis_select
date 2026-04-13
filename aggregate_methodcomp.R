#!/usr/bin/env Rscript
# ==============================================================================
# aggregate_methodcomp.R — Aggregate Section 6.3 method comparison results
# ==============================================================================
# Trio: prop_0.75pct / prop_1p25pct / prop_1p75pct
# Seeds 1-50, eps=0.6 R=5
#
# Collects: DT (MSE + NLL), Oracle/DIC/WAIC, ESIM
# Output: results_summary/methodcomp_results.RDS
# ==============================================================================

library(parallel)
library(doParallel)
library(tidyverse)

N_CORES    <- 11
EPS        <- 0.6
REPS       <- 5
ESIM_ITERS <- 100
SEEDS      <- 1:50

DESIGNS <- list(
  list(name = "prop_0.75pct", dir = "_results_prop0.75pct"),
  list(name = "prop_1p25pct", dir = "_results_prop1p25pct"),
  list(name = "prop_1p75pct", dir = "_results_prop1p75pct")
)

dt_list     <- list()
oracle_list <- list()
esim_list   <- list()

for (d in DESIGNS) {
  cat(sprintf("\n%s\n", toupper(d$name)))

  if (!dir.exists(d$dir)) {
    cat(sprintf("  Directory not found: %s — skipping\n", d$dir)); next
  }

  cl <- makeForkCluster(N_CORES)
  registerDoParallel(cl)

  # --- DT: MSE and plugin_NLL at eps=0.6 R=5 ---------------------------------
  cat("  DT...")
  dt_raw <- parLapply(cl, SEEDS, function(x, res_dir) {
    source("sim_functions/summary_dt.R")
    tryCatch(
      bind_rows(
        summary_dt(x, 1, "MSE",        res_dir, REPS, EPS),
        summary_dt(x, 1, "plugin_NLL", res_dir, REPS, EPS)
      ) %>% mutate(epsilon = EPS, n_reps_used = REPS),
      error = function(e) NULL
    )
  }, res_dir = d$dir)
  dt_list[[d$name]] <- bind_rows(dt_raw) %>% mutate(config = d$name)
  cat(sprintf(" %d rows\n", nrow(dt_list[[d$name]])))

  # --- Oracle: mse_true, DIC, WAIC --------------------------------------------
  cat("  Oracle/DIC/WAIC...")
  oracle_raw <- parLapply(cl, SEEDS, function(x, res_dir) {
    source("sim_functions/summary_oracle.R")
    tryCatch(summary_oracle(x, res_dir), error = function(e) NULL)
  }, res_dir = d$dir)
  oracle_list[[d$name]] <- bind_rows(oracle_raw) %>% mutate(config = d$name)
  cat(sprintf(" %d rows\n", nrow(oracle_list[[d$name]])))

  # --- ESIM -------------------------------------------------------------------
  cat("  ESIM...")
  esim_raw <- parLapply(cl, SEEDS, function(x, res_dir, n_iters) {
    source("sim_functions/summary_esim.R")
    tryCatch(
      summary_esim(x, results_dir = res_dir, validation = "standard",
                   n_iters = n_iters),
      error = function(e) NULL
    )
  }, res_dir = d$dir, n_iters = ESIM_ITERS)
  esim_list[[d$name]] <- bind_rows(esim_raw) %>% mutate(config = d$name)
  cat(sprintf(" %d rows\n", nrow(esim_list[[d$name]])))

  stopCluster(cl)
}

# ==============================================================================
# QUICK SUMMARY
# ==============================================================================

oracle_all <- bind_rows(oracle_list)
dt_all     <- bind_rows(dt_list)
esim_all   <- bind_rows(esim_list)

oracle_sel <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n = 1, with_ties = FALSE) %>%
  select(config, comp_no, oracle_nbasis = nbasis)

dic_sel  <- oracle_all %>% filter(!is.na(DIC)) %>%
  group_by(config, comp_no) %>% slice_min(DIC,  n=1, with_ties=FALSE) %>%
  select(config, comp_no, selected_nbasis=nbasis) %>% mutate(method="DIC")

waic_sel <- oracle_all %>% filter(!is.na(WAIC)) %>%
  group_by(config, comp_no) %>% slice_min(WAIC, n=1, with_ties=FALSE) %>%
  select(config, comp_no, selected_nbasis=nbasis) %>% mutate(method="WAIC")

nll_sel  <- dt_all %>%
  filter(loss_function=="plugin_NLL", abs(epsilon-EPS)<0.001, n_reps_used==REPS) %>%
  group_by(config, comp_no) %>% slice_min(metric, n=1, with_ties=FALSE) %>%
  select(config, comp_no, selected_nbasis=nbasis) %>% mutate(method="DT-NLL")

mse_sel  <- dt_all %>%
  filter(loss_function=="MSE", abs(epsilon-EPS)<0.001, n_reps_used==REPS) %>%
  group_by(config, comp_no) %>% slice_min(metric, n=1, with_ties=FALSE) %>%
  select(config, comp_no, selected_nbasis=nbasis) %>% mutate(method="DT-MSE")

esim_sel <- esim_all %>%
  filter(metric_type == "MSE", !is.na(nbasis)) %>%
  group_by(config, comp_no) %>% slice_min(metric, n=1, with_ties=FALSE) %>%
  select(config, comp_no, selected_nbasis=nbasis) %>% mutate(method="ESIM")

all_sel <- bind_rows(dic_sel, waic_sel, nll_sel, mse_sel, esim_sel) %>%
  left_join(oracle_sel, by=c("config","comp_no")) %>%
  mutate(deviation = selected_nbasis - oracle_nbasis)

cat("\n=== MAE and bias (S=50, eps=0.6, R=5) ===\n")
all_sel %>%
  group_by(config, method) %>%
  summarise(MAE=round(mean(abs(deviation)),2), bias=round(mean(deviation),1),
            n=n(), .groups="drop") %>%
  pivot_wider(names_from=method, values_from=c(MAE,bias), names_glue="{method}_{.value}") %>%
  print()

cat("\n=== Overall ===\n")
all_sel %>%
  group_by(method) %>%
  summarise(MAE=round(mean(abs(deviation)),2), bias=round(mean(deviation),1),
            .groups="drop") %>%
  print()

# ==============================================================================
# SAVE
# ==============================================================================

dir.create("results_summary", showWarnings=FALSE, recursive=TRUE)

results <- list(
  dt       = dt_all,
  oracle   = oracle_all,
  esim     = esim_all,
  metadata = list(
    designs        = sapply(DESIGNS, `[[`, "name"),
    dirs           = sapply(DESIGNS, `[[`, "dir"),
    seeds          = SEEDS,
    epsilon        = EPS,
    n_reps         = REPS,
    esim_iters     = ESIM_ITERS,
    loss_functions = c("MSE", "plugin_NLL"),
    created        = Sys.time()
  )
)

saveRDS(results, "results_summary/methodcomp_results.RDS")
cat("\nSaved: results_summary/methodcomp_results.RDS\n")
