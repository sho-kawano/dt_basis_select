#!/usr/bin/env Rscript
# ==============================================================================
# Quick check: S=100 results for equal_30, equal_75, equal_125
# Only aggregates eps=0.7, R=5 DT + oracle/DIC/WAIC
# ==============================================================================

library(parallel)
library(doParallel)
library(tidyverse)

DESIGNS <- list(
  list(name = "equal_30",  dir = "_results_equal30"),
  list(name = "equal_75",  dir = "_results_equal75"),
  list(name = "equal_125", dir = "_results_equal125")
)
N_CORES <- 11

source("sim_functions/summary_dt.R")
source("sim_functions/summary_oracle.R")

dt_list     <- list()
oracle_list <- list()

for (design in DESIGNS) {
  cat(sprintf("\nAggregating %s...\n", design$name))

  comp_dirs <- list.dirs(design$dir, recursive = FALSE)
  comp_nos  <- as.integer(sub(".*comparison_0*", "",
                comp_dirs[grepl("comparison_", comp_dirs)]))
  cat(sprintf("  Found %d comparisons\n", length(comp_nos)))

  cl <- makeForkCluster(N_CORES)
  registerDoParallel(cl)

  # DT eps=0.7, both loss functions
  dt_raw_mse <- parLapply(cl, comp_nos, function(x, res_dir) {
    source("sim_functions/summary_dt.R")
    tryCatch(
      summary_dt(x, n_folds=1, loss_function="MSE", results_dir=res_dir,
                 n_reps_to_use=5, eps=0.7) %>%
        mutate(epsilon=0.7, n_reps_used=5, loss="MSE"),
      error = function(e) NULL)
  }, res_dir = design$dir)

  dt_raw_nll <- parLapply(cl, comp_nos, function(x, res_dir) {
    source("sim_functions/summary_dt.R")
    tryCatch(
      summary_dt(x, n_folds=1, loss_function="plugin_NLL", results_dir=res_dir,
                 n_reps_to_use=5, eps=0.7) %>%
        mutate(epsilon=0.7, n_reps_used=5, loss="NLL"),
      error = function(e) NULL)
  }, res_dir = design$dir)

  dt_list[[design$name]] <- bind_rows(dt_raw_mse, dt_raw_nll) %>% mutate(config = design$name)

  # Oracle + DIC/WAIC
  oracle_raw <- parLapply(cl, comp_nos, function(x, res_dir) {
    source("sim_functions/summary_oracle.R")
    tryCatch(summary_oracle(x, res_dir), error = function(e) NULL)
  }, res_dir = design$dir)

  oracle_list[[design$name]] <- bind_rows(oracle_raw) %>% mutate(config = design$name)
  stopCluster(cl)
}

dt_all     <- bind_rows(dt_list)
oracle_all <- bind_rows(oracle_list)
mse_lookup <- oracle_all %>% select(config, comp_no, nbasis, mse_true) %>% filter(!is.na(mse_true))

oracle_sel <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n=1, with_ties=FALSE) %>%
  select(config, comp_no, oracle_nbasis=nbasis, mse_oracle=mse_true)

dt_sel <- dt_all %>%
  group_by(config, comp_no, loss) %>%
  slice_min(metric, n=1, with_ties=FALSE) %>%
  select(config, comp_no, loss, selected_nbasis=nbasis) %>%
  mutate(method=ifelse(loss=="MSE", "DT-MSE e=0.7", "DT-NLL e=0.7")) %>%
  select(-loss)

ic_sel <- oracle_all %>%
  pivot_longer(cols=c(DIC, WAIC), names_to="method", values_to="ic") %>%
  filter(!is.na(ic)) %>%
  group_by(config, comp_no, method) %>%
  slice_min(ic, n=1, with_ties=FALSE) %>%
  select(config, comp_no, method, selected_nbasis=nbasis)

all_sel <- bind_rows(dt_sel, ic_sel) %>%
  left_join(oracle_sel, by=c("config","comp_no")) %>%
  left_join(mse_lookup, by=c("config","comp_no","selected_nbasis"="nbasis")) %>%
  mutate(rel_excess=(mse_true - mse_oracle)/mse_oracle*100,
         deviation=selected_nbasis - oracle_nbasis,
         n_label=as.integer(gsub("equal_","",config)),
         method=factor(method, levels=c("DIC","DT-MSE e=0.7","DT-NLL e=0.7","WAIC")))

cat("\n=== S=100 results: MSE excess% (MAE) ===\n")
all_sel %>%
  group_by(n_label, method) %>%
  summarise(excess=mean(rel_excess), mae=mean(abs(deviation)),
            n=n(), .groups="drop") %>%
  mutate(cell=sprintf("%.1f%% (%.1f)", excess, mae)) %>%
  select(n_label, method, cell) %>%
  pivot_wider(names_from=method, values_from=cell) %>%
  arrange(n_label) %>% print()

cat("\n=== Oracle stats ===\n")
oracle_sel %>%
  mutate(n_label=as.integer(gsub("equal_","",config))) %>%
  group_by(n_label) %>%
  summarise(n=n(), mean_oracle=mean(oracle_nbasis), sd=sd(oracle_nbasis)) %>%
  print()
