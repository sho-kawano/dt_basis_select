#!/usr/bin/env Rscript
# Full table with 6 dimensions + method performance

library(tidyverse)

oracle_all <- readRDS('results_multi_config/oracle_all.RDS')
dt_all <- readRDS('results_multi_config/dt_all.RDS')

configs <- c('equal_40', 'equal_50', 'equal_75', 'prop_0.5pct', 'prop_1pct', 'prop_2pct')

config_to_dir <- c(
  'equal_40' = '_results_equal40_comparison',
  'equal_50' = '_results_equal50_comparison',
  'equal_75' = '_results_equal75_rerun',
  'prop_0.5pct' = '_results_prop0.5pct_comparison',
  'prop_1pct' = '_results_prop1pct_comparison',
  'prop_2pct' = '_results_prop2pct_comparison'
)

# Oracle selections
oracle_sel <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  slice_min(mse_true, n=1, with_ties=FALSE) %>%
  select(config, comp_no, oracle_nbasis=nbasis)

# 1. Signal-to-Noise (CV)
cv_stats <- map_df(configs, function(cfg) {
  results_dir <- config_to_dir[cfg]
  if (!dir.exists(results_dir)) return(tibble(config=cfg, CV=NA_real_))

  comp_dirs <- list.files(results_dir, pattern='^comparison_', full.names=TRUE)
  if (length(comp_dirs) == 0) return(tibble(config=cfg, CV=NA_real_))

  comp_data <- map_df(comp_dirs, function(comp_dir) {
    z_list <- readRDS(file.path(comp_dir, 'z.RDS'))
    d_list <- readRDS(file.path(comp_dir, 'd.RDS'))
    tibble(cv=sqrt(d_list$values)/abs(z_list$values))
  })

  tibble(config=cfg, CV=mean(comp_data$cv, na.rm=TRUE))
})

# 2. Oracle Signal (MSE variation)
mse_var <- oracle_all %>%
  filter(!is.na(mse_true)) %>%
  group_by(config, comp_no) %>%
  summarise(mse_range=max(mse_true)-min(mse_true), mse_min=min(mse_true), .groups='drop') %>%
  group_by(config) %>%
  summarise(MSE_var=mean(mse_range/mse_min)*100, .groups='drop')

# 3. Variance Heterogeneity (d_ratio)
d_ratio <- map_df(configs, function(cfg) {
  results_dir <- config_to_dir[cfg]
  if (!dir.exists(results_dir)) return(tibble(config=cfg, d_ratio=NA_real_))

  comp_dirs <- list.files(results_dir, pattern='^comparison_', full.names=TRUE)
  if (length(comp_dirs) == 0) return(tibble(config=cfg, d_ratio=NA_real_))

  ratios <- map_dbl(comp_dirs, function(comp_dir) {
    d <- readRDS(file.path(comp_dir, 'd.RDS'))$values
    max(d)/min(d)
  })
  tibble(config=cfg, d_ratio=mean(ratios))
})

# 4. Oracle Stability (oracle SD)
oracle_sd <- oracle_sel %>%
  group_by(config) %>%
  summarise(Oracle_SD=sd(oracle_nbasis), .groups='drop')

# 5. Model Complexity Support
complexity <- oracle_sel %>%
  group_by(config) %>%
  summarise(Mean_nbasis=mean(oracle_nbasis), .groups='drop')

# 6. Selection Task Difficulty
difficulty <- mse_var %>%
  left_join(oracle_sd, by='config') %>%
  mutate(Difficulty=Oracle_SD/(MSE_var/100)) %>%
  select(config, Difficulty)

# Method performance
dic <- oracle_all %>%
  filter(!is.na(DIC)) %>%
  group_by(config, comp_no) %>%
  slice_min(DIC, n=1, with_ties=FALSE) %>%
  left_join(oracle_sel, by=c('config','comp_no')) %>%
  mutate(dev=nbasis-oracle_nbasis) %>%
  group_by(config) %>%
  summarise(DIC=mean(abs(dev)), .groups='drop')

waic <- oracle_all %>%
  filter(!is.na(WAIC)) %>%
  group_by(config, comp_no) %>%
  slice_min(WAIC, n=1, with_ties=FALSE) %>%
  left_join(oracle_sel, by=c('config','comp_no')) %>%
  mutate(dev=nbasis-oracle_nbasis) %>%
  group_by(config) %>%
  summarise(WAIC=mean(abs(dev)), .groups='drop')

dt_e5n1 <- dt_all %>%
  filter(method_name=='dt_1fold', abs(epsilon-0.5)<0.01, n_reps_used==1) %>%
  group_by(config, comp_no) %>%
  slice_min(metric, n=1, with_ties=FALSE) %>%
  left_join(oracle_sel, by=c('config','comp_no')) %>%
  mutate(dev=nbasis-oracle_nbasis) %>%
  group_by(config) %>%
  summarise(DT_e5n1=mean(abs(dev)), .groups='drop')

dt_e5n5 <- dt_all %>%
  filter(method_name=='dt_1fold', abs(epsilon-0.5)<0.01, n_reps_used==5) %>%
  group_by(config, comp_no) %>%
  slice_min(metric, n=1, with_ties=FALSE) %>%
  left_join(oracle_sel, by=c('config','comp_no')) %>%
  mutate(dev=nbasis-oracle_nbasis) %>%
  group_by(config) %>%
  summarise(DT_e5n5=mean(abs(dev)), .groups='drop')

dt_e7n5 <- dt_all %>%
  filter(method_name=='dt_1fold', abs(epsilon-0.7)<0.01, n_reps_used==5) %>%
  group_by(config, comp_no) %>%
  slice_min(metric, n=1, with_ties=FALSE) %>%
  left_join(oracle_sel, by=c('config','comp_no')) %>%
  mutate(dev=nbasis-oracle_nbasis) %>%
  group_by(config) %>%
  summarise(DT_e7n5=mean(abs(dev)), .groups='drop')

# Combine all
full_table <- cv_stats %>%
  left_join(mse_var, by='config') %>%
  left_join(d_ratio, by='config') %>%
  left_join(oracle_sd, by='config') %>%
  left_join(complexity, by='config') %>%
  left_join(difficulty, by='config') %>%
  left_join(dic, by='config') %>%
  left_join(waic, by='config') %>%
  left_join(dt_e5n1, by='config') %>%
  left_join(dt_e5n5, by='config') %>%
  left_join(dt_e7n5, by='config') %>%
  arrange(Difficulty)

cat('\n=== SIX DIMENSIONS + METHOD PERFORMANCE ===\n\n')
print(full_table %>% knitr::kable(digits=2))

cat('\n\n=== SUMMARY STATISTICS ===\n\n')
cat('Method Performance (Mean MAD, Range):\n\n')
print(data.frame(
  Method = c('DIC', 'WAIC', 'DT ε=0.5 n=1', 'DT ε=0.5 n=5', 'DT ε=0.7 n=5'),
  Mean_MAD = c(mean(full_table$DIC, na.rm=T), mean(full_table$WAIC, na.rm=T),
               mean(full_table$DT_e5n1, na.rm=T), mean(full_table$DT_e5n5, na.rm=T),
               mean(full_table$DT_e7n5, na.rm=T)),
  Range = c(max(full_table$DIC, na.rm=T)-min(full_table$DIC, na.rm=T),
            max(full_table$WAIC, na.rm=T)-min(full_table$WAIC, na.rm=T),
            max(full_table$DT_e5n1, na.rm=T)-min(full_table$DT_e5n1, na.rm=T),
            max(full_table$DT_e5n5, na.rm=T)-min(full_table$DT_e5n5, na.rm=T),
            max(full_table$DT_e7n5, na.rm=T)-min(full_table$DT_e7n5, na.rm=T))
) %>% arrange(Range) %>% knitr::kable(digits=2))

cat('\n')
