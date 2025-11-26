#!/usr/bin/env Rscript
# Minimal CA oracle test - no printing, easy to review

library(parallel)
library(doParallel)
library(LaplacesDemon)
library(tidyverse)
library(Matrix)

# === CONFIGURATION ===
model_config <- list(
  population_file = "data/ca_pums_population.rds",
  adjacency_file = "data/ca_puma_adjacency.RDA",
  response_var = "employed",
  response_type = "binary",
  response_filter = NULL, # NULL = employment-to-population ratio (includes children 0-15)
  X_covariates = NULL,
  nbasis_values = seq(3, 60, by = 3),
  hyp = list(c = 0.001, d = 0.001),
  ndesired = 2000,
  nburn = 1500,
  nthin = 1,
  spatial_type = "fixed",
  equal_allocation = TRUE,
  equal_n = 60,
  n_cores = 10,
  results_dir = "_results_ca_employed_equal60_oracle"
)

n_comp <- 20

# === SETUP ===
source("sim_functions/sampling_and_setup.R")
source("sim_functions/full_data_fit.R")
source("sim_functions/summary_oracle.R")

setup_comp(ncomps = n_comp, results_dir = model_config$results_dir, model_config = model_config)
print("Setup complete!")

# === RUN ORACLE FITS ===
cl <- makeCluster(model_config$n_cores, type = "FORK")
registerDoParallel(cl)

foreach(comp_no = 1:n_comp) %dopar% {
  full_data_fit(comp_no, model_config$results_dir, model_config)
}

stopCluster(cl)

# === ANALYZE ===
oracle_results <- lapply(1:n_comp, function(comp_no) {
  summary_oracle(comp_no, model_config$results_dir)
}) %>% bind_rows()

oracle_models <- oracle_results %>%
  filter(model != "Direct" & model != "D.Est") %>%
  mutate(nbasis = as.numeric(gsub("nbasis_", "", model)))

oracle_optimal <- oracle_models %>%
  group_by(comp_no) %>%
  summarise(
    optimal_nbasis = nbasis[which.min(mse_true)],
    min_mse = min(mse_true),
    max_mse = max(mse_true),
    variation_pct = 100 * (max_mse - min_mse) / min_mse,
    .groups = "drop"
  )

oracle_avg <- oracle_models %>%
  group_by(nbasis) %>%
  summarise(avg_mse = mean(mse_true), .groups = "drop")

# === PLOT ===
plots_dir <- file.path(model_config$results_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

y_min <- min(oracle_models$mse_true) * 0.95
y_max <- max(oracle_models$mse_true) * 1.05

p <- ggplot() +
  geom_line(
    data = oracle_models, aes(x = nbasis, y = mse_true, group = comp_no),
    color = "#FF9999", linewidth = 0.8, alpha = 0.7
  ) +
  geom_line(
    data = oracle_avg, aes(x = nbasis, y = avg_mse),
    color = "#8B0000", linewidth = 2
  ) +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(
    title = sprintf("CA Employed: Equal Allocation n=%d", model_config$equal_n),
    x = "Number of basis functions",
    y = "MSE"
  ) +
  theme_minimal(base_size = 16)
print(p)
ggsave(file.path(plots_dir, "oracle_overlay.png"), p, width = 12, height = 8, dpi = 300)


# === SUMMARY ===
cat(sprintf("\nResults: %s\n", model_config$results_dir))
cat(sprintf(
  "Optimal nbasis: mean=%.1f, SD=%.1f, range=[%d,%d]\n",
  mean(oracle_optimal$optimal_nbasis),
  sd(oracle_optimal$optimal_nbasis),
  min(oracle_optimal$optimal_nbasis),
  max(oracle_optimal$optimal_nbasis)
))
cat(sprintf("MSE variation: %.1f%%\n", mean(oracle_optimal$variation_pct)))
print(oracle_avg, n = Inf)
