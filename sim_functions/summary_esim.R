summary_esim <- function(comp_no, results_dir,
                         validation = c("standard", "data_fission"),
                         n_iters = 100) {
  validation <- match.arg(validation)

  # to help calculate it for one iteration for a given comparison
  helper <- function(sim_no) {
    # folder for comparison
    comp_folder <- file.path(getwd(), results_dir, sprintf("comparison_%03d", comp_no))

    # load the mcmc_summary
    results <- readRDS(file.path(comp_folder, "emp_sim", sprintf("%03d", sim_no), "results.RDS"))

    # load z (original noisy estimate) - new format with PUMA IDs
    z_data <- readRDS(file.path(comp_folder, "z.RDS"))
    puma_ids <- z_data$puma
    z <- z_data$values

    # Choose validation target based on method
    if (validation == "standard") {
      # Standard ESIM: validate on z
      truth <- z
    } else { # data_fission
      # Data fission: validate on z - e, where e = w - z
      # This removes the added noise from validation target
      w <- readRDS(file.path(comp_folder, "emp_sim", sprintf("%03d", sim_no), "w.RDS")) %>% as.numeric()
      e <- w - z
      truth <- z - e # = 2z - w (noise-corrected estimate)
    }

    # create the results table
    res <- results %>%
      select(domain, mean, method, nbasis) %>%
      rename(fips = domain, estim = mean) %>%
      left_join(tibble(fips = puma_ids, truth = truth), by = join_by(fips))

    # calculate mse
    res %>%
      group_by(method, nbasis) %>%
      reframe(mse_sim = mean((estim - truth)^2)) %>%
      mutate(sim_no = sim_no, comp_no = comp_no) %>%
      relocate(comp_no, sim_no) %>%
      arrange(mse_sim)
  }
  #-----------------------------------------------------------------------------
  # run the helper on all n_iters iterations for a comparison, return the results
  results <- lapply(1:n_iters, helper) %>% bind_rows()

  # aggregate by comparison number/method
  results %>%
    group_by(comp_no, method, nbasis) %>%
    reframe(metric = mean(mse_sim)) %>%
    arrange(metric) %>%
    rename(model = method) %>%
    mutate(
      method = paste0("ESIM_", validation),
      metric_type = "Av.MSE"
    ) %>%
    relocate(comp_no, method, model, metric, metric_type, nbasis)
}
