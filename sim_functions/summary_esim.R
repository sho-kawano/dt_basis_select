summary_esim <- function(comp_no, results_dir,
                         validation = c("standard", "data_fission"),
                         n_iters = 100) {
  validation <- match.arg(validation)

  helper <- function(sim_no) {
    comp_folder <- file.path(getwd(), results_dir, sprintf("comparison_%03d", comp_no))
    results <- readRDS(file.path(comp_folder, "emp_sim", sprintf("%03d", sim_no), "results.RDS"))

    z_data <- readRDS(file.path(comp_folder, "z.RDS"))
    puma_ids <- z_data$puma
    z <- z_data$values

    truth <- if (validation == "standard") {
      z
    } else {
      w <- readRDS(file.path(comp_folder, "emp_sim", sprintf("%03d", sim_no), "w.RDS")) %>% as.numeric()
      z - (w - z)  # = 2z - w
    }

    truth_df <- tibble(fips = puma_ids, truth = truth)

    mse_res <- results %>%
      select(domain, mean, method, nbasis) %>%
      rename(fips = domain, estim = mean) %>%
      left_join(truth_df, by = "fips") %>%
      group_by(method, nbasis) %>%
      reframe(metric = mean((estim - truth)^2)) %>%
      mutate(metric_type = "MSE")

    # IS: alpha=0.10 (90% CI); cr_int.lower/upper from mcmc_summary(level=0.9)
    # Only computed when credible interval columns exist (Paper 2)
    if (all(c("cr_int.lower", "cr_int.upper") %in% names(results))) {
      is_res <- results %>%
        filter(method != "Direct") %>%
        select(domain, cr_int.lower, cr_int.upper, method, nbasis) %>%
        rename(fips = domain) %>%
        left_join(truth_df, by = "fips") %>%
        group_by(method, nbasis) %>%
        reframe(metric = mean(
          (cr_int.upper - cr_int.lower) +
          (2 / 0.10) * pmax(cr_int.lower - truth, 0) +
          (2 / 0.10) * pmax(truth - cr_int.upper, 0)
        )) %>%
        mutate(metric_type = "IS")
      mse_res <- bind_rows(mse_res, is_res)
    }

    mse_res %>%
      mutate(sim_no = sim_no, comp_no = comp_no) %>%
      relocate(comp_no, sim_no)
  }

  results <- lapply(1:n_iters, helper) %>% bind_rows()

  results %>%
    group_by(comp_no, method, nbasis, metric_type) %>%
    reframe(metric = mean(metric)) %>%
    rename(model = method) %>%
    mutate(method = paste0("ESIM_", validation)) %>%
    relocate(comp_no, method, model, metric, metric_type, nbasis)
}
