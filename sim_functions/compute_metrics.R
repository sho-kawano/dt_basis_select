library(dplyr)

#' Compute evaluation metrics for model selection methods
#'
#' Computes Mean Absolute Deviation (MAD), directional bias, and oracle penalty
#' for comparing model selection methods against the oracle optimal model.
#'
#' @param oracle_data Data frame with columns: config, comp_no, model, nbasis, mse_true
#'   Contains oracle MSE values for all models
#' @param method_data Data frame with method selection results.
#'   Must contain: config, comp_no, selected_model, selected_nbasis
#'   And a column identifying the method (specified by method_col)
#' @param method_col Name of the column in method_data that identifies the method
#' @param group_cols Character vector of additional grouping columns (e.g., c("epsilon", "n_reps"))
#'
#' @return Data frame with columns:
#'   - config: Survey design configuration
#'   - method: Method name (from method_col)
#'   - Any additional group_cols
#'   - mean_selected_nbasis: Average nbasis selected by method
#'   - mean_oracle_nbasis: Average oracle-optimal nbasis (for reference)
#'   - MAD: Mean absolute deviation from oracle nbasis
#'   - directional_bias: Mean signed deviation (negative = under-selection)
#'   - oracle_penalty: Mean % increase in MSE from suboptimal selection
#'   - n_comparisons: Number of comparisons in the config
#'
#' @examples
#' # For DIC/WAIC (simple case)
#' dic_selections <- oracle_all %>%
#'   group_by(config, comp_no) %>%
#'   slice_min(DIC) %>%
#'   select(config, comp_no, selected_model = model, selected_nbasis = nbasis, method = "DIC")
#'
#' metrics <- compute_metrics(oracle_all, dic_selections, method_col = "method")
#'
#' # For DT with epsilon/n_reps grouping
#' dt_selections <- dt_all %>%
#'   group_by(config, comp_no, epsilon, n_reps_used) %>%
#'   slice_min(metric) %>%
#'   select(config, comp_no, epsilon, n_reps_used,
#'          selected_model = model, selected_nbasis = nbasis) %>%
#'   mutate(method = sprintf("DT ε=%.1f R=%d", epsilon, n_reps_used))
#'
#' metrics <- compute_metrics(oracle_all, dt_selections,
#'                            method_col = "method",
#'                            group_cols = c("epsilon", "n_reps_used"))
compute_metrics <- function(oracle_data, method_data, method_col = "method", group_cols = NULL) {

  # Validate inputs
  required_oracle_cols <- c("config", "comp_no", "model", "nbasis", "mse_true")
  required_method_cols <- c("config", "comp_no", "selected_model", "selected_nbasis", method_col)

  missing_oracle <- setdiff(required_oracle_cols, names(oracle_data))
  missing_method <- setdiff(required_method_cols, names(method_data))

  if (length(missing_oracle) > 0) {
    stop("oracle_data missing columns: ", paste(missing_oracle, collapse = ", "))
  }
  if (length(missing_method) > 0) {
    stop("method_data missing columns: ", paste(missing_method, collapse = ", "))
  }

  # Get oracle selections (best model per comparison)
  oracle_selections <- oracle_data %>%
    filter(!is.na(mse_true)) %>%
    group_by(config, comp_no) %>%
    slice_min(mse_true, n = 1, with_ties = FALSE) %>%
    select(config, comp_no, oracle_model = model, oracle_nbasis = nbasis, oracle_mse = mse_true) %>%
    ungroup()

  # Join method selections with oracle selections
  joined <- method_data %>%
    left_join(oracle_selections, by = c("config", "comp_no"))

  # Get MSE at method's selected model
  joined <- joined %>%
    left_join(
      oracle_data %>% select(config, comp_no, model, mse_true),
      by = c("config", "comp_no", "selected_model" = "model")
    ) %>%
    rename(selected_mse = mse_true)

  # Compute metrics
  # Build grouping columns
  group_by_cols <- c("config", method_col, group_cols)

  metrics <- joined %>%
    mutate(
      deviation = abs(selected_nbasis - oracle_nbasis),
      signed_deviation = selected_nbasis - oracle_nbasis,
      oracle_penalty_pct = (selected_mse - oracle_mse) / oracle_mse * 100
    ) %>%
    group_by(across(all_of(group_by_cols))) %>%
    summarise(
      mean_selected_nbasis = mean(selected_nbasis),
      mean_oracle_nbasis = mean(oracle_nbasis),
      MAD = mean(deviation),
      directional_bias = mean(signed_deviation),
      oracle_penalty = mean(oracle_penalty_pct),
      n_comparisons = n(),
      .groups = "drop"
    )

  return(metrics)
}


#' Get method selections from aggregated DT data
#'
#' Helper function to extract model selections from DT results
#'
#' @param dt_data Data frame with DT results (from aggregate script)
#' @param loss_function Which loss function to use ("MSE" or "plugin_NLL")
#' @param epsilon_values Vector of epsilon values to include (default: all)
#' @param n_reps_values Vector of n_reps values to include (default: all)
#'
#' @return Data frame formatted for compute_metrics()
get_dt_selections <- function(dt_data, loss_function = "MSE",
                               epsilon_values = NULL, n_reps_values = NULL) {

  result <- dt_data %>%
    filter(loss_function == !!loss_function)

  if (!is.null(epsilon_values)) {
    result <- result %>% filter(epsilon %in% epsilon_values)
  }

  if (!is.null(n_reps_values)) {
    result <- result %>% filter(n_reps_used %in% n_reps_values)
  }

  result <- result %>%
    group_by(config, comp_no, epsilon, n_reps_used) %>%
    slice_min(metric, n = 1, with_ties = FALSE) %>%
    select(config, comp_no, epsilon, n_reps_used,
           selected_model = model, selected_nbasis = nbasis) %>%
    mutate(method = sprintf("DT-%s ε=%.1f R=%d", loss_function, epsilon, n_reps_used)) %>%
    ungroup()

  return(result)
}


#' Get method selections from oracle data (for DIC/WAIC)
#'
#' Helper function to extract model selections from DIC/WAIC
#'
#' @param oracle_data Data frame with oracle results (includes DIC/WAIC)
#' @param criterion Which criterion to use ("DIC" or "WAIC")
#'
#' @return Data frame formatted for compute_metrics()
get_ic_selections <- function(oracle_data, criterion = "DIC") {

  if (!criterion %in% c("DIC", "WAIC")) {
    stop("criterion must be 'DIC' or 'WAIC'")
  }

  result <- oracle_data %>%
    filter(!is.na(.data[[criterion]])) %>%
    group_by(config, comp_no) %>%
    slice_min(.data[[criterion]], n = 1, with_ties = FALSE) %>%
    select(config, comp_no, selected_model = model, selected_nbasis = nbasis) %>%
    mutate(method = criterion) %>%
    ungroup()

  return(result)
}
