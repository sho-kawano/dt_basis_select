library(dplyr)
library(ggplot2)

# Load results
results <- readRDS("_results_wagp/all_summaries.RDS")
results <- results %>%
  mutate(nbasis = as.numeric(gsub("nbasis_", "", model))) %>%
  filter(!is.na(nbasis))

# Prepare data for each method
method_data <- list()

# 1. Oracle (OD-Oracle MSE)
oracle <- results %>%
  filter(source == "oracle") %>%
  select(comp_no, nbasis, value = mse_true) %>%
  mutate(method = "OD-Oracle MSE")

# 2. WAIC
waic <- results %>%
  filter(source == "oracle") %>%
  select(comp_no, nbasis, value = WAIC) %>%
  mutate(method = "WAIC")

# 3. DIC
dic <- results %>%
  filter(source == "oracle") %>%
  select(comp_no, nbasis, value = DIC) %>%
  mutate(method = "DIC")

# 4. DT 5-fold
dt5 <- results %>%
  filter(source == "dt_5fold") %>%
  select(comp_no, nbasis, value = metric) %>%
  mutate(method = "DT 5-fold")

# 5. DT 1-fold (eps=0.3)
dt1_03 <- results %>%
  filter(source == "dt_1fold", eps == 0.3) %>%
  select(comp_no, nbasis, value = metric) %>%
  mutate(method = "DT 1-fold (eps=0.3)")

# 6. DT 1-fold (eps=0.5)
dt1_05 <- results %>%
  filter(source == "dt_1fold", eps == 0.5) %>%
  select(comp_no, nbasis, value = metric) %>%
  mutate(method = "DT 1-fold (eps=0.5)")

# 7. ESIM Standard
esim_std <- results %>%
  filter(source == "esim_standard") %>%
  select(comp_no, nbasis, value = metric) %>%
  mutate(method = "ESIM Standard")

# 8. ESIM Data Fission
esim_fis <- results %>%
  filter(source == "esim_fission") %>%
  select(comp_no, nbasis, value = metric) %>%
  mutate(method = "ESIM Data Fission")

# Combine all
all_methods <- bind_rows(oracle, waic, dic, dt5, dt1_03, dt1_05, esim_std, esim_fis)

# Normalize each comparison within each method to [0,1]
# For WAIC and DIC, lower is better, so we need to handle direction
all_methods_norm <- all_methods %>%
  group_by(method, comp_no) %>%
  mutate(
    value_norm = case_when(
      method %in% c("WAIC", "DIC") ~ (max(value) - value) / (max(value) - min(value)),
      TRUE ~ (value - min(value)) / (max(value) - min(value))
    )
  ) %>%
  ungroup()

# Create faceted overlay plot
p <- ggplot(all_methods_norm, aes(x = nbasis, y = value_norm)) +
  geom_line(aes(group = comp_no), alpha = 0.25, color = "gray60") +
  stat_summary(aes(group = 1), fun = mean, geom = "line", 
               color = "blue", linewidth = 1.2) +
  facet_wrap(~method, ncol = 3) +
  labs(
    title = "All Methods: Individual vs Average Curves (WAGP)",
    subtitle = "Normalized to [0,1] per comparison (0=best, 1=worst). Gray=individual, Blue=average",
    x = "Number of basis functions",
    y = "Normalized metric (0=best)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.text = element_text(face = "bold", size = 9),
    axis.text = element_text(size = 7)
  )

ggsave("/tmp/wagp_all_methods_overlay.png", p, width = 12, height = 8, dpi = 150)

cat("✓ Plot saved to: /tmp/wagp_all_methods_overlay.png\n")
