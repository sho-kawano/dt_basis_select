# Process PUMS data for simulation study
# Creates individual-level finite population from CA PUMS data
# Output: cleaned data saved as .rds for efficient reading

library(readr)
library(dplyr)

# Read household-level variables
households <- read_csv(
  file = file.path("data", "psam_h06.csv"),
  col_select = c(
    "SERIALNO",
    "TEN",         # Tenure (own/rent) - for homeownership_rate
    "VEH",         # Vehicles available - for mean_vehicles
    "BDSP"         # Bedrooms - for mean_bedrooms
  )
)

# Read person-level variables
persons <- read_csv(
  file = file.path("data", "psam_p06.csv"),
  col_select = c(
    "SERIALNO",
    "PWGTP",       # Person weight
    "PUMA",        # Geographic area
    # Demographics
    "AGEP",        # Age - for mean_age, sd_age
    "SEX",         # Sex - for pct_male
    "RAC1P",       # Race - for pct_asian, pct_black
    "HISP",        # Hispanic origin - for pct_hispanic
    "MAR",         # Marital status - for pct_married
    "CIT",         # Citizenship - for pct_citizen
    # Economic
    "WAGP",        # Wages/salary - for median_wage
    "ESR",         # Employment status - for employment_rate
    # Education
    "SCHL",        # Educational attainment - for pct_bachelor
    # Health
    "DIS",         # Disability status - for disability_rate
    # Response variables
    "POVPIP",      # Income-to-poverty ratio - for poverty_rate response
    "HICOV"        # Health insurance - for hicov_rate response
  )
)

# Join person and household data
acs_pop <- persons %>% left_join(households, by = "SERIALNO")

# Exclude group quarters (only households)
# Group quarters = institutional settings (dorms, prisons, nursing homes)
# HU = housing units (regular households)
acs_pop <- acs_pop %>%
  filter(substr(SERIALNO, start = 5, stop = 6) == "HU")

# Summary
cat("Processed CA PUMS data:\n")
cat("  Total records:", nrow(acs_pop), "\n")
cat("  Number of PUMAs:", n_distinct(acs_pop$PUMA), "\n")
cat("  Variables included:", ncol(acs_pop), "\n")

# Save as .rds for efficient reading
saveRDS(acs_pop, file = file.path("data", "ca_pums_population.rds"))

cat("\nSaved to: data/ca_pums_population.rds\n")

