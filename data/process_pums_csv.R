# Process PUMS data for simulation study
# Creates individual-level finite population from CA PUMS data
# Output: cleaned data saved as .rds for efficient reading

library(readr)
library(dplyr)

# Read household-level variables
households <- read_csv(
  file = "psam_h06.csv",
  col_select = c(
    "SERIALNO",
    "TEN", # Tenure (own/rent) - for homeownership_rate
    "VEH", # Vehicles available - for mean_vehicles
    "BDSP" # Bedrooms - for mean_bedrooms
  )
)

# Read person-level variables
persons <- read_csv(
  file = "psam_p06.csv",
  col_select = c(
    "SERIALNO",
    "PWGTP", # Person weight
    "PUMA", # Geographic area
    # Demographics
    "AGEP", # Age - for age-based indicators
    "SEX", # Sex - for sex-based indicators
    "RAC1P", # Race - for race-based indicators
    "HISP", # Hispanic origin - for pct_hispanic
    "MAR", # Marital status - for married indicator
    "CIT", # Citizenship - for pct_citizen
    # Economic
    "WAGP", # Wages/salary
    "ESR", # Employment status
    "JWMNP", # Travel time to work (commute)
    # Education
    "SCHL", # Educational attainment
    # Health
    "DIS", # Disability status
    "PUBCOV", # Public health insurance coverage
    # Response variables
    "POVPIP" # Income-to-poverty ratio
  )
)

# Join person and household data
acs_pop <- persons %>% left_join(households, by = "SERIALNO")

# Exclude group quarters (only households)
# Group quarters = institutional settings (dorms, prisons, nursing homes)
# HU = housing units (regular households)
# We'll use wage as part of the weights so filter out individuals that have NAs for that
acs_pop <- acs_pop %>%
  filter(substr(SERIALNO, start = 5, stop = 6) == "HU")
# %>%
# filter(!is.na(WAGP))
# Summary
cat("Processed CA PUMS data:\n")
cat("  Total records:", nrow(acs_pop), "\n")
cat("  Number of PUMAs:", n_distinct(acs_pop$PUMA), "\n")
cat("  Variables included:", ncol(acs_pop), "\n")

# Save as .rds for efficient reading
saveRDS(acs_pop, file = "ca_pums_population.rds")

cat("\nSaved to: data/ca_pums_population.rds\n")
