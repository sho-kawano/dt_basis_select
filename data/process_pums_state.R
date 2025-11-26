#!/usr/bin/env Rscript
# Process PUMS data for any state
# Usage: Rscript process_pums_state.R <state_code>
# Example: Rscript process_pums_state.R 48  (for Texas)

library(readr)
library(dplyr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage: Rscript process_pums_state.R <state_code>\n")
  cat("Example: Rscript process_pums_state.R 48  (for Texas)\n")
  cat("\nAvailable states:\n")
  cat("  06 - California\n")
  cat("  48 - Texas\n")
  cat("  36 - New York\n")
  cat("  17 - Illinois\n")
  cat("  22 - Louisiana\n")
  cat("  34 - New Jersey\n")
  quit(status = 1)
}

state_code <- args[1]

# State names for output
state_names <- c(
  "06" = "ca",
  "48" = "tx",
  "36" = "ny",
  "17" = "il",
  "22" = "la",
  "34" = "nj"
)

state_abbrev <- state_names[state_code]
if (is.na(state_abbrev)) {
  cat(sprintf("Error: Unknown state code '%s'\n", state_code))
  quit(status = 1)
}

cat(sprintf("\nProcessing PUMS data for state %s (%s)...\n\n", state_code, toupper(state_abbrev)))

# Read household-level variables
h_file <- sprintf("psam_h%s.csv", state_code)
if (!file.exists(h_file)) {
  cat(sprintf("Error: Household file '%s' not found\n", h_file))
  quit(status = 1)
}

households <- read_csv(
  file = h_file,
  col_select = c(
    "SERIALNO",
    "TEN",  # Tenure (own/rent)
    "VEH",  # Vehicles available
    "BDSP" # Bedrooms
  ),
  show_col_types = FALSE
)

# Read person-level variables
p_file <- sprintf("psam_p%s.csv", state_code)
if (!file.exists(p_file)) {
  cat(sprintf("Error: Person file '%s' not found\n", p_file))
  quit(status = 1)
}

persons <- read_csv(
  file = p_file,
  col_select = c(
    "SERIALNO",
    "PWGTP",   # Person weight
    "PUMA",    # Geographic area
    # Demographics
    "AGEP",    # Age
    "SEX",     # Sex
    "RAC1P",   # Race
    "HISP",    # Hispanic origin
    "MAR",     # Marital status
    "CIT",     # Citizenship
    # Economic
    "WAGP",    # Wages/salary
    "ESR",     # Employment status
    "JWMNP",   # Travel time to work
    # Education
    "SCHL",    # Educational attainment
    # Health
    "DIS",     # Disability status
    "PUBCOV",  # Public health insurance
    # Response variables
    "POVPIP"   # Income-to-poverty ratio
  ),
  show_col_types = FALSE
)

# Join person and household data
pop <- persons %>% left_join(households, by = "SERIALNO")

# Exclude group quarters (only households)
pop <- pop %>%
  filter(substr(SERIALNO, start = 5, stop = 6) == "HU")

# Summary
cat(sprintf("Processed %s PUMS data:\n", toupper(state_abbrev)))
cat(sprintf("  Total records: %d\n", nrow(pop)))
cat(sprintf("  Number of PUMAs: %d\n", n_distinct(pop$PUMA)))
cat(sprintf("  Variables included: %d\n\n", ncol(pop)))

# Save as .rds for efficient reading
output_file <- sprintf("%s_pums_population.rds", state_abbrev)
saveRDS(pop, file = output_file)

cat(sprintf("Saved to: data/%s\n", output_file))
