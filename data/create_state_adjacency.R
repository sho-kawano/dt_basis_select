#!/usr/bin/env Rscript
# Create PUMA Adjacency Matrix for any state
# Usage: Rscript create_state_adjacency.R <state_code>

library(tigris)
library(spdep)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript create_state_adjacency.R <state_code>\n")
  quit(status = 1)
}

state_code <- args[1]
state_names <- c("06"="ca", "48"="tx", "36"="ny", "17"="il", "22"="la", "34"="nj")
state_abbrev <- state_names[state_code]

if (is.na(state_abbrev)) {
  cat(sprintf("Error: Unknown state code '%s'\n", state_code))
  quit(status = 1)
}

cat(sprintf("\n=== Creating %s PUMA Adjacency Matrix ===\n\n", toupper(state_abbrev)))

# Load population data
pop_file <- sprintf("%s_pums_population.rds", state_abbrev)
if (!file.exists(pop_file)) {
  cat(sprintf("Error: Population file '%s' not found. Run process_pums_state.R first.\n", pop_file))
  quit(status = 1)
}

pop <- readRDS(pop_file)
pop_pumas <- sort(unique(pop$PUMA))
n_pumas <- length(pop_pumas)

cat(sprintf("Found %d PUMAs in population data\n\n", n_pumas))

# Download PUMA shapefiles - try multiple years
cat("Downloading PUMA shapefiles...\n")
years_to_try <- c(2022, 2021, 2020, 2019, 2018)
puma_shape <- NULL

for (year in years_to_try) {
  cat(sprintf("  Trying year %d...\n", year))
  
  shape_temp <- tryCatch({
    tigris::pumas(state = state_code, cb = FALSE, year = year)
  }, error = function(e) NULL)
  
  if (is.null(shape_temp)) next
  
  # Identify PUMA field
  if ("PUMACE10" %in% names(shape_temp)) {
    puma_field <- "PUMACE10"
  } else if ("PUMACE20" %in% names(shape_temp)) {
    puma_field <- "PUMACE20"
  } else next
  
  # Check alignment
  shape_pumas <- sort(unique(shape_temp[[puma_field]]))
  if (length(shape_pumas) == n_pumas && all(shape_pumas == pop_pumas)) {
    cat(sprintf("  SUCCESS! Year %d perfectly aligns\n\n", year))
    puma_shape <- shape_temp
    puma_shape$PUMA <- puma_shape[[puma_field]]
    break
  }
}

if (is.null(puma_shape)) {
  cat("ERROR: Could not find matching PUMA shapefile\n")
  quit(status = 1)
}

# Create adjacency matrix
cat("Creating adjacency matrix...\n")
neighbors <- poly2nb(puma_shape, queen = FALSE)
A <- nb2mat(neighbors, style = "B")
rownames(A) <- colnames(A) <- sort(unique(puma_shape$PUMA))

# Summary
n_edges <- sum(A) / 2
cat(sprintf("  %d PUMAs\n", nrow(A)))
cat(sprintf("  %d edges\n", n_edges))
cat(sprintf("  %.2f average neighbors per PUMA\n\n", mean(rowSums(A))))

# Save
output_file <- sprintf("%s_puma_adjacency.RDA", state_abbrev)
save(A, puma_shape, file = output_file)
cat(sprintf("Saved to: %s\n", output_file))
