# Create PUMA Adjacency Matrix for California
# Generates spatial adjacency matrix needed for spatial basis function model
# Output: data/ca_puma_adjacency.RDA containing A (adjacency matrix) and puma_shape

library(tigris)
library(spdep)
library(dplyr)

cat("=== Creating CA PUMA Adjacency Matrix ===\n\n")

# Load population data to get PUMA IDs
cat("1. Loading PUMA IDs from population data...\n")
pop <- readRDS("data/ca_pums_population.rds")
pop_pumas <- sort(unique(pop$PUMA))
n_pumas <- length(pop_pumas)
cat("   Found", n_pumas, "PUMAs in population data\n")
cat("   Range:", head(pop_pumas, 1), "to", tail(pop_pumas, 1), "\n\n")

# Download PUMA shapefiles - try multiple years to find matching vintage
cat("2. Downloading CA PUMA shapefiles...\n")

# Try multiple years to find one that aligns with our data
years_to_try <- c(2022, 2021, 2020, 2019, 2018)
puma_shape <- NULL
alignment_success <- FALSE

for (year in years_to_try) {
  cat("   Trying year", year, "...\n")

  # Download shapefile
  shape_temp <- tryCatch({
    tigris::pumas(state = "CA", cb = FALSE, year = year)
  }, error = function(e) {
    cat("     Error downloading:", e$message, "\n")
    NULL
  })

  if (is.null(shape_temp)) next

  # Check which field contains PUMA codes
  if ("PUMACE10" %in% names(shape_temp)) {
    puma_field <- "PUMACE10"
    cat("     Using PUMACE10 field (2010-based),", nrow(shape_temp), "PUMAs\n")
  } else if ("PUMACE20" %in% names(shape_temp)) {
    puma_field <- "PUMACE20"
    cat("     Using PUMACE20 field (2020-based),", nrow(shape_temp), "PUMAs\n")
  } else {
    cat("     Cannot find PUMA code field, skipping\n")
    next
  }

  # Check alignment
  shape_pumas <- as.character(shape_temp[[puma_field]])
  missing_in_shape <- setdiff(pop_pumas, shape_pumas)

  if (length(missing_in_shape) == 0) {
    cat("     ✓ Perfect alignment! All", n_pumas, "PUMAs found\n")
    puma_shape <- shape_temp
    alignment_success <- TRUE
    break
  } else {
    cat("     ✗ Missing", length(missing_in_shape), "PUMAs from population data\n")
  }
}

if (!alignment_success) {
  stop("Could not find a shapefile year that aligns with population data.\n",
       "Tried years: ", paste(years_to_try, collapse = ", "))
}

cat("\n3. Aligning PUMA IDs...\n")

# Extract and format PUMA codes
shape_pumas <- as.character(puma_shape[[puma_field]])
cat("   Shapefile PUMA range:", min(shape_pumas), "to", max(shape_pumas), "\n")

# Check for extra PUMAs in shapefile
missing_in_pop <- setdiff(shape_pumas, pop_pumas)

if (length(missing_in_pop) > 0) {
  cat("   Note:", length(missing_in_pop), "PUMAs in shapefile not in population (will be removed)\n")
}

# Filter shapefile to only PUMAs in population data
puma_shape <- puma_shape[shape_pumas %in% pop_pumas, ]

# Sort by PUMA code for consistent ordering
puma_shape <- puma_shape[order(puma_shape[[puma_field]]), ]
rownames(puma_shape) <- puma_shape[[puma_field]]

cat("   Final shapefile:", nrow(puma_shape), "PUMAs\n")
cat("   Alignment successful!\n\n")

# Create adjacency matrix
cat("4. Creating adjacency matrix...\n")

# Get neighbor list using rook contiguity (shared borders)
nb <- poly2nb(puma_shape, queen = FALSE, row.names = puma_shape[[puma_field]])

# Convert to binary adjacency matrix
A <- nb2mat(nb, style = 'B', zero.policy = FALSE)

# Ensure dimensions are correct
stopifnot(nrow(A) == nrow(puma_shape))
stopifnot(ncol(A) == nrow(puma_shape))

cat("   Matrix dimensions:", nrow(A), "x", ncol(A), "\n")

# Calculate diagnostics
n_neighbors <- rowSums(A)
n_edges <- sum(A) / 2  # Divide by 2 since symmetric
sparsity <- sum(A) / (nrow(A)^2) * 100

cat("   Total edges:", n_edges, "\n")
cat("   Average neighbors per PUMA:", round(mean(n_neighbors), 2), "\n")
cat("   Neighbor range:", min(n_neighbors), "to", max(n_neighbors), "\n")
cat("   Sparsity:", round(sparsity, 2), "%\n")

# Check for islands (PUMAs with no neighbors)
islands <- which(n_neighbors == 0)
if (length(islands) > 0) {
  cat("\n   WARNING:", length(islands), "island PUMAs with no neighbors:\n")
  island_codes <- rownames(A)[islands]
  cat("   ", paste(island_codes, collapse = ", "), "\n")
  cat("   (Islands may be actual geographic islands)\n")
} else {
  cat("   No island PUMAs detected\n")
}

cat("\n")

# Save results
cat("5. Saving results...\n")
save(A, puma_shape, file = "data/ca_puma_adjacency.RDA")
cat("   Saved to: data/ca_puma_adjacency.RDA\n")
cat("   Objects saved: A (adjacency matrix), puma_shape (sf object)\n\n")

cat("=== Adjacency matrix creation complete! ===\n")
