#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(logging)
  library(gaia)
  library(parallel)
  library(foreach)
  library(doParallel)
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript run_gaia.R <tree_name>")
}
tree_name <- args[1]

# Set up logging
basicConfig()
addHandler(writeToConsole)

# Test logging setup
loginfo("Logging system initialized")
loginfo(sprintf("Processing tree: %s", tree_name))

# Directory paths and configuration
trees_dir <- file.path("trees", "subsets", tree_name)
locations_dir <- file.path("sample_locations", "subsets", tree_name)
output_base <- file.path("inferred_locations", "subsets", tree_name)
max_samples <- 10000  # Maximum number of samples allowed

# Determine number of cores (total available - 2, minimum 1)
num_cores <- max(1, parallel::detectCores() - 2)
loginfo(sprintf("Using %d cores for processing", num_cores))

# Create output directories if needed
output_dirs <- c(
  file.path(output_base, "mpr"),
  file.path(output_base, "locations")
)
sapply(output_dirs, dir.create, showWarnings = FALSE, recursive = TRUE)

# Check number of samples in a locations file
check_sample_size <- function(locations_file) {
  tryCatch({
    num_samples <- nrow(fread(locations_file))
    return(list(valid = num_samples <= max_samples, count = num_samples))
  }, error = function(e) {
    return(list(valid = FALSE, error = e$message))
  })
}

# Process a single tree/location pair
process_subset <- function(subset_name) {
  tree_file <- file.path(trees_dir, paste0(subset_name, ".trees"))
  locations_file <- file.path(locations_dir, paste0(subset_name, "_locations.csv"))
  
  tryCatch({
    # Verify files exist
    if (!file.exists(tree_file)) {
      return(list(success = FALSE, subset = subset_name, 
                  error = "Tree file not found"))
    }
    if (!file.exists(locations_file)) {
      return(list(success = FALSE, subset = subset_name, 
                  error = "Locations file not found"))
    }
    
    # Check sample size
    sample_check <- check_sample_size(locations_file)
    if (!sample_check$valid) {
      if (!is.null(sample_check$error)) {
        return(list(success = FALSE, subset = subset_name, 
                    error = paste("Error checking samples:", sample_check$error)))
      } else {
        return(list(success = FALSE, subset = subset_name, 
                    error = sprintf("Too many samples: %d (limit: %d)", 
                                    sample_check$count, max_samples)))
      }
    }
    
    loginfo(sprintf("Processing subset: %s", subset_name))
    
    # Load data
    tree <- gaia::treeseq_load(tree_file)
    sample_locations <- as.matrix(fread(locations_file))
    
    # Run MPR analysis
    mpr_output <- gaia::treeseq_quadratic_mpr(tree, sample_locations, TRUE)
    
    # Save MPR results
    mpr_path <- file.path(output_base, "mpr", paste0(subset_name, "_mpr.rds"))
    saveRDS(mpr_output, mpr_path)
    
    # Minimize and get locations
    locations_output <- gaia::treeseq_quadratic_mpr_minimize(mpr_output)
    
    # Convert to data.table and add node_id column
    locations_dt <- as.data.table(locations_output)
    locations_dt[, node_id := .I - 1]  # 0-based indexing
    
    # Rename columns
    setnames(locations_dt, c("V1", "V2", "V3"), c("x", "y", "z"))
    
    # Read sample locations to get sample node IDs
    sample_dt <- fread(locations_file)
    sample_nodes <- sample_dt$node_id
    
    # Filter out rows corresponding to samples
    filtered_locations <- locations_dt[!node_id %in% sample_nodes]
    
    # Reorder columns to put node_id first
    setcolorder(filtered_locations, c("node_id", "x", "y", "z"))
    
    # Save filtered locations
    loc_path <- file.path(output_base, "locations", 
                          paste0(subset_name, "_inferred_locations.csv"))
    fwrite(filtered_locations, loc_path)
    
    loginfo(sprintf("Successfully processed subset: %s", subset_name))
    return(list(success = TRUE, subset = subset_name))
    
  }, error = function(e) {
    return(list(success = FALSE, subset = subset_name, error = e$message))
  })
}

# Main execution
main <- function() {
  loginfo(sprintf("Starting analysis for tree: %s", tree_name))
  
  # Check if directories exist
  if (!dir.exists(trees_dir)) {
    stop(sprintf("Trees directory not found: %s", trees_dir))
  }
  if (!dir.exists(locations_dir)) {
    stop(sprintf("Locations directory not found: %s", locations_dir))
  }
  
  # Get all tree files
  tree_files <- list.files(path = trees_dir, pattern = "\\.trees$")
  subset_names <- tools::file_path_sans_ext(tree_files)
  
  if (length(subset_names) == 0) {
    stop(sprintf("No .trees files found in: %s", trees_dir))
  }
  
  loginfo(sprintf("Found %d subsets to process", length(subset_names)))
  
  # Set up parallel processing
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Process in parallel
  results <- foreach(
    subset_name = subset_names,
    .packages = c("data.table", "gaia", "logging"),
    .export = c("process_subset", "check_sample_size", "max_samples",
                "trees_dir", "locations_dir", "output_base"),
    .errorhandling = "pass"
  ) %dopar% {
    process_subset(subset_name)
  }
  
  # Stop cluster
  stopCluster(cl)
  
  # Process results
  successful <- sapply(results, function(x) isTRUE(x$success))
  success_count <- sum(successful)
  total_count <- length(results)
  
  # Log results
  loginfo(sprintf("\nProcessing complete: %d/%d successful", success_count, total_count))
  
  if (!all(successful)) {
    logerror("\nFailed subsets:")
    failed_results <- results[!successful]
    for (res in failed_results) {
      if (!is.null(res$error)) {
        logerror(sprintf("%s: %s", res$subset, res$error))
      }
    }
  }
  
  return(success_count > 0)
}

# Run main function
main()