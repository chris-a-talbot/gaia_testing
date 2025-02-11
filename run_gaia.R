#!/usr/bin/env Rscript

library(data.table)
library(logging)
library(gaia)
library(parallel)
library(foreach)
library(doParallel)

# Set up logging
basicConfig()
addHandler(writeToConsole)

# Test logging setup
loginfo("Logging system initialized")
logerror("Error logging test")

# Directory paths and configuration
trees_dir <- "trees/subsets"
locations_dir <- "sample_locations"
output_dir <- "inferred_locations"
num_cores <- 4  # Set to use 4 cores
max_samples <- 10000  # Maximum number of samples allowed

# Create output directories if needed
dir.create(file.path(output_dir, "mpr"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "locations"), showWarnings = FALSE, recursive = TRUE)

# Check number of samples in a locations file
check_sample_size <- function(locations_file) {
  tryCatch({
    # Just count rows in the locations CSV
    num_samples <- nrow(fread(locations_file))
    return(list(valid = num_samples <= max_samples, count = num_samples))
  }, error = function(e) {
    return(list(valid = FALSE, error = e$message))
  })
}

# Process a single tree/location pair
process_subset <- function(tree_file, locations_file, output_prefix) {
  tryCatch({
    # First verify the files exist
    if (!file.exists(tree_file)) {
      return(list(success = FALSE, file = basename(tree_file), 
                  error = "Tree file not found"))
    }
    if (!file.exists(locations_file)) {
      return(list(success = FALSE, file = basename(tree_file), 
                  error = "Locations file not found"))
    }
    
    # Check sample size
    sample_check <- check_sample_size(locations_file)
    if (!sample_check$valid) {
      if (!is.null(sample_check$error)) {
        return(list(success = FALSE, file = basename(tree_file), 
                    error = paste("Error checking samples:", sample_check$error)))
      } else {
        return(list(success = FALSE, file = basename(tree_file), 
                    error = sprintf("Too many samples: %d (limit: %d)", 
                                    sample_check$count, max_samples)))
      }
    }
    
    # Load data
    tree <- gaia::treeseq_load(tree_file)
    sample_locations <- as.matrix(fread(locations_file))
    
    # Run MPR analysis
    mpr_output <- gaia::treeseq_quadratic_mpr(tree, sample_locations, TRUE)
    
    # Save MPR results
    mpr_path <- file.path(output_dir, "mpr", paste0(output_prefix, "_mpr.rds"))
    saveRDS(mpr_output, mpr_path)
    
    # Minimize and save locations
    locations_output <- gaia::treeseq_quadratic_mpr_minimize(mpr_output)
    loc_path <- file.path(output_dir, "locations", paste0(output_prefix, "_locations.csv"))
    fwrite(as.data.table(locations_output), loc_path)
    
    return(list(success = TRUE, file = basename(tree_file)))
  }, error = function(e) {
    return(list(success = FALSE, file = basename(tree_file), error = e$message))
  })
}

# Check if files exist and match
check_files <- function(tree_files, location_files) {
  if (length(tree_files) == 0) {
    logerror(sprintf("No .trees files found in: %s", trees_dir))
    return(FALSE)
  }
  
  if (length(location_files) == 0) {
    logerror(sprintf("No location files found in: %s", locations_dir))
    return(FALSE)
  }
  
  # Print total files found
  loginfo(sprintf("Found %d tree files and %d location files", 
                  length(tree_files), length(location_files)))
  
  # Check for matching pairs
  tree_bases <- tools::file_path_sans_ext(basename(tree_files))
  loc_bases <- tools::file_path_sans_ext(
    basename(gsub("_locations\\.csv$", "", location_files))
  )
  
  matching <- tree_bases %in% loc_bases
  
  if (!all(matching)) {
    logerror("Missing location files for some trees:")
    logerror(paste(tree_bases[!matching], collapse = ", "))
    return(FALSE)
  }
  
  return(TRUE)
}

# Main execution
main <- function() {
  loginfo("Starting parallel batch analysis...")
  
  # Get all tree files
  loginfo(sprintf("Looking for tree files in: %s", trees_dir))
  tree_files <- list.files(
    path = trees_dir,
    pattern = "\\.trees$",
    full.names = TRUE
  )
  loginfo(sprintf("Found %d tree files", length(tree_files)))
  
  # Get all location files
  loginfo(sprintf("Looking for location files in: %s", locations_dir))
  location_files <- list.files(
    path = locations_dir,
    pattern = "_locations\\.csv$",
    full.names = TRUE
  )
  loginfo(sprintf("Found %d location files", length(location_files)))
  
  # Debug: print first few files if any found
  if (length(tree_files) > 0) {
    loginfo("First few tree files:")
    loginfo(paste(head(basename(tree_files)), collapse="\n"))
  }
  if (length(location_files) > 0) {
    loginfo("First few location files:")
    loginfo(paste(head(basename(location_files)), collapse="\n"))
  }
  
  if (!check_files(tree_files, location_files)) {
    logerror("Aborting due to file matching issues")
    return(FALSE)
  }
  
  # Create processing pairs for files under the sample limit
  processing_pairs <- list()
  skip_count <- 0
  for (tree_file in tree_files) {
    tree_base <- tools::file_path_sans_ext(basename(tree_file))
    locations_file <- file.path(locations_dir, paste0(tree_base, "_locations.csv"))
    
    sample_check <- check_sample_size(locations_file)
    if (!sample_check$valid) {
      if (!is.null(sample_check$error)) {
        logerror(sprintf("Error checking samples for %s: %s", 
                         basename(tree_file), sample_check$error))
      } else {
        loginfo(sprintf("Skipping %s: %d samples (limit: %d)", 
                        basename(tree_file), sample_check$count, max_samples))
        skip_count <- skip_count + 1
      }
    } else {
      processing_pairs[[length(processing_pairs) + 1]] <- list(
        tree_file = tree_file,
        locations_file = locations_file,
        output_prefix = tree_base
      )
    }
  }
  
  loginfo(sprintf("\nTotal files to skip due to sample size: %d", skip_count))
  loginfo(sprintf("Remaining files to process: %d", length(processing_pairs)))
  
  if (length(processing_pairs) == 0) {
    loginfo("No files to process after filtering")
    return(TRUE)
  }
  
  # Set up parallel processing
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  loginfo(sprintf("Processing files using %d cores...", num_cores))
  
  # Process in parallel
  results <- foreach(
    pair = processing_pairs,
    .packages = c("data.table", "gaia", "logging"),  # Added 'logging'
    .export = c("process_subset", "check_sample_size", "max_samples", "output_dir"),  # Export variables/functions
    .errorhandling = "pass"
  ) %dopar% {
    process_subset(pair$tree_file, pair$locations_file, pair$output_prefix)
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
    logerror("\nFailed files:")
    failed_results <- results[!successful]
    for (res in failed_results) {
      if (!is.null(res$error)) {
        logerror(sprintf("%s: %s", res$file, res$error))
      }
    }
  }
  
  return(success_count > 0)
}

# Run main function
if (!interactive()) {
  main()
}