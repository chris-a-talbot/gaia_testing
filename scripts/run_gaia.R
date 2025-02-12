#!/usr/bin/env Rscript

# Function to ensure package is installed
ensure_package <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    install.packages(package_name, repos = "https://cloud.r-project.org")
  }
}

# Install required packages if missing
required_packages <- c("data.table", "logging", "gaia", "parallel", "foreach", "doParallel")
sapply(required_packages, ensure_package)

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

# Determine number of cores from SLURM allocation
num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(num_cores) || num_cores < 1) {
  # Fallback if not running under SLURM
  num_cores <- max(1, parallel::detectCores() - 2)
}
loginfo(sprintf("Using %d cores for processing", num_cores))

# Create output directories if needed
output_dirs <- c(
  file.path(output_base, "mpr"),
  file.path(output_base, "locations")
)
sapply(output_dirs, dir.create, showWarnings = FALSE, recursive = TRUE)

# Rest of your existing functions remain the same...
[Previous check_sample_size and process_subset functions]

# Modified main execution
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

  # Set up parallel processing with error handling
  cl <- tryCatch({
    makeCluster(num_cores)
  }, error = function(e) {
    logwarn(sprintf("Failed to create cluster with %d cores: %s", num_cores, e$message))
    logwarn("Falling back to sequential processing")
    return(NULL)
  })

  if (!is.null(cl)) {
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
  } else {
    # Sequential processing fallback
    results <- lapply(subset_names, process_subset)
  }

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

# Run main function with error handling
tryCatch({
  success <- main()
  if (!success) {
    quit(status = 1)
  }
}, error = function(e) {
  logerror(sprintf("Fatal error: %s", e$message))
  quit(status = 1)
})