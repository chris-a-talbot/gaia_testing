#!/usr/bin/env python3

import os
import sys
import logging
import numpy as np
import pandas as pd
import tskit
from pathlib import Path
from itertools import combinations

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def compute_max_sample_distance(ts):
    """
    Compute the maximum Euclidean distance between any pair of samples in the tree sequence.

    Args:
        ts (tskit.TreeSequence): The tree sequence to analyze

    Returns:
        float: Maximum distance between any pair of samples
    """
    # Get sample node IDs
    sample_nodes = ts.samples()

    # Get locations for all samples
    sample_locations = []
    for node_id in sample_nodes:
        individual_id = ts.nodes()[node_id].individual
        if individual_id != -1:  # Check if node has an individual
            location = ts.individuals()[individual_id].location
            sample_locations.append(location)

    # Compute maximum distance between all pairs
    max_distance = 0
    for loc1, loc2 in combinations(sample_locations, 2):
        distance = np.sqrt(sum((a - b) ** 2 for a, b in zip(loc1, loc2)))
        max_distance = max(max_distance, distance)

    return max_distance


def process_subset(tree_path, inferred_locations_path, output_path):
    """
    Process a single subset of data, computing distances between known and inferred locations.

    Args:
        tree_path (Path): Path to the tree sequence file
        inferred_locations_path (Path): Path to the inferred locations CSV
        output_path (Path): Path where the output CSV will be saved
    """
    logger.info(f"Processing subset: {tree_path.stem}")

    # Load tree sequence
    ts = tskit.load(str(tree_path))

    # Load inferred locations with node_id as integer
    inferred_df = pd.read_csv(inferred_locations_path, dtype={'node_id': int})

    # Compute maximum sample distance for normalization
    max_sample_dist = compute_max_sample_distance(ts)
    logger.info(f"Maximum sample distance: {max_sample_dist}")

    # Initialize lists to store results
    results = []

    # Process each node in the inferred locations
    for _, row in inferred_df.iterrows():
        node_id = int(row['node_id'])
        inferred_loc = [row['x'], row['y'], row['z']]

        # Get true location from tree sequence
        individual_id = ts.nodes()[node_id].individual
        if individual_id != -1:  # Check if node has an individual
            true_loc = ts.individuals()[individual_id].location

            # Compute Euclidean distance
            raw_dist = np.sqrt(sum((a - b) ** 2 for a, b in zip(true_loc, inferred_loc)))
            normalized_dist = raw_dist / max_sample_dist

            results.append({
                'node_id': node_id,
                'raw_dist': raw_dist,
                'normalized_dist': normalized_dist
            })

    # Create and save results DataFrame
    results_df = pd.DataFrame(results)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_path, index=False)
    logger.info(f"Results saved to: {output_path}")


def main(tree_name):
    """
    Main function to process all subsets for a given tree name.

    Args:
        tree_name (str): Name of the tree to process
    """
    logger.info(f"Starting processing for tree: {tree_name}")

    # Set up paths
    tree_base = Path("./trees/subsets") / tree_name
    inferred_base = Path("./inferred_locations/subsets") / tree_name / "locations"
    output_base = Path("./accuracy/subsets") / tree_name

    # Create output directory if it doesn't exist
    output_base.mkdir(parents=True, exist_ok=True)

    # Process all subset files
    for tree_file in tree_base.glob("*.trees"):
        subset_name = tree_file.stem

        # Construct paths
        inferred_path = inferred_base / f"{subset_name}_inferred_locations.csv"
        output_path = output_base / f"{subset_name}_accuracy.csv"

        if not inferred_path.exists():
            logger.warning(f"Skipping {subset_name}: No inferred locations file found")
            continue

        try:
            process_subset(tree_file, inferred_path, output_path)
        except Exception as e:
            logger.error(f"Error processing {subset_name}: {str(e)}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python analyze_gaia.py <tree_name>")
        sys.exit(1)

    tree_name = sys.argv[1]
    main(tree_name)