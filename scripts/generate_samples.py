#!/usr/bin/env python3

# Generate sample location data from tree sequences
# Path: scripts/generate_samples.py
# Run from the root directory with `python ./scripts/generate_samples.py <tree_name>`
# Use --subsets flag to process all subset trees in ./trees/subsets/{tree_name} (desired behavior)

import tskit
import pandas as pd
from pathlib import Path
import csv
import argparse


def read_ancient_ids(log_path):
    """
    Read pedigree IDs of ancient samples from log file

    Args:
        log_path (str): Path to log file

    Returns:
        set: Set of pedigree IDs for ancient samples
    """
    ancient_ids = set()
    try:
        with open(log_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Convert string of comma-separated IDs into individual IDs
                ids = row['pedigree_IDs'].strip('"').split(',')
                ancient_ids.update(map(int, ids))
        return ancient_ids
    except Exception as e:
        print(f"Error reading log file {log_path}: {str(e)}")
        return set()


def create_location_csv(tree_path, output_path, ancient_ids=None, subset_mode=False):
    """
    Create CSV with sample locations from a single tree sequence file

    Args:
        tree_path (str): Path to .trees file
        output_path (str): Full path for output CSV
        ancient_ids (set, optional): Set of pedigree IDs for ancient samples
        subset_mode (bool): If True, only output node_id, x, y, z columns
    """
    try:
        ts = tskit.load(tree_path)
        sample_locations = []
        processed_individuals = set()  # Track processed individuals to avoid duplicates

        # First process all sample nodes
        samples = ts.samples()
        for sample_id in samples:
            individual_id = ts.node(sample_id).individual
            if individual_id is not None and individual_id not in processed_individuals:
                individual = ts.individual(individual_id)
                if len(individual.location) >= 2:
                    x, y = individual.location[0], individual.location[1]
                    # Add an entry for each node associated with this individual
                    for node_id in individual.nodes:
                        if subset_mode:
                            sample_locations.append({
                                'node_id': node_id,
                                'x': x,
                                'y': y,
                                'z': 0.0
                            })
                        else:
                            sample_locations.append({
                                'node_id': node_id,
                                'individual_id': individual_id,
                                'pedigree_id': individual.metadata.get('pedigree_id', None),
                                'x': x,
                                'y': y,
                                'time': ts.individuals()[individual_id].time,
                                'is_ancient': individual.metadata.get('pedigree_id', None) in ancient_ids
                            })
                    processed_individuals.add(individual_id)

        # Then ensure all ancient samples are included (only in non-subset mode)
        if not subset_mode and ancient_ids:
            for individual_id in range(ts.num_individuals):
                if individual_id not in processed_individuals:
                    individual = ts.individual(individual_id)
                    pedigree_id = individual.metadata.get('pedigree_id', None)
                    if pedigree_id in ancient_ids:
                        if len(individual.location) >= 2:
                            x, y = individual.location[0], individual.location[1]
                            # Add an entry for each node associated with this individual
                            for node_id in individual.nodes:
                                sample_locations.append({
                                    'node_id': node_id,
                                    'individual_id': individual_id,
                                    'pedigree_id': pedigree_id,
                                    'x': x,
                                    'y': y,
                                    'z': 0.0,
                                    'is_ancient': True
                                })
                            processed_individuals.add(individual_id)

        # Create DataFrame and save
        df = pd.DataFrame(sample_locations)
        df.to_csv(output_path, index=False)

        # Print summary statistics
        total_samples = len(df)
        if subset_mode:
            return f"Successfully created {output_path} with {total_samples} nodes"
        else:
            ancient_samples = len(df[df['is_ancient']])
            unique_individuals = len(df['individual_id'].unique())
            return (f"Successfully created {output_path} with {total_samples} nodes "
                    f"from {unique_individuals} individuals ({ancient_samples} ancient nodes)")

    except Exception as e:
        return f"Error processing {tree_path}: {str(e)}"


def process_subsets(tree_name):
    """
    Process all tree sequences in the subsets directory for a given tree name

    Args:
        tree_name (str): Name of the tree to process subsets for
    """
    subsets_dir = Path(f"./trees/subsets/{tree_name}")
    output_dir = Path(f"./sample_locations/subsets/{tree_name}")
    output_dir.mkdir(parents=True, exist_ok=True)

    if not subsets_dir.exists():
        print(f"Error: Subsets directory not found: {subsets_dir}")
        return

    # Process each .trees file in the subsets directory
    for tree_file in subsets_dir.glob("*.trees"):
        subset_name = tree_file.stem
        output_path = output_dir / f"{subset_name}_locations.csv"
        result = create_location_csv(str(tree_file), str(output_path), subset_mode=True)
        print(result)


def main():
    parser = argparse.ArgumentParser(description="Generate sample location data from tree sequences")
    parser.add_argument("tree_name", help="Name of the tree to process")
    parser.add_argument("--subsets", action="store_true",
                      help="Process subset trees in ./trees/subsets/{tree_name}")
    args = parser.parse_args()

    if args.subsets:
        # Process all subsets for the given tree name
        process_subsets(args.tree_name)
    else:
        # Original processing for single tree
        input_path = f"./trees/{args.tree_name}.trees"
        log_path = Path("./logs") / f"{args.tree_name}.txt"
        output_dir = "./sample_locations"
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # Read ancient sample IDs from log file
        ancient_ids = read_ancient_ids(log_path)
        if not ancient_ids:
            print("Warning: No ancient sample IDs found in log file")
            return

        # Create output filename
        output_path = Path(output_dir) / f"{args.tree_name}_locations.csv"

        # Process the file
        result = create_location_csv(input_path, output_path, ancient_ids)
        print(result)


if __name__ == "__main__":
    main()