#!/usr/bin/env python3

# Generate sample location data from tree sequences
# Path: scripts/generate_samples.py
# Run from the root directory with `python ./scripts/generate_samples.py <tree_name>`
# Use --subsets flag to process all subset trees in ./trees/subsets/{tree_name} (desired behavior)

import tskit
import pandas as pd
import logging
from pathlib import Path
import csv
import argparse
from typing import Set, List, Dict, Optional, Union
import sys

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)


def read_ancient_ids(log_path: Union[str, Path]) -> Set[int]:
    """
    Read pedigree IDs of ancient samples from log file

    Args:
        log_path: Path to log file

    Returns:
        Set of pedigree IDs for ancient samples
    """
    ancient_ids: Set[int] = set()
    log_path = Path(log_path)

    if not log_path.exists():
        logging.warning(f"Log file not found: {log_path}")
        return ancient_ids

    try:
        with open(log_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Convert string of comma-separated IDs into individual IDs
                ids = row['pedigree_IDs'].strip('"').split(',')
                ancient_ids.update(map(int, ids))
        logging.info(f"Read {len(ancient_ids)} ancient IDs from {log_path}")
        return ancient_ids
    except Exception as e:
        logging.error(f"Error reading log file {log_path}: {str(e)}")
        return set()


def create_location_csv(
        tree_path: Union[str, Path],
        output_path: Union[str, Path],
        ancient_ids: Optional[Set[int]] = None,
        subset_mode: bool = False
) -> str:
    """
    Create CSV with sample locations from a single tree sequence file

    Args:
        tree_path: Path to .trees file
        output_path: Full path for output CSV
        ancient_ids: Set of pedigree IDs for ancient samples
        subset_mode: If True, only output node_id, x, y, z columns

    Returns:
        Status message string
    """
    tree_path = Path(tree_path)
    output_path = Path(output_path)

    if not tree_path.exists():
        return f"Error: Tree file not found: {tree_path}"

    try:
        ts = tskit.load(str(tree_path))
        sample_locations: List[Dict] = []
        processed_individuals: Set[int] = set()

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
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, index=False)

        # Generate summary statistics
        total_samples = len(df)
        if subset_mode:
            return f"Successfully created {output_path} with {total_samples} nodes"
        else:
            ancient_samples = len(df[df['is_ancient']])
            unique_individuals = len(df['individual_id'].unique())
            return (f"Successfully created {output_path} with {total_samples} nodes "
                    f"from {unique_individuals} individuals ({ancient_samples} ancient nodes)")

    except Exception as e:
        logging.error(f"Error processing {tree_path}: {str(e)}")
        return f"Error processing {tree_path}: {str(e)}"


def process_subsets(tree_name: str) -> None:
    """
    Process all tree sequences in the subsets directory for a given tree name

    Args:
        tree_name: Name of the tree to process subsets for
    """
    subsets_dir = Path(f"./trees/subsets/{tree_name}")
    output_dir = Path(f"./sample_locations/subsets/{tree_name}")

    if not subsets_dir.exists():
        logging.error(f"Subsets directory not found: {subsets_dir}")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    # Process each .trees file in the subsets directory
    tree_files = list(subsets_dir.glob("*.trees"))
    if not tree_files:
        logging.warning(f"No .trees files found in {subsets_dir}")
        return

    for tree_file in tree_files:
        subset_name = tree_file.stem
        output_path = output_dir / f"{subset_name}_locations.csv"
        result = create_location_csv(tree_file, output_path, subset_mode=True)
        logging.info(result)


def main() -> int:
    """
    Main function to generate sample location data from tree sequences.

    Returns:
        0 for success, 1 for failure
    """
    parser = argparse.ArgumentParser(description="Generate sample location data from tree sequences")
    parser.add_argument("tree_name", help="Name of the tree to process")
    parser.add_argument("--subsets", action="store_true",
                        help="Process subset trees in ./trees/subsets/{tree_name}")
    args = parser.parse_args()

    try:
        if args.subsets:
            process_subsets(args.tree_name)
        else:
            input_path = Path("./trees") / f"{args.tree_name}.trees"
            log_path = Path("./logs") / f"{args.tree_name}.txt"
            output_dir = Path("./sample_locations")
            output_dir.mkdir(parents=True, exist_ok=True)

            if not input_path.exists():
                logging.error(f"Tree file not found: {input_path}")
                return 1

            # Read ancient sample IDs from log file
            ancient_ids = read_ancient_ids(log_path)
            if not ancient_ids:
                logging.warning("No ancient sample IDs found in log file")
                return 1

            # Create output filename and process
            output_path = output_dir / f"{args.tree_name}_locations.csv"
            result = create_location_csv(input_path, output_path, ancient_ids)
            logging.info(result)

            return 0 if "Successfully" in result else 1

    except Exception as e:
        logging.error(f"Unexpected error: {str(e)}")
        return 1


if __name__ == "__main__":
    sys.exit(main())