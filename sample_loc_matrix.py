# Sample location matrix maker

import os
import tskit
import pandas as pd
from pathlib import Path


def create_location_csv(tree_path, output_path):
    """
    Create CSV with sample locations from a single tree sequence file

    Args:
        tree_path (str): Path to .trees file
        output_path (str): Full path for output CSV
    """
    try:
        ts = tskit.load(tree_path)
        sample_locations = []

        # Get all sample nodes (should be 500 nodes for 250 diploid individuals)
        samples = ts.samples()

        for sample_id in samples:
            individual_id = ts.node(sample_id).individual
            if individual_id is not None:
                individual = ts.individual(individual_id)
                if len(individual.location) >= 2:
                    x, y = individual.location[0], individual.location[1]
                    sample_locations.append({
                        'node_id': sample_id,
                        'x': x,
                        'y': y,
                        'time': ts.node(sample_id).time
                    })

        # Create DataFrame and save
        df = pd.DataFrame(sample_locations)
        df.to_csv(output_path, index=False)
        return f"Successfully created {output_path} with {len(df)} samples"

    except Exception as e:
        return f"Error processing {os.path.basename(tree_path)}: {str(e)}"


def main():
    # Input and output paths
    input_dir = "./trees/subsets"
    output_dir = "./sample_locations"

    # Create output directory if needed
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Get all .trees files in the input directory
    tree_files = [f for f in os.listdir(input_dir) if f.endswith('.trees')]

    if not tree_files:
        print(f"Error: No .trees files found in {input_dir}")
        return

    # Process each tree file
    for tree_file in tree_files:
        input_path = os.path.join(input_dir, tree_file)

        # Create output filename
        tree_name = os.path.splitext(tree_file)[0]
        output_path = os.path.join(output_dir, f"{tree_name}_locations.csv")

        # Process the file
        result = create_location_csv(input_path, output_path)
        print(result)


if __name__ == "__main__":
    main()