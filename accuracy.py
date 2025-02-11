import tskit
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import json
import re


class AccuracyAnalyzer:
    def __init__(self, base_dir="."):
        self.base_dir = Path(base_dir)
        self.trees_dir = self.base_dir / "trees" / "subsets"
        self.sample_locations_dir = self.base_dir / "sample_locations"
        self.inferred_locations_dir = self.base_dir / "inferred_locations" / "locations"
        self.accuracy_dir = self.base_dir / "accuracy"
        self.accuracy_dir.mkdir(exist_ok=True)

        # Time bins edges for analysis (log scale)
        self.time_bins = np.logspace(0, 4, 21)  # 20 bins from 1 to 10000

    def get_sampling_scheme(self, filename):
        """Extract sampling scheme from filename."""
        match = re.search(r'_(extant|temporal|combined|geographic).*?(?=\.trees)', filename)
        return match.group(1) if match else 'unknown'

    def get_sample_max_distance(self, sample_coords):
        """Calculate maximum pairwise distance between samples."""
        if len(sample_coords) < 2:
            return 1.0
        return np.max(pdist(sample_coords))

    def process_single_tree(self, tree_file):
        """Process a single tree sequence and calculate accuracy metrics."""
        try:
            # Check for corresponding files
            tree_stem = tree_file.stem
            sample_loc_file = self.sample_locations_dir / f"{tree_stem}_locations.csv"
            inferred_loc_file = self.inferred_locations_dir / f"{tree_stem}_locations.csv"

            if not sample_loc_file.exists() or not inferred_loc_file.exists():
                print(f"Missing files for {tree_stem}")
                return None

            # Load data
            tree = tskit.load(str(tree_file))
            sample_locations = pd.read_csv(sample_loc_file)
            inferred_locations = pd.read_csv(inferred_loc_file)
            inferred_locations['node_id'] = inferred_locations.index

            # Get true locations from tree sequence
            node_data = []
            for ind in tree.individuals():
                node_id = ind.nodes[0]
                node = tree.node(node_id)
                location = ind.location[:2]
                node_data.append({
                    'node_id': node_id,
                    'x': location[0],
                    'y': location[1],
                    'time': node.time
                })
            tree_df = pd.DataFrame(node_data)

            # Calculate normalization factor
            sample_coords = sample_locations[['x', 'y']].values
            max_distance = self.get_sample_max_distance(sample_coords)

            # Calculate errors for non-sampled nodes
            sampled_nodes = set(sample_locations['node_id'])
            errors = []

            for _, row in inferred_locations.iterrows():
                node_id = row['node_id']
                if node_id not in sampled_nodes and node_id in tree_df['node_id'].values:
                    true_loc = tree_df[tree_df['node_id'] == node_id].iloc[0]
                    inferred_loc = np.array([row.iloc[0], row.iloc[1]])
                    true_loc_array = np.array([true_loc['x'], true_loc['y']])

                    error = np.linalg.norm(true_loc_array - inferred_loc) / max_distance
                    time = float(true_loc['time'])
                    errors.append({'time': time, 'error': error})

            if not errors:
                print(f"No valid errors calculated for {tree_stem}")
                return None

            errors_df = pd.DataFrame(errors)

            # Calculate binned statistics
            bins = self.time_bins
            digitized = np.digitize(errors_df['time'], bins)
            binned_stats = []

            for i in range(1, len(bins)):
                mask = digitized == i
                if np.any(mask):
                    bin_errors = errors_df.loc[mask, 'error']
                    binned_stats.append({
                        'time_min': float(bins[i - 1]),
                        'time_max': float(bins[i]),
                        'error_mean': float(bin_errors.mean()),
                        'error_median': float(bin_errors.median()),
                        'error_std': float(bin_errors.std()) if len(bin_errors) > 1 else 0.0,
                        'count': int(np.sum(mask))
                    })

            # Calculate overall statistics
            results = {
                'tree_file': tree_file.name,
                'scheme': self.get_sampling_scheme(tree_file.name),
                'num_samples': len(sample_locations),
                'overall_stats': {
                    'mean_error': float(errors_df['error'].mean()),
                    'median_error': float(errors_df['error'].median()),
                    'num_nodes': len(errors_df),
                    'time_range': [float(errors_df['time'].min()), float(errors_df['time'].max())],
                    'max_sample_distance': float(max_distance)
                },
                'binned_stats': binned_stats
            }

            # Save results
            output_file = self.accuracy_dir / "individual" / f"{tree_stem}_accuracy.json"
            output_file.parent.mkdir(exist_ok=True)
            with open(output_file, 'w') as f:
                json.dump(results, f, indent=2)

            return results

        except Exception as e:
            print(f"Error processing {tree_file}: {str(e)}")
            return None

    def analyze_all_trees(self):
        """Process all tree sequences."""
        (self.accuracy_dir / "individual").mkdir(exist_ok=True)
        total_trees = len(list(self.trees_dir.glob("*.trees")))

        successful = 0
        failed = 0

        for idx, tree_file in enumerate(self.trees_dir.glob("*.trees"), 1):
            print(f"Processing tree {idx}/{total_trees}: {tree_file.name}")
            result = self.process_single_tree(tree_file)
            if result is not None:
                successful += 1
            else:
                failed += 1

        print(f"\nProcessing complete:")
        print(f"Successful: {successful}")
        print(f"Failed: {failed}")

        if successful > 0:
            self.generate_comparisons()

    def generate_comparisons(self):
        """Generate comparisons from saved results."""
        results_dir = self.accuracy_dir / "individual"
        results = []

        for result_file in results_dir.glob("*_accuracy.json"):
            with open(result_file) as f:
                results.append(json.load(f))

        if not results:
            print("No results found to generate comparisons")
            return

        # Group by sampling scheme
        scheme_results = defaultdict(list)
        for r in results:
            scheme = r['scheme']
            scheme_results[scheme].append(r)

        self.plot_comparisons(scheme_results)

    def plot_comparisons(self, scheme_results):
        """Generate comparison plots."""
        # Accuracy over time plot
        plt.figure(figsize=(12, 8))

        for scheme, results in scheme_results.items():
            times = []
            errors = []

            for r in results:
                for stat in r['binned_stats']:
                    if stat['count'] > 0:
                        mid_time = (stat['time_min'] + stat['time_max']) / 2
                        times.append(mid_time)
                        errors.append(stat['error_mean'])

            if times:
                plt.plot(times, errors, label=scheme, alpha=0.7)

        plt.xscale('log')
        plt.xlabel('Time (generations)')
        plt.ylabel('Mean Normalized Error')
        plt.title('Accuracy Over Time by Sampling Scheme')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(self.accuracy_dir / 'accuracy_over_time.png')
        plt.close()

        # Save scheme comparisons
        comparisons = {}
        for scheme, results in scheme_results.items():
            scheme_stats = {
                'mean_error': np.mean([r['overall_stats']['mean_error'] for r in results]),
                'median_error': np.median([r['overall_stats']['mean_error'] for r in results]),
                'num_trees': len(results)
            }
            comparisons[scheme] = scheme_stats

        with open(self.accuracy_dir / 'scheme_comparisons.json', 'w') as f:
            json.dump(comparisons, f, indent=2)


def main():
    analyzer = AccuracyAnalyzer()
    analyzer.analyze_all_trees()


if __name__ == "__main__":
    main()