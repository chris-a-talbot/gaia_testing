import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re
import tskit
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("plot_errors.log"),
        logging.StreamHandler()
    ]
)

# Set plotting style
sns.set(style="whitegrid", palette="colorblind")
plt.rcParams["figure.figsize"] = (12, 8)


def parse_scheme_details(filename):
    """Improved parser for scheme details with better validation"""
    details = {
        "scheme_type": "other",
        "num_extant": None,
        "num_ancient": None,
        "spacing": None
    }

    try:
        # Combined samples
        combined_match = re.search(r"combined_(\d+)_(\d+)", filename)
        if combined_match:
            details.update({
                "scheme_type": "combined",
                "num_extant": int(combined_match.group(1)),
                "num_ancient": int(combined_match.group(2))
            })
            return details

        # Extant samples
        extant_match = re.search(r"extant_(\d+)", filename)
        if extant_match:
            details.update({
                "scheme_type": "extant",
                "num_extant": int(extant_match.group(1))
            })
            return details

    except Exception as e:
        logging.error(f"Error parsing {filename}: {str(e)}")

    return details


def load_comparison_data(accuracy_dir, trees_dir):
    """Load data with extensive validation"""
    results = []
    results_dir = Path(accuracy_dir) / "individual"

    if not results_dir.exists():
        logging.error(f"Missing results directory: {results_dir}")
        return pd.DataFrame()

    json_files = list(results_dir.glob("*.json"))
    logging.info(f"Found {len(json_files)} JSON result files")

    for result_file in json_files:
        try:
            with open(result_file) as f:
                data = json.load(f)

            details = parse_scheme_details(data["tree_file"])
            if details["scheme_type"] not in ["extant", "combined"]:
                continue

            tree_path = Path(trees_dir) / data["tree_file"]
            if not tree_path.exists():
                logging.warning(f"Missing tree file: {tree_path}")
                continue

            try:
                ts = tskit.load(str(tree_path))
                max_time = max(n.time for n in ts.nodes())
            except Exception as e:
                logging.error(f"Error loading {tree_path}: {str(e)}")
                continue

            for bin_data in data["binned_stats"]:
                if bin_data["count"] < 5:  # Skip poorly populated bins
                    continue

                entry = {
                    "tree_file": data["tree_file"],
                    "num_samples": data["num_samples"],
                    "mean_error": data["overall_stats"]["mean_error"],
                    "time_mid": (bin_data["time_min"] + bin_data["time_max"]) / 2,
                    "bin_mean_error": bin_data["error_mean"],
                    "normalized_time": (bin_data["time_min"] + bin_data["time_max"]) / (2 * max_time),
                    "max_time": max_time,
                    **details
                }
                results.append(entry)

        except Exception as e:
            logging.error(f"Error processing {result_file.name}: {str(e)}")

    logging.info(f"Loaded {len(results)} valid data points")
    return pd.DataFrame(results)


def create_comparison_plots(df):
    """Generate plots with validation"""
    if df.empty:
        logging.error("No data available for plotting")
        return

    # Filter and group data
    extant_df = df[df["scheme_type"] == "extant"]
    combined_df = df[df["scheme_type"] == "combined"]

    logging.info(f"Found {len(extant_df)} extant entries and {len(combined_df)} combined entries")

    # Get unique extant sample counts
    common_counts = np.intersect1d(extant_df.num_extant.unique(),
                                   combined_df.num_extant.unique())
    logging.info(f"Common extant counts: {common_counts}")

    if not common_counts.size:
        logging.error("No matching extant counts between extant and combined schemes")
        return

    # Create output directory
    output_dir = Path("accuracy/comparison_plots")
    output_dir.mkdir(exist_ok=True, parents=True)

    # Create one plot per extant sample count
    for extant_count in common_counts:
        # Get relevant data subsets
        extant_subset = extant_df[(extant_df.num_extant == extant_count)]
        combined_subset = combined_df[combined_df.num_extant == extant_count]

        logging.info(f"\nProcessing extant count {extant_count}:")
        logging.info(f"Extant entries: {len(extant_subset)}")
        logging.info(f"Combined entries: {len(combined_subset)}")

        if combined_subset.empty or extant_subset.empty:
            logging.warning(f"Skipping count {extant_count} due to missing data")
            continue

        # Get min/max ancient counts for this extant count
        ancient_counts = combined_subset.num_ancient.unique()
        logging.info(f"Found ancient counts: {ancient_counts}")

        if len(ancient_counts) < 1:
            logging.warning(f"No ancient counts for extant count {extant_count}")
            continue

        min_ancient = min(ancient_counts)
        max_ancient = max(ancient_counts)

        # Create comparison grid
        fig, ax = plt.subplots()
        plotted_lines = 0

        # Plot error differences for min/max ancient samples
        for ancient_count, color in zip([min_ancient, max_ancient], ['#1f77b4', '#ff7f0e']):
            combined_data = combined_subset[combined_subset.num_ancient == ancient_count]

            if combined_data.empty:
                logging.warning(f"No data for ancient count {ancient_count}")
                continue

            merged = pd.merge(combined_data, extant_subset,
                              on=["time_mid", "normalized_time"],
                              suffixes=('_combined', '_extant'),
                              how="inner")

            if merged.empty:
                logging.warning(f"No overlapping time bins for ancient {ancient_count}")
                continue

            merged["error_diff"] = merged["bin_mean_error_extant"] - merged["bin_mean_error_combined"]

            if merged["error_diff"].isna().all():
                logging.warning(f"All NaN differences for ancient {ancient_count}")
                continue

            sns.lineplot(
                x="normalized_time",
                y="error_diff",
                data=merged,
                color=color,
                label=f"{ancient_count} Ancient Samples",
                ax=ax,
                estimator='median',
                errorbar='sd'
            )
            plotted_lines += 1

        if plotted_lines == 0:
            logging.warning(f"No plottable data for extant count {extant_count}")
            plt.close()
            continue

        # Format plot
        ax.set(
            xscale="log",
            xlabel="Normalized Time (generations / max tree time)",
            ylabel="Accuracy Improvement (Extant Error - Combined Error)",
            title=f"Extant Samples: {extant_count}\nAccuracy Improvement from Adding Historical Samples"
        )
        ax.axhline(0, color='gray', linestyle='--')
        ax.legend(title="Historical Samples Added")

        # Save plot
        plot_path = output_dir / f"comparison_extant_{extant_count}.png"
        plt.savefig(plot_path, bbox_inches="tight")
        plt.close()
        logging.info(f"Saved plot: {plot_path}")


def main():
    logging.info("Starting plot generation")

    # Load data with temporal context
    df = load_comparison_data("accuracy", "trees/subsets")

    if df.empty:
        logging.error("No data loaded - aborting")
        return

    # Generate plots
    create_comparison_plots(df)
    logging.info("Plot generation complete")


if __name__ == "__main__":
    main()