# Tree subsetter

import tskit
import numpy as np
import pandas as pd
from pathlib import Path
import json
from typing import List, Dict
from enum import Enum
from dataclasses import dataclass
import random
import logging
import sys
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
    ]
)
logger = logging.getLogger(__name__)


class TemporalCoverage(Enum):
    QUARTER = 0.25
    HALF = 0.50
    THREE_QUARTERS = 0.75
    FULL = 1.0


class TemporalDistribution(Enum):
    EVEN = "even"
    RANDOM = "random"
    CLUSTERED = "clustered"


@dataclass
class SamplingScheme:
    name: str
    n_samples: int
    ancient_ratio: float = 0
    temporal_coverage: TemporalCoverage = TemporalCoverage.FULL
    temporal_distribution: TemporalDistribution = TemporalDistribution.EVEN
    n_clusters: int = 3  # For clustered distribution


class AncientSampler:
    def __init__(self, tree_path: str, locations_path: str):
        """Initialize sampler with tree and location data"""
        logger.info(f"Initializing AncientSampler with tree: {tree_path}, locations: {locations_path}")

        try:
            self.ts = tskit.load(tree_path)
            logger.info(f"Loaded tree sequence with {self.ts.num_nodes} nodes")
        except Exception as e:
            logger.error(f"Failed to load tree sequence: {e}")
            raise

        try:
            self.sample_df = pd.read_csv(locations_path)
            logger.info(f"Loaded sample data with {len(self.sample_df)} entries")
        except Exception as e:
            logger.error(f"Failed to load sample locations: {e}")
            raise

        self._validate_data()
        self._compute_time_bounds()

        # Define a constant for majority ancient ratio
        self.MAJORITY_ANCIENT_RATIO = 0.9  # 90% ancient, 10% extant
        logger.info(f"Set majority ancient ratio to {self.MAJORITY_ANCIENT_RATIO}")

        # Define sampling schemes
        self.schemes = [
            # 1. EXTANT-ONLY BASELINE SCHEMES
            SamplingScheme("n25_extant_only", n_samples=25),
            SamplingScheme("n50_extant_only", n_samples=50),
            SamplingScheme("n100_extant_only", n_samples=100),
            SamplingScheme("n200_extant_only", n_samples=200),
            SamplingScheme("n400_extant_only", n_samples=400),

            # 2. ANCIENT-ONLY SCHEMES
            # 2a. Varying sample sizes with even distribution
            SamplingScheme("n50_ancient90_even", n_samples=50, ancient_ratio=self.MAJORITY_ANCIENT_RATIO,
                           temporal_distribution=TemporalDistribution.EVEN),
            SamplingScheme("n100_ancient90_even", n_samples=100, ancient_ratio=self.MAJORITY_ANCIENT_RATIO,
                           temporal_distribution=TemporalDistribution.EVEN),
            SamplingScheme("n200_ancient90_even", n_samples=200, ancient_ratio=self.MAJORITY_ANCIENT_RATIO,
                           temporal_distribution=TemporalDistribution.EVEN),

            # 2b. Varying sample sizes with random distribution
            SamplingScheme("n50_ancient90_random", n_samples=50, ancient_ratio=self.MAJORITY_ANCIENT_RATIO,
                           temporal_distribution=TemporalDistribution.RANDOM),
            SamplingScheme("n100_ancient90_random", n_samples=100, ancient_ratio=self.MAJORITY_ANCIENT_RATIO,
                           temporal_distribution=TemporalDistribution.RANDOM),
            SamplingScheme("n200_ancient90_random", n_samples=200, ancient_ratio=self.MAJORITY_ANCIENT_RATIO,
                           temporal_distribution=TemporalDistribution.RANDOM),

            # 2c. Varying sample sizes with clustered distribution
            SamplingScheme("n50_ancient90_clustered", n_samples=50, ancient_ratio=self.MAJORITY_ANCIENT_RATIO,
                           temporal_distribution=TemporalDistribution.CLUSTERED),
            SamplingScheme("n100_ancient90_clustered", n_samples=100, ancient_ratio=self.MAJORITY_ANCIENT_RATIO,
                           temporal_distribution=TemporalDistribution.CLUSTERED),
            SamplingScheme("n200_ancient90_clustered", n_samples=200, ancient_ratio=self.MAJORITY_ANCIENT_RATIO,
                           temporal_distribution=TemporalDistribution.CLUSTERED),

            # 2d. Full coverage with varying cluster numbers
            SamplingScheme("n100_ancient90_2clusters", n_samples=100, ancient_ratio=self.MAJORITY_ANCIENT_RATIO,
                           temporal_distribution=TemporalDistribution.CLUSTERED, n_clusters=2),
            SamplingScheme("n100_ancient90_5clusters", n_samples=100, ancient_ratio=self.MAJORITY_ANCIENT_RATIO,
                           temporal_distribution=TemporalDistribution.CLUSTERED, n_clusters=5),
            SamplingScheme("n100_ancient90_8clusters", n_samples=100, ancient_ratio=self.MAJORITY_ANCIENT_RATIO,
                           temporal_distribution=TemporalDistribution.CLUSTERED, n_clusters=8),

            # 3. MIXED SCHEMES - TEMPORAL COVERAGE VARIATION
            # 3a. Quarter coverage with different distributions
            SamplingScheme("n200_ancient50_quarter_even", n_samples=200, ancient_ratio=0.5,
                           temporal_coverage=TemporalCoverage.QUARTER,
                           temporal_distribution=TemporalDistribution.EVEN),
            SamplingScheme("n200_ancient50_quarter_random", n_samples=200, ancient_ratio=0.5,
                           temporal_coverage=TemporalCoverage.QUARTER,
                           temporal_distribution=TemporalDistribution.RANDOM),
            SamplingScheme("n200_ancient50_quarter_clustered", n_samples=200, ancient_ratio=0.5,
                           temporal_coverage=TemporalCoverage.QUARTER,
                           temporal_distribution=TemporalDistribution.CLUSTERED),

            # 3b. Half coverage with different distributions
            SamplingScheme("n200_ancient50_half_even", n_samples=200, ancient_ratio=0.5,
                           temporal_coverage=TemporalCoverage.HALF,
                           temporal_distribution=TemporalDistribution.EVEN),
            SamplingScheme("n200_ancient50_half_random", n_samples=200, ancient_ratio=0.5,
                           temporal_coverage=TemporalCoverage.HALF,
                           temporal_distribution=TemporalDistribution.RANDOM),
            SamplingScheme("n200_ancient50_half_clustered", n_samples=200, ancient_ratio=0.5,
                           temporal_coverage=TemporalCoverage.HALF,
                           temporal_distribution=TemporalDistribution.CLUSTERED),

            # 3c. Three-quarter coverage with different distributions
            SamplingScheme("n200_ancient50_75pct_even", n_samples=200, ancient_ratio=0.5,
                           temporal_coverage=TemporalCoverage.THREE_QUARTERS,
                           temporal_distribution=TemporalDistribution.EVEN),
            SamplingScheme("n200_ancient50_75pct_random", n_samples=200, ancient_ratio=0.5,
                           temporal_coverage=TemporalCoverage.THREE_QUARTERS,
                           temporal_distribution=TemporalDistribution.RANDOM),
            SamplingScheme("n200_ancient50_75pct_clustered", n_samples=200, ancient_ratio=0.5,
                           temporal_coverage=TemporalCoverage.THREE_QUARTERS,
                           temporal_distribution=TemporalDistribution.CLUSTERED),

            # 4. MIXED SCHEMES - RATIO VARIATION
            # 4a. Low ancient ratio with different distributions
            SamplingScheme("n200_ancient25_even", n_samples=200, ancient_ratio=0.25,
                           temporal_distribution=TemporalDistribution.EVEN),
            SamplingScheme("n200_ancient25_random", n_samples=200, ancient_ratio=0.25,
                           temporal_distribution=TemporalDistribution.RANDOM),
            SamplingScheme("n200_ancient25_clustered", n_samples=200, ancient_ratio=0.25,
                           temporal_distribution=TemporalDistribution.CLUSTERED),

            # 4b. Balanced ratio with different distributions
            SamplingScheme("n200_ancient50_even", n_samples=200, ancient_ratio=0.5,
                           temporal_distribution=TemporalDistribution.EVEN),
            SamplingScheme("n200_ancient50_random", n_samples=200, ancient_ratio=0.5,
                           temporal_distribution=TemporalDistribution.RANDOM),
            SamplingScheme("n200_ancient50_clustered", n_samples=200, ancient_ratio=0.5,
                           temporal_distribution=TemporalDistribution.CLUSTERED),

            # 4c. High ancient ratio with different distributions
            SamplingScheme("n200_ancient75_even", n_samples=200, ancient_ratio=0.75,
                           temporal_distribution=TemporalDistribution.EVEN),
            SamplingScheme("n200_ancient75_random", n_samples=200, ancient_ratio=0.75,
                           temporal_distribution=TemporalDistribution.RANDOM),
            SamplingScheme("n200_ancient75_clustered", n_samples=200, ancient_ratio=0.75,
                           temporal_distribution=TemporalDistribution.CLUSTERED),

            # 5. MIXED SCHEMES - SAMPLE SIZE VARIATION
            # 5a. Small total sample size with different ratios
            SamplingScheme("n50_ancient50", n_samples=50, ancient_ratio=0.5),
            SamplingScheme("n50_ancient75", n_samples=50, ancient_ratio=0.75),
            SamplingScheme("n50_ancient25", n_samples=50, ancient_ratio=0.25),

            # 5b. Medium total sample size with different ratios
            SamplingScheme("n100_ancient50", n_samples=100, ancient_ratio=0.5),
            SamplingScheme("n100_ancient75", n_samples=100, ancient_ratio=0.75),
            SamplingScheme("n100_ancient25", n_samples=100, ancient_ratio=0.25),

            # 5c. Large total sample size with different ratios
            SamplingScheme("n400_ancient50", n_samples=400, ancient_ratio=0.5),
            SamplingScheme("n400_ancient75", n_samples=400, ancient_ratio=0.75),
            SamplingScheme("n400_ancient25", n_samples=400, ancient_ratio=0.25),

            # 6. SPECIAL COMBINATIONS FOR INTERACTION EFFECTS
            # 6a. Size + Coverage + Distribution interactions
            SamplingScheme("n50_ancient50_quarter_clustered", n_samples=50, ancient_ratio=0.5,
                           temporal_coverage=TemporalCoverage.QUARTER,
                           temporal_distribution=TemporalDistribution.CLUSTERED),
            SamplingScheme("n400_ancient50_quarter_clustered", n_samples=400, ancient_ratio=0.5,
                           temporal_coverage=TemporalCoverage.QUARTER,
                           temporal_distribution=TemporalDistribution.CLUSTERED),
            SamplingScheme("n50_ancient50_full_clustered", n_samples=50, ancient_ratio=0.5,
                           temporal_coverage=TemporalCoverage.FULL,
                           temporal_distribution=TemporalDistribution.CLUSTERED),
            SamplingScheme("n400_ancient50_full_clustered", n_samples=400, ancient_ratio=0.5,
                           temporal_coverage=TemporalCoverage.FULL,
                           temporal_distribution=TemporalDistribution.CLUSTERED),

            # 6b. Ratio + Coverage + Distribution interactions
            SamplingScheme("n200_ancient25_quarter_even", n_samples=200, ancient_ratio=0.25,
                           temporal_coverage=TemporalCoverage.QUARTER,
                           temporal_distribution=TemporalDistribution.EVEN),
            SamplingScheme("n200_ancient75_quarter_even", n_samples=200, ancient_ratio=0.75,
                           temporal_coverage=TemporalCoverage.QUARTER,
                           temporal_distribution=TemporalDistribution.EVEN),
            SamplingScheme("n200_ancient25_full_even", n_samples=200, ancient_ratio=0.25,
                           temporal_coverage=TemporalCoverage.FULL,
                           temporal_distribution=TemporalDistribution.EVEN),
            SamplingScheme("n200_ancient75_full_even", n_samples=200, ancient_ratio=0.75,
                           temporal_coverage=TemporalCoverage.FULL,
                           temporal_distribution=TemporalDistribution.EVEN),
        ]
        logger.info(f"Initialized {len(self.schemes)} sampling schemes")

    def _validate_data(self):
        """Ensure required columns exist in sample data"""
        required_cols = ['node_id', 'individual_id', 'pedigree_id', 'time', 'is_ancient']
        missing = [col for col in required_cols if col not in self.sample_df.columns]
        if missing:
            logger.error(f"Missing required columns: {missing}")
            raise ValueError(f"Missing required columns: {missing}")
        logger.info("Data validation successful")

    def _compute_time_bounds(self):
        """Compute time boundaries for ancient samples"""
        ancient_samples = self.sample_df[self.sample_df.is_ancient]
        self.min_ancient_time = ancient_samples.time.min()
        self.max_ancient_time = ancient_samples.time.max()
        self.time_span = self.max_ancient_time - self.min_ancient_time
        logger.info(
            f"Time bounds computed: min={self.min_ancient_time}, max={self.max_ancient_time}, span={self.time_span}")

    def _get_time_threshold(self, coverage: TemporalCoverage) -> float:
        """Get time threshold based on coverage"""
        threshold = self.max_ancient_time - (self.time_span * coverage.value)
        logger.debug(f"Time threshold for coverage {coverage.value}: {threshold}")
        return threshold

    def _select_ancient_samples(self, n_samples: int,
                                coverage: TemporalCoverage,
                                distribution: TemporalDistribution,
                                n_clusters: int = 3) -> pd.DataFrame:
        """Select ancient samples based on temporal parameters

        Args:
            n_samples: Number of samples to select
            coverage: Temporal coverage enum specifying time range to sample from
            distribution: Distribution strategy for sampling
            n_clusters: Number of clusters for clustered distribution

        Returns:
            DataFrame of selected samples

        Raises:
            ValueError: If not enough samples available
        """
        logger.info(f"Selecting {n_samples} ancient samples with coverage={coverage.value}, "
                    f"distribution={distribution.value}, clusters={n_clusters}")

        ancient = self.sample_df[self.sample_df.is_ancient].copy()
        time_threshold = self._get_time_threshold(coverage)
        eligible = ancient[ancient.time >= time_threshold].copy()

        logger.debug(f"Found {len(eligible)} eligible samples for time threshold {time_threshold}")

        if len(eligible) < n_samples:
            error_msg = f"Not enough ancient samples (requested {n_samples}, available {len(eligible)})"
            logger.error(error_msg)
            raise ValueError(error_msg)

        if distribution == TemporalDistribution.EVEN:
            logger.debug("Using EVEN distribution strategy")
            unique_times = eligible.time.unique()

            # Handle case where we have fewer unique timepoints than requested samples
            if len(unique_times) < n_samples:
                logger.warning(
                    f"Fewer unique times ({len(unique_times)}) than requested samples ({n_samples}). "
                    "Will select one sample per unique timepoint and randomly select remaining samples."
                )
                n_bins = len(unique_times)
            else:
                n_bins = n_samples

            bins = pd.qcut(eligible.time, n_bins, labels=False, duplicates='drop')
            selected = []

            # Select one sample from each bin
            for bin_id in range(bins.max() + 1):
                bin_samples = eligible[bins == bin_id]
                if not bin_samples.empty:
                    selected.append(bin_samples.sample(1))
                    logger.debug(f"Selected 1 sample from bin {bin_id} (time range: "
                                 f"{bin_samples.time.min():.2f}-{bin_samples.time.max():.2f})")

            selected = pd.concat(selected, ignore_index=True)

            # If we need additional samples, randomly select from remaining samples
            if len(selected) < n_samples:
                additional_needed = n_samples - len(selected)
                logger.info(
                    f"Need {additional_needed} additional samples to reach target. "
                    "These will be randomly selected from eligible timepoints."
                )
                remaining = eligible[~eligible.index.isin(selected.index)]
                if len(remaining) < additional_needed:
                    logger.warning(
                        f"Only {len(remaining)} remaining samples available for {additional_needed} needed. "
                        "Will use all remaining samples."
                    )
                additional = remaining.sample(min(additional_needed, len(remaining)))
                selected = pd.concat([selected, additional], ignore_index=True)

        elif distribution == TemporalDistribution.RANDOM:
            logger.debug("Using RANDOM distribution strategy")
            selected = eligible.sample(n_samples)

        else:  # CLUSTERED
            logger.debug(f"Using CLUSTERED distribution strategy with {n_clusters} clusters")
            cluster_times = np.linspace(time_threshold, self.max_ancient_time, n_clusters)
            logger.debug(f"Cluster center times: {cluster_times}")

            eligible.loc[:, 'cluster'] = eligible.time.apply(
                lambda x: np.argmin(np.abs(x - cluster_times)))

            # Calculate samples needed per cluster
            base_samples = n_samples // n_clusters
            extra_samples = n_samples % n_clusters
            samples_per_cluster = [base_samples + (1 if i < extra_samples else 0)
                                   for i in range(n_clusters)]

            logger.debug(f"Target samples per cluster: {samples_per_cluster}")

            selected = []
            for cluster in range(n_clusters):
                cluster_samples = eligible[eligible.cluster == cluster]
                if cluster_samples.empty:
                    logger.warning(f"No samples available in cluster {cluster}")
                    continue

                n_to_sample = min(samples_per_cluster[cluster], len(cluster_samples))
                if n_to_sample < samples_per_cluster[cluster]:
                    logger.warning(
                        f"Cluster {cluster} has fewer samples ({len(cluster_samples)}) "
                        f"than requested ({samples_per_cluster[cluster]})"
                    )

                logger.debug(f"Selecting {n_to_sample} samples from cluster {cluster} "
                             f"(available: {len(cluster_samples)})")
                selected.append(cluster_samples.sample(n_to_sample))

            selected = pd.concat(selected, ignore_index=True)

            # If we didn't get enough samples from the clusters, randomly select remaining
            if len(selected) < n_samples:
                additional_needed = n_samples - len(selected)
                logger.warning(
                    f"Only got {len(selected)} samples from clusters, need {additional_needed} more. "
                    "Will randomly select remaining samples from eligible timepoints."
                )
                remaining = eligible[~eligible.index.isin(selected.index)]
                additional = remaining.sample(min(additional_needed, len(remaining)))
                selected = pd.concat([selected, additional], ignore_index=True)

        logger.info(f"Selected {len(selected)} ancient samples spanning time range "
                    f"{selected.time.min():.2f} to {selected.time.max():.2f}")

        return selected.head(n_samples)

    def create_sampling_scheme(self, scheme: SamplingScheme) -> Dict:
        """Create a single sampling scheme with enhanced metadata

        Args:
            scheme: SamplingScheme object defining the sampling parameters

        Returns:
            Dictionary containing nodes and metadata for the scheme
        """
        logger.info(f"Creating sampling scheme: {scheme.name}")

        n_ancient = int(scheme.n_samples * scheme.ancient_ratio)
        n_extant = scheme.n_samples - n_ancient
        logger.debug(f"Target samples - ancient: {n_ancient}, extant: {n_extant}")

        extant = self.sample_df[~self.sample_df.is_ancient]
        logger.debug(f"Available extant samples: {len(extant)}")

        if scheme.ancient_ratio >= 0.9:
            logger.debug("Using consistent extant samples for majority ancient scheme")
            rng = np.random.RandomState(seed=hash(tuple(sorted(extant.individual_id.unique()))) % 2 ** 32)

            if not hasattr(self, '_baseline_extant_samples'):
                min_scheme_size = min(s.n_samples for s in self.schemes
                                      if getattr(s, 'ancient_ratio', 0) >= 0.9)
                n_baseline = int(min_scheme_size * (1 - self.MAJORITY_ANCIENT_RATIO))
                logger.info(f"Selecting {n_baseline} baseline extant samples")
                self._baseline_extant_samples = extant.sample(n=n_baseline, random_state=rng)

            selected_extant = self._baseline_extant_samples
        else:
            logger.debug(f"Randomly selecting {n_extant} extant samples")
            selected_extant = extant.sample(n_extant) if n_extant > 0 else pd.DataFrame()

        selected_ancient = self._select_ancient_samples(
            n_ancient,
            scheme.temporal_coverage,
            scheme.temporal_distribution,
            scheme.n_clusters
        ) if n_ancient > 0 else pd.DataFrame()

        selected = pd.concat([selected_extant, selected_ancient])
        logger.debug(f"Combined {len(selected_extant)} extant and {len(selected_ancient)} ancient samples")

        nodes = []
        for ind_id in selected.individual_id.unique():
            individual = self.ts.individual(ind_id)
            nodes.extend(individual.nodes)

        logger.info(f"Found {len(nodes)} nodes for {len(selected)} individuals")

        # Calculate actual ratios based on selected samples
        total_samples = len(selected)
        actual_ancient_ratio = len(selected_ancient) / total_samples if total_samples > 0 else 0.0

        # Only include n_clusters in metadata if using clustered distribution
        clustering_info = {}
        if scheme.temporal_distribution == TemporalDistribution.CLUSTERED:
            clustering_info['n_clusters'] = scheme.n_clusters

        metadata = {
            'name': scheme.name,
            'n_samples': len(selected),
            'n_nodes': len(nodes),
            'n_ancient': len(selected_ancient),
            'n_extant': len(selected_extant),
            'temporal_coverage': scheme.temporal_coverage.value,
            'temporal_distribution': scheme.temporal_distribution.value,
            'time_range': [float(selected.time.min()), float(selected.time.max())],
            'target_ancient_ratio': float(scheme.ancient_ratio),
            'actual_ancient_ratio': float(actual_ancient_ratio),
            **clustering_info  # This will only add n_clusters if it was set above
        }

        return {
            'name': scheme.name,
            'nodes': nodes,
            'metadata': metadata
        }

    def create_subsets(self, output_dir: Path):
        """Create multiple sampling schemes"""
        logger.info(f"Creating subsets in directory: {output_dir}")

        for scheme in self.schemes:
            try:
                logger.info(f"Processing scheme: {scheme.name}")
                result = self.create_sampling_scheme(scheme)

                # Save subsetted tree sequence
                tree_file = output_dir / f"{scheme.name}.trees"
                subset_ts = self.ts.simplify(result['nodes'])
                subset_ts.dump(str(tree_file))
                logger.info(f"Saved tree subset to {tree_file}")

                # Save metadata
                meta_file = output_dir / f"{scheme.name}_meta.json"
                with open(meta_file, 'w') as f:
                    json.dump(result['metadata'], f, indent=2)
                logger.info(f"Saved metadata to {meta_file}")

            except Exception as e:
                logger.error(f"Error creating scheme {scheme.name}: {str(e)}", exc_info=True)


def main():
    logger.info("Starting ancient sampling process")

    if len(sys.argv) != 2:
        logger.error("Incorrect number of arguments")
        print("Usage: python subset_tree.py <tree_name>")
        sys.exit(1)

    tree_name = sys.argv[1]
    logger.info(f"Processing tree: {tree_name}")

    # Setup paths
    tree_path = f"./trees/{tree_name}.trees"
    locations_path = f"./sample_locations/{tree_name}_locations.csv"
    output_dir = Path(f"./trees/subsets/{tree_name}")
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Output directory created: {output_dir}")

    # Create subsets
    try:
        sampler = AncientSampler(tree_path, locations_path)
        sampler.create_subsets(output_dir)
        logger.info("Sampling process completed successfully")
    except Exception as e:
        logger.error(f"Fatal error during sampling process: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()