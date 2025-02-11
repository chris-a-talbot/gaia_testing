# Tree subsetter

import tskit
import numpy as np
from pathlib import Path
import pandas as pd
from sklearn.cluster import KMeans


def get_individual_times(tree):
    """
    Get time for each individual by looking at their nodes.
    Returns dict mapping individual ID to time.
    """
    times = {}
    for ind in tree.individuals():
        node_time = tree.nodes()[ind.nodes[0]].time
        times[ind.id] = node_time
    return times


def create_temporal_sample(tree, sample_times, num_ancient=4, temporal_spacing="even",
                           time_period=None):
    """
    Create a temporal sample of individuals.

    Parameters:
    -----------
    tree: tskit.TreeSequence
    sample_times: dict
        Mapping of individual IDs to their times
    num_ancient: int
        Number of ancient samples to include
    temporal_spacing: str
        "even", "random", or "clustered"
    time_period: float
        For clustered sampling, focus on most recent fraction of time
        (0.25, 0.5, or 0.75)

    Returns:
    --------
    list of unique individual IDs to keep
    """
    # Get modern samples (time = 1.0) - ensure unique and exist in sample_times
    modern_inds = list(set(ind.id for ind in tree.individuals()
                           if ind.id in sample_times and sample_times[ind.id] == 1.0))

    # Get ancient individuals sorted by time - ensure unique
    ancient_inds = list(set((ind_id, time) for ind_id, time in sample_times.items()
                            if time > 1.0))
    ancient_inds.sort(key=lambda x: x[1])

    if len(ancient_inds) == 0:
        return modern_inds

    if num_ancient > len(ancient_inds):
        raise ValueError(f"Requested {num_ancient} ancient samples but only {len(ancient_inds)} available")

    max_time = max(t for _, t in ancient_inds)

    if temporal_spacing == "even":
        # For even spacing, we want exactly num_ancient samples
        indices = np.linspace(0, len(ancient_inds) - 1, num_ancient).astype(int)
        selected_ancient = [ancient_inds[i][0] for i in indices]

    elif temporal_spacing == "clustered":
        if time_period not in [0.25, 0.5, 0.75]:
            raise ValueError("time_period must be 0.25, 0.5, or 0.75")

        time_cutoff = max_time * 0.01 * (1 - time_period)
        recent_inds = list(set((i, t) for i, t in ancient_inds if t <= time_cutoff))

        if len(recent_inds) < num_ancient:
            raise ValueError(f"Not enough samples in time period (have {len(recent_inds)}, need {num_ancient})")

        # Use random choice without replacement to ensure uniqueness
        selected_ancient = np.random.choice(
            [x[0] for x in recent_inds],
            size=num_ancient,
            replace=False
        )

    else:  # random
        # Use random choice without replacement to ensure uniqueness
        selected_ancient = np.random.choice(
            [x[0] for x in ancient_inds],
            size=num_ancient,
            replace=False
        )

    # Combine and ensure final uniqueness
    return list(set(modern_inds + list(selected_ancient)))


def create_combined_sample(tree, num_extant, num_ancient, temporal_spacing="even"):
    """
    Create a sample with both extant and ancient individuals.

    Parameters:
    -----------
    tree: tskit.TreeSequence
    num_extant: int
        Number of extant individuals to keep
    num_ancient: int
        Number of ancient samples to include
    temporal_spacing: str
        "even", "random", or "clustered"

    Returns:
    --------
    list of unique individual IDs to keep
    """
    times = get_individual_times(tree)

    # Get unique extant sample
    extant_inds = create_extant_sample(tree, num_extant)

    # Create modified times dict excluding selected extant individuals
    reduced_times = {k: v for k, v in times.items()
                     if k not in set(extant_inds)}

    # Get ancient sample from remaining individuals
    ancient_sample = create_temporal_sample(
        tree, reduced_times, num_ancient, temporal_spacing
    )

    # Combine and ensure final uniqueness
    return list(set(extant_inds + ancient_sample))


def create_extant_sample(tree, num_extant):
    """
    Create a sample of only extant individuals (time=1.0).

    Returns:
    --------
    list of unique individual IDs
    """
    times = get_individual_times(tree)
    # Get unique extant individuals
    extant_inds = list(set(ind_id for ind_id, time in times.items() if time == 1.0))

    if num_extant > len(extant_inds):
        raise ValueError(f"Requested {num_extant} extant individuals but only {len(extant_inds)} available")

    # Use random choice without replacement to ensure uniqueness
    return list(np.random.choice(extant_inds, size=num_extant, replace=False))


def create_geographic_sample(tree, num_regions=4, min_per_region=5):
    """
    Create a geographically stratified sample.

    Returns:
    --------
    list of unique individual IDs
    """
    # Get xy coordinates (ignore z=0)
    locations = tree.individual_locations[:, :2]

    # Cluster into regions
    kmeans = KMeans(n_clusters=num_regions)
    labels = kmeans.fit_predict(locations)

    # Sample from each region ensuring uniqueness
    kept_inds = set()
    for i in range(num_regions):
        region_inds = np.where(labels == i)[0]
        if len(region_inds) == 0:
            continue

        n_sample = max(min_per_region, len(region_inds) // 10)
        n_sample = min(n_sample, len(region_inds))  # Ensure we don't try to sample more than available

        # Use random choice without replacement for each region
        region_sample = np.random.choice(region_inds, size=n_sample, replace=False)
        kept_inds.update(region_sample)

    return list(kept_inds)


def simplify_to_individuals(tree, keep_individual_ids):
    """
    Simplify tree sequence to only include specified individuals.
    Ensures no duplicate nodes are included.
    """
    keep_nodes = []
    for ind_id in keep_individual_ids:
        keep_nodes.extend(tree.individuals()[ind_id].nodes)
    # Remove duplicates while maintaining order
    keep_nodes = list(dict.fromkeys(keep_nodes))
    return tree.simplify(keep_nodes)


def process_trees(input_pattern, output_dir):
    """
    Process multiple tree files with different subsetting schemes.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    for tree_path in Path('.').glob(input_pattern):
        tree = tskit.load(str(tree_path))
        stem = tree_path.stem
        times = get_individual_times(tree)
        num_extant = sum(1 for time in times.values() if time == 1.0)

        # Base configurations
        configs = [
            # Keep all samples
            ("all_samples", list(range(tree.num_individuals))),

            # Extant-only samples with wider range
            ("extant_10", create_extant_sample(tree, 10)),
            ("extant_25", create_extant_sample(tree, 25)),
            ("extant_50", create_extant_sample(tree, 50)),
            ("extant_100", create_extant_sample(tree, 100)),
            ("extant_250", create_extant_sample(tree, min(250, num_extant))),
            ("extant_500", create_extant_sample(tree, min(500, num_extant))),
            ("extant_1000", create_extant_sample(tree, min(1000, num_extant))),

            # Temporal samples with wider range
            ("temporal_even_4", create_temporal_sample(tree, times, 4, "even")),
            ("temporal_even_8", create_temporal_sample(tree, times, 8, "even")),
            ("temporal_even_16", create_temporal_sample(tree, times, 16, "even")),
            ("temporal_even_32", create_temporal_sample(tree, times, 32, "even")),
            ("temporal_random_4", create_temporal_sample(tree, times, 4, "random")),
            ("temporal_random_8", create_temporal_sample(tree, times, 8, "random")),
            ("temporal_random_16", create_temporal_sample(tree, times, 16, "random")),
            ("temporal_random_32", create_temporal_sample(tree, times, 32, "random")),

            # Clustered temporal samples
            ("temporal_clustered_4_25", create_temporal_sample(tree, times, 4, "clustered", 0.25)),
            ("temporal_clustered_8_25", create_temporal_sample(tree, times, 8, "clustered", 0.25)),
            ("temporal_clustered_16_25", create_temporal_sample(tree, times, 16, "clustered", 0.25)),
            ("temporal_clustered_4_50", create_temporal_sample(tree, times, 4, "clustered", 0.50)),
            ("temporal_clustered_8_50", create_temporal_sample(tree, times, 8, "clustered", 0.50)),
            ("temporal_clustered_16_50", create_temporal_sample(tree, times, 16, "clustered", 0.50)),
            ("temporal_clustered_4_75", create_temporal_sample(tree, times, 4, "clustered", 0.75)),
            ("temporal_clustered_8_75", create_temporal_sample(tree, times, 8, "clustered", 0.75)),
            ("temporal_clustered_16_75", create_temporal_sample(tree, times, 16, "clustered", 0.75)),

            # Combined samples (extant + ancient)
            ("combined_10_4", create_combined_sample(tree, 10, 4)),
            ("combined_10_8", create_combined_sample(tree, 10, 8)),
            ("combined_10_16", create_combined_sample(tree, 10, 16)),
            ("combined_10_32", create_combined_sample(tree, 10, 32)),
            ("combined_25_4", create_combined_sample(tree, 25, 4)),
            ("combined_25_8", create_combined_sample(tree, 25, 8)),
            ("combined_25_16", create_combined_sample(tree, 25, 16)),
            ("combined_25_32", create_combined_sample(tree, 25, 32)),
            ("combined_50_4", create_combined_sample(tree, 50, 4)),
            ("combined_50_8", create_combined_sample(tree, 50, 8)),
            ("combined_50_16", create_combined_sample(tree, 50, 16)),
            ("combined_50_32", create_combined_sample(tree, 50, 32)),
            ("combined_100_4", create_combined_sample(tree, 100, 4)),
            ("combined_100_8", create_combined_sample(tree, 100, 8)),
            ("combined_100_16", create_combined_sample(tree, 100, 16)),
            ("combined_100_32", create_combined_sample(tree, 100, 32)),

            # Geographic samples
            ("geographic_4", create_geographic_sample(tree, 4)),
            ("geographic_8", create_geographic_sample(tree, 8)),
            ("geographic_16", create_geographic_sample(tree, 16))
        ]

        # Process each configuration
        for name, keep_ids in configs:
            subset = simplify_to_individuals(tree, keep_ids)
            out_path = output_dir / f"{stem}_{name}.trees"
            subset.dump(str(out_path))

            # Save metadata
            metadata = {
                'original_individuals': tree.num_individuals,
                'subset_individuals': subset.num_individuals,
                'original_nodes': tree.num_nodes,
                'subset_nodes': subset.num_nodes,
                'config': name
            }
            meta_path = out_path.with_suffix('.json')
            pd.Series(metadata).to_json(meta_path)


if __name__ == "__main__":
    process_trees("trees/tree-S0.3-R*.trees", "./trees/subsets")