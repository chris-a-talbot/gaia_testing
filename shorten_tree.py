import tskit
import numpy as np


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


def create_extant_sample(tree, num_extant):
    """
    Create a sample of only extant individuals (time=1.0).
    """
    times = get_individual_times(tree)
    extant_inds = [ind_id for ind_id, time in times.items() if time == 1.0]

    if num_extant > len(extant_inds):
        raise ValueError(f"Requested {num_extant} extant individuals but only {len(extant_inds)} available")

    return list(np.random.choice(extant_inds, size=num_extant, replace=False))


def simplify_to_individuals(tree, keep_individual_ids):
    """
    Simplify tree sequence to only include specified individuals.
    """
    keep_nodes = []
    for ind_id in keep_individual_ids:
        keep_nodes.extend(tree.individuals()[ind_id].nodes)
    return tree.simplify(keep_nodes)


def main():
    # Load the first tree sequence (R1)
    tree_path = "trees/tree-S0.3-R1.trees"  # Update path if necessary
    tree = tskit.load(tree_path)

    # Create a sample of 250 extant individuals
    try:
        sample = create_extant_sample(tree, 250)
    except ValueError as e:
        print(f"Error creating sample: {e}")
        return

    # Simplify the tree to the selected individuals
    subset_tree = simplify_to_individuals(tree, sample)

    # Save the subset tree
    output_path = "trees/tree-S0.3-R1_subset_250.trees"  # Update path if necessary
    subset_tree.dump(output_path)
    print(f"Subset tree saved to {output_path}")


if __name__ == "__main__":
    main()