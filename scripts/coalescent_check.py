#!/usr/bin/env python3

# Checks a tree sequence for coalescence properties in a sexual diploid population
# Path: scripts/coalescent_check.py
# Run from the root directory with `python ./scripts/coalescent_check.py <tree_name>`

import sys
import logging
import tskit
import argparse
import math
from pathlib import Path
from typing import Tuple, List, Set

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)


def has_unary_descendants(tree: tskit.Tree, node: int) -> bool:
    """
    Check if a node has any unary descendants.

    Args:
        tree: The tree to check
        node: The node ID to start from

    Returns:
        bool: True if the node has any unary descendants
    """
    stack: List[int] = [node]
    children = tree.children(node)
    if len(children) == 1:
        stack.extend(children)
    else:
        stack.extend(children)

    while stack:
        current = stack.pop()
        children = tree.children(current)
        if len(children) == 1:
            return True
        stack.extend(children)
    return False


def check_sexual_diploid_coalescence(tree: tskit.Tree, ts: tskit.TreeSequence) -> bool:
    """
    Check if roots belong to at most two diploid individuals.

    Args:
        tree: The tree to check
        ts: The tree sequence containing the tree

    Returns:
        bool: True if the tree has proper coalescence
    """
    if tree.num_roots <= 2:
        return True

    roots = list(tree.roots)
    root_individuals = [ts.node(root).individual for root in roots]
    specified_individuals = {ind for ind in root_individuals if ind != -1}
    return len(specified_individuals) <= 2


def analyze_coalescence(ts: tskit.TreeSequence) -> Tuple[bool, float, int, float, bool, bool]:
    """
    Analyze coalescence properties of a tree sequence for a sexual diploid population.

    Args:
        ts: The tree sequence to analyze

    Returns:
        Tuple containing:
        - bool: Whether complete coalescence is achieved
        - float: Maximum TMRCA across all trees
        - int: Maximum number of roots across all trees
        - float: Fraction of trees with maximum TMRCA
        - bool: Whether any roots have unary descendants
        - bool: Whether all trees properly coalesce
    """
    has_complete_coalescence = True
    max_tmrca = float('-inf')
    max_roots = 0
    trees_with_max_tmrca = 0
    total_trees = ts.num_trees
    has_unary_root_descendants = False
    all_trees_properly_coalesced = True

    # First pass to find max TMRCA
    for tree in ts.trees():
        proper_coalescence = check_sexual_diploid_coalescence(tree, ts)
        if not proper_coalescence:
            has_complete_coalescence = False
            all_trees_properly_coalesced = False

        if tree.num_roots > 0:
            tree_tmrca = max(tree.time(root) for root in tree.roots)
            if math.isfinite(tree_tmrca):  # Check for valid TMRCA
                max_tmrca = max(max_tmrca, tree_tmrca)

        max_roots = max(max_roots, tree.num_roots)

        for root in tree.roots:
            if has_unary_descendants(tree, root):
                has_unary_root_descendants = True

    # Second pass to count trees with max TMRCA
    if math.isfinite(max_tmrca):
        for tree in ts.trees():
            if tree.num_roots > 0:
                tree_tmrca = max(tree.time(root) for root in tree.roots)
                if math.isfinite(tree_tmrca) and abs(tree_tmrca - max_tmrca) < 1e-10:
                    trees_with_max_tmrca += 1

    fraction_max_tmrca = trees_with_max_tmrca / total_trees if total_trees > 0 else 0

    return (has_complete_coalescence, max_tmrca, max_roots,
            fraction_max_tmrca, has_unary_root_descendants,
            all_trees_properly_coalesced)


def main() -> int:
    parser = argparse.ArgumentParser(description='Analyze coalescence properties of a tree sequence.')
    parser.add_argument('tree_name', type=str, help='Name of the tree sequence file (without .trees extension)')
    args = parser.parse_args()

    tree_path = Path('./trees') / f'{args.tree_name}.trees'

    if not tree_path.exists():
        logging.error(f"Tree sequence file not found at {tree_path}")
        return 1

    try:
        ts = tskit.load(str(tree_path))
    except Exception as e:
        logging.error(f"Failed to load tree sequence: {e}")
        return 1

    logging.info(f"Analyzing coalescence for {args.tree_name}")

    try:
        (coalescence_status, max_tmrca, max_roots, fraction_max_tmrca,
         has_unary, all_properly_coalesced) = analyze_coalescence(ts)
    except Exception as e:
        logging.error(f"Analysis failed: {e}")
        return 1

    print(f"\nAnalysis for {args.tree_name}:")
    print("-" * 80)
    print(f"Complete coalescence (allowing up to 2 ancestral individuals): {'Yes' if coalescence_status else 'No'}")
    print(f"Maximum TMRCA across trees: {max_tmrca:.2f}")
    print(f"Maximum number of roots across trees: {max_roots}")
    print(f"Fraction of trees with maximum TMRCA: {fraction_max_tmrca:.3f}")
    print(f"Has unary descendants from roots: {'Yes' if has_unary else 'No'}")
    print(f"All trees coalesce to ≤2 individuals: {'Yes' if all_properly_coalesced else 'No'}")

    if coalescence_status and all_properly_coalesced:
        print(0)
        return 0
    else:
        print(1)
        return 1


if __name__ == "__main__":
    sys.exit(main())