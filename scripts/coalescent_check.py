import tskit
import argparse
from pathlib import Path


def has_unary_descendants(tree, node):
    """
    Check if a node has any unary descendants.
    """
    stack = [node]
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


def check_sexual_diploid_coalescence(tree, ts):
    """
    Check if roots belong to at most two diploid individuals.
    """
    if tree.num_roots <= 2:
        return True

    roots = list(tree.roots)
    root_individuals = [ts.node(root).individual for root in roots]
    specified_individuals = set(ind for ind in root_individuals if ind != -1)
    return len(specified_individuals) <= 2


def analyze_coalescence(ts):
    """
    Analyze coalescence properties of a tree sequence for a sexual diploid population.
    """
    has_complete_coalescence = True
    max_tmrca = float('-inf')
    max_roots = 0
    trees_with_max_tmrca = 0
    total_trees = ts.num_trees
    has_unary_root_descendants = False
    all_trees_properly_coalesced = True

    for tree in ts.trees():
        proper_coalescence = check_sexual_diploid_coalescence(tree, ts)
        if not proper_coalescence:
            has_complete_coalescence = False
            all_trees_properly_coalesced = False

        if tree.num_roots > 0:
            tree_tmrca = max(tree.time(root) for root in tree.roots)
            max_tmrca = max(max_tmrca, tree_tmrca)

        max_roots = max(max_roots, tree.num_roots)

        for root in tree.roots:
            if has_unary_descendants(tree, root):
                has_unary_root_descendants = True

    for tree in ts.trees():
        if tree.num_roots > 0:
            tree_tmrca = max(tree.time(root) for root in tree.roots)
            if abs(tree_tmrca - max_tmrca) < 1e-10:
                trees_with_max_tmrca += 1

    fraction_max_tmrca = trees_with_max_tmrca / total_trees if total_trees > 0 else 0

    return (has_complete_coalescence, max_tmrca, max_roots,
            fraction_max_tmrca, has_unary_root_descendants,
            all_trees_properly_coalesced)


def main():
    parser = argparse.ArgumentParser(description='Analyze coalescence properties of a tree sequence.')
    parser.add_argument('tree_name', type=str, help='Name of the tree sequence file (without .trees extension)')
    args = parser.parse_args()

    tree_path = Path('./trees') / f'{args.tree_name}.trees'

    if not tree_path.exists():
        print(f"Error: Tree sequence file not found at {tree_path}")
        return 1

    ts = tskit.load(str(tree_path))

    (coalescence_status, max_tmrca, max_roots, fraction_max_tmrca,
     has_unary, all_properly_coalesced) = analyze_coalescence(ts)

    print(f"\nAnalysis for {args.tree_name}:")
    print("-" * 80)
    print(f"Complete coalescence (allowing up to 2 ancestral individuals): {'Yes' if coalescence_status else 'No'}")
    print(f"Maximum TMRCA across trees: {max_tmrca:.2f}")
    print(f"Maximum number of roots across trees: {max_roots}")
    print(f"Fraction of trees with maximum TMRCA: {fraction_max_tmrca:.3f}")
    print(f"Has unary descendants from roots: {'Yes' if has_unary else 'No'}")
    print(f"All trees coalesce to ≤2 individuals: {'Yes' if all_properly_coalesced else 'No'}")

    return 0


if __name__ == "__main__":
    exit(main())