#!/usr/bin/env python

import argparse
from pathlib import Path
from typing import List
import sys

import ete3
import numpy as np

# Reference: Testing for phylogenetic signal in phenotypic traits: New matrices of phylogenetic proximities


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Calculate Abouheif's C_mean measure of phylogenetic signal",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="file to read newick tree"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to write",
    )
    args = args[1:]
    return parser.parse_args(args)


def get_tree(input_path: Path) -> ete3.TreeNode:
    tree = ete3.Tree(input_path, format=4)
    return tree


def get_path_to_ancestor(ancestor: ete3.TreeNode, node: ete3.TreeNode) -> List[ete3.TreeNode]:
    path = []
    while node != ancestor:
        node = node.up
        path.append(node)
    return path


def get_branch_path(tree: ete3.TreeNode, node1: ete3.TreeNode, node2: ete3.TreeNode):
    common_ancestor = tree.get_common_ancestor(node1, node2)
    path = [common_ancestor]
    path.extend(get_path_to_ancestor(common_ancestor, node1))
    path.extend(get_path_to_ancestor(common_ancestor, node2))

    while node1 != common_ancestor:
        node1 = node1.up
        path.append(node1)

    while node2 != common_ancestor:
        node2 = node2.up
        path.append(node2)

    return path


def get_Abouheif_A_matrix(tree: ete3.TreeNode) -> np.ndarray:
    """
    Notes:
        - For a species i, aii is therefore the inverse of the product of
        the number of branches descending from each node from the species to the root.
        - For a couple (i, j), aij is the inverse of the product of
        the number of branches descending from each node in the path connecting i and j.
        - Abouheif proximity refers to the phylogenetic proximity underlying the test of Abouheif (see references).
             - Let P be the set of all the nodes in the path going from node1 to node2.
             - Let DDP be the number of direct descendants from each node in P.
             Then, the so-called ’Abouheif’ distance is the inverse of the product of all terms in DDP.
             - oriAbouheif returns a matrix with non-null diagonal elements, as formulated in Pavoine et al. (2008).
             - This matrix is bistochastic (all marginal sums equal 1),
             but this bilinear symmetric form does not give rise to a Moran’s index,
             since it requires a null diagonal.
            - Abouheif contains Abouheif’s proximities but has a null diagonal, giving rise to a Moran’s index.
    """
    # Get the leaves of the tree
    leaves = tree.get_leaves()
    n_leaves = len(leaves)

    # Initialize the matrix
    mtx = np.zeros((n_leaves, n_leaves))

    # Calculate the proximity scores for off-diagonal elements of the matrix
    for i in range(n_leaves):
        for j in range(n_leaves):
            if i != j:
                p = set(get_branch_path(tree, leaves[i], leaves[j]))
                ddp = [len(node.get_children()) for node in p]
                mtx[i, j] = 1 / np.prod(ddp)

    # Calculate the diagonal scores for originality of the species
    for k in range(n_leaves):
        mtx[k, k] = 1 - np.sum(mtx[k, :])

    return mtx


def get_Abouheif_C_mean(input_path: Path, output_path: Path):
    tree = get_tree(input_path)
    A_mtx = get_Abouheif_A_matrix(tree)


def main() -> int:
    options = parse_args(sys.argv)
    get_Abouheif_C_mean(options.input, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
