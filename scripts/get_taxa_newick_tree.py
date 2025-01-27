#!/usr/bin/env python

import argparse
from collections import defaultdict
from io import StringIO
from pathlib import Path
from typing import Dict, List, Union
import sys

from Bio import Phylo
from Bio.Phylo.Newick import Clade, Tree
import pandas as pd


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Get newick tree from taxonomic classification",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="file to read taxonomic classifications per sample",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to write newick tree (.nwk)",
    )
    args = args[1:]
    return parser.parse_args(args)


def write_taxa_newick_tree(input_path: Path, output_path: Path) -> None:
    def get_child(parent_clade: Clade, child_name: str):
        for child in parent_clade.clades:
            if child.name == child_name:
                return child
        new_child = Clade(name=child_name)
        parent_clade.clades.append(new_child)
        return new_child

    root = Clade(name="root")  # Root node for the tree
    for line in open(input_path):
        if line.startswith("Reference"):
            continue
        current_node = root
        taxa = line.rstrip().split(",")[5:]
        for rank in taxa:
            if rank == ".":
                continue
            current_node = get_child(current_node, rank)
    tree = Tree(root=root, rooted=True)
    Phylo.write(tree, output_path, "newick")


def main() -> int:
    options = parse_args(sys.argv)
    write_taxa_newick_tree(options.input, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
