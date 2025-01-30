#!/usr/bin/env python

import argparse
from pathlib import Path
import sys

from ete3 import Tree
import pandas as pd

TAXONOMIC_RANKS = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Sample']


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


# Use ete3
def write_sample_taxonomic_tree(input_path: Path, output_path: Path) -> None:
    def get_child(parent_clade, child_name: str):
        """
        Find the child with the specified name or create a new one if not found.
        """
        for child in parent_clade.get_children():
            if child.name == child_name:
                return child
        new_child = parent_clade.add_child(name=child_name)
        return new_child

    # Initialize tree
    root = Tree(name="root")

    # Import and iterate through the DToL samplesheet
    df = pd.read_csv(input_path)
    for (_idx, row) in df.iterrows():
        current_node = root
        sample_taxonomic_ranks = list(row[TAXONOMIC_RANKS])
        for rank in sample_taxonomic_ranks:
            if rank == ".":
                continue
            current_node = get_child(current_node, rank)

    # Export tree to newick format
    # format = 4 compatible with Julia Phylo package
    root.write(format=1, outfile=str(output_path))

    # Show tree for debugging
    # root.show()


def main() -> int:
    options = parse_args(sys.argv)
    write_sample_taxonomic_tree(options.input, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
