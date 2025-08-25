#!/usr/bin/env python

import argparse
from pathlib import Path
from typing import List, Optional
import sys

from ete3 import Tree
import pandas as pd

TAXONOMIC_RANKS = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    if args is None:
        args = sys.argv[1:]

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
def write_newick_tree(input_path: Path, output_path: Path) -> None:
    def get_child(parent_clade, child_name: str):
        """
        Find the child with the specified name or create a new one if not found.
        """
        for child in parent_clade.get_children():
            if child.name == child_name:
                return child
        new_child = parent_clade.add_child(name=child_name)
        return new_child
    
    # Import and iterate through the DToL samplesheet
    df = pd.read_csv(input_path)

    # Subset data frame # Remove samples with incomplete taxonomic classification
    df_subset = df[~(df.eq(".").any(axis=1))]

    # Initialize tree
    root = Tree(name="root")

    # Build newick tree
    for (_idx, row) in df_subset.iterrows():
        current_node = root
        sample_taxonomic_ranks = list(row[TAXONOMIC_RANKS])
        for rank in sample_taxonomic_ranks:
            current_node = get_child(current_node, rank)

    # Export tree to newick format
    root.write(format=1, outfile=str(output_path))


def main() -> int:
    options = parse_args()
    write_newick_tree(options.input, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
