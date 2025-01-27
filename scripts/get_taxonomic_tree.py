#!/usr/bin/env python

import argparse
from pathlib import Path
import sys

from ete3 import Tree


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

## Use Bio.Phylo and Bio.Phylo.Netwick
## def write_taxa_newick_tree(input_path: Path, output_path: Path) -> None:
##     def get_child(parent_clade: Clade, child_name: str):
##         for child in parent_clade.clades:
##             if child.name == child_name:
##                 return child
##         new_child = Clade(name=child_name)
##         parent_clade.clades.append(new_child)
##         return new_child
## 
##     root = Clade(name="root")  # Root node for the tree
##     for line in open(input_path):
##         if line.startswith("Reference"):
##             continue
##         current_node = root
##         taxa = line.rstrip().split(",")[5:]
##         for rank in taxa:
##             if rank == ".":
##                 continue
##             current_node = get_child(current_node, rank)
##     tree = Tree(root=root, rooted=True)
##     # Phylo.draw_ascii(tree)
##     Phylo.write(tree, output_path, "newick")


# Use ete3
def write_taxonomic_tree(input_path: Path, output_path: Path) -> None:
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

    # Build tree based on taxonomic classifications per sample
    with open(input_path) as infile:
        for line in infile:
            if line.startswith("Reference"):
                continue
            current_node = root
            taxonomic_ranks = line.rstrip().split(",")[5:]
            for rank in taxonomic_ranks:
                if rank == ".":
                    continue
                current_node = get_child(current_node, rank)

    # Export tree to newick format
    root.write(format=1, outfile=str(output_path))

    # Show tree for debugging
    # root.show()

def main() -> int:
    options = parse_args(sys.argv)
    write_taxonomic_tree(options.input, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
