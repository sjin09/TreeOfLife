#!/usr/bin/env python

import argparse
from pathlib import Path
import sys

# from Bio import Phylo
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

# Use Bio.Phylo and Bio.Phylo.Netwick
# def write_sample_newick_tree(input_path: Path, output_path: Path) -> None:
#     def get_newick_tree(nested_taxonomy: TaxonomicTree) -> str:
#         if isinstance(nested_taxonomy, list):
#             species = [f"'{element}'" for element in nested_taxonomy]
#             return "(" + ",".join(species) + ")"
#         if isinstance(nested_taxonomy, dict):
#             return "(" + ",".join(f"{key}:{get_newick_tree(value)}" for key, value in nested_taxonomy.items()) + ")"
#         return str(nested_taxonomy)
#
#     # Initialize
#     nested_taxonomy = defaultdict(
#         lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
#     )
#
#     # Import
#     df = pd.read_csv(input_path)
#     for _ridx, row in df.iterrows():
#         species = row["Species"]
#         kingdom = row["Kingdom"]
#         phylum = row["Phylum"]
#         subphylum = row["Subphylum"]
#         class_taxon = row["Class"]
#         order = row["Order"]
#         if species not in nested_taxonomy[kingdom][phylum][subphylum][class_taxon][order]:
#             nested_taxonomy[kingdom][phylum][subphylum][class_taxon][order].append(species)
#
#     # Reformat
#     newick_tree = get_newick_tree(nested_taxonomy)
#     tree = Phylo.read(StringIO(newick_tree), "newick")
#     ## Phylo.draw_ascii(tree)
#     Phylo.write(tree, output_path, "newick")


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

    # Build tree based on taxonomic classifications per sample
    with open(input_path) as infile:
        for line in infile:
            if line.startswith("Reference"):
                continue
            current_node = root
            fields = line.rstrip().split(",")
            taxonomic_ranks = fields[5:] + [fields[2]]
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
    write_sample_taxonomic_tree(options.input, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
