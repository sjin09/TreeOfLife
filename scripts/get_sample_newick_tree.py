#!/usr/bin/env python

import argparse
from collections import defaultdict
from io import StringIO
from pathlib import Path
from typing import Dict, List, Union
import sys

from Bio import Phylo
import pandas as pd

TaxonomicTree = Union[Dict[str, Union["TaxonomicTree", List[str]]], List[str]]


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


def write_sample_newick_tree(input_path: Path, output_path: Path) -> None:
    def get_newick_tree(nested_taxonomy: TaxonomicTree) -> str:
        if isinstance(nested_taxonomy, list):
            species = [f"'{element}'" for element in nested_taxonomy]
            return "(" + ",".join(species) + ")"
        if isinstance(nested_taxonomy, dict):
            return "(" + ",".join(f"{key}:{get_newick_tree(value)}" for key, value in nested_taxonomy.items()) + ")"
        return str(nested_taxonomy)

    # Initialize
    nested_taxonomy = defaultdict(
        lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
    )

    # Import
    df = pd.read_csv(input_path)
    for _ridx, row in df.iterrows():
        species = row["Species"]
        kingdom = row["Kingdom"]
        phylum = row["Phylum"]
        subphylum = row["Subphylum"]
        class_taxon = row["Class"]
        order = row["Order"]
        if species not in nested_taxonomy[kingdom][phylum][subphylum][class_taxon][order]:
            nested_taxonomy[kingdom][phylum][subphylum][class_taxon][order].append(species)

    # Reformat
    newick_tree = get_newick_tree(nested_taxonomy)
    tree = Phylo.read(StringIO(newick_tree), "newick")
    ## Phylo.draw_ascii(tree)
    Phylo.write(tree, output_path, "newick")


def main() -> int:
    options = parse_args(sys.argv)
    write_sample_newick_tree(options.input, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
