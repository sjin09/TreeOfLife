#!/usr/bin/env python

import argparse
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
import sys

import pandas as pd


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Summarise dataset",
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
        help="file to write",
    )
    args = args[1:]
    return parser.parse_args(args)


TAXONOMIC_RANKS = ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


@dataclass
class Taxonomy:
    hierarchy: defaultdict = field(
        default_factory=lambda: defaultdict(
            lambda: defaultdict(
                lambda: defaultdict(
                    lambda: defaultdict(
                        lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
                    )
                )
            )
        )
    )
    count_per_rank: defaultdict = field(
        default_factory=lambda: defaultdict(
            lambda: defaultdict(
                lambda: defaultdict(
                    lambda: defaultdict(
                        lambda: defaultdict(
                            lambda: defaultdict(
                                lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
                            )
                        )
                    )
                )
            )
        )
    )

    # Method to append a new taxonomy entry dynamically
    def add_taxonomy_entry(self, row):
        current_level = self.hierarchy  # Start at the root level of the hierarchy

        # Iterate over each rank and dynamically navigate through the hierarchy
        for rank in TAXONOMIC_RANKS:
            # Ensure current_level is a defaultdict (or dict) at each step
            current_level = current_level[row[rank]]  # Move to the next level

        # Update the sample count at each taxonomic level
        self.update_member_counts(row)

    def update_member_counts(self, row):
        current_level = self.count_per_rank
        for rank in TAXONOMIC_RANKS:
            current_level = current_level[row[rank]]
            if "Count" not in current_level:
                current_level["Count"] = 0
            current_level["Count"] += 1


def get_taxonomic_hierarchy(input_path: Path):
    # Initialize variables
    taxonomy = Taxonomy()

    # Read the data from CSV
    df = pd.read_csv(input_path)

    # Iterate through each row in the dataframe
    for (_op_idx, row) in df.iterrows():
        taxonomy.add_taxonomy_entry(row)

    return taxonomy


def get_taxonomic_classification_statistics(input_path: Path, output_path: Path) -> None:
    # Initialize variables
    taxonomy = get_taxonomic_hierarchy(input_path)
    # Write taxonomy and count of members belong to each taxonomic rank
    with open(output_path, "w") as outfile:
        for domain, kingdoms in taxonomy.hierarchy.items():
            outfile.write("Domain: {} (n = {})\n".format(domain, taxonomy.count_per_rank[domain]["Count"]))
            for kingdom, phylums in kingdoms.items():
                outfile.write(
                    "{}Kingdom: {} (n = {})\n".format(
                        " " * 4, kingdom, taxonomy.count_per_rank[domain][kingdom]["Count"]
                    )
                )
                for phylum, tax_classes in phylums.items():
                    outfile.write(
                        "{}Phylum: {} (n = {})\n".format(
                            " " * 8, phylum, taxonomy.count_per_rank[domain][kingdom][phylum]["Count"]
                        )
                    )
                    for tax_class, orders in tax_classes.items():
                        outfile.write(
                            "{}Class: {} (n = {})\n".format(
                                " " * 12,
                                tax_class,
                                taxonomy.count_per_rank[domain][kingdom][phylum][tax_class]["Count"],
                            )
                        )
                        for order, families in orders.items():
                            outfile.write(
                                "{}Order: {} (n = {})\n".format(
                                    " " * 16,
                                    order,
                                    taxonomy.count_per_rank[domain][kingdom][phylum][tax_class][order]["Count"],
                                )
                            )
                            for family, genera in families.items():
                                outfile.write(
                                    "{}Family: {} (n = {})\n".format(
                                        " " * 20,
                                        family,
                                        taxonomy.count_per_rank[domain][kingdom][phylum][tax_class][order][family][
                                            "Count"
                                        ],
                                    )
                                )
                                for genus, species in genera.items():
                                    outfile.write(
                                        "{}Genus: {} (n = {})\n".format(
                                            " " * 24,
                                            genus,
                                            taxonomy.count_per_rank[domain][kingdom][phylum][tax_class][order][family][
                                                genus
                                            ]["Count"],
                                        )
                                    )
                                    for species_name in species:
                                        outfile.write(
                                            "{}Species: {} (n = {})\n".format(
                                                " " * 28,
                                                species_name,
                                                taxonomy.count_per_rank[domain][kingdom][phylum][tax_class][order][
                                                    family
                                                ][genus][species_name]["Count"],
                                            )
                                        )


def main() -> int:
    options = parse_args(sys.argv)
    get_taxonomic_classification_statistics(options.input, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
