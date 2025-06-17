#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

import pandas as pd


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the CSV file containing phylogentic analysis of mutational signatures",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to a file to write"
    )
    args = args[1:]
    return parser.parse_args(args)


def write_mutational_signatures_with_phylogenetic_signal(input_path: Path, output_path: Path):
    # Load signatures
    df = pd.read_csv(input_path)

    # Ensure q_value column is numeric
    df['q_value'] = pd.to_numeric(df['q_value'], errors='coerce')

    # filter to only rows with q_value < 0.01
    df_filtered = df[df['q_value'] < 0.01]

    # Write mutational signatures with phylogentic signal analysis
    with open(output_path, "w") as outfile:
        for df_idx, row in df_filtered.iterrows():
            outfile.write("X{}\n".format(df_idx + 1))


def main():
    options = parse_args(sys.argv)
    write_mutational_signatures_with_phylogenetic_signal(options.input, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
