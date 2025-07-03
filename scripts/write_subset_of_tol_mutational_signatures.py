#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

import natsort
import pandas as pd


def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--sig",
        type=Path,
        required=True,
        help="Path to the CSV file containing HDP somatic mutational signatures",
    )
    parser.add_argument(
        "--sig-names",
        type=Path,
        required=True,
        help="Path to a file containing the list of signature names to include in the subset",
    )
    parser.add_argument(
        "--rtol-sig-names",
        type=Path,
        required=True,
        help="Path to a file containing the list of artefactual signature names to exclude from the subset",
    )
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to a file to write")
    args = args[1:]
    return parser.parse_args(args)


def write_mutational_signature_subset(
    sig_path: Path, sig_names_path: Path, rtol_sig_names_path: Path, output_path: Path
):
    # Load signatures
    df = pd.read_csv(sig_path, sep=",", index_col=0)

    # Load signature names
    sig_names = [line.rstrip() for line in open(sig_names_path)]

    # Load artefactual signature names
    rtol_sig_names = [line.rstrip() for line in open(rtol_sig_names_path)]
    stol_sig_names = list(set(sig_names).difference(set(rtol_sig_names)))
    stol_sig_names = natsort.natsorted(stol_sig_names)

    # Subset data frame
    df_subset = df[stol_sig_names]

    # Write the subsetted data frame as CSV, preserving the row index
    df_subset.to_csv(output_path, index=True)


def main():
    options = parse_args(sys.argv)
    write_mutational_signature_subset(options.sig, options.sig_names, options.rtol_sig_names, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
