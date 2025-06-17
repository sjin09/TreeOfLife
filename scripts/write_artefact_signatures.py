#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys
from typing import List

import numpy as np
from numpy import dot
from numpy.linalg import norm
import pandas as pd


NTS = ["A", "C", "G", "T"]
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SBS96_CLASSIFICATIONS = [f"{nti}[{sub}]{ntj}" for sub in SUBS for nti in NTS for ntj in NTS]


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--sigs",
        type=Path,
        required=True,
        help="Path to the CSV file containing phylogentic analysis of mutational signatures",
    )
    parser.add_argument(
        "--cord-blood",
        type=Path,
        required=True,
        help="Path to the TSV file containing corrected SBS96 counts from cord blood",
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


def calculate_cosine_similarity(va: List[float], vb: List[float]) -> float:
    a = np.asarray(va, dtype=float)
    b = np.asarray(vb, dtype=float)
    sim = dot(a, b)/(norm(a)*norm(b))
    return sim


def load_cord_blood_sbs96_counts(cord_blood_path: Path) -> List[float]:
    df = pd.read_csv(cord_blood_path, sep="\t")
    df = df.set_index('sbs96')
    cord_counts = df.reindex(SBS96_CLASSIFICATIONS)["normcounts"].values.astype(float).tolist()
    return cord_counts


def write_artefact_mutational_signatures(sigs_path: Path, cord_blood_path: Path, output_path: Path):
    # Load cord blood mutatioin counts
    cord_counts = load_cord_blood_sbs96_counts(cord_blood_path)

    # Load signatures
    sigs_df = pd.read_csv(sigs_path)
    sigs_df = sigs_df.iloc[:, 1:]
    sigs_df.index = SBS96_CLASSIFICATIONS
    sigs_df.index.name = "SBS96"

    # Calculate cosine similarity between signatures and cord blood mutation counts
    results = []
    for sig_name, col_series in sigs_df.items():
        vec = col_series.astype(float).tolist()
        sim = calculate_cosine_similarity(vec, cord_counts)
        results.append((sig_name, sim))

    # Write data frame
    df_sims = pd.DataFrame(results, columns=["Signature", "Similarity"])
    df_sims = df_sims.sort_values("Similarity", ascending=False)
    df_sims.to_csv(output_path, index=False)


def main():
    options = parse_args(sys.argv)
    write_artefact_mutational_signatures(options.sigs, options.cord_blood, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
