#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity


NTS = ["A", "C", "G", "T"]
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SBS96_CLASSIFICATIONS = [f"{nti}[{sub}]{ntj}" for sub in SUBS for nti in NTS for ntj in NTS]


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--sig1",
        type=Path,
        required=True,
        help="Path to the first HDP somatic mutational signature file",
    )
    parser.add_argument(
        "--sig2",
        type=Path,
        required=True,
        help="Path to the second HDP somatic mutational signature file",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Path to the file to plot"
    )
    args = args[1:]
    return parser.parse_args(args)


def calculate_cosine_similarity_between_signatures(sig1_path: Path, sig2_path: Path) -> pd.DataFrame:
    # Load the signatures
    sig1 = pd.read_csv(sig1_path, sep=",")
    sig2 = pd.read_csv(sig2_path, sep=",")

    # calculate cosine similarity
    sig_sim_mtx = cosine_similarity(sig1.values.T, sig2.values.T)
    return sig_sim_mtx


def calculate_cosine_similarity_between_exposures(exp1_path: Path, exp2_path: Path) -> pd.DataFrame:
    # Load the exposures
    exp1 = pd.read_csv(exp1_path, sep=",", index_col=0)
    exp2 = pd.read_csv(exp2_path, sep=",", index_col=0)

    # Get samples
    exp1_samples = exp1.index
    exp2_samples = exp2.index

    # Get the intersection of samples
    common_samples = exp1_samples.intersection(exp2_samples)

    # Filtered and ordered the exposures to only include common samples
    exp1 = exp1.loc[common_samples]
    exp2 = exp2.loc[common_samples]

    # calculate cosine similarity
    exp_sim_mtx = cosine_similarity(exp1.values.T, exp2.values.T)
    return exp_sim_mtx


def plot_sim_heatmap(sig_sim_mtx: pd.DataFrame, output_path: Path):

    # Set the row and column names
    nrows, ncols = sig_sim_mtx.shape
    row_names = [f"HDP{i}" for i in range(1, nrows+1)]
    column_names = [f"HDP{i}" for i in range(1, ncols+1)]

    # Create a heatmap using matplotlib
    fig, ax = plt.subplots(figsize=(24, 20))
    fig.patch.set_facecolor('white')  # Set the figure background color
    ax.set_facecolor("white")  # Set the axes background color

    # remove the default spines
    for spine in ax.spines.values():
        spine.set_visible(False)

    # draw the heatmap
    im = ax.imshow(sig_sim_mtx, cmap="Blues", vmin=0.4, vmax=1.0, aspect='auto')

    # —— Add white gridlines between cells ——
    ax.set_xticks(np.arange(-.5, len(column_names), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(row_names), 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=2)
    ax.tick_params(which='minor', length=0)

    # set x-ticks & labels
    ax.set_xticks(np.arange(len(column_names)))
    ax.set_xticklabels(column_names, rotation=90, ha='center')

    # set y-ticks & labels
    ax.set_yticks(np.arange(len(row_names)))
    ax.set_yticklabels(row_names)

    # Set the title and labels
    ax.set_xlabel("\nHDP somatic mutational signature set 2\n", fontsize=14, labelpad=10)
    ax.set_ylabel("\nHDP somatic mutational signature set 1\n", fontsize=14, labelpad=10)

    # Add colorbar
    # cbar = fig.colorbar(im, ax=ax, orientation='horizontal')
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.02, pad=0.08, shrink=0.7)
    cbar.set_label('Cosine similarity', fontsize=14)

    if output_path:
        plt.savefig(output_path, bbox_inches='tight', dpi=300)


def main():
    options = parse_args(sys.argv)
    sig_sim_mtx = calculate_cosine_similarity_between_signatures(
        options.sig1,
        options.sig2,
    )
    plot_sim_heatmap(sig_sim_mtx, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
