#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity


NTS = ["A", "C", "G", "T"]
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SBS96_CLASSIFICATIONS = [f"{nti}[{sub}]{ntj}" for sub in SUBS for nti in NTS for ntj in NTS]


def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--cosmic",
        type=Path,
        required=True,
        help="Path to the COSMIC somatic mutational signature file",
    )
    parser.add_argument(
        "--stol",
        type=Path,
        required=True,
        help="Path to the tree of life somatic mutational signature file",
    )
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the file to plot")
    args = args[1:]
    return parser.parse_args(args)


def write_table_mapping_hdp_to_stol(hdp_sig_names: List[str], stol_sig_names: List[str]):
    with open("hdp_to_stol_mapping_table.tsv", "w") as outfile:
        for (hdp_sig_name, stol_sig_name) in zip(hdp_sig_names, stol_sig_names):
            outfile.write(f"{hdp_sig_name}\t{stol_sig_name}\n")


def calculate_cosine_similarity_between_signatures(cosmic_path: Path, stol_path: Path) -> pd.DataFrame:
    # Load the signatures
    cosmic_sigs = pd.read_csv(cosmic_path, sep="\t")
    stol_sigs = pd.read_csv(stol_path, sep=",")

    # Reindex COSMIC signatures
    cosmic_sigs = cosmic_sigs.set_index(cosmic_sigs.columns[0])
    cosmic_sigs = cosmic_sigs.reindex(SBS96_CLASSIFICATIONS)
    cosmic_sigs.reset_index(inplace=True)
    cosmic_sigs.rename(columns={cosmic_sigs.columns[0]: "SBS96"}, inplace=True)

    # Set first column as SBS96 classifications
    stol_sigs.iloc[:, 0] = SBS96_CLASSIFICATIONS
    stol_sigs.rename(columns={stol_sigs.columns[0]: "SBS96"}, inplace=True)

    # Drop the first column
    cosmic_sigs = cosmic_sigs.iloc[:, 1:]
    stol_sigs = stol_sigs.iloc[:, 1:]

    # Get column names
    cosmic_sig_names = cosmic_sigs.columns.tolist()
    stol_sig_names = [f"sTOL{op_idx}" for op_idx, sig_name in enumerate(stol_sigs.columns.tolist(), start=1)]

    # Write a table mapping HDP signature names to sTOL signature names
    write_table_mapping_hdp_to_stol(stol_sigs.columns.tolist(), stol_sig_names)

    # calculate cosine similarity
    sig_sim_mtx = cosine_similarity(cosmic_sigs.values.T, stol_sigs.values.T)
    return sig_sim_mtx, cosmic_sig_names, stol_sig_names


def plot_sim_heatmap(
    sig_sim_mtx: pd.DataFrame, cosmic_sig_names: List[str], stol_sig_names: List[str], output_path: Path
):

    # Create a heatmap using matplotlib
    fig, ax = plt.subplots(figsize=(24, 20))
    fig.patch.set_facecolor("white")  # Set the figure background color
    ax.set_facecolor("white")  # Set the axes background color

    # remove the default spines
    for spine in ax.spines.values():
        spine.set_visible(False)

    # draw the heatmap
    im = ax.imshow(sig_sim_mtx, cmap="Greens", vmin=0.4, vmax=1.0, aspect="auto")

    # —— Add white gridlines between cells ——
    ax.set_xticks(np.arange(-0.5, len(stol_sig_names), 1), minor=True)  #  columns
    ax.set_yticks(np.arange(-0.5, len(cosmic_sig_names), 1), minor=True)  # rows
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", length=0)

    # set x-ticks & labels
    ax.set_xticks(np.arange(len(stol_sig_names)))
    ax.set_xticklabels(stol_sig_names, rotation=90, ha="center")

    # set y-ticks & labels
    ax.set_yticks(np.arange(len(cosmic_sig_names)))
    ax.set_yticklabels(cosmic_sig_names)

    # Set the title and labels
    ax.set_xlabel("\nTree of life somatic mutational signatures\n", fontsize=14, labelpad=10)
    ax.set_ylabel("\nCOSMIC somatic mutational signatures\n", fontsize=14, labelpad=10)

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.02, pad=0.08, shrink=0.7)
    cbar.set_label("Cosine similarity", fontsize=14)

    if output_path:
        plt.savefig(output_path, bbox_inches="tight", dpi=300)


def main():
    options = parse_args(sys.argv)
    sig_sim_mtx, cosmic_sig_names, stol_sig_names = calculate_cosine_similarity_between_signatures(
        options.cosmic,
        options.stol,
    )
    plot_sim_heatmap(sig_sim_mtx, cosmic_sig_names, stol_sig_names, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
