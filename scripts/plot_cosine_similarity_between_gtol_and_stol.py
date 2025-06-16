#!/usr/bin/env python3

import argparse
from collections import defaultdict
from pathlib import Path
import sys
from typing import Dict, List

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
        "--gtol-sigs",
        type=Path,
        required=True,
        help="Path to the HDP germline mutational signature file",
    )
    parser.add_argument(
        "--stol-sigs",
        type=Path,
        required=True,
        help="Path to the HDP somatic mutational signature file",
    )
    parser.add_argument(
        "--sbs52",
        type=Path,
        required=True,
        help="Path to the SBS52 table"
    )
    parser.add_argument(
        "--sbs96-to-sbs52",
        type=Path,
        required=True,
        help="Path to the SBS96 to SBS52 lookup table"
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


def get_sbs96_to_sbs52_lookup_table(sbs96_to_sbs52_path: Path) -> Dict[str, str]:
    sbs96_to_sbs52_lookup_table = {}
    for line in open(sbs96_to_sbs52_path):
        if line.startswith("SBS96"):
            continue
        fields = line.strip().split("\t")
        sbs96, sbs52 = fields[0], fields[1]
        sbs96_to_sbs52_lookup_table[sbs96] = sbs52
    return sbs96_to_sbs52_lookup_table


def load_stol_sigs(
    stol_sig_path: Path,
    sbs96_to_sbs52_path: Dict[str, str],
    sbs52_classifications: List[str]
) -> pd.DataFrame:
    # Load the signatures
    stol_sbs96_sigs = pd.read_csv(stol_sig_path, sep=",")
    stol_sbs96_sigs.iloc[:, 0] = SBS96_CLASSIFICATIONS
    stol_sbs96_sigs.rename(columns={stol_sbs96_sigs.columns[0]: 'SBS96'}, inplace=True)

    # Load the SBS96 to SBS52 lookup table
    sbs96_to_sbs52_lookup_table = get_sbs96_to_sbs52_lookup_table(sbs96_to_sbs52_path)

    # group by SBS52
    stol_sbs52_sigs_dict = {}
    for col_name in stol_sbs96_sigs.columns[1:]:
        stol_sbs52_sigs_dict[col_name] = defaultdict(lambda: 0.0)
        probs = stol_sbs96_sigs[col_name]
        for sbs96, prob in zip(SBS96_CLASSIFICATIONS, probs):
            sbs52 = sbs96_to_sbs52_lookup_table.get(sbs96)
            stol_sbs52_sigs_dict[col_name][sbs52] += prob

    # Convert the dictionary to a DataFrame
    stol_sbs52_sigs = pd.DataFrame.from_dict(
        stol_sbs52_sigs_dict,
        orient='index'
    )
    stol_sbs52_sigs = stol_sbs52_sigs.reindex(columns=sbs52_classifications)
    return stol_sbs52_sigs.T


def calculate_cosine_similarity_between_signatures(
    gtol_sig_path: Path,
    stol_sig_path: Path,
    sbs52_path: Path,
    sbs96_to_sbs52_path: Path
) -> pd.DataFrame:
    # Load SBS52 classifications
    sbs52_classifications = [line.rstrip() for line in open(sbs52_path)]

    # Load signatures
    gtol_sigs = pd.read_csv(gtol_sig_path, sep=",")
    stol_sigs = load_stol_sigs(stol_sig_path, sbs96_to_sbs52_path, sbs52_classifications)

    # calculate cosine similarity
    sig_sim_mtx = cosine_similarity(stol_sigs.values.T, gtol_sigs.values.T)
    return sig_sim_mtx


def plot_sim_heatmap(sig_sim_mtx: pd.DataFrame, output_path: Path):

    # Set the row and column names
    nrows, ncols = sig_sim_mtx.shape
    row_names = [f"sTOL{i}" for i in range(1, nrows+1)]
    column_names = [f"gTOL{i}" for i in range(1, ncols+1)]

    # Create a heatmap using matplotlib
    fig, ax = plt.subplots(figsize=(24, 20))
    fig.patch.set_facecolor('white')  # Set the figure background color
    ax.set_facecolor("white")  # Set the axes background color

    # remove the default spines
    for spine in ax.spines.values():
        spine.set_visible(False)

    # draw the heatmap
    im = ax.imshow(sig_sim_mtx, cmap="Greens", vmin=0.4, vmax=1.0, aspect='auto')

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
    ax.set_xlabel("\nHDP somatic mutational signature set\n", fontsize=14, labelpad=10)
    ax.set_ylabel("\nHDP germline mutational signature set\n", fontsize=14, labelpad=10)

    # Add colorbar
    # cbar = fig.colorbar(im, ax=ax, orientation='horizontal')
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.02, pad=0.08, shrink=0.7)
    cbar.set_label('Cosine similarity', fontsize=14)

    if output_path:
        plt.savefig(output_path, bbox_inches='tight', dpi=300)


def main():
    options = parse_args(sys.argv)
    sig_sim_mtx = calculate_cosine_similarity_between_signatures(
        gtol_sig_path=options.gtol_sigs,
        stol_sig_path=options.stol_sigs,
        sbs52_path=options.sbs52,
        sbs96_to_sbs52_path=options.sbs96_to_sbs52
    )
    plot_sim_heatmap(sig_sim_mtx, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
