#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys
from typing import List

from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NTS = ["A", "C", "G", "T"]
PYRS = ["C", "T"]
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
TRIS = [f"{nti}{ntj}{ntk}" for ntj in PYRS for nti in NTS for ntk in NTS]


def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to read cord blood granulocyte match and mismatch count by substitution type and context",
    )
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the file to plot")
    args = args[1:]
    return parser.parse_args(args)


def plot_heatmap(ax: Axes, df: pd.DataFrame, vmin: float, vmax: float):
    # Pivot the data to prepare for the heatmap
    hdf = df.pivot_table(index="substitution", columns="trinucleotide", values="phred_q_score")
    im = ax.imshow(hdf.values, cmap="BuPu", aspect="auto", vmin=vmin, vmax=vmax)

    # Annotate each cell with its value
    for i in range(hdf.shape[0]):  # rows
        for j in range(hdf.shape[1]):  # columns
            value = hdf.values[i, j]
            ax.text(
                j,
                i,
                f"Q{value:.1f}",  # format with 1 decimal place
                ha="center",
                va="center",
                color="black",  # text color
                fontsize=6,
                fontfamily="Helvetica",  # use Helvetica font
            )

    # ticks & labels
    ax.set_xticks(np.arange(hdf.shape[1]))
    ax.set_xticklabels(hdf.columns, rotation=90, ha="center", fontsize=6, fontfamily="Helvetica")
    ax.set_yticks(np.arange(hdf.shape[0]))
    ax.set_yticklabels(hdf.index, fontsize=6, fontfamily="Helvetica")

    # Add white gridlines between cells
    ax.set_xticks(np.arange(hdf.shape[1]) - 0.5, minor=True)
    ax.set_yticks(np.arange(hdf.shape[0]) - 0.5, minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)

    # x-axis labels
    ax.set_xlabel("\nTrinucleotide\n", fontsize=7, fontfamily="Helvetica")
    ax.set_ylabel("\nSubstitution\n", fontsize=7, fontfamily="Helvetica")
    return im


def plot_q_score_heatmap(input_path: Path, output_path: Path):
    # Import data frame
    df = pd.read_csv(input_path, header=0, sep="\t")

    # Create trinucleotide column using apply to combine 5_base, substitution (taking first base), and 3_base
    df["trinucleotide"] = df.apply(lambda row: row["5_base"] + row["substitution"][0] + row["3_base"], axis=1)

    # Sort by natural order of substitutions
    df["substitution"] = pd.Categorical(df["substitution"], categories=SUBS, ordered=True)

    # Subset data frame
    df_c_subs = df[df["substitution"].str.startswith(("C"))]
    df_t_subs = df[df["substitution"].str.startswith(("T"))]

    # Set minimum and maximum values
    vmin = min(df_c_subs["phred_q_score"].min(), df_t_subs["phred_q_score"].min())
    vmax = max(df_c_subs["phred_q_score"].max(), df_t_subs["phred_q_score"].max())

    # Create a heatmap using matplotlib # Set figure size
    fig, (ax_c, ax_t) = plt.subplots(2, 1, figsize=(7.09, 7.09))
    fig.patch.set_facecolor("white")  # Set the figure background color

    # Plot heatmap
    im_c = plot_heatmap(ax_c, df_c_subs, vmin, vmax)
    im_t = plot_heatmap(ax_t, df_t_subs, vmin, vmax)
    
    # Create a colorbar for the heatmap.
    cbar_ax = fig.add_axes([0.20, -0.01, 0.7, 0.03]) # left margin, bottom margin, width of colorbar, height of color bar
    cbar = fig.colorbar(im_t, ax=[ax_c, ax_t], cax=cbar_ax, orientation="horizontal")
    cbar.set_label("Empirical Phred Q score", fontsize=7, fontfamily="Helvetica")
    cbar.ax.tick_params(labelsize=6)
    
    # Explort plot
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", dpi=300)  # uncomment if you want to save


def main():
    options = parse_args(sys.argv)
    plot_q_score_heatmap(options.input, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
