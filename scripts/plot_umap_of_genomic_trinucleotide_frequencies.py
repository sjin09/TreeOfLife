#!/usr/bin/env python

import argparse
from pathlib import Path
import sys
from typing import Dict, List, Tuple

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import umap
from sklearn.preprocessing import StandardScaler

KINGDOMS = ["Viridiplantae", "Metazoa", "Fungi"]
TAXNOMIC_RANKS_OF_INTEREST = [
    "Fungi",
    "Viridiplantae",
    "Aves",
    "Actinopteri",
    "Mammalia",
    "Mollusca",
    "Lepidoptera",
    "Coleoptera",
    "Diptera",
    "Hymenoptera",
]


def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "-i",
        "--fofn",
        type=Path,
        required=True,
        help="file of file names (FOFN) containing paths to trinucleotide frequency files",
    )
    parser.add_argument(
        "--taxonomic-classification", type=Path, required=False, help="file to read taxonomic classification per sample"
    )
    parser.add_argument("-o", "--output", type=Path, required=True, help="file to write")
    args = args[1:]
    return parser.parse_args(args)


def assign_umap_label(df: pd.DataFrame) -> pd.DataFrame:
    labels = []
    for _, row in df.iterrows():
        kingdom = row.get("Kingdom")
        phylum = row.get("Phylum")
        _class = row.get("Class")
        order = row.get("Order")

        if kingdom == "Viridiplantae":
            labels.append("Viridiplantae")
        elif kingdom == "Fungi":
            labels.append("Fungi")
        else:
            if any(rank in TAXNOMIC_RANKS_OF_INTEREST for rank in (phylum, _class, order)):
                if phylum == "Arthropoda":
                    labels.append(f"{order} ({_class})")
                elif phylum == "Chordata":
                    labels.append(f"{_class} ({phylum})")
                else:
                    labels.append(phylum)
            else:
                labels.append(".")
    df["Label"] = labels


def generate_label_color_lookup(df: pd.DataFrame) -> Dict[str, Tuple[float, float, float]]:
    metazoa_color_lookup = {
       'Aves (Chordata)': "#031CA6",
       'Actinopteri (Chordata)': "#0468BF",
       'Lepidoptera (Insecta)': "#A63603",
       'Coleoptera (Insecta)': "#E6550D",
       'Diptera (Insecta)': "#FBAE17",
       'Hymenoptera (Insecta)': "#FFD92F",
       'Mammalia (Chordata)': "#04B2D9",
       'Mollusca': "#C04BF2"
    }

    # Group labels by kingdom
    kingdom_groups = {}
    for kingdom in KINGDOMS:
        kingdom_groups[kingdom] = df.loc[df["Kingdom"] == kingdom, "Label"].unique().tolist()

    # Assign colors within each kingdom group
    label_color_lookup = {}
    for (kingdom, kingdom_labels) in kingdom_groups.items():
        if kingdom == "Viridiplantae":
            label_color_lookup["Viridiplantae"] = "#115923"
        elif kingdom == "Fungi":
            label_color_lookup["Fungi"] = "#7F2704"
        else:  # Metazoa
            for lbl in kingdom_labels:
                label_color_lookup[lbl] = metazoa_color_lookup[lbl]
    return label_color_lookup


def load_taxonomic_classification_table(taxonomic_classification_path: Path) -> Tuple[pd.DataFrame, Dict[str, str]]:
    # Load data frame
    df = pd.read_csv(taxonomic_classification_path, header=0)
    df = df.drop(columns=["Sample", "Common name"])
    df = df.drop_duplicates(keep="first")

    # Subset data frame
    df_subset = df[~(df.iloc[:, 3:10].eq(".").any(axis=1))]
    # Assign UMAP labels
    assign_umap_label(df_subset)

    # Remove samples with unassigned label
    df_subset = df_subset[df_subset["Label"] != "."]

    # Generate sample to label lookup
    sample_to_label_lookup = dict(zip(df_subset["ref_sample"], df_subset["Label"]))
    return df_subset, sample_to_label_lookup


def load_genomic_trinucleotide_frequencies(fofn_path: Path, sample_to_label_lookup: Dict[str, str]) -> pd.DataFrame:
    data_dict = {}
    sample_names = []
    for line in open(fofn_path).readlines():
        file_path = Path(line.rstrip())

        # Collect sample names
        sample_name = Path(file_path).stem
        if sample_name not in sample_to_label_lookup:
            continue
        sample_names.append(sample_name)

        # Read the trinucleotide counts
        df = pd.read_csv(file_path, sep="\t", header=None, names=["tri", "count"])

        # Convert counts to frequencies
        total_count = df["count"].sum()
        frequencies = df["count"] / total_count

        # Store as dictionary with trinucleotide as key
        data_dict[sample_name] = dict(zip(df["tri"], frequencies))

    # Convert to DataFrame
    frequency_matrix = pd.DataFrame.from_dict(data_dict, orient="index")

    # Fill any missing trinucleotides with 0 (shouldn't happen with complete data)
    frequency_matrix = frequency_matrix.fillna(0)
    return frequency_matrix, sample_names


def perform_umap_analysis(
    trinucleotide_frequency_matrix: pd.DataFrame,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    n_components: int = 2,
    random_state: int = 42,
    standardize: bool = True,
) -> Tuple[np.ndarray, umap.UMAP]:
    """
    Perform UMAP dimensionality reduction on trinucleotide frequencies

    Args:
        frequency_matrix (pd.DataFrame): Matrix of trinucleotide frequencies
        n_neighbors (int): UMAP parameter for local neighborhood size
        min_dist (float): UMAP parameter for minimum distance between points
        n_components (int): Number of dimensions for output
        random_state (int): Random seed for reproducibility
        standardize (bool): Whether to standardize the data before UMAP

    Returns:
        np.ndarray: UMAP embedding coordinates
        umap.UMAP: Fitted UMAP model
    """
    X = trinucleotide_frequency_matrix.values

    if standardize:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
    else:
        X_scaled = X

    # Initialize and fit UMAP
    umap_model = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        random_state=random_state,
        metric="euclidean",
    )

    embedding = umap_model.fit_transform(X_scaled)
    return embedding, umap_model


def plot_umap(
    embedding: np.ndarray,
    sample_names: List[str],
    sample_to_label_lookup: Dict[str, str],
    label_color_lookup: Dict[str, Tuple[float, float, float]],
    figsize: Tuple[int, int] = (7.4, 7.4),
    point_size: int = 80,
    alpha: float = 0.85,
    output_path: Path = Path("umap_plot.pdf"),
):
    """
    Scatter plot of UMAP colored by the assigned Label (per-sample), using label_color_lookup.
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Colors for each point
    lbl_colors = []
    lbls_for_legend = []
    for sample_name in sample_names:
        lbl = sample_to_label_lookup[sample_name]
        lbls_for_legend.append(lbl)
        lbl_colors.append(label_color_lookup[lbl])

    # Set Helvetica globally
    mpl.rcParams['font.family'] = 'Helvetica'

    # Plot the UMAP embedding as a scatter plot
    ax.scatter(
        embedding[:, 0],  # X-coordinates (UMAP dimension 1)
        embedding[:, 1],  # Y-coordinates (UMAP dimension 2)
        c=lbl_colors,  # point colors (from your label_color_lookup)
        s=point_size,  # size of each point (area, not radius)
        alpha=alpha,  # transparency of the points (0=transparent, 1=opaque)
        edgecolors="black",  # draw a thin black outline around each point
        linewidth=0.5,  # thickness of that outline
    )

    # Add x-axis and y-axis labels
    ax.set_xlabel("UMAP 1", fontsize=10)
    ax.set_ylabel("UMAP 2", fontsize=10)

    from matplotlib.lines import Line2D
    unique_labels = sorted(set(lbls_for_legend))
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=label_color_lookup.get(lbl, (0.6, 0.6, 0.6)),
            markeredgecolor="black",
            markersize=8,
            label=lbl,
        )
        for lbl in unique_labels
        if lbl != "Unknown"
    ]
    if legend_handles:
        ax.legend(handles=legend_handles, bbox_to_anchor=(0.5, -0.1), loc="upper center", fontsize=7, ncol=4)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")


def plot_umap_of_genomic_trinucleotide_frequencies(
    fofn_path: Path, taxonomic_classification_path: Path, output_path: Path
):
    # Load taxonomic classification data frame and a lookup table mapping sample to label
    taxonomic_classification_df, sample_to_label_lookup = load_taxonomic_classification_table(
        taxonomic_classification_path
    )

    # Generate a colour for each label
    label_color_lookup = generate_label_color_lookup(taxonomic_classification_df)

    # Generate a data frame containing trinucleotide frequencies from samples with labels
    trinucleotide_frequency_df, sample_names = load_genomic_trinucleotide_frequencies(fofn_path, sample_to_label_lookup)

    # Perform dimensionality reduction
    # n_neighbors = [11, 12, 13, 14]
    # min_dists = [0.2, 0.3]
    # for i in n_neighbors:
    #     for j in min_dists:
    #         embedding, _umap_model = perform_umap_analysis(
    #             trinucleotide_frequency_df, n_neighbors=i, min_dist=j, standardize=True
    #         )
    #         # Plot UMAP visualisation
    #         plot_umap(
    #             embedding=embedding,
    #             sample_names=sample_names,
    #             sample_to_label_lookup=sample_to_label_lookup,
    #             label_color_lookup=label_color_lookup,
    #             alpha=1,
    #             point_size=20,
    #             output_path=f"{output_path.stem}_n_neighbors_{i}_min_dist_{j}.pdf",
    #         )

    # Perform dimensionality reduction
    embedding, _umap_model = perform_umap_analysis(
        trinucleotide_frequency_df, n_neighbors=13, min_dist=0.2, standardize=True
    )

    # Plot UMAP visualisation
    plot_umap(
        embedding=embedding,
        sample_names=sample_names,
        sample_to_label_lookup=sample_to_label_lookup,
        label_color_lookup=label_color_lookup,
        alpha=0.95,
        point_size=25,
        output_path=output_path,
    )


def main() -> int:
    options = parse_args(sys.argv)
    plot_umap_of_genomic_trinucleotide_frequencies(options.fofn, options.taxonomic_classification, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
