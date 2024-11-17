#!/usr/bin/env python

import argparse
from collections import defaultdict
from pathlib import Path
import sys

import natsort
import pandas as pd
import plotnine as p9


NTS = ["A", "C", "G", "T"]
SBS52_MUTSIG_FILL_COLOURS = ("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#F5ABCC")
SBS52_SUBS = ["C>A", "C>G", "C>T", "T>A", "T>G"]
SBS96_SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SBS96_CLASSIFICATION = [f"{nti}[{sub}]{ntj}" for sub in SBS96_SUBS for nti in NTS for ntj in NTS]
SBS96_TO_SBS52 = {
    "A[C>A]A": "T[T>G]T",
    "T[T>G]T": "T[T>G]T",
    "A[C>A]C": "A[C>A]C",
    "G[T>G]T": "A[C>A]C",
    "A[C>A]G": "C[T>G]T",
    "C[T>G]T": "C[T>G]T",
    "A[C>A]T": "A[C>A]T",
    "A[T>G]T": "A[C>A]T",
    "C[C>A]A": "C[C>A]A",
    "T[T>G]G": "C[C>A]A",
    "C[C>A]C": "C[C>A]C",
    "G[T>G]G": "C[C>A]C",
    "C[C>A]G": "C[C>A]G",
    "C[T>G]G": "C[C>A]G",
    "C[C>A]T": "C[C>A]T",
    "A[T>G]G": "C[C>A]T",
    "G[C>A]A": "T[T>G]C",
    "T[T>G]C": "T[T>G]C",
    "G[C>A]C": "G[C>A]C",
    "G[T>G]C": "G[C>A]C",
    "G[C>A]G": "C[T>G]C",
    "C[T>G]C": "C[T>G]C",
    "G[C>A]T": "A[T>G]C",
    "A[T>G]C": "A[T>G]C",
    "T[C>A]A": "T[C>A]A",
    "T[T>G]A": "T[C>A]A",
    "T[C>A]C": "T[C>A]C",
    "G[T>G]A": "T[C>A]C",
    "T[C>A]G": "C[T>G]A",
    "C[T>G]A": "C[T>G]A",
    "T[C>A]T": "T[C>A]T",
    "A[T>G]A": "T[C>A]T",
    "A[C>T]A": "A[C>T]A",
    "A[T>C]A": "A[C>T]A",
    "A[C>T]C": "A[C>T]C",
    "A[T>C]C": "A[C>T]C",
    "A[C>T]G": "A[C>T]G",
    "A[T>C]G": "A[C>T]G",
    "A[C>T]T": "A[C>T]T",
    "A[T>C]T": "A[C>T]T",
    "C[C>T]A": "C[C>T]A",
    "C[T>C]A": "C[C>T]A",
    "C[C>T]C": "C[C>T]C",
    "C[T>C]C": "C[C>T]C",
    "C[C>T]G": "C[C>T]G",
    "C[T>C]G": "C[C>T]G",
    "C[C>T]T": "C[C>T]T",
    "C[T>C]T": "C[C>T]T",
    "G[C>T]A": "G[C>T]A",
    "G[T>C]A": "G[C>T]A",
    "G[C>T]C": "G[C>T]C",
    "G[T>C]C": "G[C>T]C",
    "G[C>T]G": "G[C>T]G",
    "G[T>C]G": "G[C>T]G",
    "G[C>T]T": "G[C>T]T",
    "G[T>C]T": "G[C>T]T",
    "T[C>T]A": "T[C>T]A",
    "T[T>C]A": "T[C>T]A",
    "T[C>T]C": "T[C>T]C",
    "T[T>C]C": "T[C>T]C",
    "T[C>T]G": "T[C>T]G",
    "T[T>C]G": "T[C>T]G",
    "T[C>T]T": "T[C>T]T",
    "T[T>C]T": "T[C>T]T",
    "A[C>G]A": "T[C>G]T",
    "T[C>G]T": "T[C>G]T",
    "A[C>G]C": "A[C>G]C",
    "G[C>G]T": "A[C>G]C",
    "A[C>G]G": "C[C>G]T",
    "C[C>G]T": "C[C>G]T",
    "A[C>G]T": "A[C>G]T",
    "C[C>G]A": "C[C>G]A",
    "T[C>G]G": "C[C>G]A",
    "C[C>G]C": "C[C>G]C",
    "G[C>G]G": "C[C>G]C",
    "C[C>G]G": "C[C>G]G",
    "G[C>G]A": "T[C>G]C",
    "T[C>G]C": "T[C>G]C",
    "G[C>G]C": "G[C>G]C",
    "T[C>G]A": "T[C>G]A",
    "A[T>A]A": "T[T>A]T",
    "T[T>A]T": "T[T>A]T",
    "A[T>A]C": "A[T>A]C",
    "G[T>A]T": "A[T>A]C",
    "A[T>A]G": "C[T>A]T",
    "C[T>A]T": "C[T>A]T",
    "A[T>A]T": "A[T>A]T",
    "C[T>A]A": "C[T>A]A",
    "T[T>A]G": "C[T>A]A",
    "C[T>A]C": "C[T>A]C",
    "G[T>A]G": "C[T>A]C",
    "C[T>A]G": "C[T>A]G",
    "G[T>A]A": "T[T>A]C",
    "T[T>A]C": "T[T>A]C",
    "G[T>A]C": "G[T>A]C",
    "T[T>A]A": "T[T>A]A",
}

SBS52_TO_SBS96_CLASSIFICATIONS = defaultdict(list)
for sbs96, sbs52 in SBS96_TO_SBS52.items():
    SBS52_TO_SBS96_CLASSIFICATIONS[sbs52].append(sbs96)
SBS52_TO_SBS96_CLASSIFICATIONS = {k: natsort.natsorted(set(v)) for k, v in SBS52_TO_SBS96_CLASSIFICATIONS.items()}

SBS52_CLASSIFICATIONS_PER_SUB = defaultdict(set)
for sbs52 in SBS96_TO_SBS52.values():
    ref, alt = sbs52[2], sbs52[4]
    SBS52_CLASSIFICATIONS_PER_SUB[f"{ref}>{alt}"].add(sbs52)
SBS52_CLASSIFICATIONS_PER_SUB = {k: natsort.natsorted(list(v)) for k, v in SBS52_CLASSIFICATIONS_PER_SUB.items()}
SBS52_CLASSIFICATION = [sbs52 for sub in SBS52_SUBS for sbs52 in SBS52_CLASSIFICATIONS_PER_SUB[sub]]

SBS52_BARPLOT_XAXIS = [
    f"({';'.join(SBS52_TO_SBS96_CLASSIFICATIONS[sbs52])}) {sbs52[0]}{sbs52[2]}{sbs52[6]}"
    for sub in SBS52_SUBS
    for sbs52 in SBS52_CLASSIFICATIONS_PER_SUB[sub]
]


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Get SBS52 barplot",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="file to read"
    )
    parser.add_argument(
        "--sample",
        type=str,
        required=True,
        help="sample id"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="file to draw"
    )
    args = args[1:]
    return parser.parse_args(args)


def draw_sbs52_barplot(input_path: Path, sample: str, output_path: Path) -> None:
    df = pd.read_csv(input_path, sep="\t")
    df["TRINUCLEOTIDE"] = pd.Categorical(
        df["TRINUCLEOTIDE"],
        categories=SBS52_BARPLOT_XAXIS,
        ordered=True
    )
    plot = (
        p9.ggplot(df, p9.aes(x="TRINUCLEOTIDE", y="COUNT", fill="SUBSTITUTION"))
        + p9.geom_bar(stat="identity")
        + p9.theme_bw(16)
        + p9.facet_grid(". ~ SUBSTITUTION", scales="free")
        + p9.scale_fill_manual(
            values=SBS52_MUTSIG_FILL_COLOURS
        )
        + p9.labs(x="\nTrinucleotide\n", y="\nCounts\n")
        + p9.ggtitle("\n{}\n".format(sample))
        + p9.theme(
            text=p9.element_text(size=10),
            legend_title=p9.element_blank(),
            axis_text_x=p9.element_text(
                family="monospace", angle=90, ha="center"
            ),
        )
    )
    plot.save(output_path, width=22, height=12)


def main():
    options = parse_args(sys.argv)
    draw_sbs52_barplot(options.input, options.sample, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(0)
