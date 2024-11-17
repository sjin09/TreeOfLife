#!/usr/bin/env python

import argparse
from pathlib import Path
import sys

import pandas as pd
import plotnine as p9

PUR = set(["A", "G"])
PUR_TO_PYR_LOOKUP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

NTS = ["A", "C", "G", "T"]
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
MUTSIG_FILL_COLOURS = ("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#83A603", "#F5ABCC")
SBS96_CLASSIFICATION = [f"{nti}[{sub}]{ntj}" for sub in SUBS for nti in NTS for ntj in NTS]


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Get SBS96 barplot",
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
        required=False,
        help="file to draw"
    )
    args = args[1:]
    return parser.parse_args(args)


def get_sbs96_barplot(input_path: Path, sample: str, output_path: Path):
    df = pd.read_csv(input_path, sep="\t")
    plot = (
        p9.ggplot(df, p9.aes(x="TRINUCLEOTIDE", y="COUNT", fill="SUBSTITUTION"))
        + p9.geom_bar(stat="identity")
        + p9.theme_bw(16)
        + p9.facet_grid(". ~ SUBSTITUTION", scales="free")
        + p9.scale_fill_manual(
            values=MUTSIG_FILL_COLOURS
        )
        + p9.labs(x="\nTrinucleotide\n", y="\nCounts\n")
        + p9.ggtitle("\n{}\n".format(sample))
        + p9.theme(
            legend_title=p9.element_blank(),
            axis_text_x=p9.element_text(
                family="monospace", angle=90, ha="center"
            ),
        )
    )
    plot.save(output_path, width=22, height=12)


def main():
    options = parse_args(sys.argv)
    get_sbs96_barplot(options.input, options.sample, options.output)
    return 0


if __name__ == "__main__":
    main()
    sys.exit(0)
