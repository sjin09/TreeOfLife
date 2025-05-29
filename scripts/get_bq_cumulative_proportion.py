#!/usr/bin/env python

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import plotnine as p9
import pysam

MIN_BQ = 1
MAX_BQ = 93
BQ_THRESHOLDS = np.arange(MIN_BQ, MAX_BQ + 1)


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Get base quality score cumulative proportion", formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-i", "--bam", type=Path, required=True, help="BAM file to read")
    parser.add_argument("--pdf", type=Path, required=True, help="file to plot")
    parser.add_argument("-o", "--output", type=Path, required=True, help="file to write")
    args = args[1:]
    return parser.parse_args(args)


def calculate_cumulative_bq_distribution(bam_path: Path, output_path: Path):
    # Read through alignments and collect base qualities
    counter = 0
    global_base_quality_counts = np.zeros(94, dtype=float)
    alignments = pysam.AlignmentFile(bam_path, "rb")
    for alignment in alignments.fetch(until_eof=True):
        if alignment.is_secondary or alignment.is_supplementary:
            continue
        base_qualities = alignment.query_qualities
        alignment_base_quality_counts = np.bincount(base_qualities, minlength=94)
        global_base_quality_counts += alignment_base_quality_counts
        counter += 1
        if counter > 999:
            break
    alignments.close()

    # Calculate total number of bases
    total_bases = global_base_quality_counts.sum()

    # Calculate cumulative counts and proportions
    cumulative_counts = np.cumsum(global_base_quality_counts[::-1])[::-1]  # Reverse cumsum
    cumulative_proportions = np.round((cumulative_counts / total_bases) * 100, 1)

    # Create DataFrame
    df = pd.DataFrame({
        'BQ': np.arange(94),
        'Cumulative_counts': cumulative_counts,
        'Cumulative_proportion': cumulative_proportions
    })
    df.to_csv(output_path, sep='\t', index=False)
    return df


def plot_cumulative_bq_proportion(df: pd.DataFrame, pdf_path: Path):
    # Create the plot
    plot = (
        p9.ggplot(df, p9.aes(x='BQ', y='Cumulative_proportion'))
        + p9.geom_line(size=1)
        + p9.labs(
            x='\nBase Quality (BQ) Score\n',
            y='\nCumulative Proportion of Bases ≥ BQ Threshold\n'
        )
        + p9.theme_bw(32)
        + p9.scale_x_continuous(breaks=range(0, 91, 10))
        + p9.scale_y_continuous(limits=(0, None))
    )

    # Save the plot
    plot.save(filename=pdf_path, width=24, height=24)


def main() -> int:
    options = parse_args(sys.argv)
    df = calculate_cumulative_bq_distribution(options.bam, options.output)
    plot_cumulative_bq_proportion(df, options.pdf)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
