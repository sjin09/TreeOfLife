#!/usr/bin/env python

import argparse
import csv
from pathlib import Path
import sys

import pysam


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Get genomic base quality score distribution", formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-i", "--bam", type=str, required=True, help="BAM file to read")
    parser.add_argument("-o", "--output", type=str, required=True, help="file to write")
    args = args[1:]
    return parser.parse_args(args)


def write_genomic_base_quality_score_distribution(bam_path: Path, output_path: Path):

    alignments = pysam.AlignmentFile(bam_path, "rb")
    chromosomes = [str(i) for i in range(1, 23)] + ["X", "Y"]
    fieldnames = ["chrom", "pos", "bq93_base_count", "base_count", "q93_base_proportion"]
    with open(output_path, "w") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for chrom in chromosomes:
            for pileupcolumn in alignments.pileup(chrom):
                bq93_base_count = 0
                total_base_count = 0
                for pileupread in pileupcolumn.pileups:
                    if pileupread.is_refskip:
                        break
                    if pileupread.is_del:
                        continue
                    bq = pileupread.alignment.query_qualities[pileupread.query_position]
                    if bq >= 93:
                        bq93_base_count += 1
                    total_base_count += 1
                writer.writerow(
                    {
                        "chrom": chrom,
                        "pos": pileupcolumn.pos,
                        "bq93_base_count": bq93_base_count,
                        "base_count": total_base_count,
                        "q93_base_proportion": bq93_base_count / total_base_count if total_base_count > 0 else 0.0,
                    }
                )
    alignments.close()


def main() -> int:
    options = parse_args(sys.argv)
    write_genomic_base_quality_score_distribution(options.bam, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
