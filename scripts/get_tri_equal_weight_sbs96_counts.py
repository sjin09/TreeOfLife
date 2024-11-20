#!/usr/bin/env python

import argparse
from pathlib import Path
from typing import Dict
import sys

NTS = ["A", "C", "G", "T"]
PYR = ["C", "T"]
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SBS96_CLASSIFICATION = [f"{nti}[{sub}]{ntj}" for sub in SUBS for nti in NTS for ntj in NTS]
SBS96_TRINUCLEOTIDES = [f"{ntj}{nti}{ntk}" for nti in PYR for ntj in NTS for ntk in NTS]
SBS96_TRINUCLEOTIDE_COUNT = len(SBS96_TRINUCLEOTIDES)
SBS96_TRINUCLEOTIDE_WEIGHT = 1/float(SBS96_TRINUCLEOTIDE_COUNT)


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Get equally weighted SBS96 counts",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="file to read SBS96 counts"
    )
    parser.add_argument(
        "--tri",
        type=Path,
        required=True,
        help="file to read trinucleotide counts"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="file to write equally weighted SBS96 counts"
    )
    args = args[1:]
    return parser.parse_args(args)


def load_trinucleotide_counts(tri_path: Path) -> Dict[str, int]:
    count_per_tri = {
        line.rstrip().split()[0]: int(line.rstrip().split()[1])
        for line in open(tri_path).readlines()
        if not line.startswith("TRI")
    }
    return count_per_tri


def load_sbs96_counts(input_path: Path) -> Dict[str, int]:
    count_per_sbs96 = {}
    for line in open(input_path).readlines():
        if line.startswith("SUBSTITUTION"):
            continue
        fields = line.rstrip().split()
        count_per_sbs96[fields[2]] = int(fields[3])
    return count_per_sbs96


def write_tri_equal_weight_sbs96_counts(input_path: Path, tri_path: Path, output_path: Path):
    count_per_tri = load_trinucleotide_counts(tri_path)
    count_per_sbs96 = load_sbs96_counts(input_path)
    total_tri_count = sum(count_per_tri.values())
    proportion_per_tri = {
        tri: count / total_tri_count
        for tri, count in count_per_tri.items()
    }
    with open(output_path, "w") as outfile:
        outfile.write("SUBSTITUTION\ttri\tSBS96\tCOUNT\n")
        for sbs96 in SBS96_CLASSIFICATION:
            sub = f"{sbs96[2]}>{sbs96[4]}"
            tri = f"{sbs96[0]}{sbs96[2]}{sbs96[6]}"
            count = count_per_sbs96[sbs96]
            weight = 1/(proportion_per_tri[tri]/SBS96_TRINUCLEOTIDE_WEIGHT)
            weighted_count = count * weight
            outfile.write("{}\t{}\t{}\t{}\n".format(sub, tri, sbs96, weighted_count))


def main() -> int:
    options = parse_args(sys.argv)
    write_tri_equal_weight_sbs96_counts(options.input, options.tri, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
