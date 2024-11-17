#!/usr/bin/env python

import argparse
from collections import defaultdict
import multiprocessing as mp
from pathlib import Path
import sys
from typing import Dict, List

import pysam


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="FASTA file to read"
    )
    parser.add_argument(
        "--target",
        type=Path,
        required=False,
        help="target chromosomes per line"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        required=False,
        help="number of threads"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="file to write"
    )
    args = args[1:]
    return parser.parse_args(args)


PUR = set(["A", "G"])
PYR = ["C", "T"]
NTS = ["A", "C", "G", "T"]
TRINUCLEOTIDES = [f"{ntj}{nti}{ntk}" for nti in PYR for ntj in NTS for ntk in NTS]
PURINE_TO_PYRIMIDINE_LOOKUP_TABLE = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def load_chroms(target_path: Path) -> List[str]:
    chroms = [line.strip() for line in open(target_path)]
    return chroms


def get_trinucleotide_count_per_chrom(
    chrom: str,
    sequence: str,
    trinucleotide_count_per_chrom: Dict[str, Dict[str, int]]
):
    count_per_trinucleotide = defaultdict(lambda: 0)
    for i in range(len(sequence)-2):
        base = sequence[i]
        if base == "N":
            continue
        trinucleotide = sequence[i:i+3]
        if trinucleotide[1] in PUR:
            trinucleotide_pyr = "".join(
                [PURINE_TO_PYRIMIDINE_LOOKUP_TABLE.get(base, "N") for base in trinucleotide[::-1]]
            )
            count_per_trinucleotide[trinucleotide_pyr] += 1
        else:
            count_per_trinucleotide[trinucleotide] += 1
    trinucleotide_count_per_chrom[chrom] = dict(count_per_trinucleotide)


def get_trinucleotide_count(ref_fasta_path: Path, target_path, thread_count: int, output_path: Path):
    p = mp.Pool(thread_count)
    manager = mp.Manager()
    sequence_lookup = pysam.FastaFile(ref_fasta_path)
    if target_path:
        chroms = load_chroms(target_path)
    else:
        chroms = sequence_lookup.references
    trinucleotide_count_per_chrom = manager.dict()
    get_trinucleotide_count_per_chrom_arguments = [
        (
            chrom,
            sequence_lookup[chrom],
            trinucleotide_count_per_chrom
        )
        for chrom in chroms
    ]
    p.starmap(get_trinucleotide_count_per_chrom, get_trinucleotide_count_per_chrom_arguments)
    p.close()
    p.join()
    # write trinucleotide count
    write_trinucleotide_count(dict(trinucleotide_count_per_chrom), output_path)


def write_trinucleotide_count(trinucleotide_count_per_chrom: Dict[str, Dict[str, int]], output_path: Path) -> None:
    trinucleotide_count = defaultdict(lambda: 0)
    for chrom in trinucleotide_count_per_chrom:
        for tri, count in trinucleotide_count_per_chrom[chrom].items():
            trinucleotide_count[tri] += count

    with open(output_path, "w") as outfile:
        for tri in TRINUCLEOTIDES:
            outfile.write("{}\t{}\n".format(tri, trinucleotide_count[tri]))


def main() -> int:
    options = parse_args(sys.argv)
    get_trinucleotide_count(options.input, options.target, options.threads, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
