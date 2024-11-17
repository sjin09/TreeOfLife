#!/usr/bin/env python

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, List
import sys

import pysam

PUR = set(["A", "G"])
PUR_TO_PYR_LOOKUP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

NTS = ["A", "C", "G", "T"]
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
MUTSIG_FILL_COLOURS = ("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#83A603", "#F5ABCC")
SBS96_CLASSIFICATION = [f"{nti}[{sub}]{ntj}" for sub in SUBS for nti in NTS for ntj in NTS]


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Get SBS96 counts from a VCF file",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--vcf",
        type=Path,
        required=True,
        help="VCF file to read"
    )
    parser.add_argument(
        "--ref-fasta",
        type=str,
        required=True,
        help="reference FASTA file to read"
    )
    parser.add_argument(
        "--target",
        type=Path,
        required=False,
        help="target chromosome per line"
    )
    parser.add_argument(
        "-o",
        "--out",
        type=Path,
        required=True,
        help="file to write"
    )
    args = args[1:]
    return parser.parse_args(args)


def load_chroms(target_path: Path) -> List[str]:
    chroms = [line.rstrip() for line in open(target_path).readlines()]
    return chroms


def get_sbs96(record: pysam.VariantRecord, reference_sequence_lookup: pysam.FastaFile) -> str:
    trinucleotide = reference_sequence_lookup.fetch(record.chrom, record.pos - 2, record.pos + 1)
    if record.ref in PUR:
        ubase, _, dbase = trinucleotide[::-1]
        sbs96 = "{}[{}>{}]{}".format(
            PUR_TO_PYR_LOOKUP.get(ubase, "N"),
            PUR_TO_PYR_LOOKUP.get(record.ref, "N"),
            PUR_TO_PYR_LOOKUP.get(record.alts[0], "N"),
            PUR_TO_PYR_LOOKUP.get(dbase, "N"),
        )
    else:
        ubase, _, dbase = trinucleotide
        sbs96 = "{}[{}>{}]{}".format(ubase, record.ref, record.alts[0], dbase)
    return sbs96


def load_sbs96_counts(vcf_path: Path, ref_fasta_path: Path, target_path: Path) -> Dict[str, int]:
    count_per_sbs96 = defaultdict(lambda: 0)
    variant_records = pysam.VariantFile(vcf_path)
    sequence_lookup = pysam.FastaFile(ref_fasta_path)
    if target_path:
        chroms = load_chroms(target_path)
    else:
        chroms = sequence_lookup.references
    for chrom in chroms:
        for record in variant_records.fetch(chrom):
            if list(record.filter)[0] != "PASS":
                continue
            if len(record.alts) != 1:
                continue
            if len(record.ref) and len(record.alts[0]) != 1:
                continue
            sbs96 = get_sbs96(record, sequence_lookup)
            count_per_sbs96[sbs96] += 1
    return dict(count_per_sbs96)


def write_sbs96_counts(count_per_sbs96: Dict[str, int], out_path: Path):
    with open(out_path, "w") as outfile:
        outfile.write("SUBSTITUTION\tTRINUCLEOTIDE\tSBS96\tCOUNT\n")
        for sbs96 in SBS96_CLASSIFICATION:
            ubase, _, ref, _, alt, _, dbase = list(sbs96)
            sub = f"{ref}>{alt}"
            tri = f"{ubase}{ref}{dbase}"
            outfile.write("{}\t{}\t{}\t{}\n".format(sub, tri, sbs96, count_per_sbs96[sbs96]))


def get_sbs96_counts(vcf_path: Path, ref_fasta_path: Path, target_path: Path, out_path: Path):
    count_per_sbs96 = load_sbs96_counts(vcf_path, ref_fasta_path, target_path)
    write_sbs96_counts(count_per_sbs96, out_path)


def main():
    options = parse_args(sys.argv)
    get_sbs96_counts(options.vcf, options.ref_fasta, options.target, options.out)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
