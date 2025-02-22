#!/usr/bin/env python

"""
Given a VCF file with either germline or somatmic mutations and a reference FASTA file, write a table with SBS52 counts.
"""

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, List
import sys

import natsort
import pysam


NTS = ["A", "C", "G", "T"]
PUR = set(["A", "G"])
PUR_TO_PYR_LOOKUP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
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


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Get SBS52 counts from a VCF file",
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
        "--is-sample-reference-sample",
        required=False,
        action="store_true",
        help="reads and the reference FASTA file is derived from the same sample"
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


def is_biallelic_snp(record: pysam.VariantRecord) -> bool:
    if len(record.alts) == 1:
        if len(record.ref) == 1 and len(record.alts[0]) == 1:
            return True
    return False


def is_pass(record: pysam.VariantRecord) -> bool:
    return True if list(record.filter)[0] == "PASS" else False


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


def load_sbs96_counts(
    vcf_path: Path,
    ref_fasta_path: Path,
    target_path: Path,
    is_sample_reference_sample: bool
) -> Dict[str, int]:
    count_per_sbs96 = defaultdict(lambda: 0)
    variant_records = pysam.VariantFile(vcf_path)
    sample = variant_records.header.samples[0]
    sequence_lookup = pysam.FastaFile(ref_fasta_path)
    if target_path:
        chroms = load_chroms(target_path)
    else:
        chroms = sequence_lookup.references
    if is_sample_reference_sample:
        for chrom in chroms:
            for record in variant_records.fetch(chrom):
                if is_pass(record) and is_biallelic_snp(record):
                    sample_gt = record.samples[sample]["GT"]
                    if sample_gt != (0, 1):
                        continue
                    sbs96 = get_sbs96(record, sequence_lookup)
                    count_per_sbs96[sbs96] += 1
    else:
        for record in variant_records:
            if is_pass(record) and is_biallelic_snp(record):
                sbs96 = get_sbs96(record, sequence_lookup)
                count_per_sbs96[sbs96] += 1
    variant_records.close()
    return dict(count_per_sbs96)


def load_sbs52_counts(
    vcf_path: Path,
    ref_fasta_path: Path,
    target_path: Path,
    is_sample_reference_sample: bool
) -> Dict[str, int]:
    count_per_sbs52 = defaultdict(lambda: 0)
    count_per_sbs96 = load_sbs96_counts(vcf_path, ref_fasta_path, target_path, is_sample_reference_sample)
    for sbs96 in SBS96_CLASSIFICATION:
        count_per_sbs52[SBS96_TO_SBS52[sbs96]] += count_per_sbs96[sbs96]
    return dict(count_per_sbs52)


def write_sbs52_counts(count_per_sbs96: Dict[str, int], out_path: Path):
    with open(out_path, "w") as outfile:
        outfile.write("SUBSTITUTION\tTRINUCLEOTIDE\tSBS52\tCOUNT\n")
        for sbs52 in SBS52_CLASSIFICATION:
            ubase, _, ref, _, alt, _, dbase = list(sbs52)
            sub = f"{ref}>{alt}"
            tri = f"{ubase}{ref}{dbase}"
            tri = "({}) {}".format(";".join(SBS52_TO_SBS96_CLASSIFICATIONS[sbs52]), tri)
            outfile.write("{}\t{}\t{}\t{}\n".format(sub, tri, sbs52, count_per_sbs96[sbs52]))


def get_sbs52_counts(
    vcf_path: Path,
    ref_fasta_path: Path,
    target_path: Path,
    is_sample_reference_sample: bool,
    out_path: Path
):
    count_per_sbs52 = load_sbs52_counts(vcf_path, ref_fasta_path, target_path, is_sample_reference_sample)
    write_sbs52_counts(count_per_sbs52, out_path)


def main() -> 0:
    options = parse_args(sys.argv)
    get_sbs52_counts(options.vcf, options.ref_fasta, options.target, options.is_sample_reference_sample, options.out)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(0)
