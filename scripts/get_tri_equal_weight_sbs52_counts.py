#!/usr/bin/env python

import argparse
import logging
import multiprocessing as mp
import sys
## from collections import defaultdict
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List

import natsort
import pandas as pd
import plotnine as p9
import psutil
import pysam


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


PUR_SET = set(["A", "G"])
NTS = ["A", "C", "G", "T"]
SBS52_SUB_LST = ["C>A", "C>G", "C>T", "T>A", "T>G"]
SBS96_SUB_LST = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
COMPLEMENTARY_BASE_LOOKUP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

SBS96_LST = []
for sbs96_sub in SBS96_SUB_LST:
    for nti in NTS:
        for ntj in NTS:
            sbs96 = f"{nti}[{sbs96_sub}]{ntj}"
            SBS96_LST.append(sbs96)

SBS96_TRI_LST = []
for nti in NTS:
    for ntj in ["C", "G"]:
        for ntk in NTS:
            tri = "{}{}{}".format(nti, ntj, ntk)
            SBS96_TRI_LST.append(tri)

SBS52_TO_SBS96_LST = defaultdict(list)
for SBS96, SBS52 in SBS96_TO_SBS52.items():
    SBS52_TO_SBS96_LST[SBS52].append(SBS96)

SBS52_TO_SBS96_TRI_LST = defaultdict(set)
for SBS96, SBS52 in SBS96_TO_SBS52.items():
    SBS52_TO_SBS96_LST[SBS52].append(SBS96)
    ubase, _, ref, _, alt, _, dbase = list(SBS96)
    tri = f"{ubase}{ref}{dbase}"
    SBS52_TO_SBS96_TRI_LST[SBS52].add(tri)

SBS52_TRI_SET = set()
SUB_TO_SBS52_LST = defaultdict(list)
for SBS52 in set(list(SBS96_TO_SBS52.values())):
    ubase, _, ref, _, alt, _, dbase = list(SBS52)
    sub = "{}>{}".format(ref, alt)
    tri = "{}{}{}".format(ubase, ref, dbase)
    SUB_TO_SBS52_LST[sub].append(SBS52)
    SBS52_TRI_SET.add(tri)
SBS52_TRI_LST = natsort.natsorted(list(SBS52_TRI_SET))
SBS52_TRI_LST = natsort.natsorted(SBS52_TRI_LST)
SBS52_TRI_WEIGHT = len(SBS52_TRI_LST)

SBS52_LST = []
for sbs52_sub in SBS52_SUB_LST:
    SUB_TO_SBS52_LST[sbs52_sub] = natsort.natsorted(SUB_TO_SBS52_LST[sbs52_sub])
    SBS52_LST.extend(SUB_TO_SBS52_LST[sbs52_sub])

for SBS52 in SBS52_TO_SBS96_LST:
    SBS52_TO_SBS96_LST[SBS52] = natsort.natsorted(list(set(SBS52_TO_SBS96_LST[SBS52])))

TRI_WITH_ANNOTATION = []
for sbs52_sub in SBS52_SUB_LST:
    for sbs52 in SUB_TO_SBS52_LST[sbs52_sub]:
        ubase, _, ref, _, alt, _, dbase = list(sbs52)
        tri = f"{ubase}{ref}{dbase}"
        TRI_WITH_ANNOTATION.append("({}) {}".format(";".join(SBS52_TO_SBS96_LST[sbs52]), tri))

SBS52_MUTSIG_FILL_COLOURS = ("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#F5ABCC")


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--vcf",
        type=Path,
        required=True,
        help="VCF file to read"
    )
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="reference FASTA file to read"
    )
    parser.add_argument(
        "--region",
        type=Path,
        required=False,
        help="target chromosome"
    )
    parser.add_argument(
        "--region-list",
        type=Path,
        required=False,
        help="list of target chromosomes separated by new line"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        required=False,
        help="number of threads"
    )
    args = args[1:]
    return parser.parse_args(args)


def check_num_threads(thread_count: int):
    system_thread_count = psutil.cpu_count()
    if thread_count > system_thread_count:
        logging.info("Operating system does not have {} threads".format(thread_count))
        sys.exit()


def get_sample(vcf_file_path: Path):
    vcf_file = pysam.VariantFile(vcf_file_path)
    sample = vcf_file.header.samples[0]
    vcf_file.close()
    return sample


def load_loci(
    region: str,
    region_list: str
):
    chrom_lst = []
    if region is None and region_list is not None:
        for line in open(region_list).readlines():
            chrom_lst.append(line.strip())
    elif region is not None and region_list is None:
        chrom_lst.append(region)
    elif region is not None and region_list is not None:
        for line in open(region_list).readlines():
            chrom_lst.append(line.strip())
    else:
        sys.exit()
    return natsort.natsorted(chrom_lst)


def get_sbs96(
    xvariant: ExpandedVariantRecord,
    reference_sequence_lookup: pysam.FastaFile,
):
    tri = reference_sequence_lookup.fetch(xvariant.chrom, xvariant.pos - 2, xvariant.pos + 1)
    if xvariant.ref in PUR_SET:
        ubase, _, dbase = tri[::-1]
        sbs96 = "{}[{}>{}]{}".format(
            COMPLEMENTARY_BASE_LOOKUP.get(ubase, "N"),
            COMPLEMENTARY_BASE_LOOKUP.get(xvariant.ref, "N"),
            COMPLEMENTARY_BASE_LOOKUP.get(xvariant.alt, "N"),
            COMPLEMENTARY_BASE_LOOKUP.get(dbase, "N"),
        )
    else:
        ubase, _, dbase = tri
        sbs96 = "{}[{}>{}]{}".format(ubase, xvariant.ref, xvariant.alt, dbase)
    if sbs96.count("N") == 0:
        return sbs96
    return "N[N>N]N"


def load_sbs96_counts(
    chromosomes: List[str],
    vcf_file_path: Path,
    ref_file_path: Path,
) -> Dict[str, int]:
    vcf_file = pysam.VariantFile(vcf_file_path)
    reference_sequence_lookup = pysam.FastaFile(ref_file_path)
    if vcf_file_path.suffix == ".vcf":
        sbs96_per_chrom = defaultdict(lambda: {sbs96: 0 for sbs96 in SBS96_LST})
        for variant in vcf_file:
            xvariant = ExpandedVariantRecord(variant)
            if xvariant.is_pass and xvariant.is_biallelic_snp:
                sbs96 = get_sbs96(xvariant, reference_sequence_lookup)
                sbs96_per_chrom[variant.chrom][sbs96] += 1
    elif vcf_file_path.suffix == ".bgz":
        sbs96_per_chrom = {chrom: {sbs96: 0 for sbs96 in SBS96_LST} for chrom in chromosomes}
        for chrom in chromosomes:
            for variant in vcf_file.fetch(chrom):
                xvariant = ExpandedVariantRecord(variant)
                if xvariant.is_pass and xvariant.is_biallelic_snp:
                    sbs96 = get_sbs96(xvariant, reference_sequence_lookup)
                    sbs96_per_chrom[chrom][sbs96] += 1
    sbs96_counts = {sbs96: 0 for sbs96 in SBS96_LST}
    for chrom in chromosomes:
        for sbs96, count in sbs96_per_chrom[chrom].items():
            if sbs96.count("N") != 0:
                continue
            sbs96_counts[sbs96] += count
    return sbs96_counts


def load_sbs52_counts(
    chromosomes: List[str],
    vcf_file_path: Path,
    ref_file_path: Path,
):
    sbs52_counts = {sbs52: 0 for sbs52 in SBS52_LST}
    sbs96_counts = load_sbs96_counts(chromosomes, vcf_file_path, ref_file_path)
    for sbs96, count in sbs96_counts.items():
        if sbs96.count("N") != 0:
            continue
        sbs52_counts[SBS96_TO_SBS52[sbs96]] += count
    return sbs52_counts


def get_sbs96_tricounts_per_target(
    chrom: str,
    sequence: str,
    tricount_lookup_per_chrom: Dict[str, Dict[str, int]],
):
    tricount_lookup = defaultdict(lambda: 0)
    for seq_idx, _ in enumerate(sequence[:-2]):
        trinucleotide = sequence[seq_idx:seq_idx+3]
        if trinucleotide.count("N") != 0:
            continue
        if trinucleotide[1] in PUR_SET:
            trinucleotide = "".join([COMPLEMENTARY_BASE_LOOKUP[nt] for nt in trinucleotide[::-1]])
        tricount_lookup[trinucleotide] += 1
    tricount_lookup_per_chrom[chrom] = dict(tricount_lookup)


def get_target_sbs96_tricounts(
    chromosomes: List[str],
    ref_file_path: Path,
    threads: int
) -> Dict[str, Dict[str, int]]:

    p = mp.Pool(threads)
    manager = mp.Manager()
    ref_tricount_lookup_per_chrom = manager.dict()
    reference_sequence_lookup = pysam.FastaFile(ref_file_path)
    get_tricounts_per_target_arg_lst = [
        (
            chrom,
            str(reference_sequence_lookup[chrom]),
            ref_tricount_lookup_per_chrom
        )
        for chrom in chromosomes
    ]
    p.starmap(get_sbs96_tricounts_per_target, get_tricounts_per_target_arg_lst)
    p.close()
    p.join()

    ref_tricount_lookup = defaultdict(lambda: 0)
    for chrom in chromosomes:
        for tri, count in ref_tricount_lookup_per_chrom[chrom].items():
            ref_tricount_lookup[tri] += count
    return ref_tricount_lookup


def get_target_sbs52_tricounts(
    chromosomes: List[str],
    ref_file_path: Path,
    threads: int
):
    ref_sbs52_tri_counts = {tri: 0 for tri in SBS52_TRI_LST}
    ref_sbs96_tri_counts = get_target_sbs96_tricounts(chromosomes, ref_file_path, threads)
    for sbs52 in SBS52_LST:
        ubase, _, ref, _, alt, _, dbase = list(sbs52)
        sbs52_tri = f"{ubase}{ref}{dbase}"
        for sbs96_tri in SBS52_TO_SBS96_TRI_LST[sbs52]:
            ref_sbs52_tri_counts[sbs52_tri] += ref_sbs96_tri_counts[sbs96_tri]
    for sbs52, sbs96_lst in SBS52_TO_SBS96_LST.items():  # TODO
        sbs96_tri_lst = []
        for sbs96 in sbs96_lst:
            ubase, _, ref, _, alt, _, dbase = list(sbs96)
            sbs96_tri = f"{ubase}{ref}{dbase}"
            sbs96_tri_lst.append(sbs96_tri)
        print(sbs52, ",".join(sbs96_tri_lst))
    return ref_sbs52_tri_counts


def draw_sbs52_barplot(
    sbs52_file: Path,
    sample: str,
    sbs52_pdf_file: str
):

    df = pd.read_csv(sbs52_file, sep="\t")
    df["TRI"] = pd.Categorical(
        df["TRI"],
        categories=TRI_WITH_ANNOTATION,
        ordered=True
    )
    plot = (
        p9.ggplot(df, p9.aes(x="TRI", y="NORMCOUNT", fill="SUB"))
        + p9.geom_bar(stat="identity")
        + p9.theme_bw()
        + p9.facet_grid(". ~ SUB", scales="free")
        + p9.scale_fill_manual(
            values=SBS52_MUTSIG_FILL_COLOURS
        )
        + p9.labs(x="\nTrinucleotide Context\n", y="\nCounts\n")
        + p9.ggtitle("\n{}\n".format(sample))
        + p9.theme(
            text=p9.element_text(size=10),
            legend_title=p9.element_blank(),
            axis_text_x=p9.element_text(
                family="monospace", angle=90, ha="center"
            ),
        )
    )
    plot.save(sbs52_pdf_file, width=22, height=12)


def write_tri_equal_weight_sbs52_counts(
    vcf_file_path: Path,
    ref_file_path: Path,
    region: Path,
    region_list: Path,
    threads: int,
):

    sample = get_sample(vcf_file_path)
    chromosomes = load_loci(region, region_list)
    sbs52_counts = load_sbs52_counts(chromosomes, vcf_file_path, ref_file_path)
    ref_sbs52_tricount_lookup = get_target_sbs52_tricounts(chromosomes, ref_file_path, threads)
    ref_sbs52_tri_sum = sum(ref_sbs52_tricount_lookup.values())
    ref_sbs52_tri_frequencies_lookup = {
        tri: tricount/ref_sbs52_tri_sum
        for tri, tricount in ref_sbs52_tricount_lookup.items()
    }
    sbs52_file_prefix = str(vcf_file_path).replace(".vcf", "").replace(".bgz", "")
    sbs52_file_path = Path(f"{sbs52_file_prefix}.tri_equal_weight.sbs52.tsv")
    sbs52_pdf_file_path = Path(f"{sbs52_file_prefix}.tri_equal_weight.sbs52.pdf")
    with open(sbs52_file_path, "w") as outfile:
        print(
            "SUB",
            "TRI",
            "SBS52",
            "NORMCOUNT",
            "COUNT",
            "WEIGHT",
            sep="\t",
            file=outfile
        )
        for sbs52 in SBS52_LST:
            ubase, _, ref, _, alt, _, dbase = list(sbs52)
            sub = f"{ref}>{alt}"
            tri = f"{ubase}{ref}{dbase}"
            tri_with_annotation = "({}) {}".format(";".join(SBS52_TO_SBS96_LST[sbs52]), tri)
            sbs52_tri_weight = 1/(ref_sbs52_tri_frequencies_lookup[tri]/SBS52_TRI_WEIGHT)
            sbs52_count = sbs52_counts[sbs52]
            normalised_sbs52_count = sbs52_count * sbs52_tri_weight
            print(
                sub,
                tri_with_annotation,
                sbs52,
                normalised_sbs52_count,
                sbs52_count,
                sbs52_tri_weight,
                sep="\t",
                file=outfile
            )
    draw_sbs52_barplot(sbs52_file_path, sample, sbs52_pdf_file_path)


def main():
    options = parse_args(sys.argv)
    write_tri_equal_weight_sbs52_counts(options.vcf, options.ref, options.region, options.region_list, options.threads)
    sys.exit(0)


if __name__ == "__main__":
    main()
