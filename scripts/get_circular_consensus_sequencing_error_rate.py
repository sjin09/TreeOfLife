#!/usr/bin/env python

import argparse
import math
from pathlib import Path
from typing import Dict, Tuple
import sys

import pysam

AUTOSOMES = [str(i) for i in range(1, 23)]
CORD_BLOOD_MUTATION_BURDEN = 40


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Get CCS substitution error rate",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="file to read PD47269d haplotype-phased and normalised somatic mutation counts"
    )
    parser.add_argument(
        "--ref-fasta",
        type=str,
        required=True,
        help="reference FASTA file"
    )
    parser.add_argument(
        "--gold-standard",
        type=str,
        required=True,
        help="file to read CB001 somatic mutation counts"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to write"
    )
    args = args[1:]
    return parser.parse_args(args)


def get_reference_base_counts(ref_fasta_path: Path) -> Tuple[int, int]:
    # get autosome base counts
    autosome_base_count = 0
    sequence_lookup = pysam.FastaFile(ref_fasta_path)
    for chrom in AUTOSOMES:
        seq_length = sequence_lookup.get_reference_length(chrom)
        n_count = sequence_lookup[chrom].count("N")
        autosome_base_count += (seq_length - n_count)
    # get genome base counts
    genome_base_count = 0
    x_length = sequence_lookup.get_reference_length("X")
    y_length = sequence_lookup.get_reference_length("Y")
    xn_count = sequence_lookup["X"].count("N")
    yn_count = sequence_lookup["Y"].count("N")
    genome_base_count += (autosome_base_count)
    genome_base_count += (x_length - xn_count)
    genome_base_count += (y_length - yn_count)
    sequence_lookup.close()
    return autosome_base_count, genome_base_count


def get_autosome_diploid_sequence_coverage(input_path: Path) -> float:
    ref_callable_tri_count_sum = 0
    ccs_callable_tri_count_sum = 0
    for line in open(input_path).readlines():
        if line.startswith("#"):
            continue
        elif line.startswith("sub"):
            continue
        fields = line.rstrip().split()
        ref_callable_tri_count = int(fields[-2])
        ccs_callable_tri_count = int(fields[-1])
    ref_callable_tri_count_sum += ref_callable_tri_count
    ccs_callable_tri_count_sum += ccs_callable_tri_count
    dipcov = ccs_callable_tri_count_sum / (2 * ref_callable_tri_count_sum)
    return dipcov


def load_gold_standard_mutational_probabilities(gold_standard_path: Path) -> Dict[str, float]:
    count_per_sbs96 = {}
    prob_per_sbs96 = {}
    for line in open(gold_standard_path).readlines():
        fields = line.rstrip().split()
        sbs96 = fields[2]
        count = int(fields[3])
        count_per_sbs96[sbs96] = int(count)
    total_mutation_count = sum(count_per_sbs96.values())
    for k, v in count_per_sbs96.items():
        prob_per_sbs96[k] = v/total_mutation_count
    return prob_per_sbs96


def get_expected_somatic_mutation_count(dipcov: float, ref_fasta_path: Path) -> float:
    autosome_base_count, genome_base_count = get_reference_base_counts(ref_fasta_path)
    expected_somatic_mutation_count = dipcov * (
        CORD_BLOOD_MUTATION_BURDEN * autosome_base_count/float(genome_base_count)
    )
    return expected_somatic_mutation_count


def get_phred_scaled_q_score(mismatch_rate: float):
    phred_q_score = -10 * math.log10(mismatch_rate)
    return phred_q_score


def get_circular_consensus_sequencing_error_rate(
    input_path: Path,
    ref_fasta_path: Path,
    gold_standard_path: Path,
    output_path: Path
):
    dipcov = get_autosome_diploid_sequence_coverage(input_path)
    expected_total_somatic_mutaton_count = get_expected_somatic_mutation_count(dipcov, ref_fasta_path)
    gold_standard_mutational_probabilities = load_gold_standard_mutational_probabilities(gold_standard_path)
    with open(output_path, "w") as outfile:
        outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            "substitution",
            "5_base",
            "3_base",
            "sbs96",
            "error_rate",
            "error_count",
            "callable_ccs_tri_count",
            "phred_q_score"
        ))
        for line in open(input_path).readlines():
            if line.startswith("#"):
                continue
            elif line.startswith("sub"):
                continue
            fields = line.rstrip().split()
            sub = fields[0]
            tri = fields[1]
            ubase = tri[0]
            dbase = tri[2]
            sbs96 = fields[2]
            observed_count = float(fields[4])
            callable_ccs_tri_count = int(fields[-1])
            prob = gold_standard_mutational_probabilities[sbs96]
            expected_somatic_mutation_count = expected_total_somatic_mutaton_count * prob
            ccs_error_count = observed_count - expected_somatic_mutation_count
            if ccs_error_count > 0:
                ccs_error_rate = ccs_error_count/callable_ccs_tri_count
            else:
                ccs_error_rate = 1/callable_ccs_tri_count
            phred_q_score = get_phred_scaled_q_score(ccs_error_rate)
            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                sub,
                ubase,
                dbase,
                sbs96,
                ccs_error_rate,
                ccs_error_count,
                callable_ccs_tri_count,
                phred_q_score
            ))


def main() -> int:
    options = parse_args(sys.argv)
    get_circular_consensus_sequencing_error_rate(
        options.input,
        options.ref_fasta,
        options.gold_standard,
        options.output
    )
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
