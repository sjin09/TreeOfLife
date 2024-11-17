import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict
import sys

import natsort

SBS52_SUBS = ["C>A", "C>G", "C>T", "T>A", "T>G"]


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="file to read SBS96 counts",
    )
    parser.add_argument(
        "--sbs96-to-sbs52",
        type=str,
        required=True,
        help="SBS96 to SBS52 lookup table",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return SBS52 counts",
    )
    args = args[1:]
    return parser.parse_args(args)


def get_sbs96_mapped_to_sbs52(input_path: Path, sbs96_to_sbs52_lookup: Dict[str, str]) -> Dict[str, int]:
    count_per_sbs52 = defaultdict(lambda: 0)
    for line in open(input_path).readlines():
        if line.startswith("SUB"):
            continue
        fields = line.rstrip().split("\t")
        sbs96 = fields[2]
        count = fields[3]
        count_per_sbs52[sbs96_to_sbs52_lookup[sbs96]] += float(count)
    return count_per_sbs52


def write_sbs96_mapped_to_sbs52(input_path: Path, sbs96_to_sbs52_path: Path, output_path: Path):
    sbs52_classifications_per_sub = defaultdict(set)
    sbs52_to_sbs96_classifications = defaultdict(list)
    sbs96_to_sbs52_lookup = dict(line.rstrip().split("\t") for line in open(sbs96_to_sbs52_path))
    for (sbs96, sbs52) in sbs96_to_sbs52_lookup.items():
        ref, alt = sbs52[2], sbs52[4]
        sbs52_classifications_per_sub[f"{ref}>{alt}"].add(sbs52)
        sbs52_to_sbs96_classifications[sbs52].append(sbs96)
    sbs52_classifications_per_sub = {k: natsort.natsorted(list(v)) for k, v in sbs52_classifications_per_sub.items()}
    sbs52_classification = [sbs52 for sub in SBS52_SUBS for sbs52 in sbs52_classifications_per_sub[sub]]
    sbs52_to_sbs96_classifications = {k: natsort.natsorted(set(v)) for k, v in sbs52_to_sbs96_classifications.items()}
    count_per_sbs52 = get_sbs96_mapped_to_sbs52(input_path, sbs96_to_sbs52_lookup)
    with open(output_path, "w") as outfile:
        outfile.write("SUBSTITUTION\tTRINUCLEOTIDE\tSBS52\tCOUNT\n")
        for sbs52 in sbs52_classification:
            ubase, _, ref, _, alt, _, dbase = list(sbs52)
            sub = f"{ref}>{alt}"
            tri = f"{ubase}{ref}{dbase}"
            tri = "({}) {}".format(";".join(sbs52_to_sbs96_classifications[sbs52]), tri)
            count = count_per_sbs52[sbs52]
            if count.is_integer():
                outfile.write("{}\t{}\t{}\t{}\n".format(sub, tri, sbs52, int(count)))
            else:
                outfile.write("{}\t{}\t{}\t{}\n".format(sub, tri, sbs52, count))


def main():
    options = parse_args(sys.argv)
    write_sbs96_mapped_to_sbs52(options.input, options.sbs96_to_sbs52, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
