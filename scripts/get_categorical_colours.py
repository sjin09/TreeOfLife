#!/usr/bin/env python

import argparse
from pathlib import Path
import sys
from typing import List, Tuple

import colorcet as cc


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Return categorical colours",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-n",
        "--number",
        type=int,
        required=True,
        help="number of categorical colours"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="file to return"
    )
    args = args[1:]
    return parser.parse_args(args)


def rgb_to_hex(rgbs: List[Tuple[float, float, float]]) -> List[str]:
    hex_codes = []
    for rgb in rgbs:
        r, g, b = [int(round(c * 255)) for c in rgb]
        hex_codes.append(f"#{r:02X}{g:02X}{b:02X}")
    return hex_codes


def write_categorical_colours(number_of_colours: int, output_path: Path) -> None:
    palette = cc.glasbey_hv
    rgbs = palette[:number_of_colours]

    # RGB to hex codes
    hex_codes = rgb_to_hex(rgbs)

    # Return hex codes
    with open(output_path, "w") as outfile:
        for hex_code in hex_codes:
            outfile.write("{}\n".format(hex_code))


def main() -> int:
    options = parse_args(sys.argv)
    write_categorical_colours(options.number, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
