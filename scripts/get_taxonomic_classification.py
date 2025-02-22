#!/usr/bin/env python

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict
import sys

import pandas as pd
import requests

STATUS_SUMMARY = set(["1 submitted", "2 curated", "3 curation"])


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Get taxonomic classification per sample",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="file to read DToL samplesheet"
    )
    parser.add_argument(
        "--samples",
        type=Path,
        required=True,
        help="a list of samples, separated by new line"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="file to write taxonomic classifications per sample"
    )
    args = args[1:]
    return parser.parse_args(args)


@dataclass
class TaxonomicClassification:
    domain: str = "Eukaryota"
    kingdom: str = "."
    phylum: str = "."
    taxonomic_class: str = "."
    order: str = "."
    family: str = "."
    genus: str = "."
    species: str = "."
    common_name: str = "."


def get_reference_sample_lookup(samples_path: Path) -> Dict[str, str]:
    reference_sample_lookup = {}
    for line in open(samples_path).readlines():
        ref_sample, sample = line.rstrip().split()
        reference_sample_lookup[sample] = ref_sample
    return reference_sample_lookup


def make_api_request(url: str, timeout=30):
    session = requests.Session()
    try:
        response = session.get(url, timeout=timeout)
        response.raise_for_status()
        return response
    except requests.exceptions.Timeout:
        print(f"Request timed out after {timeout} seconds")
        raise
    except requests.exceptions.ConnectionError as e:
        print(f"Connection error occurred: {e}")
        raise
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")
        raise
    finally:
        session.close()


def get_taxonomy_data(family_taxon: str):
    response = make_api_request(
        f"https://goat.genomehubs.org/api/v2/record?recordId={family_taxon}&result=taxon&taxonomy=ncbi"
    )
    return response.json()


def get_taxonomic_classification_per_sample(
    input_path: Path, reference_sample_lookup: Dict[str, str]
) -> Dict[str, TaxonomicClassification]:
    # Initialize variables
    fieldnames = [
        "Sample",
        "Common name",
        "Domain",
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
    ]
    reference_samples = set(reference_sample_lookup.values())
    seen = set()
    taxonomic_classification_per_sample = {}

    # Write taxonomic classification per reference sample
    with open("dtol_reference_samples.taxonomic_classification.csv", "w") as outfile:
        # Create a DictWriter object
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)

        # Write header
        writer.writeheader()

        # Import and iterate through the DToL samplesheet
        df = pd.read_csv(input_path)
        for (_idx, row) in df.iterrows():
            sample = row["sample"]
            if sample in seen:
                continue
            if sample not in reference_samples:
                continue
            if not sample.startswith("ilYpoPade"):  # TODO
                continue
            status_summary = row["statussummary"]
            if status_summary not in STATUS_SUMMARY:
                continue
            family = row["family"]
            genus = row["genus"]
            species = row["species"]
            family_taxon = int(row["family_taxon"])
            common_name = row["common name"]
            if pd.isna(common_name):
                common_name = "."
            else:
                common_name = common_name.capitalize()
            ncbi_response = get_taxonomy_data(family_taxon)
            try:
                seen.add(sample)
                ncbi_taxonomic_records = ncbi_response["records"][0]["record"]["lineage"]
                for taxonomy_record in ncbi_taxonomic_records:
                    taxonomic_rank = taxonomy_record["taxon_rank"]
                    if "kingdom" == taxonomic_rank:
                        kingdom = taxonomy_record["scientific_name"]
                    elif "phylum" == taxonomic_rank:
                        phylum = taxonomy_record["scientific_name"]
                    elif "class" == taxonomic_rank:
                        taxonomic_class = taxonomy_record["scientific_name"]
                    elif "order" == taxonomic_rank:
                        order = taxonomy_record["scientific_name"]
                # Addresing a bug in the original samplesheet
                if sample.startswith("ilYpoPade"):
                    genus = "Yponomeuta"
                    species = "Yponomeuta padella"
                sample_taxanomic_classification = TaxonomicClassification(
                    kingdom=kingdom,
                    phylum=phylum,
                    taxonomic_class=taxonomic_class,
                    order=order,
                    family=family,
                    genus=genus,
                    species=species,
                    common_name=common_name
                )
                taxonomic_classification_per_sample[sample] = sample_taxanomic_classification
                writer.writerow({
                    "Sample": sample,
                    "Common name": sample_taxanomic_classification.common_name,
                    "Domain": sample_taxanomic_classification.domain,
                    "Kingdom": sample_taxanomic_classification.kingdom,
                    "Phylum": sample_taxanomic_classification.phylum,
                    "Class": sample_taxanomic_classification.taxonomic_class,
                    "Order": sample_taxanomic_classification.order,
                    "Family": sample_taxanomic_classification.family,
                    "Genus": sample_taxanomic_classification.genus,
                    "Species": sample_taxanomic_classification.species
                })
            except IndexError:
                raise ValueError(f"API request for sample '{sample}' failed")

    return taxonomic_classification_per_sample


def write_taxonomic_classification(input_path: Path, samples_path: Path, output_path: Path) -> None:
    # Initialize variables
    fieldnames = [
        "Reference sample",
        "Sample",
        "Common name",
        "Domain",
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
    ]
    reference_sample_lookup = get_reference_sample_lookup(samples_path)
    taxonomic_classification_per_sample = get_taxonomic_classification_per_sample(input_path, reference_sample_lookup)

    with open(output_path, "w") as outfile:
        # Create a DictWriter object
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)

        # Write header
        writer.writeheader()

        # Write taxonomic classification per sample
        for sample, reference_sample in reference_sample_lookup.items():
            taxonomic_classification = taxonomic_classification_per_sample.get(reference_sample)
            if taxonomic_classification is None:
                continue
            writer.writerow({
                "Reference sample": reference_sample,
                "Sample": sample,
                "Common name": taxonomic_classification.common_name,
                "Domain": taxonomic_classification.domain,
                "Kingdom": taxonomic_classification.kingdom,
                "Phylum": taxonomic_classification.phylum,
                "Class": taxonomic_classification.taxonomic_class,
                "Order": taxonomic_classification.order,
                "Family": taxonomic_classification.family,
                "Genus": taxonomic_classification.genus,
                "Species": taxonomic_classification.species
            })


def main() -> int:
    options = parse_args(sys.argv)
    write_taxonomic_classification(options.input, options.samples, options.output)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
