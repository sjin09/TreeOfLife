import argparse
from collections import defaultdict
from dataclasses import dataclass
import json
from pathlib import Path
import sys
from typing import Dict, List

from natsort import natsorted
import pandas as pd
import patchworklib as pw
import plotnine as p9
from plotnine.qplot import ggplot

NTS = ["A", "C", "G", "T"]
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SBS52_FILL_COLOURS = ("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#F5ABCC")
SBS96_FILL_COLOURS = ("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#83A603", "#F5ABCC")
SBS96_CLASSIFICATIONS = [f"{nti}[{sub}]{ntj}" for sub in SUBS for nti in NTS for ntj in NTS]


def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--gtol-sigs",
        type=Path,
        required=True,
        help="file to read germline mutational signatures",
    )
    parser.add_argument(
        "--gtol-exps",
        type=Path,
        required=True,
        help="file to read germilne mutational signature exposures",
    )
    parser.add_argument(
        "--gtol-cmean",
        type=Path,
        required=True,
        help="file to read Abouheif's Cmean for germilne mutational signatures",
    )
    parser.add_argument("--sbs52", type=Path, required=True, help="file to read SBS52 classification")
    parser.add_argument(
        "--stol-sigs",
        type=Path,
        required=True,
        help="file to read somatic mutational signatures",
    )
    parser.add_argument(
        "--stol-exps",
        type=Path,
        required=True,
        help="file to read somatic mutational signature exposures",
    )
    parser.add_argument(
        "--stol-cmean", type=Path, required=True, help="file to read Abouheif's Cmean for somatic mutational signatures"
    )
    parser.add_argument(
        "--rtol-sig-names",
        type=Path,
        required=True,
        help="file to read artefactual mutational signature names",
    )
    parser.add_argument(
        "-n",
        "--number",
        type=int,
        default=1,
        required=False,
        help="Starting number for Supplementary Figures"
    )
    parser.add_argument(
        "--dtol-id-to-species-lookup", type=Path, required=True, help="file to read DToL id to species lookup table"
    )
    parser.add_argument(
        "--dtol-id-to-taxonomic-classification-lookup",
        type=Path,
        required=True,
        help="file to read DToL id to species lookup table",
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


def get_substitution(sbs96: str) -> str:
    ubase, _, ref, _, alt, _, dbase = list(sbs96)
    return f"{ref}>{alt}"


def get_trinucleotide(sbs96: str) -> str:
    ubase, _, ref, _, alt, _, dbase = list(sbs96)
    return f"{ubase}{ref}{dbase}"


def load_gtol_signatures(gtol_signatures_path: Path, sbs52_classification_path: Path) -> pd.DataFrame:
    # Load SBS52 classification
    sbs52_classifications = load_sbs52_classification(sbs52_classification_path)

    # Load germline mutational signatures
    sigs_df = pd.read_csv(gtol_signatures_path, sep=",", index_col=0, header=0)
    sigs_df.index = sbs52_classifications  # Direct index assignment
    sigs_df.index.name = "SBS52"

    # Convert to long format
    long_sigs_df = sigs_df.reset_index().melt(id_vars=["SBS52"], var_name="Sig", value_name="Probability")

    # Convert fraction to percentage
    long_sigs_df["Probability"] = long_sigs_df["Probability"] * 100

    # Drop rows where Sample is 'X0'
    long_sigs_df = long_sigs_df[long_sigs_df["Sig"] != "X0"]

    # Add Substitution and Trinucleotide columns
    long_sigs_df["Substitution"] = long_sigs_df["SBS52"].apply(get_substitution)
    long_sigs_df["Trinucleotide"] = long_sigs_df["SBS52"].apply(get_trinucleotide)

    # Replace 'X' with 'gtol' in the 'Sig' column
    long_sigs_df["Sig"] = long_sigs_df["Sig"].str.replace("X", "gToL")
    return long_sigs_df


def load_gtol_signature_exposures(gtol_exposures_path: Path) -> pd.DataFrame:
    exps_df = pd.read_csv(gtol_exposures_path, sep=",", index_col=0, header=0)
    exps_df.index.name = "Sample"
    long_exps_df = exps_df.reset_index().melt(id_vars=["Sample"], var_name="Sig", value_name="Attribution")

    # Convert proportion to percentage
    long_exps_df["Attribution"] = long_exps_df["Attribution"] * 100

    # Add X prefix to Sig column
    long_exps_df["Sig"] = "gToL" + long_exps_df["Sig"].astype(str)

    # Drop rows where Sig is 'gToL0'
    long_exps_df = long_exps_df[long_exps_df["Sig"] != "gToL0"]
    return long_exps_df


def load_stol_signatures(stol_signatures_path: Path, rtol_names: List[str]) -> pd.DataFrame:
    # Load somatic mutational signatures
    sigs_df = pd.read_csv(stol_signatures_path, sep=",", index_col=0, header=0)
    sigs_df.index = SBS96_CLASSIFICATIONS  # Direct index assignment
    sigs_df.index.name = "SBS96"

    # Convert to long format
    sigs_long_df = sigs_df.reset_index().melt(id_vars=["SBS96"], var_name="Sig", value_name="Probability")

    # Convert fraction to percentage
    sigs_long_df["Probability"] = sigs_long_df["Probability"] * 100

    # Drop rows where Sample is 'sToL0'
    sigs_long_df = sigs_long_df[sigs_long_df["Sig"] != "X0"]

    # Add Substitution and Trinucleotide columns
    sigs_long_df["Substitution"] = sigs_long_df["SBS96"].apply(get_substitution)
    sigs_long_df["Trinucleotide"] = sigs_long_df["SBS96"].apply(get_trinucleotide)

    # Subset into two groups
    rtol_sigs_long_df = sigs_long_df[sigs_long_df["Sig"].isin(rtol_names)]
    stol_sigs_long_df = sigs_long_df[~sigs_long_df["Sig"].isin(rtol_names)]

    # Get current signature names
    hdp_rtol_names = natsorted(rtol_sigs_long_df["Sig"].unique())
    hdp_stol_names = natsorted(stol_sigs_long_df["Sig"].unique())

    # Build lookup dictionaries
    rtol_name_lookup = {rtol_name: f"rToL{idx}" for idx, rtol_name in enumerate(hdp_rtol_names, start=1)}
    stol_name_lookup = {stol_name: f"sToL{idx}" for idx, stol_name in enumerate(hdp_stol_names, start=1)}

    # Change signature names using lookup
    rtol_sigs_long_df["Sig"] = rtol_sigs_long_df["Sig"].map(rtol_name_lookup)
    stol_sigs_long_df["Sig"] = stol_sigs_long_df["Sig"].map(stol_name_lookup)
    return stol_sigs_long_df, rtol_sigs_long_df


def load_stol_signature_exposures(stol_exposures_path: Path, rtol_names: List[str]) -> pd.DataFrame:
    exps_df = pd.read_csv(stol_exposures_path, sep=",", index_col=0, header=0)
    exps_df.index.name = "Sample"
    exps_long_df = exps_df.reset_index().melt(id_vars=["Sample"], var_name="Sig", value_name="Attribution")

    # Convert proportion to percentage
    exps_long_df["Attribution"] = exps_long_df["Attribution"] * 100

    # Add X prefix to Sig column
    exps_long_df["Sig"] = "X" + exps_long_df["Sig"].astype(str)

    # Drop rows where Sig is 'X0'
    exps_long_df = exps_long_df[exps_long_df["Sig"] != "X0"]

    # Subset into two groups
    rtol_exps_long_df = exps_long_df[exps_long_df["Sig"].isin(rtol_names)]
    stol_exps_long_df = exps_long_df[~exps_long_df["Sig"].isin(rtol_names)]

    # Get current signature names
    hdp_rtol_names = natsorted(rtol_exps_long_df["Sig"].unique())
    hdp_stol_names = natsorted(stol_exps_long_df["Sig"].unique())

    # Build lookup dictionaries
    rtol_name_lookup = {rtol_name: f"rToL{idx}" for idx, rtol_name in enumerate(hdp_rtol_names, start=1)}
    stol_name_lookup = {stol_name: f"sToL{idx}" for idx, stol_name in enumerate(hdp_stol_names, start=1)}

    # Change signature names using lookup
    rtol_exps_long_df["Sig"] = rtol_exps_long_df["Sig"].map(rtol_name_lookup)
    stol_exps_long_df["Sig"] = stol_exps_long_df["Sig"].map(stol_name_lookup)
    return stol_exps_long_df, rtol_exps_long_df


def load_abouheif_cmean(cmean_path: Path) -> pd.DataFrame:
    cmean_df = pd.read_csv(cmean_path, sep=",", header=0)
    return cmean_df


def load_id_to_species_lookup(dtol_id_to_species_lookup_path: Path) -> Dict[str, str]:
    with open(dtol_id_to_species_lookup_path, "r") as f:
        id_to_species_lookup = json.load(f)
    return id_to_species_lookup


def load_id_to_taxonomic_classification_lookup(
    dtol_id_to_taxonomic_classification_lookup_path: Path,
) -> Dict[str, TaxonomicClassification]:
    df = pd.read_csv(dtol_id_to_taxonomic_classification_lookup_path, sep=",")
    id_to_taxonomic_classification_lookup = defaultdict(TaxonomicClassification)
    for _, row in df.iterrows():
        id_to_taxonomic_classification_lookup[row["Sample"]] = TaxonomicClassification(
            kingdom=row["Kingdom"],
            phylum=row["Phylum"],
            taxonomic_class=row["Class"],
            order=row["Order"],
            family=row["Family"],
            genus=row["Genus"],
            species=row["Species"],
        )
    return id_to_taxonomic_classification_lookup


def load_rtol_names(rtol_names_path: Path) -> List[str]:
    rtol_names = [line.rstrip() for line in open(rtol_names_path)]
    return rtol_names


def load_sample_count_per_reference_sample(
    dtol_id_to_taxonomic_classification_lookup_path: Path,
) -> Dict[str, int]:
    # Load data frame
    df = pd.read_csv(dtol_id_to_taxonomic_classification_lookup_path, sep=",")

    # Count number of samples per reference sample
    sample_count_per_reference_sample = defaultdict(lambda: 0)
    for _, row in df.iterrows():
        sample_count_per_reference_sample[row["ref_sample"]] += 1
    return sample_count_per_reference_sample


def load_sample_count_per_species_name(
    dtol_id_to_taxonomic_classification_lookup_path: Path,
) -> Dict[str, int]:
    # Load data frame
    df = pd.read_csv(dtol_id_to_taxonomic_classification_lookup_path, sep=",")

    # Count number of samples per reference sample
    sample_count_per_species_name = defaultdict(lambda: 0)
    for _, row in df.iterrows():
        sample_count_per_species_name[row["Species"]] += 1
    return sample_count_per_species_name


def load_sbs52_classification(sbs52_path: Path) -> List[str]:
    sbs52_classifications = [line.rstrip() for line in open(sbs52_path)]
    return sbs52_classifications


def plot_caption(sig_name: str, number: int, is_rtol_signature: bool = False):
    if is_rtol_signature:
        caption = (
            f"Supplementary Figure {number}.\n"
            "a) Mutational signature spectrum.\n"
            "b) Bar plot displaying the 30 samples that exhibit the greatest exposure to the mutational signature."
        )
    else:
        caption = (
            f"Supplementary Figure {number}. {sig_name} mutational signature\n"
            "a) Mutational signature spectrum.\n"
            "b) Bar plot displaying the 30 samples that exhibit the greatest exposure to the mutational signature.\n"
            "c) Phylogenetic signal of the mutational signature."
        )
    plot = (
        p9.ggplot()
        + p9.theme_void()  # no axes, no grid
        + p9.xlim(0, 1)
        + p9.ylim(0, 1)
        + p9.annotate("text", x=0, y=0, label=caption, ha="left", va="bottom", size=12, family="Helvetica")
    )
    return plot


def plot_mutational_signature(df: pd.DataFrame, sig_name: str, is_gtol: bool) -> p9.ggplot:
    # Subset data frame
    df_subset = df[df["Sig"] == sig_name]

    # Plot mutational signature
    plot = (
        p9.ggplot(df_subset, p9.aes(x="Trinucleotide", y="Probability", fill="Substitution"))
        + p9.geom_bar(stat="identity")
        + p9.theme_bw()
        + p9.facet_grid(". ~ Substitution", scales="free")
        + p9.labs(
            x="Trinucleotide",
            y="Proportion of Mutations (%)"
        )
        + p9.theme(
            text=p9.element_text(family="Courier", size=8),
            axis_text_x=p9.element_text(angle=90, ha="center", size=6),
            axis_text_y=p9.element_text(family="Helvetica"),
            axis_title_x=p9.element_text(family="Helvetica"),  # x-axis label
            axis_title_y=p9.element_text(family="Helvetica"),  # y-axis label
            legend_position="none",
        )
    )
    if is_gtol:
        plot += p9.scale_fill_manual(values=SBS52_FILL_COLOURS)
    else:
        plot += p9.scale_fill_manual(values=SBS96_FILL_COLOURS)
    return plot


def plot_mutational_signature_attribution(
    df: pd.DataFrame,
    sig_name: str,
    id_to_species_lookup: Dict[str, str],
    sample_count_per_species_name: Dict[str, int]
) -> p9.ggplot:
    # Subset data frame
    df_subset = df[df["Sig"] == sig_name].copy()
    df_subset = df_subset.nlargest(30, "Attribution").copy()

    # Generate lookup table
    sample_name_to_species_name_lookup = {}
    sample_names = df_subset["Sample"].unique()
    for sample_name in sample_names:
        if "." not in sample_name:
            species_name = id_to_species_lookup[sample_name].replace(" ", r"\ ")
            if sample_count_per_species_name[id_to_species_lookup[sample_name]] == 1:
                sample_name_to_species_name_lookup[sample_name] = fr"$\it{{{species_name}}}$"
            else:
                sample_name_to_species_name_lookup[sample_name] = fr"$\it{{{species_name}}}$ ({sample_name})"
        else:
            ref_sample, _sample_name = sample_name.split(".")
            species_name = id_to_species_lookup[ref_sample].replace(" ", r"\ ")
            sample_name_to_species_name_lookup[
                sample_name
            ] = fr"$\it{{{species_name}}}$ ({_sample_name})"

    # Add a column for the x-axis
    df_subset["xaxis"] = df_subset["Sample"].map(sample_name_to_species_name_lookup)

    # Order x-axis by decreasing Attribution
    df_order = df_subset.sort_values("Attribution", ascending=False)["xaxis"].tolist()
    df_subset["xaxis"] = pd.Categorical(df_subset["xaxis"], categories=df_order, ordered=True)

    # Generate plot
    plot = (
        ggplot(df_subset, p9.aes(x="xaxis", y="Attribution"))
        + p9.geom_bar(stat="identity")
        + p9.theme_bw()
        + p9.labs(x="Sample", y="Mutational Signature Exposure (%)")
        + p9.theme(
            text=p9.element_text(family="Helvetica", size=8),
            axis_text_x=p9.element_text(angle=90, ha="center"),
            legend_position="none",
        ) +
        p9.scale_y_continuous(limits=(0, 100))
    )
    return plot


def plot_mutational_signature_phylogenetic_signal(df: pd.DataFrame, sig_name: str) -> p9.ggplot:
    # Subset data frame
    df_copy = df.copy()
    df_copy["Group"] = df_copy["Signature"].eq(sig_name)
    df_copy["q_value"] = df_copy["q_value"].clip(lower=1e-4)

    # Generate scatterplot
    plot = (
        p9.ggplot(df_copy, p9.aes(x="C_mean", y="q_value", fill="Group", size="Group"))
        + p9.theme_bw(8)
        + p9.geom_point(alpha=0.9)
        + p9.scale_fill_manual(values={True: "#E41A1C", False: "#D3D3D3"})
        + p9.scale_size_manual(values={True: 5, False: 1})
        + p9.scale_y_log10(
            limits=(1, 1e-4),
            breaks=[1, 1e-1, 1e-2, 1e-3, 1e-4],
            labels=["1", "0.1", "0.01", "0.001", "<0.0001"]
        )
        + p9.labs(
            x="Abouheif's Cmean",
            y="q value",
        )
        + p9.theme(
            legend_position="none",  # hide legend
            text=p9.element_text(family="Helvetica", size=8),
        )
    )
    return plot


def plot_somatic_mutational_signatures(
    figure_number: int,
    stol_signatures_path: Path,
    stol_exposures_path: Path,
    stol_cmean_path: Path,
    rtol_names: List[str],
    id_to_species_lookup: Dict[str, str],
    sample_count_per_species_name: Dict[str, int],
):
    # Change margin
    pw.param["margin"] = 0.1

    # Load data frame
    stol_cmean_df = load_abouheif_cmean(stol_cmean_path)
    stol_sigs_df, rtol_sigs_df = load_stol_signatures(stol_signatures_path, rtol_names)
    stol_exps_df, rtol_exps_df = load_stol_signature_exposures(stol_exposures_path, rtol_names)

    # # Plot somatic mutational signatures
    stol_names = natsorted(stol_sigs_df["Sig"].unique())
    for stol_name in stol_names:
        pw1 = pw.load_ggplot(
            plot_mutational_signature(stol_sigs_df, stol_name, False),
            figsize=(8.27, 3.33)
        )
        pw2 = pw.load_ggplot(
            plot_mutational_signature_attribution(
                stol_exps_df,
                stol_name,
                id_to_species_lookup,
                sample_count_per_species_name=sample_count_per_species_name
            ),
            figsize=(8.27, 3.33)
        )
        pw3 = pw.load_ggplot(
            plot_mutational_signature_phylogenetic_signal(stol_cmean_df, stol_name),
            figsize=(8.27, 3.33)
        )
        pw4 = pw.load_ggplot(plot_caption(stol_name, figure_number, is_rtol_signature=False), figsize=(8.27, 0.8))

        # Add subplot panel labels
        pw1.case.set_title('a)', x=0, y=0.95, loc="right")
        pw2.case.set_title('b)', x=0, y=0.95, loc="right")
        pw3.case.set_title('c)', x=0, y=0.95, loc="right")

        # combine plots
        plot = (pw1 / pw2 / pw3) / pw4

        # return plot
        plot.savefig(f"SF{figure_number}.pdf", dpi=300)

        # increment counter
        figure_number += 1

    # Plot artefactual mutational signatures
    rtol_names = natsorted(rtol_sigs_df["Sig"].unique())
    for rtol_name in rtol_names:
        pw1 = pw.load_ggplot(
            plot_mutational_signature(rtol_sigs_df, rtol_name, False),
            figsize=(8.27, 5)
        )
        pw2 = pw.load_ggplot(
            plot_mutational_signature_attribution(
                rtol_exps_df,
                rtol_name,
                id_to_species_lookup,
                sample_count_per_species_name=sample_count_per_species_name
            ),
            figsize=(8.27, 5)
        )
        pw4 = pw.load_ggplot(plot_caption(rtol_name, figure_number, is_rtol_signature=True), figsize=(8.27, 0.8))

        # Add subplot panel labels
        pw1.case.set_title('a)', x=0, y=0.95, loc="right")
        pw2.case.set_title('b)', x=0, y=0.95, loc="right")

        # combine plots
        plot = (pw1 / pw2) / pw4

        # return plot
        plot.savefig(f"SF{figure_number}.pdf", dpi=300)

        # increment counter
        figure_number += 1
    return figure_number


def plot_germline_mutational_signatures(
    figure_number: int,
    gtol_signatures_path: Path,
    gtol_exposures_path: Path,
    gtol_cmean_path: Path,
    sbs52_path: Path,
    id_to_species_lookup: Dict[str, str],
    sample_count_per_species_name: Dict[str, int]
):
    # Load data frames
    gtol_cmean_df = load_abouheif_cmean(gtol_cmean_path)
    gtol_sigs_df = load_gtol_signatures(gtol_signatures_path, sbs52_path)
    gtol_exps_df = load_gtol_signature_exposures(gtol_exposures_path)

    # Change plot margins
    pw.param["margin"] = 0.1

    # Plot germline mutational signatures
    gtol_names = natsorted(gtol_sigs_df["Sig"].unique())
    for gtol_name in gtol_names:
        pw1 = pw.load_ggplot(
            plot_mutational_signature(gtol_sigs_df, gtol_name, True),
            figsize=(8.27, 3.33)
        )
        pw2 = pw.load_ggplot(
            plot_mutational_signature_attribution(
                gtol_exps_df,
                gtol_name,
                id_to_species_lookup,
                sample_count_per_species_name=sample_count_per_species_name
            ),
            figsize=(8.27, 3.33)
        )
        pw3 = pw.load_ggplot(
            plot_mutational_signature_phylogenetic_signal(gtol_cmean_df, gtol_name),
            figsize=(8.27, 3.33)
        )
        pw4 = pw.load_ggplot(plot_caption(gtol_name, figure_number, is_rtol_signature=False), figsize=(8.27, 0.8))

        # Add subplot panel labels
        pw1.case.set_title('a)', x=0, y=0.95, loc="right")
        pw2.case.set_title('b)', x=0, y=0.95, loc="right")
        pw3.case.set_title('c)', x=0, y=0.95, loc="right")

        # combine plots
        plot = (pw1 / pw2 / pw3) / pw4

        # return plot
        plot.savefig(f"SF{figure_number}.pdf", dpi=300)

        # increment counter
        figure_number += 1


def plot_supplementary_figures(
    gtol_signatures_path: Path,
    gtol_exposures_path: Path,
    gtol_cmean_path: Path,
    sbs52_path: Path,
    rtol_names_path: Path,
    stol_signatures_path: Path,
    stol_exposures_path: Path,
    stol_cmean_path: Path,
    dtol_id_to_species_lookup_path: Path,
    dtol_id_to_taxonomic_classification_lookup_path: Path,
    figure_number: int,
):
    # Load data
    id_to_species_lookup = load_id_to_species_lookup(dtol_id_to_species_lookup_path)
    sample_count_per_species_name = load_sample_count_per_species_name(
        dtol_id_to_taxonomic_classification_lookup_path
    )
    rtol_names = load_rtol_names(rtol_names_path)

    # Plot somatic mutational signatures
    figure_number = plot_somatic_mutational_signatures(
        figure_number=figure_number,
        stol_signatures_path=stol_signatures_path,
        stol_exposures_path=stol_exposures_path,
        stol_cmean_path=stol_cmean_path,
        rtol_names=rtol_names,
        id_to_species_lookup=id_to_species_lookup,
        sample_count_per_species_name=sample_count_per_species_name,
    )

    # Plot germline mutational signatures
    plot_germline_mutational_signatures(
        figure_number=figure_number,
        gtol_signatures_path=gtol_signatures_path,
        gtol_exposures_path=gtol_exposures_path,
        gtol_cmean_path=gtol_cmean_path,
        sbs52_path=sbs52_path,
        id_to_species_lookup=id_to_species_lookup,
        sample_count_per_species_name=sample_count_per_species_name
    )


def main():
    options = parse_args(sys.argv)
    plot_supplementary_figures(
        gtol_signatures_path=options.gtol_sigs,
        gtol_exposures_path=options.gtol_exps,
        gtol_cmean_path=options.gtol_cmean,
        sbs52_path=options.sbs52,
        rtol_names_path=options.rtol_sig_names,
        stol_signatures_path=options.stol_sigs,
        stol_exposures_path=options.stol_exps,
        stol_cmean_path=options.stol_cmean,
        dtol_id_to_species_lookup_path=options.dtol_id_to_species_lookup,
        dtol_id_to_taxonomic_classification_lookup_path=options.dtol_id_to_taxonomic_classification_lookup,
        figure_number=options.number,
    )
    sys.exit(0)


if __name__ == "__main__":
    main()

