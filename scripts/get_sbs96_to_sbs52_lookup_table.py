from dataclasses import dataclass
from typing import Dict, List

from natsort import natsorted

BASE_COMPLEMENT_LOOKUP = str.maketrans("ACGT", "TGCA")
NTS = ["A", "C", "G", "T"]
PYRS = set(["C", "T"])
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SBS96_CLASSIFICATIONS = [f"{nti}[{sub}]{ntj}" for sub in SUBS for nti in NTS for ntj in NTS]


@dataclass
class SBS96:
    ref: str = ""
    alt: str = ""
    ref_tri: str = ""
    alt_tri: str = ""
    hdp_sbs96: str = ""
    sigprofiler_sbs96: str = ""


def get_sbs96_reverse_complement(sbs96: SBS96) -> SBS96:
    sbs96_rc = SBS96()
    sbs96_rc.ref_tri = get_reverse_complementary_sequence(sbs96.ref_tri)
    sbs96_rc.alt_tri = get_reverse_complementary_sequence(sbs96.alt_tri)
    sbs96_rc.ref = sbs96_rc.ref_tri[1]
    sbs96_rc.alt = sbs96_rc.alt_tri[1]
    sbs96_rc.sigprofiler_sbs96 = f"{sbs96_rc.ref_tri[0]}[{sbs96_rc.ref}>{sbs96_rc.alt}]{sbs96_rc.ref_tri[2]}"
    return sbs96_rc


def get_reference_perspective_sbs96(sbs96: SBS96) -> SBS96:
    ref_sbs96 = SBS96()
    ref_sbs96.ref = sbs96.alt_tri[1]
    ref_sbs96.ref_tri = sbs96.alt_tri
    ref_sbs96.sigprofiler_sbs96 = f"{sbs96.ref_tri[0]}[{sbs96.alt}>{sbs96.ref}]{sbs96.ref_tri[2]}"
    return ref_sbs96


def get_complementary_sequence(seq: str) -> str:
    seq_c = seq.translate(BASE_COMPLEMENT_LOOKUP)
    return seq_c


def get_reverse_complementary_sequence(seq: str) -> str:
    seq_rc = seq[::-1].translate(BASE_COMPLEMENT_LOOKUP)
    return seq_rc


def get_sbs96_pyrimidine_count(sbs96: SBS96) -> int:
    pyr_count = sum(1 for nt in sbs96.ref_tri if nt in PYRS)
    return pyr_count


def get_max_pyrimidine_sbs96(sbs96_classifications: List[SBS96]) -> List[SBS96]:
    """
    Get the SBS96 classification with the highest pyrimidine base count.
    """
    max_pyr_count = -1
    max_pyr_sbs96_classifications = []
    for sbs96 in sbs96_classifications:
        sbs96_pyr_count = get_sbs96_pyrimidine_count(sbs96)
        if sbs96_pyr_count > max_pyr_count:
            max_pyr_count = sbs96_pyr_count
            max_pyr_sbs96_classifications = [sbs96]
        elif sbs96_pyr_count == max_pyr_count:
            max_pyr_sbs96_classifications.append(sbs96)
    return max_pyr_sbs96_classifications


def get_lexicographically_sorted_sbs52(sbs96_classifications: List[str]) -> List[str]:
    """
    Get alphabetically sorted SBS96 classifications
    """
    return natsorted(sbs96_classifications)


def get_sbs52(sbs96: SBS96, sbs96_rc: SBS96, ref_sbs96: SBS96, ref_sbs96_rc: SBS96) -> str:
    """
    Determine the SBS52 classification from a list of SBS96 classifications.

    Algorithm:
        - The middle base needs to be a pyrimidine base
        - SBS52 needs to have the highest number of pyrimidine bases
        - If there is a tie, the SBS52 needs to be lexicographically first
    """
    sbs96_classifications = [sbs96, sbs96_rc, ref_sbs96, ref_sbs96_rc]

    # Get the pyrimidine SBS96 classifications
    pyr_sbs96_classifications = []
    for sigprofiler_sbs96 in sbs96_classifications:
        if is_pyr_sbs96(sigprofiler_sbs96):
            print(sigprofiler_sbs96)
            pyr_sbs96_classifications.append(sigprofiler_sbs96)

    # Get the SBS96 classification with the highest pyrimidine base count
    max_pyr_sbs96_classifications = get_max_pyrimidine_sbs96(pyr_sbs96_classifications)
    max_pyr_sbs96_classifications = [
        sbs96.sigprofiler_sbs96
        for sbs96 in get_max_pyrimidine_sbs96(pyr_sbs96_classifications)
    ]
    if len(max_pyr_sbs96_classifications) == 1:
        return max_pyr_sbs96_classifications[0]
    return get_lexicographically_sorted_sbs52(max_pyr_sbs96_classifications)[0]


def get_sbs96_to_sbs52_lookup() -> Dict[str, str]:
    sbs96_to_sbs52_lookup = {}
    for sigprofiler_sbs96 in SBS96_CLASSIFICATIONS:
        # Initialize
        sbs96 = SBS96()
        ubase = sigprofiler_sbs96[0]
        dbase = sigprofiler_sbs96[6]
        sbs96.ref = sigprofiler_sbs96[2]
        sbs96.alt = sigprofiler_sbs96[4]
        sbs96.ref_tri = f"{ubase}{sbs96.ref}{dbase}"
        sbs96.alt_tri = f"{ubase}{sbs96.alt}{dbase}"
        sbs96.sigprofiler_sbs96 = sigprofiler_sbs96
        sbs96_rc = get_sbs96_reverse_complement(sbs96)
        ref_sbs96 = get_reference_perspective_sbs96(sbs96)
        ref_sbs96_rc = get_reference_perspective_sbs96(sbs96_rc)
        sbs52 = get_sbs52(sbs96, sbs96_rc, ref_sbs96, ref_sbs96_rc)
        sbs96_to_sbs52_lookup[sigprofiler_sbs96] = sbs52
    return sbs96_to_sbs52_lookup


def is_pyr_sbs96(sbs96: SBS96) -> bool:
    if sbs96.ref in PYRS:
        return True
    return False


def write_sbs96_to_sbs52_lookup_table() -> None:
    # get_sbs96_to_sbs52_lookup()
    sbs96_to_sbs52_lookup = get_sbs96_to_sbs52_lookup()
    with open("sbs96_to_sbs52_lookup_table.tsv", "w") as outfile:
        outfile.write("SBS96\tSBS52\n")
        for sbs96 in SBS96_CLASSIFICATIONS:
            outfile.write(f"{sbs96}\t{sbs96_to_sbs52_lookup[sbs96]}\n")


write_sbs96_to_sbs52_lookup_table()
