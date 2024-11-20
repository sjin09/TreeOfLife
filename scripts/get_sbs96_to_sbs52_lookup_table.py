from typing import Dict, List

from natsort import natsorted

BASE_COMPLEMENT_LOOKUP = str.maketrans("ACGT", "TGCA")
NTS = ["A", "C", "G", "T"]
PYR = set(["T", "C"])
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SBS96_CLASSIFICATION = [f"{nti}[{sub}]{ntj}" for sub in SUBS for nti in NTS for ntj in NTS]


def get_sbs96_reverse_complement(sbs96: str):
    ubase = sbs96[0]
    ref = sbs96[2]
    alt = sbs96[4]
    dbase = sbs96[6]
    fwd_ref_tri = f"{ubase}{ref}{dbase}"
    fwd_alt_tri = f"{ubase}{alt}{dbase}"
    rev_ref_tri = get_reverse_complementary_sequence(fwd_alt_tri)
    rev_alt_tri = get_reverse_complementary_sequence(fwd_ref_tri)
    rc_ubase, rc_ref, rc_dbase = list(rev_ref_tri)
    rc_alt = rev_alt_tri[1]
    sbs96_rc = f"{rc_ubase}[{rc_ref}>{rc_alt}]{rc_dbase}"
    return sbs96_rc


def get_complementary_base(base: str):
    base_complement = base.translate(BASE_COMPLEMENT_LOOKUP)
    return base_complement


def get_reverse_complementary_sequence(seq: str) -> str:
    rc_seq = seq[::-1].translate(BASE_COMPLEMENT_LOOKUP)
    return rc_seq


def get_sbs96_pyrimidine_count(sbs96: str):
    tri = f"{sbs96[0]}{sbs96[2]}{sbs96[6]}"
    pyr_count = sum(1 for nti in tri if nti in PYR)
    return pyr_count


def get_max_pyrimidine_sbs96(sbs96_classifications: List[str]) -> List[str]:
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


def get_sbs52(sbs96_classifications: List[str]) -> str:
    """
    Determine the SBS52 classification from a list of SBS96 classifications.

    This function picks the SBS96 category with the most pyrimidine bases.
    If there is a tie, it chooses the one that comes first alphabetically.
    """
    max_pyr_sbs96_classifications = get_max_pyrimidine_sbs96(sbs96_classifications)
    if len(max_pyr_sbs96_classifications) == 1:
        sbs52 = max_pyr_sbs96_classifications[0]
    else:
        sorted_sbs52_classifications = get_lexicographically_sorted_sbs52(sbs96_classifications)
        sbs52 = sorted_sbs52_classifications[0]
    return sbs52


def get_sbs96_to_sbs52_classifications() -> Dict[str, str]:
    sbs96_to_potential_sbs52_classifications = {}
    for sbs96 in SBS96_CLASSIFICATION:
        sbs96_rc = get_sbs96_reverse_complement(sbs96)
        sbs96_to_potential_sbs52_classifications[sbs96] = [sbs96, sbs96_rc]

    sbs96_to_sbs52 = {}
    for (sbs96, sbs52_classifications) in sbs96_to_potential_sbs52_classifications.items():
        sbs52 = get_sbs52(sbs52_classifications)
        sbs96_to_sbs52[sbs96] = sbs52
    return sbs96_to_sbs52


def write_sbs96_to_sbs52_lookup_table() -> None:
    sbs96_to_sbs52 = get_sbs96_to_sbs52_classifications()
    with open("sbs96_to_sbs52_lookup_table.tsv", "w") as outfile:
        outfile.write("SBS96\tSBS52\n")
        for (sbs96, sbs52) in sbs96_to_sbs52.items():
            outfile.write(f"{sbs96}\t{sbs52}\n")


write_sbs96_to_sbs52_lookup_table()
