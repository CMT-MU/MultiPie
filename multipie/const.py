__def_dict__ = {
    "subgroup": {
        "triclinic": "cubic",
        "monoclinic": "cubic",
        "orthorhombic": "cubic",
        "tetragonal": "cubic",
        "trigonal": "hexagonal",
        "hexagonal": "hexagonal",
        "cubic": "cubic",
    },  # subgroup.
    "parentgroup": {
        "cubic": ["triclinic", "monoclinic", "orthorhombic", "tetragonal", "cubic"],
        "hexagonal": ["trigonal", "hexagonal"],
    },  # parent group.
    "crystal": ["triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic"],  # crystal system.
    "head": ["Q", "G", "T", "M"],  # multipole type.
    "inversion": ["polar", "axial"],  # inversion type.
    "time_reversal": ["electric", "magnetic"],  # time-reversal type.
    "head_i": {"Q": "polar", "T": "polar", "M": "axial", "G": "axial"},  # parity of head.
    "head_t": {"Q": "electric", "T": "magnetic", "M": "magnetic", "G": "electric"},  # time-reversal property of head.
    "it_head": {("polar", "electric"): "Q", ("axial", "electric"): "G", ("polar", "magnetic"): "T", ("axial", "magnetic"): "M"},
    "head_tr": {"Q": "T", "T": "Q", "M": "G", "G": "M"},  # toggle time-reversal property.
    "head_ir": {"Q": "G", "T": "M", "M": "T", "G": "Q"},  # toggle parity.
    "head_harm": {"Q": "Q", "T": "Q", "M": "G", "G": "G"},  # conversion to harmonics head.
    "to_polar": {"Q": "Q", "T": "T", "M": "T", "G": "Q"},  # conversion to polar head.
    "irrep_dim": {"A": 1, "B": 1, "E": 2, "T": 3},  # dimension of irrep. (head letter).
    "response_head": {  # (rank,type) to response tensor head.
        (0, "s"): "C",
        (1, "s"): "C",
        (2, "s"): "S",  # symmetric part.
        (2, "a"): "A",  # anti-symmetric part.
        (3, "s"): "S",  # symmetric part.
        (3, "a"): "A",  # anti-symmetric part.
        (4, "sss"): "S",  # symmetric-symmetric (symmetric) part.
        (4, "ssa"): "Sb",  # symmetric-symmetric (anti-symmetric) part.
        (4, "aas"): "A",  # anti-symmetric-anti-symmetric (symmetric) part.
        (4, "aaa"): "Ab",  # anti-symmetric-anti-symmetric (anti-symmetric) part.
        (4, "sa"): "M",  # symmetric-anti-symmetric part.
        (4, "as"): "Mb",  # anti-symmetric-symmetric part.
    },
}


def is_electric(tag):
    return tag == __def_dict__["time_reversal"][0] or tag in ["Q", "G"]


def is_magnetic(tag):
    return tag == __def_dict__["time_reversal"][1] or tag in ["T", "M"]


def is_polar(tag):
    return tag == __def_dict__["inversion"][0] or tag in ["Q", "T"]


def is_axial(tag):
    return tag == __def_dict__["inversion"][1] or tag in ["G", "M"]
