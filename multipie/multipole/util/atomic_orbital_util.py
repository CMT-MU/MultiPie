"""
This file provides utility functions for atomic orbital.
"""
import sympy as sp
from collections import Counter

from gcoreutils.string_util import remove_space
from multipie.multipole.util.spin_orbital_basis import (
    _standard_basis,
    _abbrev_basis,
)
from multipie.const import __def_dict__


# ==================================================
def _remove_space(orb_list):
    """
    remove space

    Args:
        orb_list (str/[str]): orbital.

    Returns:
        list: orbital list ?
    """
    if type(orb_list) == str:
        return [o for o in orb_list.split(" ") if o != ""]
    elif type(orb_list) == list:
        return [remove_space(str(o)) for o in orb_list]
    else:
        raise Exception(f"invalid type of orbital list = {type(orb_list)} is given.")


# ==================================================
def is_spinful(orb_list):
    """
    spinful ?

    Args:
        orb_list (str/[str]): orbital.

    Returns:
        bool: spinful ?
    """
    orb_list = _remove_space(orb_list)
    if all([len(o) > 2 and o[2] == "/" for o in orb_list]):  # jm basis
        return True
    elif all(["U" in o or "D" in o for o in orb_list]):  # spinful basis
        return True
    else:
        return False


# ==================================================
def to_spinful(orb_list):
    """
    convert from spinless orbital list to spinful orbital list.

    Args:
        orb_list (str/[str]): orbital list.

    Returns:
        list: spinful orbital list.
    """
    orb_list = _remove_space(orb_list)

    if all([is_spinful(o) for o in orb_list]):  # spinful basis
        return orb_list
    else:
        sf_orb_list = []
        if all([any([k in o for k in ("s", "p", "d", "f")]) for o in orb_list]):
            for o in orb_list:
                sf_orb_list.append(f"({o},U)")
                sf_orb_list.append(f"({o},D)")
        else:  # lm
            for o in orb_list:
                sf_orb_list.append(f"({o[1:-1]},U)")
                sf_orb_list.append(f"({o[1:-1]},D)")

        return sf_orb_list


# ==================================================
def to_spinless(orb_list):
    """
    convert from spinful orbital list to spinless orbital list.

    Args:
        orb_list (str/[str]): orbital list.

    Returns:
        list: spinless orbital list.
    """
    orb_list = _remove_space(orb_list)

    if all([len(o) > 2 and o[2] == "/" for o in orb_list]):  # jm
        return orb_list
    elif not all(["U" in o or "D" in o for o in orb_list]):  # spinless basis
        return orb_list
    else:
        sl_orb_list = []
        for o in orb_list:
            for b_type, sb in _standard_basis[1].items():
                if o not in sb:
                    continue
                idx = sb.index(o)
                idx = sp.S(idx) / 2 if idx % 2 == 0 else sp.S(idx - 1) / 2
                sl_orb_list.append(_standard_basis[0][b_type][idx])

        sl_orb_list = list(sorted(sl_orb_list, key=sl_orb_list.index))
        return sl_orb_list


# ==================================================
def basis_type(orb_list, crystal):
    """
    identify the type of basis.

    Args:
        orb_list (str/[str]): orbital list.
        crystal (str): seven crystal systems, triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic.

    Returns:
        (tuple/NoneType): type of basis (jm/lm/cubic/hexagonal).
    """
    orb_list = _remove_space(orb_list)
    orb_list = to_spinless(orb_list)

    if all(["/" in o for o in orb_list]):
        return "jm"
    elif all([any([k in o for k in ("s", "p", "d", "f")]) for o in orb_list]):  # cubic/hexagonal
        if crystal in ("triclinic", "monoclinic", "orthorhombic", "tetragonal", "cubic"):
            return "cubic"
        elif crystal in ("trigonal", "hexagonal"):
            return "hexagonal"
        else:
            raise Exception("invalid crystal = {crystal} is given.")
    else:
        return "lm"


# ==================================================
def to_standard_form(orb_list, crystal):
    """
    convert from an abbreviation form to a standard form.

    Args:
        orb_list (str/[str]): orbital list.
        b_type (str): type of basis, "lm/jm/cubic/hexagonal".

    Returns:
        list: orbital list in standard format.

    Note:
        ex) ["p"] -> ["px", "py", "pz"].
    """
    orb_list = _remove_space(orb_list)

    b_type = basis_type(orb_list, crystal)
    spinful = 1 if b_type == "jm" else 0

    standard_orb_list = []
    for o in orb_list:
        if o in _abbrev_basis[spinful][b_type]:
            if b_type in ("cubic", "hexagonal"):
                orbs = list(filter(lambda b: o in b, _standard_basis[0][b_type]))
            elif b_type == "lm":
                orbs = list(filter(lambda b: o == b[1], _standard_basis[0]["lm"]))
            elif b_type == "jm":
                orbs = list(filter(lambda b: o == "(" + b[1:4] + "," + b[-2] + ")", _standard_basis[1]["jm"]))
            else:
                raise Exception(f"invalid b_type = {b_type} is given.")
            standard_orb_list.extend(orbs)
        else:
            standard_orb_list.append(o)

    return standard_orb_list


# ==================================================
def invalid_orb_list(orb_list):
    """
    list of orbital list with the invalid forms included in given orbital list.

    Args:
        orb_list (str/[str]): orbital list.

    Returns:
        [str]: invalid orbital list.
    """
    orb_list = _remove_space(orb_list)

    orb_list_dataset = (
        sum(list(_standard_basis[0].values()), [])
        + sum(list(_standard_basis[1].values()), [])
        + sum(list(_abbrev_basis[0].values()), [])
        + sum(list(_abbrev_basis[1].values()), [])
    )

    orb_list_ = set(orb_list) & set(orb_list_dataset)
    invalid_orbs = set(orb_list) - orb_list_

    return invalid_orbs


# ==================================================
def invalid_f_orb_list(orb_list, crystal):
    """
    list of invalid f orbital list included in given orbital list.

    Args:
        orb_list (str/[str]): orbital list.
        crystal (str): seven crystal systems, triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic.

    Returns:
        [str]: invalid f orbital list.
    """
    orb_list = _remove_space(orb_list)

    # f orb_list specific to cubic and hexagonal basis.
    f_orb_dict = {
        "cubic": ("fax", "fay", "fbx", "fby"),
        "hexagonal": ("f1", "f2", "fx", "fy"),
    }
    f_orb_list = f_orb_dict["cubic"] + f_orb_dict["hexagonal"]

    # correspondence between the seven crystal systems and harmonics types.
    invalid_f_orbs = []
    for o in orb_list:
        if o in f_orb_list:
            b_type = "cubic" if o in f_orb_dict["cubic"] else "hexagonal"
            is_f_orb = o in f_orb_dict[b_type]
            is_parent_crystal = crystal in __def_dict__["parentgroup"][b_type]
            if is_f_orb and is_parent_crystal:
                continue

            invalid_f_orbs.append(o)

    return invalid_f_orbs


# ==================================================
def is_insufficient(orb_list, crystal):
    """
    check if same irrep pair is not included in given orbital list.

    Args:
        orb_list (str/[str]): orbital list.
        crystal (str): seven crystal systems, triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic.

    Returns:
        bool, str: same irrep pair is not included in given orbital list ?, orbit list that must be included at the same time.
    """
    orb_list = _remove_space(orb_list)

    pdf_orb_list = {
        "p": {o for o in orb_list if o[0] == "p"},
        "d": {o for o in orb_list if o[0] == "d"},
        "f": {o for o in orb_list if o[0] == "f"},
    }

    pdf_orb_dict = {
        "cubic": [
            ("px", "py", "pz"),
            ("du", "dv"),
            ("dyz", "dzx", "dxy"),
            ("fax", "fay", "faz"),
            ("fbx", "fby", "fbz"),
        ],
        "tetragonal": [("px", "py"), ("dyz", "dzx"), ("fax", "fay"), ("fbx", "fby")],
        "hexagonal": [("px", "py"), ("dxy", "dv"), ("fx", "fy"), ("fbz", "fxyz"), ("fax", "fay"), ("fbx", "fby")],
        "trigonal": [("px", "py"), ("dxy", "dv"), ("fx", "fy"), ("fbz", "fxyz"), ("fax", "fay"), ("fbx", "fby")],
        "triclinic": [],
        "monoclinic": [],
        "orthorhombic": [],
    }

    orbs_list = pdf_orb_dict[crystal]
    for orbs in orbs_list:
        o = pdf_orb_list[orbs[0][0]]
        if not o.isdisjoint(set(orbs)) and not o.issuperset(set(orbs)):
            return True, tuple(orbs)

    return False, ""


# ==================================================
def rank(o, spinful, crystal):
    """
    rank of given orbital.

    Args:
        o (str): orbital.
        spinful (bool): spinful ?
        crystal (str): seven crystal systems, triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic.

    Returns:
        int: rank.
    """
    b_type = basis_type(o, crystal)
    if b_type == "jm":
        return sp.sympify(o)[2]
    elif b_type == "lm":
        return sp.sympify(o)[0]
    elif b_type in ("cubic", "hexagonal"):
        spdf_rank = {"s": 0, "p": 1, "d": 2, "f": 3}
        return spdf_rank[o[1]] if spinful else spdf_rank[o[0]]
    else:
        raise Exception(f"invalid type of basis = {b_type} is given.")


# ==================================================
def sort_orb_list(orb_list, spinful, crystal):
    """
    sort given orbital list.

    Args:
        orb_list (str/[str]): orbital list.
        spinful (bool): spinful ?
        crystal (str): seven crystal systems, triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic.

    Returns:
        ([str]): sorted orbital list.
                ["px", "py", "s"] => ["s", "px", "py"]
    """
    if type(orb_list) == str:  # abbreviation form
        orb_list = [o for o in orb_list.split(" ") if o != ""]
    elif type(orb_list) == list:  # standard form
        orb_list = [remove_space(str(o)) for o in orb_list]
    else:
        raise Exception(f"invalid type of orb_list = {type(orb_list)} is given.")

    orb_list = to_spinless(orb_list)

    b_type = basis_type(orb_list, crystal)
    spinful_ = True if b_type == "jm" else False

    orb_list_data = _standard_basis[spinful_][b_type]

    max_l = 3
    orb_list_sorted = [o for l in range(max_l + 1) for o in orb_list if l == rank(o, spinful_, crystal)]
    orb_list_sorted = [o for o in orb_list_data if o in orb_list_sorted]

    if spinful:
        orb_list_sorted = to_spinful(orb_list_sorted)

    return orb_list_sorted


# ==================================================
def parse_orb_list(orb_list, spinful, crystal):
    """
    parse given orbital list.
    convert orbital list from an abbreviation form to a standard form.

    Args:
        orb_list (str/[str]): orbital list.
        spinful (bool): spinful ?
        crystal (str): seven crystal systems, triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic.

    Returns:
        ([str]): orbital list.
    """
    if type(orb_list) == str:  # abbreviation form
        orb_list = [o for o in orb_list.split(" ") if o != ""]
        is_abbrev = True
    elif type(orb_list) == list:  # standard form
        orb_list = [remove_space(str(o)) for o in orb_list]
        is_abbrev = False
    else:
        raise Exception(f"invalid type of orb_list = {type(orb_list)} is given.")

    # check if invalid form of orbital is included in the given orb_list.
    d = Counter(orb_list)
    is_duplicated = not all(list(d.values()))
    if is_duplicated:
        raise Exception(f"duplication is included in given orb_list, {d}")

    # check if invalid orb_list are included.
    invalid_orbs = invalid_orb_list(orb_list)
    is_invalid_orbs = len(invalid_orbs) != 0
    if is_invalid_orbs:
        raise Exception(f"invalid orb_list = {invalid_orbs} is given.")

    # if one of ("fax", "fay", "fbx", "fby", "f1", "f2", "fx", "fy") is given explicitly, check consistency with crystal.
    invalid_f_orbs = invalid_f_orb_list(orb_list, crystal)
    is_invalid_f_orbs = len(invalid_f_orbs) != 0
    if is_invalid_f_orbs:
        raise Exception(f"invalid f orb_list = {invalid_f_orbs} is given.")

    # check if same irrep pair is included in given orb_list.
    if crystal in ("triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic"):
        flag, orbs = is_insufficient(orb_list, crystal)
        if flag:
            raise Exception(f"orb_list must include same irrep pairs, {orbs}, but {orb_list} was fiven.")

    # convert from an abbreviation form to a standard form.
    if is_abbrev:
        orb_list = to_standard_form(orb_list, crystal)

    # convert from spinless to spinful basis set.
    if spinful:
        orb_list = to_spinful(orb_list)

    # sort orbital list.
    orb_list = sort_orb_list(orb_list, spinful, crystal)

    return orb_list


# ==================================================
def split_orb_list_rank_block(orb_list, spinful, crystal):
    """
    split given orbital list into rank block.

    Args:
        orb_list (str/[str]): orbital list.
        spinful (bool): spinful ?
        crystal (str): seven crystal systems, triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic.

    Returns:
        [[str]]: splitted orbital list.
                ["s","px","py"] => [["s"], ["px","py"]].
    """
    max_l = 3
    splitted_orb_list = [[o for o in orb_list if l == rank(o, spinful, crystal)] for l in range(max_l + 1)]
    splitted_orb_list = [lst for lst in splitted_orb_list if len(lst) > 0]

    return splitted_orb_list


# ==================================================
def to_latex(orb_list, spinful, crystal):
    """
    convert string to latex format.

    Args:
        orb_list (str/[str]): orbital list.
        spinful (bool): spinful ?
        crystal (str): seven crystal systems, triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic.

    Returns:
        list: orbital in latex format.
    """

    def f(orb):
        if spinful:
            orb, UD = orb[1:-1].split(",")
            if len(orb) > 1:
                head, sub = orb[0], orb[1:]
                orb = head + "_{" + sub + "}"
            orb = f"({orb},{UD})"
        else:
            if len(orb) > 1:
                head, sub = orb[0], orb[1:]
                orb = head + "_{" + sub + "}"
        return orb

    b_type = basis_type(orb_list, crystal)

    if b_type in "jm":
        return orb_list
    elif b_type in "lm":
        return orb_list
    else:
        if type(orb_list) == str:
            orb_list = f(orb_list)
        elif type(orb_list) == list:
            orb_list = [f(orb) for orb in orb_list]
        return orb_list
