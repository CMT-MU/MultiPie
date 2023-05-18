from multipie.multipole.util.atomic_orbital_util import (
    parse_orb_list,
    sort_orb_list,
    split_orb_list_rank_block,
    to_spinless,
    to_latex,
)


# ==================================================
def test_atomic_orbital_util():
    input = [  # [orbitals, spinful, crystal]
        # === LM basis, [(L,M)] ===#
        (["(0,0)", "(1,1)", "(1,0)", "(1,-1)"], False, "triclinic"),
        # === LMS basis, [(L,M,Sz)] ===#
        (["(0,0)", "(1,1)", "(1,0)", "(1,-1)"], True, "monoclinic"),
        # === abbreviation of LM basis, "L1 L2 ..." ===#
        ("0 1", False, "orthorhombic"),  # ["(0, 0)", "(1, 1)", "(1, 0)", "(1, -1)"]
        # === abbreviation of LMS basis, "L1 L2 ...", spinful=True ===#
        ("0 1", True, "tetragonal"),  # ["(0, 0)", "(1, 1)", "(1, 0)", "(1, -1)"]
        # === JM basis, [(J,M,L)] ===#
        (["(1/2,1/2,0)", "(1/2,-1/2,0)"], True, "trigonal"),
        ("(1/2,0) (7/2,3)", True, "monoclinic"),
        # === abbreviation of JM basis,  "(J1,L1) (J2,L2) ..." ===#
        (
            "(1/2,0) (3/2,1)",
            True,
            "hexagonal",
        ),  # ["(1/2,1/2,0)","(1/2,-1/2,0)", "(3/2,3/2,1)","(3/2,1/2,1)", "(3/2,-1/2,1)","(3/2,-3/2,1)"]
        # === cubic === #
        (["s", "dyz", "dzx", "dxy", "fax", "fay", "faz"], False, "cubic"),
        (["s", "dyz", "dzx", "dxy", "fax", "fay", "faz"], True, "tetragonal"),
        # === hexagonal === #
        (["s", "px", "py", "dyz", "dzx", "f1", "f2"], False, "hexagonal"),
        (["s", "px", "py", "dyz", "dzx", "f1", "f2"], True, "trigonal"),
        # === monoclinic === #
        (["px", "py", "pz"], False, "monoclinic"),
        # Error
        # (["pz", "px", "py"], False, "monoclinic"),
        # === tetragonal === #
        (["px", "py"], False, "tetragonal"),
        # === trigonal === #
        (["dv", "dxy"], False, "trigonal"),
        # === abbreviation of cubic/hexagonal, "s p d f"=== #
        ("s p", False, "tetragonal"),
        ("s p d f", False, "trigonal"),
        ("s p", True, "tetragonal"),
        ("d s f f", True, "trigonal"),
    ]

    for orb_list, spinful, crystal in input:
        print(f"input orb_list = {orb_list}, spinful = {spinful}, crystal = {crystal}")
        orb_list = parse_orb_list(orb_list, spinful, crystal)
        print(f"orb_list = {orb_list}")
        print(f"to_spinless = {to_spinless(orb_list)}")
        orb_list = sort_orb_list(orb_list, spinful, crystal)
        print(f"orb_list (sorted) = {orb_list}")
        print(f"orb_list (latex) = {to_latex(orb_list, spinful, crystal)}")
        orb_list = split_orb_list_rank_block(orb_list, spinful, crystal)
        print(f"orb_list (splitted) = {orb_list}")
        print()


# ================================================== test
test_atomic_orbital_util()
