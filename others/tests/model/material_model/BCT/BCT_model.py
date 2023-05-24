"""
=== Molecule or Crystal Model (* only for crystal) ===
- info
    - model : model name.
    - molecule : molecule or crystal ?
    - group : (tag, detailed str)
    - crystal : crystal class
    - cell* : {a, b, c, alpha, beta, gamma}
    - volume* : unit cell volume
    - a1* : unit cell a1 vector
    - a2* : unit cell a2 vector
    - a3* : unit cell a3 vector
    - option :
        - view : view index
        - view_mode : QtDraw mode, standard/arrow/debug
        - output : output directory.
        - minimal_samb : output minimal SAMB ?
    - generate :
        - model_type : tight_binding/phonon
        - time_reversal_type : electric/magnetic/both
        - irrep : irrep list
    - k_point* : representative k points
    - k_path* : high-symmetry line in k space
    - dimension : dimension of full matrix
    - spinful : spinful or not
    - orbital : list of orbitals (U,D: up/down spin)
    - ket : ket basis list, orbital@site
    - ket_site : list of sites
    - site : input for "site" { name: (position, orbitals) }
    - rep_site : representative site { name: (position, wp, orbitals) }
    - cell_site : { name_idx(pset): (position, SOs) }
    - bond : input for "bond" [ (tail, head, neighbors) ]
    - rep_bond : representative bond { name: (vector@center, wp, directional, neighbor) }
    - cell_bond : { name_idx(pset): (vector@center, SOs) }

- name
    - alias : { cluster_tag: name or name: cluster_tag }
    - site : { site_tag: (name, no) }
    - site_name : { position : (site_tag, pset) }
    - bond : { bond_tag: (tail:head:n:multiplicity, no) }
    - bond_name : { vector@center : (bond_tag, pset) }

- data
    - plus_set* : [ plus_set list ]
    - cluster_site : { cluster_tag: site_list }
    - cluster_bond : { cluster_tag: bond_list }
    - site : { site_tag: (position, SO, (bra_site_no, ket_site_no)) }
    - bond : { bond_tag: (vector@center, SO, (bra_site_no, ket_site_no), vector, tail;head) }
    - cluster_atomic : { (bra_site_no, ket_site_no): [(bra_no, ket_no, matrix_tag)] }
    - atomic_braket : { matrix_tag : (bra_list, ket_list) }

- detail
    - rep_bond_all : { tail_head: [rep_bond in 0-9th neighbors] }
    - cell_range* : search range for bonds
    - max_neighbor : max. of neighbors to search
    - A* : transform matrix, [a1,a2,a3]
    - version : MultiPie version
"""
BCT = {
    "info": {
        "model": "BCT",
        "molecule": False,
        "group": ("D4h^17", "space group No. 139 : D4h^17 / I4/mmm : PG D4h"),
        "crystal": "tetragonal",
        "cell": {"a": 1.0, "b": 1.0, "c": 2.32, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
        "volume": 2.32,
        "a1": "[1.0, 0.0, 0.0]",
        "a2": "[0.0, 1.0, 0.0]",
        "a3": "[0.0, 0.0, 2.32]",
        "option": {"view": None, "view_mode": "standard", "output": "BCT", "minimal_samb": True},
        "generate": {"model_type": "tight_binding", "time_reversal_type": "electric", "irrep": ["A1g"]},
        "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
        "k_path": "Γ-X",
        "dimension": 6,
        "spinful": True,
        "orbital": ["(s,U)", "(s,D)", "(px,U)", "(px,D)", "(py,U)", "(py,D)"],
        "ket": ["(s,U)@A_1", "(s,D)@A_1", "(px,U)@A_1", "(px,D)@A_1", "(py,U)@A_1", "(py,D)@A_1"],
        "ket_site": ["A_1"],
        "site": {"A": ("[0,0,0]", ["s", "px", "py"])},
        "rep_site": {"A": ("[0, 0, 0]", "2a", [["(s,U)", "(s,D)"], ["(px,U)", "(px,D)", "(py,U)", "(py,D)"]])},
        "cell_site": {
            "A_1(1)": ("[0, 0, 0]", "[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]"),
            "A_1(2)": ("[1/2, 1/2, 1/2]", "[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]"),
        },
        "bond": [("A", "A", [1, 2, 7])],
        "rep_bond": {
            "A:A:1:1": ("[-1/2, 1/2, 1/2];[1/2, 1/2, 1/2]", "4c", "ND", 1),
            "A:A:2:1": ("[0, 0, 0];[1/2, 1/2, 1/2]", "8f", "ND", 2),
            "A:A:7:1": ("[0, 0, 0];[0, 0, 1]", "2b", "ND", 7),
        },
        "cell_bond": {
            "A:A:1:1_1(1)": ("[1, 0, 0]@[1/2, 0, 0]", "[1,-2,3,-4,-9,10,-11,12]"),
            "A:A:1:1_1(2)": ("[1, 0, 0]@[0, 1/2, 1/2]", "[1,-2,3,-4,-9,10,-11,12]"),
            "A:A:1:1_2(1)": ("[0, 1, 0]@[0, 1/2, 0]", "[5,-6,7,-8,-13,14,-15,16]"),
            "A:A:1:1_2(2)": ("[0, 1, 0]@[1/2, 0, 1/2]", "[5,-6,7,-8,-13,14,-15,16]"),
            "A:A:2:1_1(1)": ("[1/2, 1/2, 1/2]@[1/4, 1/4, 1/4]", "[1,-6,-9,14]"),
            "A:A:2:1_1(2)": ("[1/2, 1/2, 1/2]@[3/4, 3/4, 3/4]", "[1,-6,-9,14]"),
            "A:A:2:1_2(1)": ("[1/2, 1/2, -1/2]@[1/4, 1/4, 3/4]", "[-2,5,10,-13]"),
            "A:A:2:1_2(2)": ("[1/2, 1/2, -1/2]@[3/4, 3/4, 1/4]", "[-2,5,10,-13]"),
            "A:A:2:1_3(1)": ("[-1/2, 1/2, 1/2]@[3/4, 1/4, 1/4]", "[-3,7,11,-15]"),
            "A:A:2:1_3(2)": ("[-1/2, 1/2, 1/2]@[1/4, 3/4, 3/4]", "[-3,7,11,-15]"),
            "A:A:2:1_4(1)": ("[1/2, -1/2, 1/2]@[1/4, 3/4, 1/4]", "[-4,8,12,-16]"),
            "A:A:2:1_4(2)": ("[1/2, -1/2, 1/2]@[3/4, 1/4, 3/4]", "[-4,8,12,-16]"),
            "A:A:7:1_1(1)": ("[0, 0, 1]@[0, 0, 1/2]", "[1,2,-3,-4,-5,-6,7,8,-9,-10,11,12,13,14,-15,-16]"),
            "A:A:7:1_1(2)": ("[0, 0, 1]@[1/2, 1/2, 0]", "[1,2,-3,-4,-5,-6,7,8,-9,-10,11,12,13,14,-15,-16]"),
        },
    },
    "name": {
        "alias": {
            "S_001": "A",
            "A": "S_001",
            "B_001": "A:A:1:1",
            "A:A:1:1": "B_001",
            "B_002": "A:A:2:1",
            "A:A:2:1": "B_002",
            "B_003": "A:A:7:1",
            "A:A:7:1": "B_003",
        },
        "site": {"site_001": ("A", 1)},
        "site_name": {"[0, 0, 0]": ("site_001", 1), "[1/2, 1/2, 1/2]": ("site_001", 2)},
        "bond": {
            "bond_001": ("A:A:1:1", 1),
            "bond_002": ("A:A:1:1", 2),
            "bond_003": ("A:A:2:1", 1),
            "bond_004": ("A:A:2:1", 2),
            "bond_005": ("A:A:2:1", 3),
            "bond_006": ("A:A:2:1", 4),
            "bond_007": ("A:A:7:1", 1),
        },
        "bond_name": {
            "[1, 0, 0]@[1/2, 0, 0]": ("bond_001", 1),
            "[1, 0, 0]@[0, 1/2, 1/2]": ("bond_001", 2),
            "[0, 1, 0]@[0, 1/2, 0]": ("bond_002", 1),
            "[0, 1, 0]@[1/2, 0, 1/2]": ("bond_002", 2),
            "[1/2, 1/2, 1/2]@[1/4, 1/4, 1/4]": ("bond_003", 1),
            "[1/2, 1/2, 1/2]@[3/4, 3/4, 3/4]": ("bond_003", 2),
            "[1/2, 1/2, -1/2]@[1/4, 1/4, 3/4]": ("bond_004", 1),
            "[1/2, 1/2, -1/2]@[3/4, 3/4, 1/4]": ("bond_004", 2),
            "[-1/2, 1/2, 1/2]@[3/4, 1/4, 1/4]": ("bond_005", 1),
            "[-1/2, 1/2, 1/2]@[1/4, 3/4, 3/4]": ("bond_005", 2),
            "[1/2, -1/2, 1/2]@[1/4, 3/4, 1/4]": ("bond_006", 1),
            "[1/2, -1/2, 1/2]@[3/4, 1/4, 3/4]": ("bond_006", 2),
            "[0, 0, 1]@[0, 0, 1/2]": ("bond_007", 1),
            "[0, 0, 1]@[1/2, 1/2, 0]": ("bond_007", 2),
        },
    },
    "data": {
        "plus_set": ["[0, 0, 0]", "[1/2, 1/2, 1/2]"],
        "cluster_site": {"S_001": ["site_001"]},
        "cluster_bond": {
            "B_001": ["bond_001", "bond_002"],
            "B_002": ["bond_003", "bond_004", "bond_005", "bond_006"],
            "B_003": ["bond_007"],
        },
        "site": {"site_001": ("[0, 0, 0]", [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], (0, 0))},
        "bond": {
            "bond_001": ("[1, 0, 0]@[1/2, 0, 0]", [0, -1, 2, -3, -8, 9, -10, 11], (0, 0), "[1, 0, 0]", "[0, 0, 0];[1, 0, 0]"),
            "bond_002": ("[0, 1, 0]@[0, 1/2, 0]", [4, -5, 6, -7, -12, 13, -14, 15], (0, 0), "[0, 1, 0]", "[0, 0, 0];[0, 1, 0]"),
            "bond_003": (
                "[1/2, 1/2, 1/2]@[1/4, 1/4, 1/4]",
                [0, -5, -8, 13],
                (0, 0),
                "[1/2, 1/2, 1/2]",
                "[0, 0, 0];[1/2, 1/2, 1/2]",
            ),
            "bond_004": (
                "[1/2, 1/2, -1/2]@[1/4, 1/4, 3/4]",
                [-1, 4, 9, -12],
                (0, 0),
                "[1/2, 1/2, -1/2]",
                "[0, 0, 1];[1/2, 1/2, 1/2]",
            ),
            "bond_005": (
                "[-1/2, 1/2, 1/2]@[3/4, 1/4, 1/4]",
                [-2, 6, 10, -14],
                (0, 0),
                "[-1/2, 1/2, 1/2]",
                "[1, 0, 0];[1/2, 1/2, 1/2]",
            ),
            "bond_006": (
                "[1/2, -1/2, 1/2]@[1/4, 3/4, 1/4]",
                [-3, 7, 11, -15],
                (0, 0),
                "[1/2, -1/2, 1/2]",
                "[0, 1, 0];[1/2, 1/2, 1/2]",
            ),
            "bond_007": (
                "[0, 0, 1]@[0, 0, 1/2]",
                [0, 1, -2, -3, -4, -5, 6, 7, -8, -9, 10, 11, 12, 13, -14, -15],
                (0, 0),
                "[0, 0, 1]",
                "[0, 0, 0];[0, 0, 1]",
            ),
        },
        "cluster_atomic": {(0, 0): [(0, 0, "M_001"), (0, 2, "M_002"), (2, 2, "M_003")]},
        "atomic_braket": {
            "M_001": (["(s,U)", "(s,D)"], ["(s,U)", "(s,D)"]),
            "M_002": (["(s,U)", "(s,D)"], ["(px,U)", "(px,D)", "(py,U)", "(py,D)"]),
            "M_003": (["(px,U)", "(px,D)", "(py,U)", "(py,D)"], ["(px,U)", "(px,D)", "(py,U)", "(py,D)"]),
        },
    },
    "detail": {
        "rep_bond_all": {
            "A_A": [
                {},
                {"A:A:1:1": ("[-1/2, 1/2, 1/2];[1/2, 1/2, 1/2]", "4c", "ND", 1)},
                {"A:A:2:1": ("[0, 0, 0];[1/2, 1/2, 1/2]", "8f", "ND", 2)},
                {"A:A:3:1": ("[-1/2, -1/2, 1/2];[1/2, 1/2, 1/2]", "2b", "ND", 3)},
                {"A:A:4:1": ("[-1/2, 1/2, 1/2];[1, 0, 1]", "8f", "ND", 4)},
                {"A:A:5:1": ("[-1, 0, 0];[1, 0, 0]", "2a", "ND", 5)},
                {"A:A:6:1": ("[-1, 0, 0];[1, 1, 0]", "4c", "ND", 6)},
                {"A:A:7:1": ("[0, 0, 0];[0, 0, 1]", "2b", "ND", 7)},
                {"A:A:8:1": ("[-1/2, -1/2, 1/2];[1, 1, 0]", "8f", "ND", 8)},
                {"A:A:9:1": ("[-1/2, 1/2, -1/2];[1/2, 1/2, 1/2]", "4c", "ND", 9)},
            ]
        },
        "cell_range": (-2, 3, -2, 3, -2, 3),
        "max_neighbor": 10,
        "A": "[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 2.32]]",
        "version": "1.1.1",
    },
}
