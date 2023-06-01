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
Oh1 = {
    "info": {
        "model": "Oh1",
        "molecule": False,
        "group": ("Oh^1", "space group No. 221 : Oh^1 / Pm-3m : PG Oh"),
        "crystal": "cubic",
        "cell": {"a": 1.0, "b": 1.0, "c": 1.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
        "volume": 1.0,
        "a1": "[1.0, 0.0, 0.0]",
        "a2": "[0.0, 1.0, 0.0]",
        "a3": "[0.0, 0.0, 1.0]",
        "option": {"view": None, "view_mode": "standard", "output": "Oh1", "minimal_samb": True},
        "generate": {"model_type": "tight_binding", "time_reversal_type": "electric", "irrep": ["A1g"]},
        "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
        "k_path": "Γ-X",
        "dimension": 8,
        "spinful": True,
        "orbital": ["(s,U)", "(s,D)", "(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"],
        "ket": ["(s,U)@A_1", "(s,D)@A_1", "(px,U)@A_1", "(px,D)@A_1", "(py,U)@A_1", "(py,D)@A_1", "(pz,U)@A_1", "(pz,D)@A_1"],
        "ket_site": ["A_1"],
        "site": {"A": ("[0,0,0]", "s p")},
        "rep_site": {
            "A": ("[0, 0, 0]", "1a", [["(s,U)", "(s,D)"], ["(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"]])
        },
        "cell_site": {
            "A_1": (
                "[0, 0, 0]",
                "[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48]",
            )
        },
        "bond": [("A", "A", [1, 2])],
        "rep_bond": {"A:A:1:1": ("[0, 0, 0];[0, 0, 1]", "3d", "ND", 1), "A:A:2:1": ("[0, 0, 0];[0, 1, 1]", "3c", "ND", 2)},
        "cell_bond": {
            "A:A:1:1_1": ("[0, 0, 1]@[0, 0, 1/2]", "[1,2,-3,-4,-5,-8,19,22,-25,-26,27,28,29,32,-43,-46]"),
            "A:A:1:1_2": ("[1, 0, 0]@[1/2, 0, 0]", "[6,-9,11,-12,13,-14,21,-24,-30,33,-35,36,-37,38,-45,48]"),
            "A:A:1:1_3": ("[0, 1, 0]@[0, 1/2, 0]", "[7,-10,15,16,-17,-18,-20,23,-31,34,-39,-40,41,42,44,-47]"),
            "A:A:2:1_1": ("[0, 1, 1]@[0, 1/2, 1/2]", "[1,-3,7,-10,-25,27,-31,34]"),
            "A:A:2:1_2": ("[0, 1, -1]@[0, 1/2, 1/2]", "[-2,4,-20,23,26,-28,44,-47]"),
            "A:A:2:1_3": ("[1, 0, -1]@[1/2, 0, 1/2]", "[5,-12,13,-19,-29,36,-37,43]"),
            "A:A:2:1_4": ("[1, -1, 0]@[1/2, 1/2, 0]", "[6,-16,18,-24,-30,40,-42,48]"),
            "A:A:2:1_5": ("[1, 0, 1]@[1/2, 0, 1/2]", "[-8,11,-14,22,32,-35,38,-46]"),
            "A:A:2:1_6": ("[1, 1, 0]@[1/2, 1/2, 0]", "[-9,15,-17,21,33,-39,41,-45]"),
        },
    },
    "name": {
        "alias": {"S_001": "A", "A": "S_001", "B_001": "A:A:1:1", "A:A:1:1": "B_001", "B_002": "A:A:2:1", "A:A:2:1": "B_002"},
        "site": {"site_001": ("A", 1)},
        "site_name": {"[0, 0, 0]": ("site_001", 1)},
        "bond": {
            "bond_001": ("A:A:1:1", 1),
            "bond_002": ("A:A:1:1", 2),
            "bond_003": ("A:A:1:1", 3),
            "bond_004": ("A:A:2:1", 1),
            "bond_005": ("A:A:2:1", 2),
            "bond_006": ("A:A:2:1", 3),
            "bond_007": ("A:A:2:1", 4),
            "bond_008": ("A:A:2:1", 5),
            "bond_009": ("A:A:2:1", 6),
        },
        "bond_name": {
            "[0, 0, 1]@[0, 0, 1/2]": ("bond_001", 1),
            "[1, 0, 0]@[1/2, 0, 0]": ("bond_002", 1),
            "[0, 1, 0]@[0, 1/2, 0]": ("bond_003", 1),
            "[0, 1, 1]@[0, 1/2, 1/2]": ("bond_004", 1),
            "[0, 1, -1]@[0, 1/2, 1/2]": ("bond_005", 1),
            "[1, 0, -1]@[1/2, 0, 1/2]": ("bond_006", 1),
            "[1, -1, 0]@[1/2, 1/2, 0]": ("bond_007", 1),
            "[1, 0, 1]@[1/2, 0, 1/2]": ("bond_008", 1),
            "[1, 1, 0]@[1/2, 1/2, 0]": ("bond_009", 1),
        },
    },
    "data": {
        "plus_set": ["[0, 0, 0]"],
        "cluster_site": {"S_001": ["site_001"]},
        "cluster_bond": {
            "B_001": ["bond_001", "bond_002", "bond_003"],
            "B_002": ["bond_004", "bond_005", "bond_006", "bond_007", "bond_008", "bond_009"],
        },
        "site": {
            "site_001": (
                "[0, 0, 0]",
                [
                    0,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    16,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                    25,
                    26,
                    27,
                    28,
                    29,
                    30,
                    31,
                    32,
                    33,
                    34,
                    35,
                    36,
                    37,
                    38,
                    39,
                    40,
                    41,
                    42,
                    43,
                    44,
                    45,
                    46,
                    47,
                ],
                (0, 0),
            )
        },
        "bond": {
            "bond_001": (
                "[0, 0, 1]@[0, 0, 1/2]",
                [0, 1, -2, -3, -4, -7, 18, 21, -24, -25, 26, 27, 28, 31, -42, -45],
                (0, 0),
                "[0, 0, 1]",
                "[0, 0, 0];[0, 0, 1]",
            ),
            "bond_002": (
                "[1, 0, 0]@[1/2, 0, 0]",
                [5, -8, 10, -11, 12, -13, 20, -23, -29, 32, -34, 35, -36, 37, -44, 47],
                (0, 0),
                "[1, 0, 0]",
                "[0, 0, 0];[1, 0, 0]",
            ),
            "bond_003": (
                "[0, 1, 0]@[0, 1/2, 0]",
                [6, -9, 14, 15, -16, -17, -19, 22, -30, 33, -38, -39, 40, 41, 43, -46],
                (0, 0),
                "[0, 1, 0]",
                "[0, 0, 0];[0, 1, 0]",
            ),
            "bond_004": ("[0, 1, 1]@[0, 1/2, 1/2]", [0, -2, 6, -9, -24, 26, -30, 33], (0, 0), "[0, 1, 1]", "[0, 0, 0];[0, 1, 1]"),
            "bond_005": (
                "[0, 1, -1]@[0, 1/2, 1/2]",
                [-1, 3, -19, 22, 25, -27, 43, -46],
                (0, 0),
                "[0, 1, -1]",
                "[0, 0, 1];[0, 1, 0]",
            ),
            "bond_006": (
                "[1, 0, -1]@[1/2, 0, 1/2]",
                [4, -11, 12, -18, -28, 35, -36, 42],
                (0, 0),
                "[1, 0, -1]",
                "[0, 0, 1];[1, 0, 0]",
            ),
            "bond_007": (
                "[1, -1, 0]@[1/2, 1/2, 0]",
                [5, -15, 17, -23, -29, 39, -41, 47],
                (0, 0),
                "[1, -1, 0]",
                "[0, 1, 0];[1, 0, 0]",
            ),
            "bond_008": (
                "[1, 0, 1]@[1/2, 0, 1/2]",
                [-7, 10, -13, 21, 31, -34, 37, -45],
                (0, 0),
                "[1, 0, 1]",
                "[0, 0, 0];[1, 0, 1]",
            ),
            "bond_009": (
                "[1, 1, 0]@[1/2, 1/2, 0]",
                [-8, 14, -16, 20, 32, -38, 40, -44],
                (0, 0),
                "[1, 1, 0]",
                "[0, 0, 0];[1, 1, 0]",
            ),
        },
        "cluster_atomic": {(0, 0): [(0, 0, "M_001"), (0, 2, "M_002"), (2, 2, "M_003")]},
        "atomic_braket": {
            "M_001": (["(s,U)", "(s,D)"], ["(s,U)", "(s,D)"]),
            "M_002": (["(s,U)", "(s,D)"], ["(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"]),
            "M_003": (
                ["(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"],
                ["(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"],
            ),
        },
    },
    "detail": {
        "rep_bond_all": {
            "A_A": [
                {},
                {"A:A:1:1": ("[0, 0, 0];[0, 0, 1]", "3d", "ND", 1)},
                {"A:A:2:1": ("[0, 0, 0];[0, 1, 1]", "3c", "ND", 2)},
                {"A:A:3:1": ("[0, 0, 0];[1, 1, 1]", "1b", "ND", 3)},
                {"A:A:4:1": ("[-1, 0, 0];[1, 0, 0]", "1a", "ND", 4)},
                {"A:A:5:1": ("[-1, 0, 0];[1, 0, 1]", "3d", "ND", 5)},
                {"A:A:6:1": ("[-1, 0, 0];[1, 1, 1]", "3c", "ND", 6)},
                {"A:A:7:1": ("[-1, -1, 0];[1, 1, 0]", "1a", "ND", 7)},
                {"A:A:8:1": ("[-1, -1, 0];[1, 1, 1]", "3d", "ND", 8)},
                {"A:A:9:1": ("[-1, -1, -1];[1, 1, 1]", "1a", "ND", 9)},
            ]
        },
        "cell_range": (-2, 3, -2, 3, -2, 3),
        "max_neighbor": 10,
        "A": "[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]",
        "version": "1.1.2",
    },
}
