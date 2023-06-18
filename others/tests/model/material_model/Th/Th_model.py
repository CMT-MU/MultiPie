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
    - rep_site : representative site { name: (position, wp, orbitals, site-symmetry) }
    - cell_site : { name_idx(pset): (position, SOs) }
    - bond : input for "bond" [ (tail, head, neighbors) ]
    - rep_bond : representative bond { name: (vector@center, wp, directional, neighbor, site-symmetry) }
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
Th = {
    "info": {
        "model": "Th",
        "molecule": True,
        "group": ("Th", "point group No. 29 : Th / m-3"),
        "crystal": "cubic",
        "option": {"view": None, "view_mode": "standard", "output": "Th", "minimal_samb": True},
        "generate": {"model_type": "tight_binding", "time_reversal_type": "electric", "irrep": ["Ag"]},
        "dimension": 64,
        "spinful": True,
        "orbital": ["(s,U)", "(s,D)", "(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"],
        "ket": [
            "(s,U)@A_1",
            "(s,D)@A_1",
            "(px,U)@A_1",
            "(px,D)@A_1",
            "(py,U)@A_1",
            "(py,D)@A_1",
            "(pz,U)@A_1",
            "(pz,D)@A_1",
            "(s,U)@A_2",
            "(s,D)@A_2",
            "(px,U)@A_2",
            "(px,D)@A_2",
            "(py,U)@A_2",
            "(py,D)@A_2",
            "(pz,U)@A_2",
            "(pz,D)@A_2",
            "(s,U)@A_3",
            "(s,D)@A_3",
            "(px,U)@A_3",
            "(px,D)@A_3",
            "(py,U)@A_3",
            "(py,D)@A_3",
            "(pz,U)@A_3",
            "(pz,D)@A_3",
            "(s,U)@A_4",
            "(s,D)@A_4",
            "(px,U)@A_4",
            "(px,D)@A_4",
            "(py,U)@A_4",
            "(py,D)@A_4",
            "(pz,U)@A_4",
            "(pz,D)@A_4",
            "(s,U)@A_5",
            "(s,D)@A_5",
            "(px,U)@A_5",
            "(px,D)@A_5",
            "(py,U)@A_5",
            "(py,D)@A_5",
            "(pz,U)@A_5",
            "(pz,D)@A_5",
            "(s,U)@A_6",
            "(s,D)@A_6",
            "(px,U)@A_6",
            "(px,D)@A_6",
            "(py,U)@A_6",
            "(py,D)@A_6",
            "(pz,U)@A_6",
            "(pz,D)@A_6",
            "(s,U)@A_7",
            "(s,D)@A_7",
            "(px,U)@A_7",
            "(px,D)@A_7",
            "(py,U)@A_7",
            "(py,D)@A_7",
            "(pz,U)@A_7",
            "(pz,D)@A_7",
            "(s,U)@A_8",
            "(s,D)@A_8",
            "(px,U)@A_8",
            "(px,D)@A_8",
            "(py,U)@A_8",
            "(py,D)@A_8",
            "(pz,U)@A_8",
            "(pz,D)@A_8",
        ],
        "ket_site": ["A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8"],
        "site": {"A": ("[1,1,1]", "s p")},
        "rep_site": {
            "A": ("[1, 1, 1]", "8b", [["(s,U)", "(s,D)"], ["(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"]], ".3.")
        },
        "cell_site": {
            "A_1": ("[1, 1, 1]", "[1,5,9]"),
            "A_2": ("[-1, -1, 1]", "[2,6,11]"),
            "A_3": ("[1, -1, -1]", "[3,7,12]"),
            "A_4": ("[-1, 1, -1]", "[4,8,10]"),
            "A_5": ("[-1, -1, -1]", "[13,17,21]"),
            "A_6": ("[1, 1, -1]", "[14,18,23]"),
            "A_7": ("[-1, 1, 1]", "[15,19,24]"),
            "A_8": ("[1, -1, 1]", "[16,20,22]"),
        },
        "bond": [("A", "A", 1)],
        "rep_bond": {"A:A:1:1": ("[1, 1, 1];[1, 1, -1]", "12c", "ND", 1, "m..")},
        "cell_bond": {
            "A:A:1:1_1": ("[0, 0, -2]@[1, 1, 0]", "[1,-14]"),
            "A:A:1:1_2": ("[0, 0, -2]@[-1, -1, 0]", "[2,-13]"),
            "A:A:1:1_3": ("[0, 0, 2]@[1, -1, 0]", "[3,-16]"),
            "A:A:1:1_4": ("[0, 0, 2]@[-1, 1, 0]", "[4,-15]"),
            "A:A:1:1_5": ("[-2, 0, 0]@[0, 1, 1]", "[5,-19]"),
            "A:A:1:1_6": ("[2, 0, 0]@[0, -1, 1]", "[6,-20]"),
            "A:A:1:1_7": ("[-2, 0, 0]@[0, -1, -1]", "[7,-17]"),
            "A:A:1:1_8": ("[2, 0, 0]@[0, 1, -1]", "[8,-18]"),
            "A:A:1:1_9": ("[0, -2, 0]@[1, 0, 1]", "[9,-22]"),
            "A:A:1:1_10": ("[0, -2, 0]@[-1, 0, -1]", "[10,-21]"),
            "A:A:1:1_11": ("[0, 2, 0]@[-1, 0, 1]", "[11,-24]"),
            "A:A:1:1_12": ("[0, 2, 0]@[1, 0, -1]", "[12,-23]"),
        },
    },
    "name": {
        "alias": {"S_001": "A", "A": "S_001", "B_001": "A:A:1:1", "A:A:1:1": "B_001"},
        "site": {
            "site_001": ("A", 1),
            "site_002": ("A", 2),
            "site_003": ("A", 3),
            "site_004": ("A", 4),
            "site_005": ("A", 5),
            "site_006": ("A", 6),
            "site_007": ("A", 7),
            "site_008": ("A", 8),
        },
        "site_name": {
            "[1, 1, 1]": ("site_001", 1),
            "[-1, -1, 1]": ("site_002", 1),
            "[1, -1, -1]": ("site_003", 1),
            "[-1, 1, -1]": ("site_004", 1),
            "[-1, -1, -1]": ("site_005", 1),
            "[1, 1, -1]": ("site_006", 1),
            "[-1, 1, 1]": ("site_007", 1),
            "[1, -1, 1]": ("site_008", 1),
        },
        "bond": {
            "bond_001": ("A:A:1:1", 1),
            "bond_002": ("A:A:1:1", 2),
            "bond_003": ("A:A:1:1", 3),
            "bond_004": ("A:A:1:1", 4),
            "bond_005": ("A:A:1:1", 5),
            "bond_006": ("A:A:1:1", 6),
            "bond_007": ("A:A:1:1", 7),
            "bond_008": ("A:A:1:1", 8),
            "bond_009": ("A:A:1:1", 9),
            "bond_010": ("A:A:1:1", 10),
            "bond_011": ("A:A:1:1", 11),
            "bond_012": ("A:A:1:1", 12),
        },
        "bond_name": {
            "[0, 0, -2]@[1, 1, 0]": ("bond_001", 1),
            "[0, 0, -2]@[-1, -1, 0]": ("bond_002", 1),
            "[0, 0, 2]@[1, -1, 0]": ("bond_003", 1),
            "[0, 0, 2]@[-1, 1, 0]": ("bond_004", 1),
            "[-2, 0, 0]@[0, 1, 1]": ("bond_005", 1),
            "[2, 0, 0]@[0, -1, 1]": ("bond_006", 1),
            "[-2, 0, 0]@[0, -1, -1]": ("bond_007", 1),
            "[2, 0, 0]@[0, 1, -1]": ("bond_008", 1),
            "[0, -2, 0]@[1, 0, 1]": ("bond_009", 1),
            "[0, -2, 0]@[-1, 0, -1]": ("bond_010", 1),
            "[0, 2, 0]@[-1, 0, 1]": ("bond_011", 1),
            "[0, 2, 0]@[1, 0, -1]": ("bond_012", 1),
        },
    },
    "data": {
        "cluster_site": {
            "S_001": ["site_001", "site_002", "site_003", "site_004", "site_005", "site_006", "site_007", "site_008"]
        },
        "cluster_bond": {
            "B_001": [
                "bond_001",
                "bond_002",
                "bond_003",
                "bond_004",
                "bond_005",
                "bond_006",
                "bond_007",
                "bond_008",
                "bond_009",
                "bond_010",
                "bond_011",
                "bond_012",
            ]
        },
        "site": {
            "site_001": ("[1, 1, 1]", [0, 4, 8], (0, 0)),
            "site_002": ("[-1, -1, 1]", [1, 5, 10], (1, 1)),
            "site_003": ("[1, -1, -1]", [2, 6, 11], (2, 2)),
            "site_004": ("[-1, 1, -1]", [3, 7, 9], (3, 3)),
            "site_005": ("[-1, -1, -1]", [12, 16, 20], (4, 4)),
            "site_006": ("[1, 1, -1]", [13, 17, 22], (5, 5)),
            "site_007": ("[-1, 1, 1]", [14, 18, 23], (6, 6)),
            "site_008": ("[1, -1, 1]", [15, 19, 21], (7, 7)),
        },
        "bond": {
            "bond_001": ("[0, 0, -2]@[1, 1, 0]", [0, -13], (5, 0), "[0, 0, -2]", "[1, 1, 1];[1, 1, -1]"),
            "bond_002": ("[0, 0, -2]@[-1, -1, 0]", [1, -12], (4, 1), "bond_001", "[-1, -1, 1];[-1, -1, -1]"),
            "bond_003": ("[0, 0, 2]@[1, -1, 0]", [2, -15], (7, 2), "-bond_001", "[1, -1, -1];[1, -1, 1]"),
            "bond_004": ("[0, 0, 2]@[-1, 1, 0]", [3, -14], (6, 3), "-bond_001", "[-1, 1, -1];[-1, 1, 1]"),
            "bond_005": ("[-2, 0, 0]@[0, 1, 1]", [4, -18], (6, 0), "[-2, 0, 0]", "[1, 1, 1];[-1, 1, 1]"),
            "bond_006": ("[2, 0, 0]@[0, -1, 1]", [5, -19], (7, 1), "-bond_005", "[-1, -1, 1];[1, -1, 1]"),
            "bond_007": ("[-2, 0, 0]@[0, -1, -1]", [6, -16], (4, 2), "bond_005", "[1, -1, -1];[-1, -1, -1]"),
            "bond_008": ("[2, 0, 0]@[0, 1, -1]", [7, -17], (5, 3), "-bond_005", "[-1, 1, -1];[1, 1, -1]"),
            "bond_009": ("[0, -2, 0]@[1, 0, 1]", [8, -21], (7, 0), "[0, -2, 0]", "[1, 1, 1];[1, -1, 1]"),
            "bond_010": ("[0, -2, 0]@[-1, 0, -1]", [9, -20], (4, 3), "bond_009", "[-1, 1, -1];[-1, -1, -1]"),
            "bond_011": ("[0, 2, 0]@[-1, 0, 1]", [10, -23], (6, 1), "-bond_009", "[-1, -1, 1];[-1, 1, 1]"),
            "bond_012": ("[0, 2, 0]@[1, 0, -1]", [11, -22], (5, 2), "-bond_009", "[1, -1, -1];[1, 1, -1]"),
        },
        "cluster_atomic": {
            (0, 0): [(0, 0, "M_001"), (0, 2, "M_002"), (2, 2, "M_003")],
            (1, 1): [(8, 8, "M_001"), (8, 10, "M_002"), (10, 10, "M_003")],
            (2, 2): [(16, 16, "M_001"), (16, 18, "M_002"), (18, 18, "M_003")],
            (3, 3): [(24, 24, "M_001"), (24, 26, "M_002"), (26, 26, "M_003")],
            (4, 4): [(32, 32, "M_001"), (32, 34, "M_002"), (34, 34, "M_003")],
            (5, 5): [(40, 40, "M_001"), (40, 42, "M_002"), (42, 42, "M_003")],
            (6, 6): [(48, 48, "M_001"), (48, 50, "M_002"), (50, 50, "M_003")],
            (7, 7): [(56, 56, "M_001"), (56, 58, "M_002"), (58, 58, "M_003")],
            (5, 0): [(40, 0, "M_001"), (40, 2, "M_002"), (42, 0, "M_004"), (42, 2, "M_003")],
            (4, 1): [(32, 8, "M_001"), (32, 10, "M_002"), (34, 8, "M_004"), (34, 10, "M_003")],
            (7, 2): [(56, 16, "M_001"), (56, 18, "M_002"), (58, 16, "M_004"), (58, 18, "M_003")],
            (6, 3): [(48, 24, "M_001"), (48, 26, "M_002"), (50, 24, "M_004"), (50, 26, "M_003")],
            (6, 0): [(48, 0, "M_001"), (48, 2, "M_002"), (50, 0, "M_004"), (50, 2, "M_003")],
            (7, 1): [(56, 8, "M_001"), (56, 10, "M_002"), (58, 8, "M_004"), (58, 10, "M_003")],
            (4, 2): [(32, 16, "M_001"), (32, 18, "M_002"), (34, 16, "M_004"), (34, 18, "M_003")],
            (5, 3): [(40, 24, "M_001"), (40, 26, "M_002"), (42, 24, "M_004"), (42, 26, "M_003")],
            (7, 0): [(56, 0, "M_001"), (56, 2, "M_002"), (58, 0, "M_004"), (58, 2, "M_003")],
            (4, 3): [(32, 24, "M_001"), (32, 26, "M_002"), (34, 24, "M_004"), (34, 26, "M_003")],
            (6, 1): [(48, 8, "M_001"), (48, 10, "M_002"), (50, 8, "M_004"), (50, 10, "M_003")],
            (5, 2): [(40, 16, "M_001"), (40, 18, "M_002"), (42, 16, "M_004"), (42, 18, "M_003")],
        },
        "atomic_braket": {
            "M_001": (["(s,U)", "(s,D)"], ["(s,U)", "(s,D)"]),
            "M_002": (["(s,U)", "(s,D)"], ["(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"]),
            "M_003": (
                ["(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"],
                ["(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"],
            ),
            "M_004": (["(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"], ["(s,U)", "(s,D)"]),
        },
    },
    "detail": {
        "rep_bond_all": {
            "A_A": [
                {},
                {"A:A:1:1": ("[1, 1, 1];[1, 1, -1]", "12c", "ND", 1, "m..")},
                {"A:A:2:1": ("[1, 1, 1];[-1, -1, 1]", "6a", "ND", 2, "2mm..")},
                {"A:A:3:1": ("[1, 1, 1];[-1, -1, -1]", "1o", "ND", 3, "m-3")},
            ]
        },
        "max_neighbor": 10,
        "version": "1.1.10",
    },
}
