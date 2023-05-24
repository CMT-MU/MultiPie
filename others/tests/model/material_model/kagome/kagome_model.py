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
kagome = {
    "info": {
        "model": "kagome",
        "molecule": False,
        "group": ("C3i^1", "space group No. 147 : C3i^1 / P-3 : PG C3i"),
        "crystal": "trigonal",
        "cell": {"a": 1.0, "b": 1.0, "c": 1.0, "alpha": 90.0, "beta": 90.0, "gamma": 120.0},
        "volume": 0.8660254037844388,
        "a1": "[1.0, 0.0, 0.0]",
        "a2": "[-0.5, 0.86602540378444, 0.0]",
        "a3": "[0.0, 0.0, 1.0]",
        "option": {"view": [0, 0, 1], "view_mode": "standard", "output": "kagome", "minimal_samb": True},
        "generate": {"model_type": "tight_binding", "time_reversal_type": "electric", "irrep": ["Ag"]},
        "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
        "k_path": "Γ-X",
        "dimension": 24,
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
        ],
        "ket_site": ["A_1", "A_2", "A_3"],
        "site": {"A": ("[1/2,0,0]", "s p")},
        "rep_site": {
            "A": ("[1/2, 0, 0]", "3e", [["(s,U)", "(s,D)"], ["(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"]])
        },
        "cell_site": {"A_1": ("[1/2, 0, 0]", "[1,4]"), "A_2": ("[0, 1/2, 0]", "[2,5]"), "A_3": ("[1/2, 1/2, 0]", "[3,6]")},
        "bond": [("A", "A", 1)],
        "rep_bond": {"A:A:1:1": ("[1/2, 0, 0];[1, 1/2, 0]", "6g", "D", 1)},
        "cell_bond": {
            "A:A:1:1_1": ("[1/2, 1/2, 0]@[3/4, 1/4, 0]", "[1]"),
            "A:A:1:1_2": ("[-1/2, 0, 0]@[3/4, 1/2, 0]", "[2]"),
            "A:A:1:1_3": ("[0, 1/2, 0]@[1/2, 1/4, 0]", "[-3]"),
            "A:A:1:1_4": ("[-1/2, -1/2, 0]@[1/4, 3/4, 0]", "[4]"),
            "A:A:1:1_5": ("[1/2, 0, 0]@[1/4, 1/2, 0]", "[5]"),
            "A:A:1:1_6": ("[0, -1/2, 0]@[1/2, 3/4, 0]", "[-6]"),
        },
    },
    "name": {
        "alias": {"S_001": "A", "A": "S_001", "B_001": "A:A:1:1", "A:A:1:1": "B_001"},
        "site": {"site_001": ("A", 1), "site_002": ("A", 2), "site_003": ("A", 3)},
        "site_name": {"[1/2, 0, 0]": ("site_001", 1), "[0, 1/2, 0]": ("site_002", 1), "[1/2, 1/2, 0]": ("site_003", 1)},
        "bond": {
            "bond_001": ("A:A:1:1", 1),
            "bond_002": ("A:A:1:1", 2),
            "bond_003": ("A:A:1:1", 3),
            "bond_004": ("A:A:1:1", 4),
            "bond_005": ("A:A:1:1", 5),
            "bond_006": ("A:A:1:1", 6),
        },
        "bond_name": {
            "[1/2, 1/2, 0]@[3/4, 1/4, 0]": ("bond_001", 1),
            "[-1/2, 0, 0]@[3/4, 1/2, 0]": ("bond_002", 1),
            "[0, 1/2, 0]@[1/2, 1/4, 0]": ("bond_003", 1),
            "[-1/2, -1/2, 0]@[1/4, 3/4, 0]": ("bond_004", 1),
            "[1/2, 0, 0]@[1/4, 1/2, 0]": ("bond_005", 1),
            "[0, -1/2, 0]@[1/2, 3/4, 0]": ("bond_006", 1),
        },
    },
    "data": {
        "plus_set": ["[0, 0, 0]"],
        "cluster_site": {"S_001": ["site_001", "site_002", "site_003"]},
        "cluster_bond": {"B_001": ["bond_001", "bond_002", "bond_003", "bond_004", "bond_005", "bond_006"]},
        "site": {
            "site_001": ("[1/2, 0, 0]", [0, 3], (0, 0)),
            "site_002": ("[0, 1/2, 0]", [1, 4], (1, 1)),
            "site_003": ("[1/2, 1/2, 0]", [2, 5], (2, 2)),
        },
        "bond": {
            "bond_001": ("[1/2, 1/2, 0]@[3/4, 1/4, 0]", [0], (0, 1), "[1/2, 1/2, 0]", "[1/2, 0, 0];[1, 1/2, 0]"),
            "bond_002": ("[-1/2, 0, 0]@[3/4, 1/2, 0]", [1], (1, 2), "[-1/2, 0, 0]", "[1, 1/2, 0];[1/2, 1/2, 0]"),
            "bond_003": ("[0, 1/2, 0]@[1/2, 1/4, 0]", [-2], (0, 2), "[0, 1/2, 0]", "[1/2, 0, 0];[1/2, 1/2, 0]"),
            "bond_004": ("[-1/2, -1/2, 0]@[1/4, 3/4, 0]", [3], (0, 1), "-bond_001", "[1/2, 1, 0];[0, 1/2, 0]"),
            "bond_005": ("[1/2, 0, 0]@[1/4, 1/2, 0]", [4], (1, 2), "-bond_002", "[0, 1/2, 0];[1/2, 1/2, 0]"),
            "bond_006": ("[0, -1/2, 0]@[1/2, 3/4, 0]", [-5], (0, 2), "-bond_003", "[1/2, 1, 0];[1/2, 1/2, 0]"),
        },
        "cluster_atomic": {
            (0, 0): [(0, 0, "M_001"), (0, 2, "M_002"), (2, 2, "M_003")],
            (1, 1): [(8, 8, "M_001"), (8, 10, "M_002"), (10, 10, "M_003")],
            (2, 2): [(16, 16, "M_001"), (16, 18, "M_002"), (18, 18, "M_003")],
            (0, 1): [(0, 8, "M_001"), (0, 10, "M_002"), (2, 8, "M_004"), (2, 10, "M_003")],
            (1, 2): [(8, 16, "M_001"), (8, 18, "M_002"), (10, 16, "M_004"), (10, 18, "M_003")],
            (0, 2): [(0, 16, "M_001"), (0, 18, "M_002"), (2, 16, "M_004"), (2, 18, "M_003")],
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
                {"A:A:1:1": ("[1/2, 0, 0];[1, 1/2, 0]", "6g", "D", 1)},
                {"A:A:2:1": ("[1/2, 0, 0];[0, 1/2, 0]", "6g", "D", 2)},
                {
                    "A:A:3:1": ("[-1/2, 0, 0];[1/2, 0, 0]", "1a", "ND", 3),
                    "A:A:3:2": ("[-1/2, 0, 0];[1/2, 1, 0]", "3e", "ND", 3),
                    "A:A:3:3": ("[1/2, 0, 0];[1/2, 1, 0]", "3e", "ND", 3),
                    "A:A:3:4": ("[1/2, 0, 0];[1/2, 0, 1]", "3f", "ND", 3),
                },
                {"A:A:4:1": ("[1/2, 0, 1];[1, 1/2, 0]", "6g", "D", 4), "A:A:4:2": ("[1/2, 0, 0];[1, 1/2, 1]", "6g", "D", 4)},
                {
                    "A:A:5:1": ("[1/2, 0, 0];[1, 3/2, 0]", "6g", "D", 5),
                    "A:A:5:2": ("[1/2, 0, 1];[0, 1/2, 0]", "6g", "D", 5),
                    "A:A:5:3": ("[-1/2, 0, 0];[1, 1/2, 0]", "6g", "D", 5),
                    "A:A:5:4": ("[1/2, 0, 0];[0, 1/2, 1]", "6g", "D", 5),
                },
                {
                    "A:A:6:1": ("[-1/2, 0, 1];[1/2, 0, 0]", "1b", "ND", 6),
                    "A:A:6:2": ("[1/2, 0, 0];[1/2, 1, 1]", "3f", "ND", 6),
                    "A:A:6:3": ("[-1/2, 0, 1];[1/2, 1, 0]", "3f", "ND", 6),
                    "A:A:6:4": ("[-1/2, 0, 0];[1/2, 1, 1]", "3f", "ND", 6),
                    "A:A:6:5": ("[1/2, 0, 1];[1/2, 1, 0]", "3f", "ND", 6),
                    "A:A:6:6": ("[-1/2, 0, 0];[1/2, 0, 1]", "1b", "ND", 6),
                },
                {"A:A:7:1": ("[-1/2, 0, 0];[1, 3/2, 0]", "6g", "D", 7)},
                {
                    "A:A:8:1": ("[-1/2, 0, 1];[1, 1/2, 0]", "6g", "D", 8),
                    "A:A:8:2": ("[-1/2, 0, 0];[1, 1/2, 1]", "6g", "D", 8),
                    "A:A:8:3": ("[1/2, 0, 1];[1, 3/2, 0]", "6g", "D", 8),
                    "A:A:8:4": ("[1/2, 0, 0];[1, 3/2, 1]", "6g", "D", 8),
                },
                {
                    "A:A:9:1": ("[-1/2, 0, 0];[3/2, 1, 0]", "3e", "ND", 9),
                    "A:A:9:2": ("[-1/2, 1, 0];[1/2, 0, 0]", "3e", "ND", 9),
                    "A:A:9:3": ("[-1/2, -1, 0];[1/2, 1, 0]", "1a", "ND", 9),
                },
            ]
        },
        "cell_range": (-2, 3, -2, 3, -2, 3),
        "max_neighbor": 10,
        "A": "[[1.0, -0.5, 0.0], [0.0, 0.86602540378444, 0.0], [0.0, 0.0, 1.0]]",
        "version": "1.1.1",
    },
}
