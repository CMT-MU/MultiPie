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
D4h1 = {
    "info": {
        "model": "D4h1",
        "molecule": False,
        "group": ("D4h^1", "space group No. 123 : D4h^1 / P4/mmm : PG D4h"),
        "crystal": "tetragonal",
        "cell": {"a": 1.0, "b": 1.0, "c": 1.5, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
        "volume": 1.5,
        "a1": "[1.0, 0.0, 0.0]",
        "a2": "[0.0, 1.0, 0.0]",
        "a3": "[0.0, 0.0, 1.5]",
        "option": {"view": None, "view_mode": "standard", "output": "D4h1", "minimal_samb": True},
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
        "cell_site": {"A_1": ("[0, 0, 0]", "[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]")},
        "bond": [("A", "A", [1, 2, 3, 4])],
        "rep_bond": {
            "A:A:1:1": ("[0, 0, 0];[0, 1, 0]", "2f", "ND", 1),
            "A:A:2:1": ("[0, 0, 0];[1, 1, 0]", "1c", "ND", 2),
            "A:A:3:1": ("[0, 0, 0];[0, 0, 1]", "1b", "ND", 3),
            "A:A:4:1": ("[0, 0, 0];[0, 1, 1]", "2e", "ND", 4),
        },
        "cell_bond": {
            "A:A:1:1_1": ("[0, 1, 0]@[0, 1/2, 0]", "[1,-2,-3,4,-9,10,11,-12]"),
            "A:A:1:1_2": ("[1, 0, 0]@[1/2, 0, 0]", "[5,-6,-7,8,-13,14,15,-16]"),
            "A:A:2:1_1": ("[1, 1, 0]@[1/2, 1/2, 0]", "[1,-2,5,-6,-9,10,-13,14]"),
            "A:A:2:1_2": ("[1, -1, 0]@[1/2, 1/2, 0]", "[3,-4,-7,8,-11,12,15,-16]"),
            "A:A:3:1_1": ("[0, 0, 1]@[0, 0, 1/2]", "[1,2,-3,-4,-5,-6,7,8,-9,-10,11,12,13,14,-15,-16]"),
            "A:A:4:1_1": ("[0, 1, 1]@[0, 1/2, 1/2]", "[1,-3,-9,11]"),
            "A:A:4:1_2": ("[0, 1, -1]@[0, 1/2, 1/2]", "[-2,4,10,-12]"),
            "A:A:4:1_3": ("[1, 0, -1]@[1/2, 0, 1/2]", "[5,-7,-13,15]"),
            "A:A:4:1_4": ("[1, 0, 1]@[1/2, 0, 1/2]", "[-6,8,14,-16]"),
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
            "B_003": "A:A:3:1",
            "A:A:3:1": "B_003",
            "B_004": "A:A:4:1",
            "A:A:4:1": "B_004",
        },
        "site": {"site_001": ("A", 1)},
        "site_name": {"[0, 0, 0]": ("site_001", 1)},
        "bond": {
            "bond_001": ("A:A:1:1", 1),
            "bond_002": ("A:A:1:1", 2),
            "bond_003": ("A:A:2:1", 1),
            "bond_004": ("A:A:2:1", 2),
            "bond_005": ("A:A:3:1", 1),
            "bond_006": ("A:A:4:1", 1),
            "bond_007": ("A:A:4:1", 2),
            "bond_008": ("A:A:4:1", 3),
            "bond_009": ("A:A:4:1", 4),
        },
        "bond_name": {
            "[0, 1, 0]@[0, 1/2, 0]": ("bond_001", 1),
            "[1, 0, 0]@[1/2, 0, 0]": ("bond_002", 1),
            "[1, 1, 0]@[1/2, 1/2, 0]": ("bond_003", 1),
            "[1, -1, 0]@[1/2, 1/2, 0]": ("bond_004", 1),
            "[0, 0, 1]@[0, 0, 1/2]": ("bond_005", 1),
            "[0, 1, 1]@[0, 1/2, 1/2]": ("bond_006", 1),
            "[0, 1, -1]@[0, 1/2, 1/2]": ("bond_007", 1),
            "[1, 0, -1]@[1/2, 0, 1/2]": ("bond_008", 1),
            "[1, 0, 1]@[1/2, 0, 1/2]": ("bond_009", 1),
        },
    },
    "data": {
        "plus_set": ["[0, 0, 0]"],
        "cluster_site": {"S_001": ["site_001"]},
        "cluster_bond": {
            "B_001": ["bond_001", "bond_002"],
            "B_002": ["bond_003", "bond_004"],
            "B_003": ["bond_005"],
            "B_004": ["bond_006", "bond_007", "bond_008", "bond_009"],
        },
        "site": {"site_001": ("[0, 0, 0]", [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], (0, 0))},
        "bond": {
            "bond_001": ("[0, 1, 0]@[0, 1/2, 0]", [0, -1, -2, 3, -8, 9, 10, -11], (0, 0), "[0, 1, 0]", "[0, 0, 0];[0, 1, 0]"),
            "bond_002": ("[1, 0, 0]@[1/2, 0, 0]", [4, -5, -6, 7, -12, 13, 14, -15], (0, 0), "[1, 0, 0]", "[0, 0, 0];[1, 0, 0]"),
            "bond_003": ("[1, 1, 0]@[1/2, 1/2, 0]", [0, -1, 4, -5, -8, 9, -12, 13], (0, 0), "[1, 1, 0]", "[0, 0, 0];[1, 1, 0]"),
            "bond_004": (
                "[1, -1, 0]@[1/2, 1/2, 0]",
                [2, -3, -6, 7, -10, 11, 14, -15],
                (0, 0),
                "[1, -1, 0]",
                "[0, 1, 0];[1, 0, 0]",
            ),
            "bond_005": (
                "[0, 0, 1]@[0, 0, 1/2]",
                [0, 1, -2, -3, -4, -5, 6, 7, -8, -9, 10, 11, 12, 13, -14, -15],
                (0, 0),
                "[0, 0, 1]",
                "[0, 0, 0];[0, 0, 1]",
            ),
            "bond_006": ("[0, 1, 1]@[0, 1/2, 1/2]", [0, -2, -8, 10], (0, 0), "[0, 1, 1]", "[0, 0, 0];[0, 1, 1]"),
            "bond_007": ("[0, 1, -1]@[0, 1/2, 1/2]", [-1, 3, 9, -11], (0, 0), "[0, 1, -1]", "[0, 0, 1];[0, 1, 0]"),
            "bond_008": ("[1, 0, -1]@[1/2, 0, 1/2]", [4, -6, -12, 14], (0, 0), "[1, 0, -1]", "[0, 0, 1];[1, 0, 0]"),
            "bond_009": ("[1, 0, 1]@[1/2, 0, 1/2]", [-5, 7, 13, -15], (0, 0), "[1, 0, 1]", "[0, 0, 0];[1, 0, 1]"),
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
                {"A:A:1:1": ("[0, 0, 0];[0, 1, 0]", "2f", "ND", 1)},
                {"A:A:2:1": ("[0, 0, 0];[1, 1, 0]", "1c", "ND", 2)},
                {"A:A:3:1": ("[0, 0, 0];[0, 0, 1]", "1b", "ND", 3)},
                {"A:A:4:1": ("[0, 0, 0];[0, 1, 1]", "2e", "ND", 4)},
                {"A:A:5:1": ("[-1, 0, 0];[1, 0, 0]", "1a", "ND", 5)},
                {"A:A:6:1": ("[0, 0, 0];[1, 1, 1]", "1d", "ND", 6)},
                {"A:A:7:1": ("[-1, 0, 0];[1, 1, 0]", "2f", "ND", 7)},
                {"A:A:8:1": ("[-1, 0, 0];[1, 0, 1]", "1b", "ND", 8)},
                {"A:A:9:1": ("[-1, 0, 0];[1, 1, 1]", "2e", "ND", 9)},
            ]
        },
        "cell_range": (-2, 3, -2, 3, -2, 3),
        "max_neighbor": 10,
        "A": "[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.5]]",
        "version": "1.1.1",
    },
}
