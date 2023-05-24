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
kappaET = {
    "info": {
        "model": "kappaET",
        "molecule": False,
        "group": ("C2v^8", "space group No. 32 : C2v^8 / Pba2 : PG C2v"),
        "crystal": "orthorhombic",
        "cell": {"a": 1.0, "b": 1.2, "c": 1.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
        "volume": 1.2,
        "a1": "[1.0, 0.0, 0.0]",
        "a2": "[0.0, 1.2, 0.0]",
        "a3": "[0.0, 0.0, 1.0]",
        "option": {"view": [0, 0, 1], "view_mode": "standard", "output": "kappaET", "minimal_samb": True},
        "generate": {"model_type": "tight_binding", "time_reversal_type": "electric", "irrep": ["A1"]},
        "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
        "k_path": "Γ-X",
        "dimension": 8,
        "spinful": True,
        "orbital": ["(s,U)", "(s,D)"],
        "ket": ["(s,U)@A_1", "(s,D)@A_1", "(s,U)@A_2", "(s,D)@A_2", "(s,U)@A_3", "(s,D)@A_3", "(s,U)@A_4", "(s,D)@A_4"],
        "ket_site": ["A_1", "A_2", "A_3", "A_4"],
        "site": {"A": ("[9/10,1/20,0]", "s")},
        "rep_site": {"A": ("[9/10, 1/20, 0]", "4c", [["(s,U)", "(s,D)"]])},
        "cell_site": {
            "A_1": ("[9/10, 1/20, 0]", "[1]"),
            "A_2": ("[1/10, 19/20, 0]", "[2]"),
            "A_3": ("[2/5, 9/20, 0]", "[3]"),
            "A_4": ("[3/5, 11/20, 0]", "[4]"),
        },
        "bond": [("A", "A", [1, 2, 3])],
        "rep_bond": {
            "A:A:1:1": ("[-1/10, 1/20, 0];[1/10, -1/20, 0]", "2a", "ND", 1),
            "A:A:2:1": ("[9/10, 1/20, 0];[3/5, 11/20, 0]", "4c", "D", 2),
            "A:A:3:1": ("[-1/10, 1/20, 0];[2/5, 9/20, 0]", "4c", "D", 3),
        },
        "cell_bond": {
            "A:A:1:1_1": ("[1/5, -1/10, 0]@[0, 0, 0]", "[1,-2]"),
            "A:A:1:1_2": ("[1/5, 1/10, 0]@[1/2, 1/2, 0]", "[3,-4]"),
            "A:A:2:1_1": ("[-3/10, 1/2, 0]@[3/4, 3/10, 0]", "[1]"),
            "A:A:2:1_2": ("[3/10, -1/2, 0]@[1/4, 7/10, 0]", "[2]"),
            "A:A:2:1_3": ("[3/10, 1/2, 0]@[1/4, 1/5, 0]", "[-3]"),
            "A:A:2:1_4": ("[-3/10, -1/2, 0]@[3/4, 4/5, 0]", "[-4]"),
            "A:A:3:1_1": ("[1/2, 2/5, 0]@[3/20, 1/4, 0]", "[1]"),
            "A:A:3:1_2": ("[-1/2, -2/5, 0]@[17/20, 3/4, 0]", "[2]"),
            "A:A:3:1_3": ("[-1/2, 2/5, 0]@[13/20, 1/4, 0]", "[-3]"),
            "A:A:3:1_4": ("[1/2, -2/5, 0]@[7/20, 3/4, 0]", "[-4]"),
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
        },
        "site": {"site_001": ("A", 1), "site_002": ("A", 2), "site_003": ("A", 3), "site_004": ("A", 4)},
        "site_name": {
            "[9/10, 1/20, 0]": ("site_001", 1),
            "[1/10, 19/20, 0]": ("site_002", 1),
            "[2/5, 9/20, 0]": ("site_003", 1),
            "[3/5, 11/20, 0]": ("site_004", 1),
        },
        "bond": {
            "bond_001": ("A:A:1:1", 1),
            "bond_002": ("A:A:1:1", 2),
            "bond_003": ("A:A:2:1", 1),
            "bond_004": ("A:A:2:1", 2),
            "bond_005": ("A:A:2:1", 3),
            "bond_006": ("A:A:2:1", 4),
            "bond_007": ("A:A:3:1", 1),
            "bond_008": ("A:A:3:1", 2),
            "bond_009": ("A:A:3:1", 3),
            "bond_010": ("A:A:3:1", 4),
        },
        "bond_name": {
            "[1/5, -1/10, 0]@[0, 0, 0]": ("bond_001", 1),
            "[1/5, 1/10, 0]@[1/2, 1/2, 0]": ("bond_002", 1),
            "[-3/10, 1/2, 0]@[3/4, 3/10, 0]": ("bond_003", 1),
            "[3/10, -1/2, 0]@[1/4, 7/10, 0]": ("bond_004", 1),
            "[3/10, 1/2, 0]@[1/4, 1/5, 0]": ("bond_005", 1),
            "[-3/10, -1/2, 0]@[3/4, 4/5, 0]": ("bond_006", 1),
            "[1/2, 2/5, 0]@[3/20, 1/4, 0]": ("bond_007", 1),
            "[-1/2, -2/5, 0]@[17/20, 3/4, 0]": ("bond_008", 1),
            "[-1/2, 2/5, 0]@[13/20, 1/4, 0]": ("bond_009", 1),
            "[1/2, -2/5, 0]@[7/20, 3/4, 0]": ("bond_010", 1),
        },
    },
    "data": {
        "plus_set": ["[0, 0, 0]"],
        "cluster_site": {"S_001": ["site_001", "site_002", "site_003", "site_004"]},
        "cluster_bond": {
            "B_001": ["bond_001", "bond_002"],
            "B_002": ["bond_003", "bond_004", "bond_005", "bond_006"],
            "B_003": ["bond_007", "bond_008", "bond_009", "bond_010"],
        },
        "site": {
            "site_001": ("[9/10, 1/20, 0]", [0], (0, 0)),
            "site_002": ("[1/10, 19/20, 0]", [1], (1, 1)),
            "site_003": ("[2/5, 9/20, 0]", [2], (2, 2)),
            "site_004": ("[3/5, 11/20, 0]", [3], (3, 3)),
        },
        "bond": {
            "bond_001": ("[1/5, -1/10, 0]@[0, 0, 0]", [0, -1], (0, 1), "[1/5, -1/10, 0]", "[-1/10, 1/20, 0];[1/10, -1/20, 0]"),
            "bond_002": ("[1/5, 1/10, 0]@[1/2, 1/2, 0]", [2, -3], (2, 3), "[1/5, 1/10, 0]", "[2/5, 9/20, 0];[3/5, 11/20, 0]"),
            "bond_003": ("[-3/10, 1/2, 0]@[3/4, 3/10, 0]", [0], (0, 3), "[-3/10, 1/2, 0]", "[9/10, 1/20, 0];[3/5, 11/20, 0]"),
            "bond_004": ("[3/10, -1/2, 0]@[1/4, 7/10, 0]", [1], (1, 2), "-bond_003", "[1/10, 19/20, 0];[2/5, 9/20, 0]"),
            "bond_005": ("[3/10, 1/2, 0]@[1/4, 1/5, 0]", [-2], (1, 2), "[3/10, 1/2, 0]", "[1/10, -1/20, 0];[2/5, 9/20, 0]"),
            "bond_006": ("[-3/10, -1/2, 0]@[3/4, 4/5, 0]", [-3], (0, 3), "-bond_005", "[9/10, 21/20, 0];[3/5, 11/20, 0]"),
            "bond_007": ("[1/2, 2/5, 0]@[3/20, 1/4, 0]", [0], (0, 2), "[1/2, 2/5, 0]", "[-1/10, 1/20, 0];[2/5, 9/20, 0]"),
            "bond_008": ("[-1/2, -2/5, 0]@[17/20, 3/4, 0]", [1], (1, 3), "-bond_007", "[11/10, 19/20, 0];[3/5, 11/20, 0]"),
            "bond_009": ("[-1/2, 2/5, 0]@[13/20, 1/4, 0]", [-2], (0, 2), "[-1/2, 2/5, 0]", "[9/10, 1/20, 0];[2/5, 9/20, 0]"),
            "bond_010": ("[1/2, -2/5, 0]@[7/20, 3/4, 0]", [-3], (1, 3), "-bond_009", "[1/10, 19/20, 0];[3/5, 11/20, 0]"),
        },
        "cluster_atomic": {
            (0, 0): [(0, 0, "M_001")],
            (1, 1): [(2, 2, "M_001")],
            (2, 2): [(4, 4, "M_001")],
            (3, 3): [(6, 6, "M_001")],
            (0, 1): [(0, 2, "M_001")],
            (2, 3): [(4, 6, "M_001")],
            (0, 3): [(0, 6, "M_001")],
            (1, 2): [(2, 4, "M_001")],
            (0, 2): [(0, 4, "M_001")],
            (1, 3): [(2, 6, "M_001")],
        },
        "atomic_braket": {"M_001": (["(s,U)", "(s,D)"], ["(s,U)", "(s,D)"])},
    },
    "detail": {
        "rep_bond_all": {
            "A_A": [
                {},
                {"A:A:1:1": ("[-1/10, 1/20, 0];[1/10, -1/20, 0]", "2a", "ND", 1)},
                {"A:A:2:1": ("[9/10, 1/20, 0];[3/5, 11/20, 0]", "4c", "D", 2)},
                {"A:A:3:1": ("[-1/10, 1/20, 0];[2/5, 9/20, 0]", "4c", "D", 3)},
                {"A:A:4:1": ("[9/10, 1/20, 0];[1/10, -1/20, 0]", "2b", "ND", 4)},
                {"A:A:5:1": ("[-1/10, 21/20, 0];[2/5, 9/20, 0]", "4c", "D", 5)},
                {"A:A:6:1": ("[-1/10, 1/20, 0];[3/5, 11/20, 0]", "4c", "D", 6)},
                {
                    "A:A:7:1": ("[9/10, 1/20, 0];[-1/10, 1/20, 0]", "4c", "D", 7),
                    "A:A:7:2": ("[9/10, 1/20, 1];[9/10, 1/20, 0]", "4c", "D", 7),
                },
                {"A:A:8:1": ("[-1/10, 1/20, 0];[1/10, -1/20, 1]", "2a", "D", 8)},
                {"A:A:9:1": ("[-1/10, 1/20, 0];[1/10, 19/20, 0]", "2b", "ND", 9)},
            ]
        },
        "cell_range": (-2, 3, -2, 3, -2, 3),
        "max_neighbor": 10,
        "A": "[[1.0, 0.0, 0.0], [0.0, 1.2, 0.0], [0.0, 0.0, 1.0]]",
        "version": "1.1.1",
    },
}
