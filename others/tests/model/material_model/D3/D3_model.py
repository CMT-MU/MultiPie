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
        - fourier_transform* : create fourier transformed SAMB ?
        - toroidal_priority : create toroidal multipoles (G,T) in priority ?
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
D3 = {
    "info": {
        "model": "D3",
        "molecule": True,
        "group": ("D3-1", "point group No. 18 : D3-1 / 321 (321)"),
        "crystal": "trigonal",
        "option": {"view": None, "view_mode": "standard", "output": "D3", "minimal_samb": True},
        "generate": {
            "model_type": "tight_binding",
            "time_reversal_type": "electric",
            "irrep": ["A1"],
            "toroidal_priority": False,
        },
        "dimension": 48,
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
        ],
        "ket_site": ["A_1", "A_2", "A_3", "A_4", "A_5", "A_6"],
        "site": {"A": ("[1,0,1]", "s p")},
        "rep_site": {
            "A": ("[1, 0, 1]", "6c", [["(s,U)", "(s,D)"], ["(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"]], "1")
        },
        "cell_site": {
            "A_1": ("[1, 0, 1]", "[1]"),
            "A_2": ("[1, 0, -1]", "[2]"),
            "A_3": ("[-1, -1, -1]", "[3]"),
            "A_4": ("[0, 1, -1]", "[4]"),
            "A_5": ("[0, 1, 1]", "[5]"),
            "A_6": ("[-1, -1, 1]", "[6]"),
        },
        "bond": [("A", "A", [1, 2])],
        "rep_bond": {
            "A:A:1:1": ("[1, 0, 1];[0, 1, 1]", "6c", "D", 1, "1"),
            "A:A:2:1": ("[1, 0, 1];[1, 0, -1]", "3b", "ND", 2, ".2."),
        },
        "cell_bond": {
            "A:A:1:1_1": ("[-1, 1, 0]@[1/2, 1/2, 1]", "[1]"),
            "A:A:1:1_2": ("[-2, -1, 0]@[0, -1/2, -1]", "[2]"),
            "A:A:1:1_3": ("[1, 2, 0]@[-1/2, 0, -1]", "[3]"),
            "A:A:1:1_4": ("[-1, 1, 0]@[1/2, 1/2, -1]", "[-4]"),
            "A:A:1:1_5": ("[-1, -2, 0]@[-1/2, 0, 1]", "[5]"),
            "A:A:1:1_6": ("[-2, -1, 0]@[0, -1/2, 1]", "[-6]"),
            "A:A:2:1_1": ("[0, 0, -2]@[1, 0, 0]", "[1,-2]"),
            "A:A:2:1_2": ("[0, 0, 2]@[-1, -1, 0]", "[3,-6]"),
            "A:A:2:1_3": ("[0, 0, 2]@[0, 1, 0]", "[4,-5]"),
        },
    },
    "name": {
        "alias": {"S_001": "A", "A": "S_001", "B_001": "A:A:1:1", "A:A:1:1": "B_001", "B_002": "A:A:2:1", "A:A:2:1": "B_002"},
        "site": {
            "site_001": ("A", 1),
            "site_002": ("A", 2),
            "site_003": ("A", 3),
            "site_004": ("A", 4),
            "site_005": ("A", 5),
            "site_006": ("A", 6),
        },
        "site_name": {
            "[1, 0, 1]": ("site_001", 1),
            "[1, 0, -1]": ("site_002", 1),
            "[-1, -1, -1]": ("site_003", 1),
            "[0, 1, -1]": ("site_004", 1),
            "[0, 1, 1]": ("site_005", 1),
            "[-1, -1, 1]": ("site_006", 1),
        },
        "bond": {
            "bond_001": ("A:A:1:1", 1),
            "bond_002": ("A:A:1:1", 2),
            "bond_003": ("A:A:1:1", 3),
            "bond_004": ("A:A:1:1", 4),
            "bond_005": ("A:A:1:1", 5),
            "bond_006": ("A:A:1:1", 6),
            "bond_007": ("A:A:2:1", 1),
            "bond_008": ("A:A:2:1", 2),
            "bond_009": ("A:A:2:1", 3),
        },
        "bond_name": {
            "[-1, 1, 0]@[1/2, 1/2, 1]": ("bond_001", 1),
            "[-2, -1, 0]@[0, -1/2, -1]": ("bond_002", 1),
            "[1, 2, 0]@[-1/2, 0, -1]": ("bond_003", 1),
            "[-1, 1, 0]@[1/2, 1/2, -1]": ("bond_004", 1),
            "[-1, -2, 0]@[-1/2, 0, 1]": ("bond_005", 1),
            "[-2, -1, 0]@[0, -1/2, 1]": ("bond_006", 1),
            "[0, 0, -2]@[1, 0, 0]": ("bond_007", 1),
            "[0, 0, 2]@[-1, -1, 0]": ("bond_008", 1),
            "[0, 0, 2]@[0, 1, 0]": ("bond_009", 1),
        },
    },
    "data": {
        "cluster_site": {"S_001": ["site_001", "site_002", "site_003", "site_004", "site_005", "site_006"]},
        "cluster_bond": {
            "B_001": ["bond_001", "bond_002", "bond_003", "bond_004", "bond_005", "bond_006"],
            "B_002": ["bond_007", "bond_008", "bond_009"],
        },
        "site": {
            "site_001": ("[1, 0, 1]", [0], (0, 0)),
            "site_002": ("[1, 0, -1]", [1], (1, 1)),
            "site_003": ("[-1, -1, -1]", [2], (2, 2)),
            "site_004": ("[0, 1, -1]", [3], (3, 3)),
            "site_005": ("[0, 1, 1]", [4], (4, 4)),
            "site_006": ("[-1, -1, 1]", [5], (5, 5)),
        },
        "bond": {
            "bond_001": ("[-1, 1, 0]@[1/2, 1/2, 1]", [0], (4, 0), "[-1, 1, 0]", "[1, 0, 1];[0, 1, 1]"),
            "bond_002": ("[-2, -1, 0]@[0, -1/2, -1]", [1], (2, 1), "[-2, -1, 0]", "[1, 0, -1];[-1, -1, -1]"),
            "bond_003": ("[1, 2, 0]@[-1/2, 0, -1]", [2], (3, 2), "[1, 2, 0]", "[-1, -1, -1];[0, 1, -1]"),
            "bond_004": ("[-1, 1, 0]@[1/2, 1/2, -1]", [-3], (3, 1), "bond_001", "[1, 0, -1];[0, 1, -1]"),
            "bond_005": ("[-1, -2, 0]@[-1/2, 0, 1]", [4], (5, 4), "-bond_003", "[0, 1, 1];[-1, -1, 1]"),
            "bond_006": ("[-2, -1, 0]@[0, -1/2, 1]", [-5], (5, 0), "bond_002", "[1, 0, 1];[-1, -1, 1]"),
            "bond_007": ("[0, 0, -2]@[1, 0, 0]", [0, -1], (1, 0), "[0, 0, -2]", "[1, 0, 1];[1, 0, -1]"),
            "bond_008": ("[0, 0, 2]@[-1, -1, 0]", [2, -5], (5, 2), "-bond_007", "[-1, -1, -1];[-1, -1, 1]"),
            "bond_009": ("[0, 0, 2]@[0, 1, 0]", [3, -4], (4, 3), "-bond_007", "[0, 1, -1];[0, 1, 1]"),
        },
        "cluster_atomic": {
            (0, 0): [(0, 0, "M_001"), (0, 2, "M_002"), (2, 2, "M_003")],
            (1, 1): [(8, 8, "M_001"), (8, 10, "M_002"), (10, 10, "M_003")],
            (2, 2): [(16, 16, "M_001"), (16, 18, "M_002"), (18, 18, "M_003")],
            (3, 3): [(24, 24, "M_001"), (24, 26, "M_002"), (26, 26, "M_003")],
            (4, 4): [(32, 32, "M_001"), (32, 34, "M_002"), (34, 34, "M_003")],
            (5, 5): [(40, 40, "M_001"), (40, 42, "M_002"), (42, 42, "M_003")],
            (4, 0): [(32, 0, "M_001"), (32, 2, "M_002"), (34, 2, "M_003")],
            (2, 1): [(16, 8, "M_001"), (16, 10, "M_002"), (18, 10, "M_003")],
            (3, 2): [(24, 16, "M_001"), (24, 18, "M_002"), (26, 18, "M_003")],
            (3, 1): [(24, 8, "M_001"), (24, 10, "M_002"), (26, 10, "M_003")],
            (5, 4): [(40, 32, "M_001"), (40, 34, "M_002"), (42, 34, "M_003")],
            (5, 0): [(40, 0, "M_001"), (40, 2, "M_002"), (42, 2, "M_003")],
            (1, 0): [(8, 0, "M_001"), (8, 2, "M_002"), (10, 2, "M_003")],
            (5, 2): [(40, 16, "M_001"), (40, 18, "M_002"), (42, 18, "M_003")],
            (4, 3): [(32, 24, "M_001"), (32, 26, "M_002"), (34, 26, "M_003")],
        },
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
                {"A:A:1:1": ("[1, 0, 1];[0, 1, 1]", "6c", "D", 1, "1")},
                {"A:A:2:1": ("[1, 0, 1];[1, 0, -1]", "3b", "ND", 2, ".2.")},
                {
                    "A:A:3:1": ("[1, 0, 1];[0, 1, -1]", "3b", "ND", 3, ".2."),
                    "A:A:3:2": ("[1, 0, 1];[-1, -1, -1]", "3b", "ND", 3, ".2."),
                },
            ]
        },
        "max_neighbor": 10,
        "version": "1.1.14",
        "cell_range": (-2, 3, -2, 3, -2, 3),
    },
}
