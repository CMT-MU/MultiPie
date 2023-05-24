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
C3h = {
    "info": {
        "model": "C3h",
        "molecule": True,
        "group": ("C3h", "point group No. 22 : C3h / -6"),
        "crystal": "hexagonal",
        "option": {"view": [0, 0, 1], "view_mode": "standard", "output": "C3h", "minimal_samb": True},
        "generate": {"model_type": "tight_binding", "time_reversal_type": "electric", "irrep": ["A'"]},
        "dimension": 32,
        "spinful": True,
        "orbital": ["(s,U)", "(s,D)", "(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"],
        "ket": [
            "(s,U)@H1_1",
            "(s,D)@H1_1",
            "(s,U)@O_1",
            "(s,D)@O_1",
            "(px,U)@O_1",
            "(px,D)@O_1",
            "(py,U)@O_1",
            "(py,D)@O_1",
            "(pz,U)@O_1",
            "(pz,D)@O_1",
            "(s,U)@O_2",
            "(s,D)@O_2",
            "(px,U)@O_2",
            "(px,D)@O_2",
            "(py,U)@O_2",
            "(py,D)@O_2",
            "(pz,U)@O_2",
            "(pz,D)@O_2",
            "(s,U)@O_3",
            "(s,D)@O_3",
            "(px,U)@O_3",
            "(px,D)@O_3",
            "(py,U)@O_3",
            "(py,D)@O_3",
            "(pz,U)@O_3",
            "(pz,D)@O_3",
            "(s,U)@H2_1",
            "(s,D)@H2_1",
            "(s,U)@H2_2",
            "(s,D)@H2_2",
            "(s,U)@H2_3",
            "(s,D)@H2_3",
        ],
        "ket_site": ["H1_1", "O_1", "O_2", "O_3", "H2_1", "H2_2", "H2_3"],
        "site": {"H1": ("[0,0,0]", "s"), "O": ("[1/3,0,0]", "s p"), "H2": ("[1/2,1/6,0]", "s")},
        "rep_site": {
            "H1": ("[0, 0, 0]", "1o", [["(s,U)", "(s,D)"]]),
            "O": ("[1/3, 0, 0]", "3b", [["(s,U)", "(s,D)"], ["(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"]]),
            "H2": ("[1/2, 1/6, 0]", "3b", [["(s,U)", "(s,D)"]]),
        },
        "cell_site": {
            "H1_1": ("[0, 0, 0]", "[1,2,3,4,5,6]"),
            "O_1": ("[1/3, 0, 0]", "[1,4]"),
            "O_2": ("[0, 1/3, 0]", "[2,6]"),
            "O_3": ("[-1/3, -1/3, 0]", "[3,5]"),
            "H2_1": ("[1/2, 1/6, 0]", "[1,4]"),
            "H2_2": ("[-1/6, 1/3, 0]", "[2,6]"),
            "H2_3": ("[-1/3, -1/2, 0]", "[3,5]"),
        },
        "bond": [("H1", "O", 1), ("O", "H2", 1)],
        "rep_bond": {
            "H1:O:1:1": ("[0, 0, 0];[1/3, 0, 0]", "3b", "D", 1),
            "O:H2:1:1": ("[1/3, 0, 0];[1/2, 1/6, 0]", "3b", "D", 1),
        },
        "cell_bond": {
            "H1:O:1:1_1": ("[1/3, 0, 0]@[1/6, 0, 0]", "[1,4]"),
            "H1:O:1:1_2": ("[0, 1/3, 0]@[0, 1/6, 0]", "[2,6]"),
            "H1:O:1:1_3": ("[-1/3, -1/3, 0]@[-1/6, -1/6, 0]", "[3,5]"),
            "O:H2:1:1_1": ("[1/6, 1/6, 0]@[5/12, 1/12, 0]", "[1,4]"),
            "O:H2:1:1_2": ("[-1/6, 0, 0]@[-1/12, 1/3, 0]", "[2,6]"),
            "O:H2:1:1_3": ("[0, -1/6, 0]@[-1/3, -5/12, 0]", "[3,5]"),
        },
    },
    "name": {
        "alias": {
            "S_001": "H1",
            "H1": "S_001",
            "S_002": "O",
            "O": "S_002",
            "S_003": "H2",
            "H2": "S_003",
            "B_001": "H1:O:1:1",
            "H1:O:1:1": "B_001",
            "B_002": "O:H2:1:1",
            "O:H2:1:1": "B_002",
        },
        "site": {
            "site_001": ("H1", 1),
            "site_002": ("O", 1),
            "site_003": ("O", 2),
            "site_004": ("O", 3),
            "site_005": ("H2", 1),
            "site_006": ("H2", 2),
            "site_007": ("H2", 3),
        },
        "site_name": {
            "[0, 0, 0]": ("site_001", 1),
            "[1/3, 0, 0]": ("site_002", 1),
            "[0, 1/3, 0]": ("site_003", 1),
            "[-1/3, -1/3, 0]": ("site_004", 1),
            "[1/2, 1/6, 0]": ("site_005", 1),
            "[-1/6, 1/3, 0]": ("site_006", 1),
            "[-1/3, -1/2, 0]": ("site_007", 1),
        },
        "bond": {
            "bond_001": ("H1:O:1:1", 1),
            "bond_002": ("H1:O:1:1", 2),
            "bond_003": ("H1:O:1:1", 3),
            "bond_004": ("O:H2:1:1", 1),
            "bond_005": ("O:H2:1:1", 2),
            "bond_006": ("O:H2:1:1", 3),
        },
        "bond_name": {
            "[1/3, 0, 0]@[1/6, 0, 0]": ("bond_001", 1),
            "[0, 1/3, 0]@[0, 1/6, 0]": ("bond_002", 1),
            "[-1/3, -1/3, 0]@[-1/6, -1/6, 0]": ("bond_003", 1),
            "[1/6, 1/6, 0]@[5/12, 1/12, 0]": ("bond_004", 1),
            "[-1/6, 0, 0]@[-1/12, 1/3, 0]": ("bond_005", 1),
            "[0, -1/6, 0]@[-1/3, -5/12, 0]": ("bond_006", 1),
        },
    },
    "data": {
        "cluster_site": {
            "S_001": ["site_001"],
            "S_002": ["site_002", "site_003", "site_004"],
            "S_003": ["site_005", "site_006", "site_007"],
        },
        "cluster_bond": {"B_001": ["bond_001", "bond_002", "bond_003"], "B_002": ["bond_004", "bond_005", "bond_006"]},
        "site": {
            "site_001": ("[0, 0, 0]", [0, 1, 2, 3, 4, 5], (0, 0)),
            "site_002": ("[1/3, 0, 0]", [0, 3], (1, 1)),
            "site_003": ("[0, 1/3, 0]", [1, 5], (2, 2)),
            "site_004": ("[-1/3, -1/3, 0]", [2, 4], (3, 3)),
            "site_005": ("[1/2, 1/6, 0]", [0, 3], (4, 4)),
            "site_006": ("[-1/6, 1/3, 0]", [1, 5], (5, 5)),
            "site_007": ("[-1/3, -1/2, 0]", [2, 4], (6, 6)),
        },
        "bond": {
            "bond_001": ("[1/3, 0, 0]@[1/6, 0, 0]", [0, 3], (0, 1), "[1/3, 0, 0]", "[0, 0, 0];[1/3, 0, 0]"),
            "bond_002": ("[0, 1/3, 0]@[0, 1/6, 0]", [1, 5], (0, 2), "[0, 1/3, 0]", "[0, 0, 0];[0, 1/3, 0]"),
            "bond_003": ("[-1/3, -1/3, 0]@[-1/6, -1/6, 0]", [2, 4], (0, 3), "[-1/3, -1/3, 0]", "[0, 0, 0];[-1/3, -1/3, 0]"),
            "bond_004": ("[1/6, 1/6, 0]@[5/12, 1/12, 0]", [0, 3], (1, 4), "[1/6, 1/6, 0]", "[1/3, 0, 0];[1/2, 1/6, 0]"),
            "bond_005": ("[-1/6, 0, 0]@[-1/12, 1/3, 0]", [1, 5], (2, 5), "[-1/6, 0, 0]", "[0, 1/3, 0];[-1/6, 1/3, 0]"),
            "bond_006": ("[0, -1/6, 0]@[-1/3, -5/12, 0]", [2, 4], (3, 6), "[0, -1/6, 0]", "[-1/3, -1/3, 0];[-1/3, -1/2, 0]"),
        },
        "cluster_atomic": {
            (0, 0): [(0, 0, "M_001")],
            (1, 1): [(2, 2, "M_001"), (2, 4, "M_002"), (4, 4, "M_003")],
            (2, 2): [(10, 10, "M_001"), (10, 12, "M_002"), (12, 12, "M_003")],
            (3, 3): [(18, 18, "M_001"), (18, 20, "M_002"), (20, 20, "M_003")],
            (4, 4): [(26, 26, "M_001")],
            (5, 5): [(28, 28, "M_001")],
            (6, 6): [(30, 30, "M_001")],
            (0, 1): [(0, 2, "M_001"), (0, 4, "M_002")],
            (0, 2): [(0, 10, "M_001"), (0, 12, "M_002")],
            (0, 3): [(0, 18, "M_001"), (0, 20, "M_002")],
            (1, 4): [(2, 26, "M_001"), (4, 26, "M_004")],
            (2, 5): [(10, 28, "M_001"), (12, 28, "M_004")],
            (3, 6): [(18, 30, "M_001"), (20, 30, "M_004")],
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
            "H1_O": [{}, {"H1:O:1:1": ("[0, 0, 0];[1/3, 0, 0]", "3b", "D", 1)}],
            "O_H2": [
                {},
                {"O:H2:1:1": ("[1/3, 0, 0];[1/2, 1/6, 0]", "3b", "D", 1)},
                {"O:H2:2:1": ("[1/3, 0, 0];[-1/3, -1/2, 0]", "3b", "D", 2)},
                {"O:H2:3:1": ("[1/3, 0, 0];[-1/6, 1/3, 0]", "3b", "D", 3)},
            ],
        },
        "max_neighbor": 10,
        "version": "1.1.1",
    },
}
