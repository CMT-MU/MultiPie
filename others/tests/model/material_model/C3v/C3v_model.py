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
C3v = {
    "info": {
        "model": "C3v",
        "molecule": True,
        "group": ("C3v-1", "point group No. 19 : C3v-1 / 31m (31m)"),
        "crystal": "trigonal",
        "option": {"view": [0, 0, 1], "view_mode": "standard", "output": "C3v", "minimal_samb": True},
        "generate": {"model_type": "tight_binding", "time_reversal_type": "electric", "irrep": ["A1"]},
        "dimension": 12,
        "spinful": False,
        "orbital": ["s", "px", "py", "pz"],
        "ket": [
            "s@A_1",
            "s@A_2",
            "s@A_3",
            "px@B_1",
            "py@B_1",
            "pz@B_1",
            "px@B_2",
            "py@B_2",
            "pz@B_2",
            "px@B_3",
            "py@B_3",
            "pz@B_3",
        ],
        "ket_site": ["A_1", "A_2", "A_3", "B_1", "B_2", "B_3"],
        "site": {"A": ("[-1/6,-1/6,0]", "s"), "B": ("[-2/3,0,0]", ["px", "py", "pz"])},
        "rep_site": {"A": ("[-1/6, -1/6, 0]", "3b", [["s"]]), "B": ("[-2/3, 0, 0]", "3b", [["px", "py", "pz"]])},
        "cell_site": {
            "A_1": ("[-1/6, -1/6, 0]", "[1,6]"),
            "A_2": ("[1/6, 0, 0]", "[2,5]"),
            "A_3": ("[0, 1/6, 0]", "[3,4]"),
            "B_1": ("[-2/3, 0, 0]", "[1,4]"),
            "B_2": ("[0, -2/3, 0]", "[2,6]"),
            "B_3": ("[2/3, 2/3, 0]", "[3,5]"),
        },
        "bond": [("A", "A", 1), ("A", "B", 1)],
        "rep_bond": {
            "A:A:1:1": ("[-1/6, -1/6, 0];[1/6, 0, 0]", "3b", "ND", 1),
            "A:B:1:1": ("[-1/6, -1/6, 0];[-2/3, 0, 0]", "6c", "D", 1),
        },
        "cell_bond": {
            "A:A:1:1_1": ("[1/3, 1/6, 0]@[0, -1/12, 0]", "[1,-5]"),
            "A:A:1:1_2": ("[-1/6, 1/6, 0]@[1/12, 1/12, 0]", "[2,-4]"),
            "A:A:1:1_3": ("[1/6, 1/3, 0]@[-1/12, 0, 0]", "[-3,6]"),
            "A:B:1:1_1": ("[-1/2, 1/6, 0]@[-5/12, -1/12, 0]", "[1]"),
            "A:B:1:1_2": ("[-1/6, -2/3, 0]@[1/12, -1/3, 0]", "[2]"),
            "A:B:1:1_3": ("[2/3, 1/2, 0]@[1/3, 5/12, 0]", "[3]"),
            "A:B:1:1_4": ("[-2/3, -1/6, 0]@[-1/3, 1/12, 0]", "[4]"),
            "A:B:1:1_5": ("[1/2, 2/3, 0]@[5/12, 1/3, 0]", "[5]"),
            "A:B:1:1_6": ("[1/6, -1/2, 0]@[-1/12, -5/12, 0]", "[6]"),
        },
    },
    "name": {
        "alias": {
            "S_001": "A",
            "A": "S_001",
            "S_002": "B",
            "B": "S_002",
            "B_001": "A:A:1:1",
            "A:A:1:1": "B_001",
            "B_002": "A:B:1:1",
            "A:B:1:1": "B_002",
        },
        "site": {
            "site_001": ("A", 1),
            "site_002": ("A", 2),
            "site_003": ("A", 3),
            "site_004": ("B", 1),
            "site_005": ("B", 2),
            "site_006": ("B", 3),
        },
        "site_name": {
            "[-1/6, -1/6, 0]": ("site_001", 1),
            "[1/6, 0, 0]": ("site_002", 1),
            "[0, 1/6, 0]": ("site_003", 1),
            "[-2/3, 0, 0]": ("site_004", 1),
            "[0, -2/3, 0]": ("site_005", 1),
            "[2/3, 2/3, 0]": ("site_006", 1),
        },
        "bond": {
            "bond_001": ("A:A:1:1", 1),
            "bond_002": ("A:A:1:1", 2),
            "bond_003": ("A:A:1:1", 3),
            "bond_004": ("A:B:1:1", 1),
            "bond_005": ("A:B:1:1", 2),
            "bond_006": ("A:B:1:1", 3),
            "bond_007": ("A:B:1:1", 4),
            "bond_008": ("A:B:1:1", 5),
            "bond_009": ("A:B:1:1", 6),
        },
        "bond_name": {
            "[1/3, 1/6, 0]@[0, -1/12, 0]": ("bond_001", 1),
            "[-1/6, 1/6, 0]@[1/12, 1/12, 0]": ("bond_002", 1),
            "[1/6, 1/3, 0]@[-1/12, 0, 0]": ("bond_003", 1),
            "[-1/2, 1/6, 0]@[-5/12, -1/12, 0]": ("bond_004", 1),
            "[-1/6, -2/3, 0]@[1/12, -1/3, 0]": ("bond_005", 1),
            "[2/3, 1/2, 0]@[1/3, 5/12, 0]": ("bond_006", 1),
            "[-2/3, -1/6, 0]@[-1/3, 1/12, 0]": ("bond_007", 1),
            "[1/2, 2/3, 0]@[5/12, 1/3, 0]": ("bond_008", 1),
            "[1/6, -1/2, 0]@[-1/12, -5/12, 0]": ("bond_009", 1),
        },
    },
    "data": {
        "cluster_site": {"S_001": ["site_001", "site_002", "site_003"], "S_002": ["site_004", "site_005", "site_006"]},
        "cluster_bond": {
            "B_001": ["bond_001", "bond_002", "bond_003"],
            "B_002": ["bond_004", "bond_005", "bond_006", "bond_007", "bond_008", "bond_009"],
        },
        "site": {
            "site_001": ("[-1/6, -1/6, 0]", [0, 5], (0, 0)),
            "site_002": ("[1/6, 0, 0]", [1, 4], (1, 1)),
            "site_003": ("[0, 1/6, 0]", [2, 3], (2, 2)),
            "site_004": ("[-2/3, 0, 0]", [0, 3], (3, 3)),
            "site_005": ("[0, -2/3, 0]", [1, 5], (4, 4)),
            "site_006": ("[2/3, 2/3, 0]", [2, 4], (5, 5)),
        },
        "bond": {
            "bond_001": ("[1/3, 1/6, 0]@[0, -1/12, 0]", [0, -4], (1, 0), "[1/3, 1/6, 0]", "[-1/6, -1/6, 0];[1/6, 0, 0]"),
            "bond_002": ("[-1/6, 1/6, 0]@[1/12, 1/12, 0]", [1, -3], (2, 1), "[-1/6, 1/6, 0]", "[1/6, 0, 0];[0, 1/6, 0]"),
            "bond_003": ("[1/6, 1/3, 0]@[-1/12, 0, 0]", [-2, 5], (2, 0), "[1/6, 1/3, 0]", "[-1/6, -1/6, 0];[0, 1/6, 0]"),
            "bond_004": ("[-1/2, 1/6, 0]@[-5/12, -1/12, 0]", [0], (3, 0), "[-1/2, 1/6, 0]", "[-1/6, -1/6, 0];[-2/3, 0, 0]"),
            "bond_005": ("[-1/6, -2/3, 0]@[1/12, -1/3, 0]", [1], (4, 1), "[-1/6, -2/3, 0]", "[1/6, 0, 0];[0, -2/3, 0]"),
            "bond_006": ("[2/3, 1/2, 0]@[1/3, 5/12, 0]", [2], (5, 2), "[2/3, 1/2, 0]", "[0, 1/6, 0];[2/3, 2/3, 0]"),
            "bond_007": ("[-2/3, -1/6, 0]@[-1/3, 1/12, 0]", [3], (3, 2), "[-2/3, -1/6, 0]", "[0, 1/6, 0];[-2/3, 0, 0]"),
            "bond_008": ("[1/2, 2/3, 0]@[5/12, 1/3, 0]", [4], (5, 1), "[1/2, 2/3, 0]", "[1/6, 0, 0];[2/3, 2/3, 0]"),
            "bond_009": ("[1/6, -1/2, 0]@[-1/12, -5/12, 0]", [5], (4, 0), "[1/6, -1/2, 0]", "[-1/6, -1/6, 0];[0, -2/3, 0]"),
        },
        "cluster_atomic": {
            (0, 0): [(0, 0, "M_001")],
            (1, 1): [(1, 1, "M_001")],
            (2, 2): [(2, 2, "M_001")],
            (3, 3): [(3, 3, "M_002")],
            (4, 4): [(6, 6, "M_002")],
            (5, 5): [(9, 9, "M_002")],
            (1, 0): [(1, 0, "M_001")],
            (2, 1): [(2, 1, "M_001")],
            (2, 0): [(2, 0, "M_001")],
            (3, 0): [(3, 0, "M_003")],
            (4, 1): [(6, 1, "M_003")],
            (5, 2): [(9, 2, "M_003")],
            (3, 2): [(3, 2, "M_003")],
            (5, 1): [(9, 1, "M_003")],
            (4, 0): [(6, 0, "M_003")],
        },
        "atomic_braket": {
            "M_001": (["s"], ["s"]),
            "M_002": (["px", "py", "pz"], ["px", "py", "pz"]),
            "M_003": (["px", "py", "pz"], ["s"]),
        },
    },
    "detail": {
        "rep_bond_all": {
            "A_A": [{}, {"A:A:1:1": ("[-1/6, -1/6, 0];[1/6, 0, 0]", "3b", "ND", 1)}],
            "A_B": [
                {},
                {"A:B:1:1": ("[-1/6, -1/6, 0];[-2/3, 0, 0]", "6c", "D", 1)},
                {"A:B:2:1": ("[-1/6, -1/6, 0];[2/3, 2/3, 0]", "3b", "D", 2)},
            ],
        },
        "max_neighbor": 10,
        "version": "1.1.7",
    },
}
