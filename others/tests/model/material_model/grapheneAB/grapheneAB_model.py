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
grapheneAB = {
    "info": {
        "model": "grapheneAB",
        "molecule": False,
        "group": ("D3h^1", "space group No. 187 : D3h^1 / P-6m2 : PG D3h"),
        "crystal": "hexagonal",
        "cell": {"a": 2.435, "b": 2.435, "c": 10.0, "alpha": 90.0, "beta": 90.0, "gamma": 120.0},
        "volume": 51.3485947475379,
        "a1": "[2.435, 0.0, 0.0]",
        "a2": "[-1.2175, 2.10877185821511, 0.0]",
        "a3": "[0.0, 0.0, 10.0]",
        "option": {"view": [0, 0, 1], "view_mode": "standard", "output": "grapheneAB", "minimal_samb": True},
        "generate": {
            "fourier_transform": False,
            "model_type": "tight_binding",
            "time_reversal_type": "electric",
            "irrep": ["A1'"],
            "toroidal_priority": False,
        },
        "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
        "k_path": "Γ-X",
        "dimension": 3,
        "spinful": False,
        "orbital": ["s", "px", "py"],
        "ket": ["s@A_1", "px@B_1", "py@B_1"],
        "ket_site": ["A_1", "B_1"],
        "site": {"A": ("[1/3,2/3,0]", "s"), "B": ("[2/3,1/3,0]", ["px", "py"])},
        "rep_site": {"A": ("[1/3, 2/3, 0]", "1c", [["s"]], "-6m2"), "B": ("[2/3, 1/3, 0]", "1e", [["px", "py"]], "-6m2")},
        "cell_site": {
            "A_1": ("[1/3, 2/3, 0]", "[1,2,3,4,5,6,7,8,9,10,11,12]"),
            "B_1": ("[2/3, 1/3, 0]", "[1,2,3,4,5,6,7,8,9,10,11,12]"),
        },
        "bond": [("A", "B", 1), ("A", "A", 1), ("B", "B", 1)],
        "rep_bond": {
            "A:B:1:1": ("[1/3, -1/3, 0];[2/3, 1/3, 0]", "3j", "D", 1, "mm2"),
            "A:A:1:1": ("[1/3, -1/3, 0];[1/3, 2/3, 0]", "3j", "ND", 1, "mm2"),
            "B:B:1:1": ("[-1/3, 1/3, 0];[2/3, 1/3, 0]", "3j", "ND", 1, "mm2"),
        },
        "cell_bond": {
            "A:B:1:1_1": ("[1/3, 2/3, 0]@[1/2, 0, 0]", "[1,2,7,10]"),
            "A:B:1:1_2": ("[1/3, -1/3, 0]@[1/2, 1/2, 0]", "[3,6,8,11]"),
            "A:B:1:1_3": ("[-2/3, -1/3, 0]@[0, 1/2, 0]", "[4,5,9,12]"),
            "A:A:1:1_1": ("[0, 1, 0]@[1/3, 1/6, 0]", "[1,-3,-8,10]"),
            "A:A:1:1_2": ("[1, 1, 0]@[5/6, 1/6, 0]", "[2,-5,7,-12]"),
            "A:A:1:1_3": ("[1, 0, 0]@[5/6, 2/3, 0]", "[-4,6,-9,11]"),
            "B:B:1:1_1": ("[1, 0, 0]@[1/6, 1/3, 0]", "[1,-2,-7,10]"),
            "B:B:1:1_2": ("[1, 1, 0]@[1/6, 5/6, 0]", "[3,-6,8,-11]"),
            "B:B:1:1_3": ("[0, 1, 0]@[2/3, 5/6, 0]", "[-4,5,-9,12]"),
        },
    },
    "name": {
        "alias": {
            "S_001": "A",
            "A": "S_001",
            "S_002": "B",
            "B": "S_002",
            "B_001": "A:B:1:1",
            "A:B:1:1": "B_001",
            "B_002": "A:A:1:1",
            "A:A:1:1": "B_002",
            "B_003": "B:B:1:1",
            "B:B:1:1": "B_003",
        },
        "site": {"site_001": ("A", 1), "site_002": ("B", 1)},
        "site_name": {"[1/3, 2/3, 0]": ("site_001", 1), "[2/3, 1/3, 0]": ("site_002", 1)},
        "bond": {
            "bond_001": ("A:B:1:1", 1),
            "bond_002": ("A:B:1:1", 2),
            "bond_003": ("A:B:1:1", 3),
            "bond_004": ("A:A:1:1", 1),
            "bond_005": ("A:A:1:1", 2),
            "bond_006": ("A:A:1:1", 3),
            "bond_007": ("B:B:1:1", 1),
            "bond_008": ("B:B:1:1", 2),
            "bond_009": ("B:B:1:1", 3),
        },
        "bond_name": {
            "[1/3, 2/3, 0]@[1/2, 0, 0]": ("bond_001", 1),
            "[1/3, -1/3, 0]@[1/2, 1/2, 0]": ("bond_002", 1),
            "[-2/3, -1/3, 0]@[0, 1/2, 0]": ("bond_003", 1),
            "[0, 1, 0]@[1/3, 1/6, 0]": ("bond_004", 1),
            "[1, 1, 0]@[5/6, 1/6, 0]": ("bond_005", 1),
            "[1, 0, 0]@[5/6, 2/3, 0]": ("bond_006", 1),
            "[1, 0, 0]@[1/6, 1/3, 0]": ("bond_007", 1),
            "[1, 1, 0]@[1/6, 5/6, 0]": ("bond_008", 1),
            "[0, 1, 0]@[2/3, 5/6, 0]": ("bond_009", 1),
        },
    },
    "data": {
        "plus_set": ["[0, 0, 0]"],
        "cluster_site": {"S_001": ["site_001"], "S_002": ["site_002"]},
        "cluster_bond": {
            "B_001": ["bond_001", "bond_002", "bond_003"],
            "B_002": ["bond_004", "bond_005", "bond_006"],
            "B_003": ["bond_007", "bond_008", "bond_009"],
        },
        "site": {
            "site_001": ("[1/3, 2/3, 0]", [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], (0, 0)),
            "site_002": ("[2/3, 1/3, 0]", [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], (1, 1)),
        },
        "bond": {
            "bond_001": ("[1/3, 2/3, 0]@[1/2, 0, 0]", [0, 1, 6, 9], (1, 0), "[1/3, 2/3, 0]", "[1/3, -1/3, 0];[2/3, 1/3, 0]"),
            "bond_002": ("[1/3, -1/3, 0]@[1/2, 1/2, 0]", [2, 5, 7, 10], (1, 0), "[1/3, -1/3, 0]", "[1/3, 2/3, 0];[2/3, 1/3, 0]"),
            "bond_003": ("[-2/3, -1/3, 0]@[0, 1/2, 0]", [3, 4, 8, 11], (1, 0), "[-2/3, -1/3, 0]", "[1/3, 2/3, 0];[-1/3, 1/3, 0]"),
            "bond_004": ("[0, 1, 0]@[1/3, 1/6, 0]", [0, -2, -7, 9], (0, 0), "[0, 1, 0]", "[1/3, -1/3, 0];[1/3, 2/3, 0]"),
            "bond_005": ("[1, 1, 0]@[5/6, 1/6, 0]", [1, -4, 6, -11], (0, 0), "[1, 1, 0]", "[1/3, -1/3, 0];[4/3, 2/3, 0]"),
            "bond_006": ("[1, 0, 0]@[5/6, 2/3, 0]", [-3, 5, -8, 10], (0, 0), "[1, 0, 0]", "[1/3, 2/3, 0];[4/3, 2/3, 0]"),
            "bond_007": ("[1, 0, 0]@[1/6, 1/3, 0]", [0, -1, -6, 9], (1, 1), "[1, 0, 0]", "[-1/3, 1/3, 0];[2/3, 1/3, 0]"),
            "bond_008": ("[1, 1, 0]@[1/6, 5/6, 0]", [2, -5, 7, -10], (1, 1), "[1, 1, 0]", "[-1/3, 1/3, 0];[2/3, 4/3, 0]"),
            "bond_009": ("[0, 1, 0]@[2/3, 5/6, 0]", [-3, 4, -8, 11], (1, 1), "[0, 1, 0]", "[2/3, 1/3, 0];[2/3, 4/3, 0]"),
        },
        "cluster_atomic": {(0, 0): [(0, 0, "M_001")], (1, 1): [(1, 1, "M_002")], (1, 0): [(1, 0, "M_003")]},
        "atomic_braket": {"M_001": (["s"], ["s"]), "M_002": (["px", "py"], ["px", "py"]), "M_003": (["px", "py"], ["s"])},
    },
    "detail": {
        "rep_bond_all": {
            "A_B": [
                {},
                {"A:B:1:1": ("[1/3, -1/3, 0];[2/3, 1/3, 0]", "3j", "D", 1, "mm2")},
                {"A:B:2:1": ("[-2/3, -1/3, 0];[2/3, 1/3, 0]", "1a", "D", 2, "-6m2")},
                {"A:B:3:1": ("[-2/3, -1/3, 0];[2/3, 4/3, 0]", "3j", "D", 3, "mm2")},
                {"A:B:4:1": ("[-2/3, -1/3, 0];[5/3, 1/3, 0]", "3j", "D", 4, "mm2")},
                {"A:B:5:1": ("[-2/3, -4/3, 0];[2/3, 4/3, 0]", "1a", "D", 5, "-6m2")},
                {"A:B:6:1": ("[-2/3, -4/3, 0];[5/3, 4/3, 0]", "3j", "D", 6, "mm2")},
                {"A:B:7:1": ("[-5/3, -1/3, 0];[5/3, 4/3, 0]", "3j", "D", 7, "mm2")},
                {"A:B:8:1": ("[-2/3, -4/3, 0];[2/3, 7/3, 0]", "3j", "D", 8, "mm2")},
                {"A:B:9:1": ("[-2/3, -7/3, 0];[5/3, 7/3, 0]", "3j", "D", 9, "mm2")},
            ],
            "A_A": [
                {},
                {"A:A:1:1": ("[1/3, -1/3, 0];[1/3, 2/3, 0]", "3j", "ND", 1, "mm2")},
                {"A:A:2:1": ("[-2/3, -1/3, 0];[4/3, 2/3, 0]", "3j", "D", 2, "mm2")},
                {"A:A:3:1": ("[-2/3, -1/3, 0];[4/3, 5/3, 0]", "1c", "ND", 3, "-6m2")},
                {"A:A:4:1": ("[-2/3, -4/3, 0];[4/3, 5/3, 0]", "3j", "D", 4, "mm2")},
                {"A:A:5:1": ("[-2/3, -4/3, 0];[4/3, 8/3, 0]", "1c", "D", 5, "-6m2")},
                {"A:A:6:1": ("[1/3, 2/3, 0];[1/3, 2/3, 1]", "1d", "ND", 6, "-6m2")},
                {"A:A:7:1": ("[1/3, -1/3, 0];[1/3, 2/3, 1]", "3k", "ND", 7, "mm2")},
                {"A:A:8:1": ("[1/3, -1/3, 0];[4/3, 5/3, 1]", "3k", "D", 8, "mm2")},
                {"A:A:9:1": ("[-2/3, -1/3, 0];[4/3, 5/3, 1]", "1d", "ND", 9, "-6m2")},
            ],
            "B_B": [
                {},
                {"B:B:1:1": ("[-1/3, 1/3, 0];[2/3, 1/3, 0]", "3j", "ND", 1, "mm2")},
                {"B:B:2:1": ("[-1/3, -2/3, 0];[2/3, 4/3, 0]", "3j", "D", 2, "mm2")},
                {"B:B:3:1": ("[-1/3, -2/3, 0];[5/3, 4/3, 0]", "1e", "ND", 3, "-6m2")},
                {"B:B:4:1": ("[-1/3, -2/3, 0];[2/3, 7/3, 0]", "3j", "D", 4, "mm2")},
                {"B:B:5:1": ("[-4/3, -2/3, 0];[8/3, 4/3, 0]", "1e", "D", 5, "-6m2")},
                {"B:B:6:1": ("[2/3, 1/3, 0];[2/3, 1/3, 1]", "1f", "ND", 6, "-6m2")},
                {"B:B:7:1": ("[-1/3, 1/3, 0];[2/3, 1/3, 1]", "3k", "ND", 7, "mm2")},
                {"B:B:8:1": ("[-1/3, 1/3, 0];[5/3, 4/3, 1]", "3k", "D", 8, "mm2")},
                {"B:B:9:1": ("[-1/3, -2/3, 0];[5/3, 4/3, 1]", "1f", "ND", 9, "-6m2")},
            ],
        },
        "cell_range": (-2, 3, -2, 3, -2, 3),
        "max_neighbor": 10,
        "A": "[[2.435, -1.2175, 0.0], [0.0, 2.10877185821511, 0.0], [0.0, 0.0, 10.0]]",
        "version": "1.1.14",
    },
}
