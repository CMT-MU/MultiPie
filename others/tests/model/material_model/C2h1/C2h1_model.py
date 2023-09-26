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
C2h1 = {
    "info": {
        "model": "C2h1",
        "molecule": False,
        "group": ("C2h^1", "space group No. 10 : C2h^1 / P2/m (b-axis) : PG C2h"),
        "crystal": "monoclinic",
        "cell": {"a": 1.0, "b": 1.2, "c": 1.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
        "volume": 1.2,
        "a1": "[1.0, 0.0, 0.0]",
        "a2": "[0.0, 1.2, 0.0]",
        "a3": "[0.0, 0.0, 1.0]",
        "option": {"view": [0, 0, 1], "view_mode": "standard", "output": "C2h1", "minimal_samb": True},
        "generate": {
            "fourier_transform": False,
            "model_type": "tight_binding",
            "time_reversal_type": "electric",
            "irrep": ["Ag"],
            "toroidal_priority": False,
        },
        "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
        "k_path": "Γ-X",
        "dimension": 4,
        "spinful": True,
        "orbital": ["(s,U)", "(s,D)"],
        "ket": ["(s,U)@A_1", "(s,D)@A_1", "(s,U)@B_1", "(s,D)@B_1"],
        "ket_site": ["A_1", "B_1"],
        "site": {"A": ("[0,0,0]", "s"), "B": ("[1/2,1/2,0]", "s")},
        "rep_site": {
            "A": ("[0, 0, 0]", "1a", [["(s,U)", "(s,D)"]], "2/m"),
            "B": ("[1/2, 1/2, 0]", "1e", [["(s,U)", "(s,D)"]], "2/m"),
        },
        "cell_site": {"A_1": ("[0, 0, 0]", "[1,2,3,4]"), "B_1": ("[1/2, 1/2, 0]", "[1,2,3,4]")},
        "bond": [("A", "A", 1), ("B", "B", 1), ("A", "B", 1)],
        "rep_bond": {
            "A:A:1:1": ("[0, 0, 0];[1, 0, 0]", "1d", "ND", 1, "2/m"),
            "A:A:1:2": ("[0, 0, 0];[0, 0, 1]", "1c", "ND", 1, "2/m"),
            "B:B:1:1": ("[-1/2, 1/2, 0];[1/2, 1/2, 0]", "1b", "ND", 1, "2/m"),
            "B:B:1:2": ("[1/2, 1/2, 0];[1/2, 1/2, 1]", "1h", "ND", 1, "2/m"),
            "A:B:1:1": ("[0, 0, 0];[1/2, 1/2, 0]", "4o", "D", 1, "1"),
        },
        "cell_bond": {
            "A:A:1:1_1": ("[1, 0, 0]@[1/2, 0, 0]", "[1,-2,-3,4]"),
            "A:A:1:2_1": ("[0, 0, 1]@[0, 0, 1/2]", "[1,-2,-3,4]"),
            "B:B:1:1_1": ("[1, 0, 0]@[0, 1/2, 0]", "[1,-2,-3,4]"),
            "B:B:1:2_1": ("[0, 0, 1]@[1/2, 1/2, 1/2]", "[1,-2,-3,4]"),
            "A:B:1:1_1": ("[1/2, 1/2, 0]@[1/4, 1/4, 0]", "[1]"),
            "A:B:1:1_2": ("[-1/2, 1/2, 0]@[3/4, 1/4, 0]", "[2]"),
            "A:B:1:1_3": ("[-1/2, -1/2, 0]@[3/4, 3/4, 0]", "[3]"),
            "A:B:1:1_4": ("[1/2, -1/2, 0]@[1/4, 3/4, 0]", "[4]"),
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
            "B_002": "A:A:1:2",
            "A:A:1:2": "B_002",
            "B_003": "B:B:1:1",
            "B:B:1:1": "B_003",
            "B_004": "B:B:1:2",
            "B:B:1:2": "B_004",
            "B_005": "A:B:1:1",
            "A:B:1:1": "B_005",
        },
        "site": {"site_001": ("A", 1), "site_002": ("B", 1)},
        "site_name": {"[0, 0, 0]": ("site_001", 1), "[1/2, 1/2, 0]": ("site_002", 1)},
        "bond": {
            "bond_001": ("A:A:1:1", 1),
            "bond_002": ("A:A:1:2", 1),
            "bond_003": ("B:B:1:1", 1),
            "bond_004": ("B:B:1:2", 1),
            "bond_005": ("A:B:1:1", 1),
            "bond_006": ("A:B:1:1", 2),
            "bond_007": ("A:B:1:1", 3),
            "bond_008": ("A:B:1:1", 4),
        },
        "bond_name": {
            "[1, 0, 0]@[1/2, 0, 0]": ("bond_001", 1),
            "[0, 0, 1]@[0, 0, 1/2]": ("bond_002", 1),
            "[1, 0, 0]@[0, 1/2, 0]": ("bond_003", 1),
            "[0, 0, 1]@[1/2, 1/2, 1/2]": ("bond_004", 1),
            "[1/2, 1/2, 0]@[1/4, 1/4, 0]": ("bond_005", 1),
            "[-1/2, 1/2, 0]@[3/4, 1/4, 0]": ("bond_006", 1),
            "[-1/2, -1/2, 0]@[3/4, 3/4, 0]": ("bond_007", 1),
            "[1/2, -1/2, 0]@[1/4, 3/4, 0]": ("bond_008", 1),
        },
    },
    "data": {
        "plus_set": ["[0, 0, 0]"],
        "cluster_site": {"S_001": ["site_001"], "S_002": ["site_002"]},
        "cluster_bond": {
            "B_001": ["bond_001"],
            "B_002": ["bond_002"],
            "B_003": ["bond_003"],
            "B_004": ["bond_004"],
            "B_005": ["bond_005", "bond_006", "bond_007", "bond_008"],
        },
        "site": {"site_001": ("[0, 0, 0]", [0, 1, 2, 3], (0, 0)), "site_002": ("[1/2, 1/2, 0]", [0, 1, 2, 3], (1, 1))},
        "bond": {
            "bond_001": ("[1, 0, 0]@[1/2, 0, 0]", [0, -1, -2, 3], (0, 0), "[1, 0, 0]", "[0, 0, 0];[1, 0, 0]"),
            "bond_002": ("[0, 0, 1]@[0, 0, 1/2]", [0, -1, -2, 3], (0, 0), "[0, 0, 1]", "[0, 0, 0];[0, 0, 1]"),
            "bond_003": ("[1, 0, 0]@[0, 1/2, 0]", [0, -1, -2, 3], (1, 1), "[1, 0, 0]", "[-1/2, 1/2, 0];[1/2, 1/2, 0]"),
            "bond_004": ("[0, 0, 1]@[1/2, 1/2, 1/2]", [0, -1, -2, 3], (1, 1), "[0, 0, 1]", "[1/2, 1/2, 0];[1/2, 1/2, 1]"),
            "bond_005": ("[1/2, 1/2, 0]@[1/4, 1/4, 0]", [0], (1, 0), "[1/2, 1/2, 0]", "[0, 0, 0];[1/2, 1/2, 0]"),
            "bond_006": ("[-1/2, 1/2, 0]@[3/4, 1/4, 0]", [1], (1, 0), "[-1/2, 1/2, 0]", "[1, 0, 0];[1/2, 1/2, 0]"),
            "bond_007": ("[-1/2, -1/2, 0]@[3/4, 3/4, 0]", [2], (1, 0), "-bond_005", "[1, 1, 0];[1/2, 1/2, 0]"),
            "bond_008": ("[1/2, -1/2, 0]@[1/4, 3/4, 0]", [3], (1, 0), "-bond_006", "[0, 1, 0];[1/2, 1/2, 0]"),
        },
        "cluster_atomic": {(0, 0): [(0, 0, "M_001")], (1, 1): [(2, 2, "M_001")], (1, 0): [(2, 0, "M_001")]},
        "atomic_braket": {"M_001": (["(s,U)", "(s,D)"], ["(s,U)", "(s,D)"])},
    },
    "detail": {
        "rep_bond_all": {
            "A_A": [
                {},
                {
                    "A:A:1:1": ("[0, 0, 0];[1, 0, 0]", "1d", "ND", 1, "2/m"),
                    "A:A:1:2": ("[0, 0, 0];[0, 0, 1]", "1c", "ND", 1, "2/m"),
                },
                {"A:A:2:1": ("[0, 0, 0];[0, 1, 0]", "1b", "ND", 2, "2/m")},
                {
                    "A:A:3:1": ("[0, 0, 0];[1, 0, 1]", "1g", "ND", 3, "2/m"),
                    "A:A:3:2": ("[0, 0, 1];[1, 0, 0]", "1g", "ND", 3, "2/m"),
                },
                {
                    "A:A:4:1": ("[0, 0, 0];[1, 1, 0]", "1e", "ND", 4, "2/m"),
                    "A:A:4:2": ("[0, 0, 0];[0, 1, 1]", "1f", "ND", 4, "2/m"),
                },
                {
                    "A:A:5:1": ("[0, 0, 0];[1, 1, 1]", "1h", "ND", 5, "2/m"),
                    "A:A:5:2": ("[0, 0, 1];[1, 1, 0]", "1h", "ND", 5, "2/m"),
                },
                {
                    "A:A:6:1": ("[-1, 0, 0];[1, 0, 0]", "1a", "ND", 6, "2/m"),
                    "A:A:6:2": ("[0, 0, -1];[0, 0, 1]", "1a", "ND", 6, "2/m"),
                },
                {
                    "A:A:7:1": ("[-1, 0, 0];[1, 0, 1]", "1c", "ND", 7, "2/m"),
                    "A:A:7:2": ("[0, 0, -1];[1, 0, 1]", "1d", "ND", 7, "2/m"),
                    "A:A:7:3": ("[-1, 0, 1];[1, 0, 0]", "1c", "ND", 7, "2/m"),
                    "A:A:7:4": ("[0, 0, 1];[1, 0, -1]", "1d", "ND", 7, "2/m"),
                },
                {
                    "A:A:8:1": ("[-1, 0, 0];[1, 1, 0]", "1b", "ND", 8, "2/m"),
                    "A:A:8:2": ("[0, 0, -1];[0, 1, 1]", "1b", "ND", 8, "2/m"),
                },
                {"A:A:9:1": ("[0, -1, 0];[0, 1, 0]", "1a", "ND", 9, "2/m")},
            ],
            "B_B": [
                {},
                {
                    "B:B:1:1": ("[-1/2, 1/2, 0];[1/2, 1/2, 0]", "1b", "ND", 1, "2/m"),
                    "B:B:1:2": ("[1/2, 1/2, 0];[1/2, 1/2, 1]", "1h", "ND", 1, "2/m"),
                },
                {"B:B:2:1": ("[1/2, -1/2, 0];[1/2, 1/2, 0]", "1d", "ND", 2, "2/m")},
                {
                    "B:B:3:1": ("[-1/2, 1/2, 1];[1/2, 1/2, 0]", "1f", "ND", 3, "2/m"),
                    "B:B:3:2": ("[-1/2, 1/2, 0];[1/2, 1/2, 1]", "1f", "ND", 3, "2/m"),
                },
                {
                    "B:B:4:1": ("[1/2, -1/2, 0];[1/2, 1/2, 1]", "1g", "ND", 4, "2/m"),
                    "B:B:4:2": ("[-1/2, -1/2, 0];[1/2, 1/2, 0]", "1a", "ND", 4, "2/m"),
                },
                {
                    "B:B:5:1": ("[-1/2, -1/2, 1];[1/2, 1/2, 0]", "1c", "ND", 5, "2/m"),
                    "B:B:5:2": ("[-1/2, -1/2, 0];[1/2, 1/2, 1]", "1c", "ND", 5, "2/m"),
                },
                {
                    "B:B:6:1": ("[1/2, 1/2, -1];[1/2, 1/2, 1]", "1e", "ND", 6, "2/m"),
                    "B:B:6:2": ("[-1/2, 1/2, 0];[3/2, 1/2, 0]", "1e", "ND", 6, "2/m"),
                },
                {
                    "B:B:7:1": ("[-1/2, 1/2, 1];[3/2, 1/2, 0]", "1h", "ND", 7, "2/m"),
                    "B:B:7:2": ("[-1/2, 1/2, -1];[1/2, 1/2, 1]", "1b", "ND", 7, "2/m"),
                    "B:B:7:3": ("[-1/2, 1/2, 1];[1/2, 1/2, -1]", "1b", "ND", 7, "2/m"),
                    "B:B:7:4": ("[-1/2, 1/2, 0];[3/2, 1/2, 1]", "1h", "ND", 7, "2/m"),
                },
                {
                    "B:B:8:1": ("[-1/2, -1/2, 0];[3/2, 1/2, 0]", "1d", "ND", 8, "2/m"),
                    "B:B:8:2": ("[1/2, -1/2, -1];[1/2, 1/2, 1]", "1d", "ND", 8, "2/m"),
                },
                {"B:B:9:1": ("[1/2, -1/2, 0];[1/2, 3/2, 0]", "1e", "ND", 9, "2/m")},
            ],
            "A_B": [
                {},
                {"A:B:1:1": ("[0, 0, 0];[1/2, 1/2, 0]", "4o", "D", 1, "1")},
                {
                    "A:B:2:1": ("[0, 0, 0];[1/2, 1/2, 1]", "4o", "D", 2, "1"),
                    "A:B:2:2": ("[0, 0, 1];[1/2, 1/2, 0]", "4o", "D", 2, "1"),
                },
                {"A:B:3:1": ("[0, 0, 0];[3/2, 1/2, 0]", "4o", "D", 3, "1")},
                {"A:B:4:1": ("[0, 0, 0];[1/2, 3/2, 0]", "4o", "D", 4, "1")},
                {
                    "A:B:5:1": ("[0, 0, 1];[3/2, 1/2, 0]", "4o", "D", 5, "1"),
                    "A:B:5:2": ("[0, 0, 0];[3/2, 1/2, 1]", "4o", "D", 5, "1"),
                },
                {
                    "A:B:6:1": ("[0, 0, 0];[1/2, 3/2, 1]", "4o", "D", 6, "1"),
                    "A:B:6:2": ("[0, 0, 1];[1/2, 3/2, 0]", "4o", "D", 6, "1"),
                },
                {
                    "A:B:7:1": ("[0, 0, -1];[1/2, 1/2, 1]", "4o", "D", 7, "1"),
                    "A:B:7:2": ("[0, 0, 1];[1/2, 1/2, -1]", "4o", "D", 7, "1"),
                },
                {"A:B:8:1": ("[0, 0, 0];[3/2, 3/2, 0]", "4o", "D", 8, "1")},
                {
                    "A:B:9:1": ("[0, 0, 0];[3/2, 3/2, 1]", "4o", "D", 9, "1"),
                    "A:B:9:2": ("[0, 0, 1];[3/2, 3/2, 0]", "4o", "D", 9, "1"),
                },
            ],
        },
        "cell_range": (-2, 3, -2, 3, -2, 3),
        "max_neighbor": 10,
        "A": "[[1.0, 0.0, 0.0], [0.0, 1.2, 0.0], [0.0, 0.0, 1.0]]",
        "version": "1.1.14",
    },
}
