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
GaAs = {
    "info": {
        "model": "GaAs",
        "molecule": False,
        "group": ("Td^2", "space group No. 216 : Td^2 / F-43m : PG Td"),
        "crystal": "cubic",
        "cell": {"a": 1.0, "b": 1.0, "c": 1.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
        "volume": 1.0,
        "a1": "[1.0, 0.0, 0.0]",
        "a2": "[0.0, 1.0, 0.0]",
        "a3": "[0.0, 0.0, 1.0]",
        "option": {"view": None, "view_mode": "standard", "output": "GaAs", "minimal_samb": True},
        "generate": {
            "fourier_transform": False,
            "model_type": "tight_binding",
            "time_reversal_type": "electric",
            "irrep": ["A1"],
            "toroidal_priority": False,
        },
        "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
        "k_path": "Γ-X",
        "dimension": 6,
        "spinful": False,
        "orbital": ["px", "py", "pz"],
        "ket": ["px@Ga_1", "py@Ga_1", "pz@Ga_1", "px@As_1", "py@As_1", "pz@As_1"],
        "ket_site": ["Ga_1", "As_1"],
        "site": {"Ga": ("[0,0,0]", ["px", "py", "pz"]), "As": ("[1/4,1/4,1/4]", ["px", "py", "pz"])},
        "rep_site": {
            "Ga": ("[0, 0, 0]", "4a", [["px", "py", "pz"]], "-43m"),
            "As": ("[1/4, 1/4, 1/4]", "4c", [["px", "py", "pz"]], "-43m"),
        },
        "cell_site": {
            "Ga_1(1)": ("[0, 0, 0]", "[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]"),
            "Ga_1(2)": ("[0, 1/2, 1/2]", "[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]"),
            "Ga_1(3)": ("[1/2, 0, 1/2]", "[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]"),
            "Ga_1(4)": ("[1/2, 1/2, 0]", "[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]"),
            "As_1(1)": ("[1/4, 1/4, 1/4]", "[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]"),
            "As_1(2)": ("[1/4, 3/4, 3/4]", "[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]"),
            "As_1(3)": ("[3/4, 1/4, 3/4]", "[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]"),
            "As_1(4)": ("[3/4, 3/4, 1/4]", "[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]"),
        },
        "bond": [("Ga", "As", [1])],
        "rep_bond": {"Ga:As:1:1": ("[0, 0, 0];[1/4, 1/4, 1/4]", "16e", "D", 1, ".3m")},
        "cell_bond": {
            "Ga:As:1:1_1(1)": ("[1/4, 1/4, 1/4]@[1/8, 1/8, 1/8]", "[1,5,9,16,17,18]"),
            "Ga:As:1:1_1(2)": ("[1/4, 1/4, 1/4]@[1/8, 5/8, 5/8]", "[1,5,9,16,17,18]"),
            "Ga:As:1:1_1(3)": ("[1/4, 1/4, 1/4]@[5/8, 1/8, 5/8]", "[1,5,9,16,17,18]"),
            "Ga:As:1:1_1(4)": ("[1/4, 1/4, 1/4]@[5/8, 5/8, 1/8]", "[1,5,9,16,17,18]"),
            "Ga:As:1:1_2(1)": ("[-1/4, -1/4, 1/4]@[3/8, 3/8, 1/8]", "[2,6,11,13,21,23]"),
            "Ga:As:1:1_2(2)": ("[-1/4, -1/4, 1/4]@[3/8, 7/8, 5/8]", "[2,6,11,13,21,23]"),
            "Ga:As:1:1_2(3)": ("[-1/4, -1/4, 1/4]@[7/8, 3/8, 5/8]", "[2,6,11,13,21,23]"),
            "Ga:As:1:1_2(4)": ("[-1/4, -1/4, 1/4]@[7/8, 7/8, 1/8]", "[2,6,11,13,21,23]"),
            "Ga:As:1:1_3(1)": ("[1/4, -1/4, -1/4]@[1/8, 3/8, 3/8]", "[3,7,12,15,19,24]"),
            "Ga:As:1:1_3(2)": ("[1/4, -1/4, -1/4]@[1/8, 7/8, 7/8]", "[3,7,12,15,19,24]"),
            "Ga:As:1:1_3(3)": ("[1/4, -1/4, -1/4]@[5/8, 3/8, 7/8]", "[3,7,12,15,19,24]"),
            "Ga:As:1:1_3(4)": ("[1/4, -1/4, -1/4]@[5/8, 7/8, 3/8]", "[3,7,12,15,19,24]"),
            "Ga:As:1:1_4(1)": ("[-1/4, 1/4, -1/4]@[3/8, 1/8, 3/8]", "[4,8,10,14,20,22]"),
            "Ga:As:1:1_4(2)": ("[-1/4, 1/4, -1/4]@[3/8, 5/8, 7/8]", "[4,8,10,14,20,22]"),
            "Ga:As:1:1_4(3)": ("[-1/4, 1/4, -1/4]@[7/8, 1/8, 7/8]", "[4,8,10,14,20,22]"),
            "Ga:As:1:1_4(4)": ("[-1/4, 1/4, -1/4]@[7/8, 5/8, 3/8]", "[4,8,10,14,20,22]"),
        },
    },
    "name": {
        "alias": {"S_001": "Ga", "Ga": "S_001", "S_002": "As", "As": "S_002", "B_001": "Ga:As:1:1", "Ga:As:1:1": "B_001"},
        "site": {"site_001": ("Ga", 1), "site_002": ("As", 1)},
        "site_name": {
            "[0, 0, 0]": ("site_001", 1),
            "[0, 1/2, 1/2]": ("site_001", 2),
            "[1/2, 0, 1/2]": ("site_001", 3),
            "[1/2, 1/2, 0]": ("site_001", 4),
            "[1/4, 1/4, 1/4]": ("site_002", 1),
            "[1/4, 3/4, 3/4]": ("site_002", 2),
            "[3/4, 1/4, 3/4]": ("site_002", 3),
            "[3/4, 3/4, 1/4]": ("site_002", 4),
        },
        "bond": {
            "bond_001": ("Ga:As:1:1", 1),
            "bond_002": ("Ga:As:1:1", 2),
            "bond_003": ("Ga:As:1:1", 3),
            "bond_004": ("Ga:As:1:1", 4),
        },
        "bond_name": {
            "[1/4, 1/4, 1/4]@[1/8, 1/8, 1/8]": ("bond_001", 1),
            "[1/4, 1/4, 1/4]@[1/8, 5/8, 5/8]": ("bond_001", 2),
            "[1/4, 1/4, 1/4]@[5/8, 1/8, 5/8]": ("bond_001", 3),
            "[1/4, 1/4, 1/4]@[5/8, 5/8, 1/8]": ("bond_001", 4),
            "[-1/4, -1/4, 1/4]@[3/8, 3/8, 1/8]": ("bond_002", 1),
            "[-1/4, -1/4, 1/4]@[3/8, 7/8, 5/8]": ("bond_002", 2),
            "[-1/4, -1/4, 1/4]@[7/8, 3/8, 5/8]": ("bond_002", 3),
            "[-1/4, -1/4, 1/4]@[7/8, 7/8, 1/8]": ("bond_002", 4),
            "[1/4, -1/4, -1/4]@[1/8, 3/8, 3/8]": ("bond_003", 1),
            "[1/4, -1/4, -1/4]@[1/8, 7/8, 7/8]": ("bond_003", 2),
            "[1/4, -1/4, -1/4]@[5/8, 3/8, 7/8]": ("bond_003", 3),
            "[1/4, -1/4, -1/4]@[5/8, 7/8, 3/8]": ("bond_003", 4),
            "[-1/4, 1/4, -1/4]@[3/8, 1/8, 3/8]": ("bond_004", 1),
            "[-1/4, 1/4, -1/4]@[3/8, 5/8, 7/8]": ("bond_004", 2),
            "[-1/4, 1/4, -1/4]@[7/8, 1/8, 7/8]": ("bond_004", 3),
            "[-1/4, 1/4, -1/4]@[7/8, 5/8, 3/8]": ("bond_004", 4),
        },
    },
    "data": {
        "plus_set": ["[0, 0, 0]", "[0, 1/2, 1/2]", "[1/2, 0, 1/2]", "[1/2, 1/2, 0]"],
        "cluster_site": {"S_001": ["site_001"], "S_002": ["site_002"]},
        "cluster_bond": {"B_001": ["bond_001", "bond_002", "bond_003", "bond_004"]},
        "site": {
            "site_001": (
                "[0, 0, 0]",
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23],
                (0, 0),
            ),
            "site_002": (
                "[1/4, 1/4, 1/4]",
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23],
                (1, 1),
            ),
        },
        "bond": {
            "bond_001": (
                "[1/4, 1/4, 1/4]@[1/8, 1/8, 1/8]",
                [0, 4, 8, 15, 16, 17],
                (1, 0),
                "[1/4, 1/4, 1/4]",
                "[0, 0, 0];[1/4, 1/4, 1/4]",
            ),
            "bond_002": (
                "[-1/4, -1/4, 1/4]@[3/8, 3/8, 1/8]",
                [1, 5, 10, 12, 20, 22],
                (1, 0),
                "[-1/4, -1/4, 1/4]",
                "[1/2, 1/2, 0];[1/4, 1/4, 1/4]",
            ),
            "bond_003": (
                "[1/4, -1/4, -1/4]@[1/8, 3/8, 3/8]",
                [2, 6, 11, 14, 18, 23],
                (1, 0),
                "[1/4, -1/4, -1/4]",
                "[0, 1/2, 1/2];[1/4, 1/4, 1/4]",
            ),
            "bond_004": (
                "[-1/4, 1/4, -1/4]@[3/8, 1/8, 3/8]",
                [3, 7, 9, 13, 19, 21],
                (1, 0),
                "[-1/4, 1/4, -1/4]",
                "[1/2, 0, 1/2];[1/4, 1/4, 1/4]",
            ),
        },
        "cluster_atomic": {(0, 0): [(0, 0, "M_001")], (1, 1): [(3, 3, "M_001")], (1, 0): [(3, 0, "M_001")]},
        "atomic_braket": {"M_001": (["px", "py", "pz"], ["px", "py", "pz"])},
    },
    "detail": {
        "rep_bond_all": {
            "Ga_As": [
                {},
                {"Ga:As:1:1": ("[0, 0, 0];[1/4, 1/4, 1/4]", "16e", "D", 1, ".3m")},
                {"Ga:As:2:1": ("[0, 0, 1];[1/4, 1/4, 1/4]", "16e", "D", 2, ".3m")},
                {"Ga:As:3:1": ("[0, 0, 0];[1/4, 3/4, 3/4]", "16e", "D", 3, ".3m")},
                {
                    "Ga:As:4:1": ("[0, 0, 1];[3/4, 3/4, 1/4]", "16e", "D", 4, ".3m"),
                    "Ga:As:4:2": ("[-1/2, 0, 1/2];[3/4, 1/4, 3/4]", "16e", "D", 4, ".3m"),
                },
                {"Ga:As:5:1": ("[-1/2, 0, 1/2];[3/4, 1/4, -1/4]", "16e", "D", 5, ".3m")},
                {"Ga:As:6:1": ("[-1/2, 0, 1/2];[3/4, 3/4, 5/4]", "16e", "D", 6, ".3m")},
                {
                    "Ga:As:7:1": ("[-1/2, -1/2, 0];[3/4, 3/4, 1/4]", "16e", "D", 7, ".3m"),
                    "Ga:As:7:2": ("[-1/2, 0, 1/2];[5/4, 1/4, 1/4]", "16e", "D", 7, ".3m"),
                },
                {
                    "Ga:As:8:1": ("[-1/2, -1/2, 1];[3/4, 3/4, 1/4]", "16e", "D", 8, ".3m"),
                    "Ga:As:8:2": ("[-1/2, 0, 1/2];[5/4, 1/4, 5/4]", "16e", "D", 8, ".3m"),
                },
                {"Ga:As:9:1": ("[-1/2, 0, 1/2];[5/4, 3/4, -1/4]", "16e", "D", 9, ".3m")},
            ]
        },
        "cell_range": (-2, 3, -2, 3, -2, 3),
        "max_neighbor": 10,
        "A": "[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]",
        "version": "1.1.14",
    },
}
