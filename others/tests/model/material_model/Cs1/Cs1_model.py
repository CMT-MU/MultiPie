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
Cs1 = {
    "info": {
        "model": "Cs1",
        "molecule": False,
        "group": ("Cs^1", "space group No. 6 : Cs^1 / Pm (b-axis) : PG Cs"),
        "crystal": "monoclinic",
        "cell": {"a": 1.0, "b": 1.0, "c": 1.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
        "volume": 1.0,
        "a1": "[1.0, 0.0, 0.0]",
        "a2": "[0.0, 1.0, 0.0]",
        "a3": "[0.0, 0.0, 1.0]",
        "option": {"view": [0, 0, 1], "view_mode": "standard", "output": "Cs1", "minimal_samb": True},
        "generate": {"model_type": "tight_binding", "time_reversal_type": "electric", "irrep": ["A'"]},
        "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
        "k_path": "Γ-X",
        "dimension": 4,
        "spinful": True,
        "orbital": ["(px,U)", "(px,D)", "(py,U)", "(py,D)"],
        "ket": ["(px,U)@A_1", "(px,D)@A_1", "(py,U)@A_1", "(py,D)@A_1"],
        "ket_site": ["A_1"],
        "site": {"A": ("[0,0,0]", ["px", "py"])},
        "rep_site": {"A": ("[0, 0, 0]", "1a", [["(px,U)", "(px,D)", "(py,U)", "(py,D)"]], "m")},
        "cell_site": {"A_1": ("[0, 0, 0]", "[1,2]")},
        "bond": [("A", "A", [1, 2])],
        "rep_bond": {
            "A:A:1:1": ("[0, 0, 0];[0, 1, 0]", "1b", "ND", 1, "m"),
            "A:A:1:2": ("[0, 0, 1];[0, 0, 0]", "1a", "D", 1, "m"),
            "A:A:1:3": ("[0, 0, 0];[1, 0, 0]", "1a", "D", 1, "m"),
            "A:A:2:1": ("[1, 0, 1];[0, 0, 0]", "1a", "D", 2, "m"),
            "A:A:2:2": ("[0, 0, 0];[1, 1, 0]", "1b", "D", 2, "m"),
            "A:A:2:3": ("[1, 0, 0];[0, 0, 1]", "1a", "D", 2, "m"),
            "A:A:2:4": ("[0, 0, 1];[0, 1, 0]", "1b", "D", 2, "m"),
        },
        "cell_bond": {
            "A:A:1:1_1": ("[0, 1, 0]@[0, 1/2, 0]", "[1,-2]"),
            "A:A:1:2_1": ("[0, 0, -1]@[0, 0, 1/2]", "[1,2]"),
            "A:A:1:3_1": ("[1, 0, 0]@[1/2, 0, 0]", "[1,2]"),
            "A:A:2:1_1": ("[-1, 0, -1]@[1/2, 0, 1/2]", "[1,2]"),
            "A:A:2:2_1": ("[1, 1, 0]@[1/2, 1/2, 0]", "[1]"),
            "A:A:2:2_2": ("[1, -1, 0]@[1/2, 1/2, 0]", "[2]"),
            "A:A:2:3_1": ("[-1, 0, 1]@[1/2, 0, 1/2]", "[1,2]"),
            "A:A:2:4_1": ("[0, 1, -1]@[0, 1/2, 1/2]", "[1]"),
            "A:A:2:4_2": ("[0, -1, -1]@[0, 1/2, 1/2]", "[2]"),
        },
    },
    "name": {
        "alias": {
            "S_001": "A",
            "A": "S_001",
            "B_001": "A:A:1:1",
            "A:A:1:1": "B_001",
            "B_002": "A:A:1:2",
            "A:A:1:2": "B_002",
            "B_003": "A:A:1:3",
            "A:A:1:3": "B_003",
            "B_004": "A:A:2:1",
            "A:A:2:1": "B_004",
            "B_005": "A:A:2:2",
            "A:A:2:2": "B_005",
            "B_006": "A:A:2:3",
            "A:A:2:3": "B_006",
            "B_007": "A:A:2:4",
            "A:A:2:4": "B_007",
        },
        "site": {"site_001": ("A", 1)},
        "site_name": {"[0, 0, 0]": ("site_001", 1)},
        "bond": {
            "bond_001": ("A:A:1:1", 1),
            "bond_002": ("A:A:1:2", 1),
            "bond_003": ("A:A:1:3", 1),
            "bond_004": ("A:A:2:1", 1),
            "bond_005": ("A:A:2:2", 1),
            "bond_006": ("A:A:2:2", 2),
            "bond_007": ("A:A:2:3", 1),
            "bond_008": ("A:A:2:4", 1),
            "bond_009": ("A:A:2:4", 2),
        },
        "bond_name": {
            "[0, 1, 0]@[0, 1/2, 0]": ("bond_001", 1),
            "[0, 0, -1]@[0, 0, 1/2]": ("bond_002", 1),
            "[1, 0, 0]@[1/2, 0, 0]": ("bond_003", 1),
            "[-1, 0, -1]@[1/2, 0, 1/2]": ("bond_004", 1),
            "[1, 1, 0]@[1/2, 1/2, 0]": ("bond_005", 1),
            "[1, -1, 0]@[1/2, 1/2, 0]": ("bond_006", 1),
            "[-1, 0, 1]@[1/2, 0, 1/2]": ("bond_007", 1),
            "[0, 1, -1]@[0, 1/2, 1/2]": ("bond_008", 1),
            "[0, -1, -1]@[0, 1/2, 1/2]": ("bond_009", 1),
        },
    },
    "data": {
        "plus_set": ["[0, 0, 0]"],
        "cluster_site": {"S_001": ["site_001"]},
        "cluster_bond": {
            "B_001": ["bond_001"],
            "B_002": ["bond_002"],
            "B_003": ["bond_003"],
            "B_004": ["bond_004"],
            "B_005": ["bond_005", "bond_006"],
            "B_006": ["bond_007"],
            "B_007": ["bond_008", "bond_009"],
        },
        "site": {"site_001": ("[0, 0, 0]", [0, 1], (0, 0))},
        "bond": {
            "bond_001": ("[0, 1, 0]@[0, 1/2, 0]", [0, -1], (0, 0), "[0, 1, 0]", "[0, 0, 0];[0, 1, 0]"),
            "bond_002": ("[0, 0, -1]@[0, 0, 1/2]", [0, 1], (0, 0), "[0, 0, -1]", "[0, 0, 1];[0, 0, 0]"),
            "bond_003": ("[1, 0, 0]@[1/2, 0, 0]", [0, 1], (0, 0), "[1, 0, 0]", "[0, 0, 0];[1, 0, 0]"),
            "bond_004": ("[-1, 0, -1]@[1/2, 0, 1/2]", [0, 1], (0, 0), "[-1, 0, -1]", "[1, 0, 1];[0, 0, 0]"),
            "bond_005": ("[1, 1, 0]@[1/2, 1/2, 0]", [0], (0, 0), "[1, 1, 0]", "[0, 0, 0];[1, 1, 0]"),
            "bond_006": ("[1, -1, 0]@[1/2, 1/2, 0]", [1], (0, 0), "[1, -1, 0]", "[0, 1, 0];[1, 0, 0]"),
            "bond_007": ("[-1, 0, 1]@[1/2, 0, 1/2]", [0, 1], (0, 0), "[-1, 0, 1]", "[1, 0, 0];[0, 0, 1]"),
            "bond_008": ("[0, 1, -1]@[0, 1/2, 1/2]", [0], (0, 0), "[0, 1, -1]", "[0, 0, 1];[0, 1, 0]"),
            "bond_009": ("[0, -1, -1]@[0, 1/2, 1/2]", [1], (0, 0), "[0, -1, -1]", "[0, 1, 1];[0, 0, 0]"),
        },
        "cluster_atomic": {(0, 0): [(0, 0, "M_001")]},
        "atomic_braket": {"M_001": (["(px,U)", "(px,D)", "(py,U)", "(py,D)"], ["(px,U)", "(px,D)", "(py,U)", "(py,D)"])},
    },
    "detail": {
        "rep_bond_all": {
            "A_A": [
                {},
                {
                    "A:A:1:1": ("[0, 0, 0];[0, 1, 0]", "1b", "ND", 1, "m"),
                    "A:A:1:2": ("[0, 0, 1];[0, 0, 0]", "1a", "D", 1, "m"),
                    "A:A:1:3": ("[0, 0, 0];[1, 0, 0]", "1a", "D", 1, "m"),
                },
                {
                    "A:A:2:1": ("[1, 0, 1];[0, 0, 0]", "1a", "D", 2, "m"),
                    "A:A:2:2": ("[0, 0, 0];[1, 1, 0]", "1b", "D", 2, "m"),
                    "A:A:2:3": ("[1, 0, 0];[0, 0, 1]", "1a", "D", 2, "m"),
                    "A:A:2:4": ("[0, 0, 1];[0, 1, 0]", "1b", "D", 2, "m"),
                },
                {"A:A:3:1": ("[1, 0, 0];[0, 1, 1]", "1b", "D", 3, "m"), "A:A:3:2": ("[0, 0, 0];[1, 1, 1]", "1b", "D", 3, "m")},
                {
                    "A:A:4:1": ("[-1, 0, 0];[1, 0, 0]", "1a", "D", 4, "m"),
                    "A:A:4:2": ("[0, 0, -1];[0, 0, 1]", "1a", "D", 4, "m"),
                    "A:A:4:3": ("[0, -1, 0];[0, 1, 0]", "1a", "ND", 4, "m"),
                },
                {
                    "A:A:5:1": ("[0, 0, -1];[0, 1, 1]", "1b", "D", 5, "m"),
                    "A:A:5:2": ("[-1, 0, 1];[1, 0, 0]", "1a", "D", 5, "m"),
                    "A:A:5:3": ("[-1, 0, 0];[1, 1, 0]", "1b", "D", 5, "m"),
                    "A:A:5:4": ("[0, 0, -1];[1, 0, 1]", "1a", "D", 5, "m"),
                    "A:A:5:5": ("[1, 0, 1];[-1, 0, 0]", "1a", "D", 5, "m"),
                    "A:A:5:6": ("[0, -1, 0];[1, 1, 0]", "1a", "D", 5, "m"),
                    "A:A:5:7": ("[1, 0, -1];[0, 0, 1]", "1a", "D", 5, "m"),
                    "A:A:5:8": ("[0, -1, 1];[0, 1, 0]", "1a", "D", 5, "m"),
                },
                {
                    "A:A:6:1": ("[0, -1, 1];[1, 1, 0]", "1a", "D", 6, "m"),
                    "A:A:6:2": ("[1, 0, 1];[-1, 1, 0]", "1b", "D", 6, "m"),
                    "A:A:6:3": ("[0, 0, -1];[1, 1, 1]", "1b", "D", 6, "m"),
                    "A:A:6:4": ("[-1, 0, 1];[1, 1, 0]", "1b", "D", 6, "m"),
                    "A:A:6:5": ("[1, -1, 1];[0, 1, 0]", "1a", "D", 6, "m"),
                    "A:A:6:6": ("[1, 0, -1];[0, 1, 1]", "1b", "D", 6, "m"),
                },
                {
                    "A:A:7:1": ("[0, -1, -1];[0, 1, 1]", "1a", "D", 7, "m"),
                    "A:A:7:2": ("[-1, 0, -1];[1, 0, 1]", "1a", "D", 7, "m"),
                    "A:A:7:3": ("[1, -1, 0];[-1, 1, 0]", "1a", "D", 7, "m"),
                    "A:A:7:4": ("[1, 0, -1];[-1, 0, 1]", "1a", "D", 7, "m"),
                },
                {
                    "A:A:8:1": ("[0, -1, 1];[1, 1, -1]", "1a", "D", 8, "m"),
                    "A:A:8:2": ("[0, -1, -1];[1, 1, 1]", "1a", "D", 8, "m"),
                    "A:A:8:3": ("[-1, -1, 1];[1, 1, 0]", "1a", "D", 8, "m"),
                    "A:A:8:4": ("[-1, -1, 0];[1, 1, 1]", "1a", "D", 8, "m"),
                    "A:A:8:5": ("[-1, 0, -1];[1, 1, 1]", "1b", "D", 8, "m"),
                    "A:A:8:6": ("[1, 0, -1];[-1, 1, 1]", "1b", "D", 8, "m"),
                },
                {
                    "A:A:9:1": ("[1, -1, -1];[-1, 1, 1]", "1a", "D", 9, "m"),
                    "A:A:9:2": ("[-1, -1, -1];[1, 1, 1]", "1a", "D", 9, "m"),
                },
            ]
        },
        "cell_range": (-2, 3, -2, 3, -2, 3),
        "max_neighbor": 10,
        "A": "[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]",
        "version": "1.1.10",
    },
}
