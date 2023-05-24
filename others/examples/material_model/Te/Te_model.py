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
Te = {
    "info": {
        "model": "Te",
        "molecule": False,
        "group": ("D3^4", "space group No. 152 : D3^4 / P3_121 : PG D3-1"),
        "crystal": "trigonal",
        "cell": {"a": 4.458, "b": 4.458, "c": 5.925, "alpha": 90.0, "beta": 90.0, "gamma": 120.0},
        "volume": 101.97626811993862,
        "a1": "[4.458, 0.0, 0.0]",
        "a2": "[-2.229, 3.86074125007103, 0.0]",
        "a3": "[0.0, 0.0, 5.925]",
        "option": {"view": None, "view_mode": "standard", "output": "Te", "minimal_samb": True},
        "generate": {"model_type": "phonon", "time_reversal_type": "electric", "irrep": ["A1"]},
        "k_point": {
            "Γ": "[0, 0, 0]",
            "A": "[0, 0, 1/2]",
            "M": "[1/2, 0, 0]",
            "K": "[1/3, 1/3, 0]",
            "H": "[1/3, 1/3, 1/2]",
            "L": "[1/2, 0, 1/2]",
        },
        "k_path": "A-Γ-H-A-L-H-K-Γ-M-K",
        "dimension": 9,
        "spinful": False,
        "orbital": ["px", "py", "pz"],
        "ket": ["px@A_1", "py@A_1", "pz@A_1", "px@A_2", "py@A_2", "pz@A_2", "px@A_3", "py@A_3", "pz@A_3"],
        "ket_site": ["A_1", "A_2", "A_3"],
        "site": {"A": ("[0.274,0,1/3]", ["px", "py", "pz"])},
        "rep_site": {"A": ("[137/500, 0, 1/3]", "3a", [["px", "py", "pz"]])},
        "cell_site": {
            "A_1": ("[137/500, 0, 1/3]", "[1,2]"),
            "A_2": ("[363/500, 363/500, 0]", "[3,6]"),
            "A_3": ("[0, 137/500, 2/3]", "[4,5]"),
        },
        "bond": [("A", "A", [1])],
        "rep_bond": {"A:A:1:1": ("[137/500, 1, 1/3];[-137/500, 363/500, 0]", "3b", "ND", 1)},
        "cell_bond": {
            "A:A:1:1_1": ("[-137/250, -137/500, -1/3]@[0, 863/1000, 1/6]", "[1,-3]"),
            "A:A:1:1_2": ("[-137/500, 137/500, 1/3]@[137/1000, 137/1000, 1/2]", "[2,-5]"),
            "A:A:1:1_3": ("[137/500, 137/250, -1/3]@[863/1000, 0, 5/6]", "[-4,6]"),
        },
    },
    "name": {
        "alias": {"S_001": "A", "A": "S_001", "B_001": "A:A:1:1", "A:A:1:1": "B_001"},
        "site": {"site_001": ("A", 1), "site_002": ("A", 2), "site_003": ("A", 3)},
        "site_name": {
            "[137/500, 0, 1/3]": ("site_001", 1),
            "[363/500, 363/500, 0]": ("site_002", 1),
            "[0, 137/500, 2/3]": ("site_003", 1),
        },
        "bond": {"bond_001": ("A:A:1:1", 1), "bond_002": ("A:A:1:1", 2), "bond_003": ("A:A:1:1", 3)},
        "bond_name": {
            "[-137/250, -137/500, -1/3]@[0, 863/1000, 1/6]": ("bond_001", 1),
            "[-137/500, 137/500, 1/3]@[137/1000, 137/1000, 1/2]": ("bond_002", 1),
            "[137/500, 137/250, -1/3]@[863/1000, 0, 5/6]": ("bond_003", 1),
        },
    },
    "data": {
        "plus_set": ["[0, 0, 0]"],
        "cluster_site": {"S_001": ["site_001", "site_002", "site_003"]},
        "cluster_bond": {"B_001": ["bond_001", "bond_002", "bond_003"]},
        "site": {
            "site_001": ("[137/500, 0, 1/3]", [0, 1], (0, 0)),
            "site_002": ("[363/500, 363/500, 0]", [2, 5], (1, 1)),
            "site_003": ("[0, 137/500, 2/3]", [3, 4], (2, 2)),
        },
        "bond": {
            "bond_001": (
                "[-137/250, -137/500, -1/3]@[0, 863/1000, 1/6]",
                [0, -2],
                (0, 1),
                "[-137/250, -137/500, -1/3]",
                "[137/500, 1, 1/3];[-137/500, 363/500, 0]",
            ),
            "bond_002": (
                "[-137/500, 137/500, 1/3]@[137/1000, 137/1000, 1/2]",
                [1, -4],
                (0, 2),
                "[-137/500, 137/500, 1/3]",
                "[137/500, 0, 1/3];[0, 137/500, 2/3]",
            ),
            "bond_003": (
                "[137/500, 137/250, -1/3]@[863/1000, 0, 5/6]",
                [-3, 5],
                (1, 2),
                "[137/500, 137/250, -1/3]",
                "[363/500, -137/500, 1];[1, 137/500, 2/3]",
            ),
        },
        "cluster_atomic": {
            (0, 0): [(0, 0, "M_001")],
            (1, 1): [(3, 3, "M_001")],
            (2, 2): [(6, 6, "M_001")],
            (0, 1): [(0, 3, "M_001")],
            (0, 2): [(0, 6, "M_001")],
            (1, 2): [(3, 6, "M_001")],
        },
        "atomic_braket": {"M_001": (["px", "py", "pz"], ["px", "py", "pz"])},
    },
    "detail": {
        "rep_bond_all": {
            "A_A": [
                {},
                {"A:A:1:1": ("[137/500, 1, 1/3];[-137/500, 363/500, 0]", "3b", "ND", 1)},
                {"A:A:2:1": ("[137/500, 0, 1/3];[363/500, 363/500, 0]", "6c", "D", 2)},
                {
                    "A:A:3:1": ("[137/500, 0, 1/3];[637/500, 1, 1/3]", "6c", "D", 3),
                    "A:A:3:2": ("[637/500, 0, 1/3];[137/500, 0, 1/3]", "3a", "D", 3),
                },
                {"A:A:4:1": ("[137/500, 1, 1/3];[-137/500, 363/500, 1]", "3a", "ND", 4)},
                {"A:A:5:1": ("[137/500, 0, 1/3];[363/500, 363/500, 1]", "6c", "D", 5)},
                {"A:A:6:1": ("[137/500, 0, 1/3];[-137/500, 363/500, 0]", "3b", "D", 6)},
                {"A:A:7:1": ("[137/500, 0, 1/3];[137/500, 0, 4/3]", "3b", "ND", 7)},
                {"A:A:8:1": ("[-363/500, 0, 1/3];[363/500, 363/500, 0]", "3b", "ND", 8)},
                {"A:A:9:1": ("[137/500, 0, 1/3];[-137/500, 363/500, 1]", "3a", "D", 9)},
            ]
        },
        "cell_range": (-2, 3, -2, 3, -2, 3),
        "max_neighbor": 10,
        "A": "[[4.458, -2.229, 0.0], [0.0, 3.86074125007103, 0.0], [0.0, 0.0, 5.925]]",
        "version": "1.1.1",
    },
}
