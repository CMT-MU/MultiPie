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
graphene = {
    "info": {
        "model": "graphene",
        "molecule": False,
        "group": ("D6h^1", "space group No. 191 : D6h^1 / P6/mmm : PG D6h"),
        "crystal": "hexagonal",
        "cell": {"a": 1.0, "b": 1.0, "c": 1.0, "alpha": 90.0, "beta": 90.0, "gamma": 120.0},
        "volume": 0.8660254037844388,
        "a1": "[1.0, 0.0, 0.0]",
        "a2": "[-0.5, 0.86602540378444, 0.0]",
        "a3": "[0.0, 0.0, 1.0]",
        "option": {"view": [0, 0, 1], "view_mode": "standard", "output": "graphene", "minimal_samb": True},
        "generate": {"model_type": "tight_binding", "time_reversal_type": "electric", "irrep": ["A1g"]},
        "k_point": {"Γ": "[0, 0, 0]", "M": "[1/2, 0, 0]", "K": "[1/3, 1/3, 0]", "K'": "[-1/3, -1/3, 0]"},
        "k_path": "Γ-M-K-Γ-K'",
        "dimension": 4,
        "spinful": True,
        "orbital": ["(pz,U)", "(pz,D)"],
        "ket": ["(pz,U)@A_1", "(pz,D)@A_1", "(pz,U)@A_2", "(pz,D)@A_2"],
        "ket_site": ["A_1", "A_2"],
        "site": {"A": ("[1/3,2/3,0]", "pz")},
        "rep_site": {"A": ("[1/3, 2/3, 0]", "2c", [["(pz,U)", "(pz,D)"]])},
        "cell_site": {
            "A_1": ("[1/3, 2/3, 0]", "[1,6,7,8,9,10,14,15,16,17,23,24]"),
            "A_2": ("[2/3, 1/3, 0]", "[2,3,4,5,11,12,13,18,19,20,21,22]"),
        },
        "bond": [("A", "A", [1])],
        "rep_bond": {"A:A:1:1": ("[1/3, -1/3, 0];[2/3, 1/3, 0]", "3f", "ND", 1)},
        "cell_bond": {
            "A:A:1:1_1": ("[1/3, 2/3, 0]@[1/2, 0, 0]", "[1,-2,-3,6,-13,14,17,-18]"),
            "A:A:1:1_2": ("[1/3, -1/3, 0]@[1/2, 1/2, 0]", "[-4,7,10,-11,15,-19,-22,23]"),
            "A:A:1:1_3": ("[-2/3, -1/3, 0]@[0, 1/2, 0]", "[-5,8,9,-12,16,-20,-21,24]"),
        },
    },
    "name": {
        "alias": {"S_001": "A", "A": "S_001", "B_001": "A:A:1:1", "A:A:1:1": "B_001"},
        "site": {"site_001": ("A", 1), "site_002": ("A", 2)},
        "site_name": {"[1/3, 2/3, 0]": ("site_001", 1), "[2/3, 1/3, 0]": ("site_002", 1)},
        "bond": {"bond_001": ("A:A:1:1", 1), "bond_002": ("A:A:1:1", 2), "bond_003": ("A:A:1:1", 3)},
        "bond_name": {
            "[1/3, 2/3, 0]@[1/2, 0, 0]": ("bond_001", 1),
            "[1/3, -1/3, 0]@[1/2, 1/2, 0]": ("bond_002", 1),
            "[-2/3, -1/3, 0]@[0, 1/2, 0]": ("bond_003", 1),
        },
    },
    "data": {
        "plus_set": ["[0, 0, 0]"],
        "cluster_site": {"S_001": ["site_001", "site_002"]},
        "cluster_bond": {"B_001": ["bond_001", "bond_002", "bond_003"]},
        "site": {
            "site_001": ("[1/3, 2/3, 0]", [0, 5, 6, 7, 8, 9, 13, 14, 15, 16, 22, 23], (0, 0)),
            "site_002": ("[2/3, 1/3, 0]", [1, 2, 3, 4, 10, 11, 12, 17, 18, 19, 20, 21], (1, 1)),
        },
        "bond": {
            "bond_001": (
                "[1/3, 2/3, 0]@[1/2, 0, 0]",
                [0, -1, -2, 5, -12, 13, 16, -17],
                (0, 1),
                "[1/3, 2/3, 0]",
                "[1/3, -1/3, 0];[2/3, 1/3, 0]",
            ),
            "bond_002": (
                "[1/3, -1/3, 0]@[1/2, 1/2, 0]",
                [-3, 6, 9, -10, 14, -18, -21, 22],
                (0, 1),
                "[1/3, -1/3, 0]",
                "[1/3, 2/3, 0];[2/3, 1/3, 0]",
            ),
            "bond_003": (
                "[-2/3, -1/3, 0]@[0, 1/2, 0]",
                [-4, 7, 8, -11, 15, -19, -20, 23],
                (0, 1),
                "[-2/3, -1/3, 0]",
                "[1/3, 2/3, 0];[-1/3, 1/3, 0]",
            ),
        },
        "cluster_atomic": {(0, 0): [(0, 0, "M_001")], (1, 1): [(2, 2, "M_001")], (0, 1): [(0, 2, "M_001")]},
        "atomic_braket": {"M_001": (["(pz,U)", "(pz,D)"], ["(pz,U)", "(pz,D)"])},
    },
    "detail": {
        "rep_bond_all": {
            "A_A": [
                {},
                {"A:A:1:1": ("[1/3, -1/3, 0];[2/3, 1/3, 0]", "3f", "ND", 1)},
                {
                    "A:A:2:1": ("[1/3, -1/3, 0];[1/3, 2/3, 0]", "6l", "ND", 2),
                    "A:A:2:2": ("[1/3, 2/3, 0];[1/3, 2/3, 1]", "2d", "ND", 2),
                },
                {
                    "A:A:3:1": ("[-2/3, -1/3, 0];[2/3, 1/3, 0]", "1a", "ND", 3),
                    "A:A:3:2": ("[1/3, -1/3, 0];[2/3, 1/3, 1]", "3g", "ND", 3),
                },
                {"A:A:4:1": ("[1/3, -1/3, 0];[1/3, 2/3, 1]", "6m", "ND", 4)},
                {
                    "A:A:5:1": ("[-2/3, -1/3, 0];[2/3, 1/3, 1]", "1b", "ND", 5),
                    "A:A:5:2": ("[-2/3, -1/3, 0];[2/3, 4/3, 0]", "3f", "ND", 5),
                },
                {"A:A:6:1": ("[-2/3, -1/3, 0];[4/3, 2/3, 0]", "6l", "D", 6)},
                {"A:A:7:1": ("[-2/3, -1/3, 0];[2/3, 4/3, 1]", "3g", "ND", 7)},
                {
                    "A:A:8:1": ("[1/3, 2/3, -1];[1/3, 2/3, 1]", "2c", "ND", 8),
                    "A:A:8:2": ("[-2/3, -1/3, 0];[4/3, 2/3, 1]", "6m", "D", 8),
                    "A:A:8:3": ("[-2/3, -1/3, 0];[4/3, 5/3, 0]", "2c", "ND", 8),
                },
                {
                    "A:A:9:1": ("[1/3, -1/3, -1];[2/3, 1/3, 1]", "3f", "ND", 9),
                    "A:A:9:2": ("[-2/3, -1/3, 0];[5/3, 1/3, 0]", "3f", "ND", 9),
                },
            ]
        },
        "cell_range": (-2, 3, -2, 3, -2, 3),
        "max_neighbor": 10,
        "A": "[[1.0, -0.5, 0.0], [0.0, 0.86602540378444, 0.0], [0.0, 0.0, 1.0]]",
        "version": "1.1.2",
    },
}
