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
C4v1 = {
    "info": {
        "model": "C4v1",
        "molecule": False,
        "group": ("C4v^1", "space group No. 99 : C4v^1 / P4mm : PG C4v"),
        "crystal": "tetragonal",
        "cell": {"a": 1.0, "b": 1.0, "c": 1.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
        "volume": 1.0,
        "a1": "[1.0, 0.0, 0.0]",
        "a2": "[0.0, 1.0, 0.0]",
        "a3": "[0.0, 0.0, 1.0]",
        "option": {"view": None, "view_mode": "standard", "output": "C4v1", "minimal_samb": True},
        "generate": {
            "fourier_transform": False,
            "model_type": "tight_binding",
            "time_reversal_type": "electric",
            "irrep": ["A1"],
            "toroidal_priority": False,
        },
        "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
        "k_path": "Γ-X",
        "dimension": 2,
        "spinful": True,
        "orbital": ["(s,U)", "(s,D)"],
        "ket": ["(s,U)@A_1", "(s,D)@A_1"],
        "ket_site": ["A_1"],
        "site": {"A": ("[0,0,0]", "s")},
        "rep_site": {"A": ("[0, 0, 0]", "1a", [["(s,U)", "(s,D)"]], "4mm")},
        "cell_site": {"A_1": ("[0, 0, 0]", "[1,2,3,4,5,6,7,8]")},
        "bond": [("A", "A", 1)],
        "rep_bond": {
            "A:A:1:1": ("[0, 0, 0];[0, 1, 0]", "2c", "ND", 1, "2mm."),
            "A:A:1:2": ("[0, 0, 0];[0, 0, 1]", "1a", "D", 1, "4mm"),
        },
        "cell_bond": {
            "A:A:1:1_1": ("[0, 1, 0]@[0, 1/2, 0]", "[1,-2,5,-6]"),
            "A:A:1:1_2": ("[1, 0, 0]@[1/2, 0, 0]", "[-3,4,-7,8]"),
            "A:A:1:2_1": ("[0, 0, 1]@[0, 0, 1/2]", "[1,2,3,4,5,6,7,8]"),
        },
    },
    "name": {
        "alias": {"S_001": "A", "A": "S_001", "B_001": "A:A:1:1", "A:A:1:1": "B_001", "B_002": "A:A:1:2", "A:A:1:2": "B_002"},
        "site": {"site_001": ("A", 1)},
        "site_name": {"[0, 0, 0]": ("site_001", 1)},
        "bond": {"bond_001": ("A:A:1:1", 1), "bond_002": ("A:A:1:1", 2), "bond_003": ("A:A:1:2", 1)},
        "bond_name": {
            "[0, 1, 0]@[0, 1/2, 0]": ("bond_001", 1),
            "[1, 0, 0]@[1/2, 0, 0]": ("bond_002", 1),
            "[0, 0, 1]@[0, 0, 1/2]": ("bond_003", 1),
        },
    },
    "data": {
        "plus_set": ["[0, 0, 0]"],
        "cluster_site": {"S_001": ["site_001"]},
        "cluster_bond": {"B_001": ["bond_001", "bond_002"], "B_002": ["bond_003"]},
        "site": {"site_001": ("[0, 0, 0]", [0, 1, 2, 3, 4, 5, 6, 7], (0, 0))},
        "bond": {
            "bond_001": ("[0, 1, 0]@[0, 1/2, 0]", [0, -1, 4, -5], (0, 0), "[0, 1, 0]", "[0, 0, 0];[0, 1, 0]"),
            "bond_002": ("[1, 0, 0]@[1/2, 0, 0]", [-2, 3, -6, 7], (0, 0), "[1, 0, 0]", "[0, 0, 0];[1, 0, 0]"),
            "bond_003": ("[0, 0, 1]@[0, 0, 1/2]", [0, 1, 2, 3, 4, 5, 6, 7], (0, 0), "[0, 0, 1]", "[0, 0, 0];[0, 0, 1]"),
        },
        "cluster_atomic": {(0, 0): [(0, 0, "M_001")]},
        "atomic_braket": {"M_001": (["(s,U)", "(s,D)"], ["(s,U)", "(s,D)"])},
    },
    "detail": {
        "rep_bond_all": {
            "A_A": [
                {},
                {
                    "A:A:1:1": ("[0, 0, 0];[0, 1, 0]", "2c", "ND", 1, "2mm."),
                    "A:A:1:2": ("[0, 0, 0];[0, 0, 1]", "1a", "D", 1, "4mm"),
                },
                {
                    "A:A:2:1": ("[0, 0, 0];[1, 1, 0]", "1b", "ND", 2, "4mm"),
                    "A:A:2:2": ("[0, 0, 1];[0, 1, 0]", "2c", "D", 2, "2mm."),
                },
                {"A:A:3:1": ("[0, 0, 1];[1, 1, 0]", "1b", "D", 3, "4mm")},
                {
                    "A:A:4:1": ("[0, 0, -1];[0, 0, 1]", "1a", "D", 4, "4mm"),
                    "A:A:4:2": ("[-1, 0, 0];[1, 0, 0]", "1a", "ND", 4, "4mm"),
                },
                {
                    "A:A:5:1": ("[-1, 0, 0];[1, 1, 0]", "2c", "ND", 5, "2mm."),
                    "A:A:5:2": ("[0, 0, -1];[0, 1, 1]", "2c", "D", 5, "2mm."),
                    "A:A:5:3": ("[-1, 0, 0];[1, 0, 1]", "1a", "D", 5, "4mm"),
                },
                {
                    "A:A:6:1": ("[0, 0, 1];[1, 1, -1]", "1b", "D", 6, "4mm"),
                    "A:A:6:2": ("[-1, 0, 0];[1, 1, 1]", "2c", "D", 6, "2mm."),
                },
                {
                    "A:A:7:1": ("[-1, -1, 0];[1, 1, 0]", "1a", "ND", 7, "4mm"),
                    "A:A:7:2": ("[-1, 0, -1];[1, 0, 1]", "1a", "D", 7, "4mm"),
                },
                {
                    "A:A:8:1": ("[-1, 0, 1];[1, 1, -1]", "2c", "D", 8, "2mm."),
                    "A:A:8:2": ("[-1, -1, 0];[1, 1, 1]", "1a", "D", 8, "4mm"),
                },
                {"A:A:9:1": ("[-1, -1, -1];[1, 1, 1]", "1a", "D", 9, "4mm")},
            ]
        },
        "cell_range": (-2, 3, -2, 3, -2, 3),
        "max_neighbor": 10,
        "A": "[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]",
        "version": "1.1.14",
    },
}
