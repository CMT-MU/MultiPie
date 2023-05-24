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
C3v1 = {
    "info": {
        "model": "C3v1",
        "molecule": False,
        "group": ("C3v^1", "space group No. 156 : C3v^1 / P3m1 : PG C3v"),
        "crystal": "trigonal",
        "cell": {"a": 1.0, "b": 1.0, "c": 1.0, "alpha": 90.0, "beta": 90.0, "gamma": 120.0},
        "volume": 0.8660254037844388,
        "a1": "[1.0, 0.0, 0.0]",
        "a2": "[-0.5, 0.86602540378444, 0.0]",
        "a3": "[0.0, 0.0, 1.0]",
        "option": {"view": [0, 0, 1], "view_mode": "standard", "output": "C3v1", "minimal_samb": True},
        "generate": {"model_type": "tight_binding", "time_reversal_type": "electric", "irrep": ["A1"]},
        "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
        "k_path": "Γ-X",
        "dimension": 8,
        "spinful": True,
        "orbital": ["(px,U)", "(px,D)", "(py,U)", "(py,D)"],
        "ket": ["(px,U)@A_1", "(px,D)@A_1", "(py,U)@A_1", "(py,D)@A_1", "(px,U)@B_1", "(px,D)@B_1", "(py,U)@B_1", "(py,D)@B_1"],
        "ket_site": ["A_1", "B_1"],
        "site": {"A": ("[1/3,2/3,0]", ["px", "py"]), "B": ("[2/3,1/3,0]", ["px", "py"])},
        "rep_site": {
            "A": ("[1/3, 2/3, 0]", "1b", [["(px,U)", "(px,D)", "(py,U)", "(py,D)"]]),
            "B": ("[2/3, 1/3, 0]", "1c", [["(px,U)", "(px,D)", "(py,U)", "(py,D)"]]),
        },
        "cell_site": {"A_1": ("[1/3, 2/3, 0]", "[1,2,3,4,5,6]"), "B_1": ("[2/3, 1/3, 0]", "[1,2,3,4,5,6]")},
        "bond": [("A", "B", 1)],
        "rep_bond": {"A:B:1:1": ("[1/3, -1/3, 0];[2/3, 1/3, 0]", "3d", "D", 1)},
        "cell_bond": {
            "A:B:1:1_1": ("[1/3, 2/3, 0]@[1/2, 0, 0]", "[1,4]"),
            "A:B:1:1_2": ("[-2/3, -1/3, 0]@[0, 1/2, 0]", "[2,6]"),
            "A:B:1:1_3": ("[1/3, -1/3, 0]@[1/2, 1/2, 0]", "[3,5]"),
        },
    },
    "name": {
        "alias": {"S_001": "A", "A": "S_001", "S_002": "B", "B": "S_002", "B_001": "A:B:1:1", "A:B:1:1": "B_001"},
        "site": {"site_001": ("A", 1), "site_002": ("B", 1)},
        "site_name": {"[1/3, 2/3, 0]": ("site_001", 1), "[2/3, 1/3, 0]": ("site_002", 1)},
        "bond": {"bond_001": ("A:B:1:1", 1), "bond_002": ("A:B:1:1", 2), "bond_003": ("A:B:1:1", 3)},
        "bond_name": {
            "[1/3, 2/3, 0]@[1/2, 0, 0]": ("bond_001", 1),
            "[-2/3, -1/3, 0]@[0, 1/2, 0]": ("bond_002", 1),
            "[1/3, -1/3, 0]@[1/2, 1/2, 0]": ("bond_003", 1),
        },
    },
    "data": {
        "plus_set": ["[0, 0, 0]"],
        "cluster_site": {"S_001": ["site_001"], "S_002": ["site_002"]},
        "cluster_bond": {"B_001": ["bond_001", "bond_002", "bond_003"]},
        "site": {
            "site_001": ("[1/3, 2/3, 0]", [0, 1, 2, 3, 4, 5], (0, 0)),
            "site_002": ("[2/3, 1/3, 0]", [0, 1, 2, 3, 4, 5], (1, 1)),
        },
        "bond": {
            "bond_001": ("[1/3, 2/3, 0]@[1/2, 0, 0]", [0, 3], (0, 1), "[1/3, 2/3, 0]", "[1/3, -1/3, 0];[2/3, 1/3, 0]"),
            "bond_002": ("[-2/3, -1/3, 0]@[0, 1/2, 0]", [1, 5], (0, 1), "[-2/3, -1/3, 0]", "[1/3, 2/3, 0];[-1/3, 1/3, 0]"),
            "bond_003": ("[1/3, -1/3, 0]@[1/2, 1/2, 0]", [2, 4], (0, 1), "[1/3, -1/3, 0]", "[1/3, 2/3, 0];[2/3, 1/3, 0]"),
        },
        "cluster_atomic": {(0, 0): [(0, 0, "M_001")], (1, 1): [(4, 4, "M_001")], (0, 1): [(0, 4, "M_001")]},
        "atomic_braket": {"M_001": (["(px,U)", "(px,D)", "(py,U)", "(py,D)"], ["(px,U)", "(px,D)", "(py,U)", "(py,D)"])},
    },
    "detail": {
        "rep_bond_all": {
            "A_B": [
                {},
                {"A:B:1:1": ("[1/3, -1/3, 0];[2/3, 1/3, 0]", "3d", "D", 1)},
                {
                    "A:B:2:1": ("[-2/3, -1/3, 0];[2/3, 1/3, 0]", "1a", "D", 2),
                    "A:B:2:2": ("[1/3, -1/3, 0];[2/3, 1/3, 1]", "3d", "D", 2),
                    "A:B:2:3": ("[1/3, -1/3, 1];[2/3, 1/3, 0]", "3d", "D", 2),
                },
                {
                    "A:B:3:1": ("[-2/3, -1/3, 0];[2/3, 4/3, 0]", "3d", "D", 3),
                    "A:B:3:2": ("[-2/3, -1/3, 0];[2/3, 1/3, 1]", "1a", "D", 3),
                    "A:B:3:3": ("[-2/3, -1/3, 1];[2/3, 1/3, 0]", "1a", "D", 3),
                },
                {
                    "A:B:4:1": ("[-2/3, -1/3, 0];[2/3, 4/3, 1]", "3d", "D", 4),
                    "A:B:4:2": ("[-2/3, -1/3, 1];[2/3, 4/3, 0]", "3d", "D", 4),
                },
                {
                    "A:B:5:1": ("[-2/3, -1/3, 0];[5/3, 1/3, 0]", "3d", "D", 5),
                    "A:B:5:2": ("[1/3, -1/3, 1];[2/3, 1/3, -1]", "3d", "D", 5),
                    "A:B:5:3": ("[1/3, -1/3, -1];[2/3, 1/3, 1]", "3d", "D", 5),
                },
                {
                    "A:B:6:1": ("[-2/3, -1/3, 0];[5/3, 1/3, 1]", "3d", "D", 6),
                    "A:B:6:2": ("[-2/3, -1/3, 1];[5/3, 1/3, 0]", "3d", "D", 6),
                    "A:B:6:3": ("[-2/3, -1/3, -1];[2/3, 1/3, 1]", "1a", "D", 6),
                    "A:B:6:4": ("[-2/3, -4/3, 0];[2/3, 4/3, 0]", "1a", "D", 6),
                    "A:B:6:5": ("[-2/3, -1/3, 1];[2/3, 1/3, -1]", "1a", "D", 6),
                },
                {
                    "A:B:7:1": ("[-2/3, -1/3, 1];[2/3, 4/3, -1]", "3d", "D", 7),
                    "A:B:7:2": ("[-2/3, -4/3, 0];[2/3, 4/3, 1]", "1a", "D", 7),
                    "A:B:7:3": ("[-2/3, -1/3, -1];[2/3, 4/3, 1]", "3d", "D", 7),
                    "A:B:7:4": ("[-2/3, -4/3, 0];[5/3, 4/3, 0]", "3d", "D", 7),
                    "A:B:7:5": ("[-2/3, -4/3, 1];[2/3, 4/3, 0]", "1a", "D", 7),
                },
                {
                    "A:B:8:1": ("[-2/3, -4/3, 1];[5/3, 4/3, 0]", "3d", "D", 8),
                    "A:B:8:2": ("[-2/3, -4/3, 0];[5/3, 4/3, 1]", "3d", "D", 8),
                },
                {
                    "A:B:9:1": ("[-2/3, -1/3, -1];[5/3, 1/3, 1]", "3d", "D", 9),
                    "A:B:9:2": ("[-2/3, -1/3, 1];[5/3, 1/3, -1]", "3d", "D", 9),
                    "A:B:9:3": ("[-5/3, -1/3, 0];[5/3, 4/3, 0]", "3d", "D", 9),
                },
            ]
        },
        "cell_range": (-2, 3, -2, 3, -2, 3),
        "max_neighbor": 10,
        "A": "[[1.0, -0.5, 0.0], [0.0, 0.86602540378444, 0.0], [0.0, 0.0, 1.0]]",
        "version": "1.1.1",
    },
}
