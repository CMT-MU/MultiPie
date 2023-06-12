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
CH4 = {
    "info": {
        "model": "CH4",
        "molecule": True,
        "group": ("Td", "point group No. 31 : Td / -43m"),
        "crystal": "cubic",
        "option": {"view": None, "view_mode": "standard", "output": "CH4", "minimal_samb": True},
        "generate": {"irrep": ["A1", "A2"], "model_type": "tight_binding", "time_reversal_type": "electric"},
        "dimension": 16,
        "spinful": True,
        "orbital": ["(s,U)", "(s,D)", "(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"],
        "ket": [
            "(s,U)@C_1",
            "(s,D)@C_1",
            "(px,U)@C_1",
            "(px,D)@C_1",
            "(py,U)@C_1",
            "(py,D)@C_1",
            "(pz,U)@C_1",
            "(pz,D)@C_1",
            "(s,U)@H_1",
            "(s,D)@H_1",
            "(s,U)@H_2",
            "(s,D)@H_2",
            "(s,U)@H_3",
            "(s,D)@H_3",
            "(s,U)@H_4",
            "(s,D)@H_4",
        ],
        "ket_site": ["C_1", "H_1", "H_2", "H_3", "H_4"],
        "site": {"C": ("[0,0,0]", "s p"), "H": ("[1/3,1/3,1/3]", "s")},
        "rep_site": {
            "C": ("[0, 0, 0]", "1o", [["(s,U)", "(s,D)"], ["(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"]]),
            "H": ("[1/3, 1/3, 1/3]", "4a", [["(s,U)", "(s,D)"]]),
        },
        "cell_site": {
            "C_1": ("[0, 0, 0]", "[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]"),
            "H_1": ("[1/3, 1/3, 1/3]", "[1,5,9,16,17,18]"),
            "H_2": ("[-1/3, -1/3, 1/3]", "[2,6,11,13,21,23]"),
            "H_3": ("[1/3, -1/3, -1/3]", "[3,7,12,15,19,24]"),
            "H_4": ("[-1/3, 1/3, -1/3]", "[4,8,10,14,20,22]"),
        },
        "bond": [("C", "H", 1)],
        "rep_bond": {"C:H:1:1": ("[0, 0, 0];[1/3, 1/3, 1/3]", "4a", "D", 1)},
        "cell_bond": {
            "C:H:1:1_1": ("[1/3, 1/3, 1/3]@[1/6, 1/6, 1/6]", "[1,5,9,16,17,18]"),
            "C:H:1:1_2": ("[-1/3, -1/3, 1/3]@[-1/6, -1/6, 1/6]", "[2,6,11,13,21,23]"),
            "C:H:1:1_3": ("[1/3, -1/3, -1/3]@[1/6, -1/6, -1/6]", "[3,7,12,15,19,24]"),
            "C:H:1:1_4": ("[-1/3, 1/3, -1/3]@[-1/6, 1/6, -1/6]", "[4,8,10,14,20,22]"),
        },
    },
    "name": {
        "alias": {"S_001": "C", "C": "S_001", "S_002": "H", "H": "S_002", "B_001": "C:H:1:1", "C:H:1:1": "B_001"},
        "site": {"site_001": ("C", 1), "site_002": ("H", 1), "site_003": ("H", 2), "site_004": ("H", 3), "site_005": ("H", 4)},
        "site_name": {
            "[0, 0, 0]": ("site_001", 1),
            "[1/3, 1/3, 1/3]": ("site_002", 1),
            "[-1/3, -1/3, 1/3]": ("site_003", 1),
            "[1/3, -1/3, -1/3]": ("site_004", 1),
            "[-1/3, 1/3, -1/3]": ("site_005", 1),
        },
        "bond": {"bond_001": ("C:H:1:1", 1), "bond_002": ("C:H:1:1", 2), "bond_003": ("C:H:1:1", 3), "bond_004": ("C:H:1:1", 4)},
        "bond_name": {
            "[1/3, 1/3, 1/3]@[1/6, 1/6, 1/6]": ("bond_001", 1),
            "[-1/3, -1/3, 1/3]@[-1/6, -1/6, 1/6]": ("bond_002", 1),
            "[1/3, -1/3, -1/3]@[1/6, -1/6, -1/6]": ("bond_003", 1),
            "[-1/3, 1/3, -1/3]@[-1/6, 1/6, -1/6]": ("bond_004", 1),
        },
    },
    "data": {
        "cluster_site": {"S_001": ["site_001"], "S_002": ["site_002", "site_003", "site_004", "site_005"]},
        "cluster_bond": {"B_001": ["bond_001", "bond_002", "bond_003", "bond_004"]},
        "site": {
            "site_001": (
                "[0, 0, 0]",
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23],
                (0, 0),
            ),
            "site_002": ("[1/3, 1/3, 1/3]", [0, 4, 8, 15, 16, 17], (1, 1)),
            "site_003": ("[-1/3, -1/3, 1/3]", [1, 5, 10, 12, 20, 22], (2, 2)),
            "site_004": ("[1/3, -1/3, -1/3]", [2, 6, 11, 14, 18, 23], (3, 3)),
            "site_005": ("[-1/3, 1/3, -1/3]", [3, 7, 9, 13, 19, 21], (4, 4)),
        },
        "bond": {
            "bond_001": (
                "[1/3, 1/3, 1/3]@[1/6, 1/6, 1/6]",
                [0, 4, 8, 15, 16, 17],
                (1, 0),
                "[1/3, 1/3, 1/3]",
                "[0, 0, 0];[1/3, 1/3, 1/3]",
            ),
            "bond_002": (
                "[-1/3, -1/3, 1/3]@[-1/6, -1/6, 1/6]",
                [1, 5, 10, 12, 20, 22],
                (2, 0),
                "[-1/3, -1/3, 1/3]",
                "[0, 0, 0];[-1/3, -1/3, 1/3]",
            ),
            "bond_003": (
                "[1/3, -1/3, -1/3]@[1/6, -1/6, -1/6]",
                [2, 6, 11, 14, 18, 23],
                (3, 0),
                "[1/3, -1/3, -1/3]",
                "[0, 0, 0];[1/3, -1/3, -1/3]",
            ),
            "bond_004": (
                "[-1/3, 1/3, -1/3]@[-1/6, 1/6, -1/6]",
                [3, 7, 9, 13, 19, 21],
                (4, 0),
                "[-1/3, 1/3, -1/3]",
                "[0, 0, 0];[-1/3, 1/3, -1/3]",
            ),
        },
        "cluster_atomic": {
            (0, 0): [(0, 0, "M_001"), (0, 2, "M_002"), (2, 2, "M_003")],
            (1, 1): [(8, 8, "M_001")],
            (2, 2): [(10, 10, "M_001")],
            (3, 3): [(12, 12, "M_001")],
            (4, 4): [(14, 14, "M_001")],
            (1, 0): [(8, 0, "M_001"), (8, 2, "M_002")],
            (2, 0): [(10, 0, "M_001"), (10, 2, "M_002")],
            (3, 0): [(12, 0, "M_001"), (12, 2, "M_002")],
            (4, 0): [(14, 0, "M_001"), (14, 2, "M_002")],
        },
        "atomic_braket": {
            "M_001": (["(s,U)", "(s,D)"], ["(s,U)", "(s,D)"]),
            "M_002": (["(s,U)", "(s,D)"], ["(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"]),
            "M_003": (
                ["(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"],
                ["(px,U)", "(px,D)", "(py,U)", "(py,D)", "(pz,U)", "(pz,D)"],
            ),
        },
    },
    "detail": {
        "rep_bond_all": {"C_H": [{}, {"C:H:1:1": ("[0, 0, 0];[1/3, 1/3, 1/3]", "4a", "D", 1)}]},
        "max_neighbor": 10,
        "version": "1.1.7",
    },
}
