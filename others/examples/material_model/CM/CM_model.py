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
CM = {
    "info": {
        "model": "CM",
        "molecule": True,
        "group": ("Td", "point group No. 31 : Td / -43m"),
        "crystal": "cubic",
        "option": {"view_mode": "debug", "view": None, "output": "CM", "minimal_samb": True},
        "generate": {"model_type": "tight_binding", "time_reversal_type": "electric", "irrep": ["A1"]},
        "dimension": 8,
        "spinful": False,
        "orbital": ["s", "px", "py", "pz"],
        "ket": ["s@C_1", "px@C_1", "py@C_1", "pz@C_1", "s@H_1", "s@H_2", "s@H_3", "s@H_4"],
        "ket_site": ["C_1", "H_1", "H_2", "H_3", "H_4"],
        "site": {"C": ("[0,0,0]", "s p"), "H": ("[1/3,1/3,1/3]", "s")},
        "rep_site": {"C": ("[0, 0, 0]", "1o", [["s"], ["px", "py", "pz"]]), "H": ("[1/3, 1/3, 1/3]", "4a", [["s"]])},
        "cell_site": {
            "C_1": ("[0, 0, 0]", "[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]"),
            "H_1": ("[1/3, 1/3, 1/3]", "[1,5,9,16,17,18]"),
            "H_2": ("[-1/3, -1/3, 1/3]", "[2,6,11,13,21,23]"),
            "H_3": ("[1/3, -1/3, -1/3]", "[3,7,12,15,19,24]"),
            "H_4": ("[-1/3, 1/3, -1/3]", "[4,8,10,14,20,22]"),
        },
        "bond": [("C", "H", 1), ("H", "H", 1)],
        "rep_bond": {
            "C:H:1:1": ("[0, 0, 0];[1/3, 1/3, 1/3]", "4a", "D", 1),
            "H:H:1:1": ("[1/3, 1/3, 1/3];[-1/3, -1/3, 1/3]", "6b", "ND", 1),
        },
        "cell_bond": {
            "C:H:1:1_1": ("[1/3, 1/3, 1/3]@[1/6, 1/6, 1/6]", "[1,5,9,16,17,18]"),
            "C:H:1:1_2": ("[-1/3, -1/3, 1/3]@[-1/6, -1/6, 1/6]", "[2,6,11,13,21,23]"),
            "C:H:1:1_3": ("[1/3, -1/3, -1/3]@[1/6, -1/6, -1/6]", "[3,7,12,15,19,24]"),
            "C:H:1:1_4": ("[-1/3, 1/3, -1/3]@[-1/6, 1/6, -1/6]", "[4,8,10,14,20,22]"),
            "H:H:1:1_1": ("[-2/3, -2/3, 0]@[0, 0, 1/3]", "[1,-2,-13,16]"),
            "H:H:1:1_2": ("[-2/3, 2/3, 0]@[0, 0, -1/3]", "[3,-4,19,-22]"),
            "H:H:1:1_3": ("[0, -2/3, -2/3]@[1/3, 0, 0]", "[5,-7,17,-24]"),
            "H:H:1:1_4": ("[0, 2/3, -2/3]@[-1/3, 0, 0]", "[6,-8,-14,21]"),
            "H:H:1:1_5": ("[-2/3, 0, -2/3]@[0, 1/3, 0]", "[9,-10,18,-20]"),
            "H:H:1:1_6": ("[2/3, 0, -2/3]@[0, -1/3, 0]", "[11,-12,-15,23]"),
        },
    },
    "name": {
        "alias": {
            "S_001": "C",
            "C": "S_001",
            "S_002": "H",
            "H": "S_002",
            "B_001": "C:H:1:1",
            "C:H:1:1": "B_001",
            "B_002": "H:H:1:1",
            "H:H:1:1": "B_002",
        },
        "site": {"site_001": ("C", 1), "site_002": ("H", 1), "site_003": ("H", 2), "site_004": ("H", 3), "site_005": ("H", 4)},
        "site_name": {
            "[0, 0, 0]": ("site_001", 1),
            "[1/3, 1/3, 1/3]": ("site_002", 1),
            "[-1/3, -1/3, 1/3]": ("site_003", 1),
            "[1/3, -1/3, -1/3]": ("site_004", 1),
            "[-1/3, 1/3, -1/3]": ("site_005", 1),
        },
        "bond": {
            "bond_001": ("C:H:1:1", 1),
            "bond_002": ("C:H:1:1", 2),
            "bond_003": ("C:H:1:1", 3),
            "bond_004": ("C:H:1:1", 4),
            "bond_005": ("H:H:1:1", 1),
            "bond_006": ("H:H:1:1", 2),
            "bond_007": ("H:H:1:1", 3),
            "bond_008": ("H:H:1:1", 4),
            "bond_009": ("H:H:1:1", 5),
            "bond_010": ("H:H:1:1", 6),
        },
        "bond_name": {
            "[1/3, 1/3, 1/3]@[1/6, 1/6, 1/6]": ("bond_001", 1),
            "[-1/3, -1/3, 1/3]@[-1/6, -1/6, 1/6]": ("bond_002", 1),
            "[1/3, -1/3, -1/3]@[1/6, -1/6, -1/6]": ("bond_003", 1),
            "[-1/3, 1/3, -1/3]@[-1/6, 1/6, -1/6]": ("bond_004", 1),
            "[-2/3, -2/3, 0]@[0, 0, 1/3]": ("bond_005", 1),
            "[-2/3, 2/3, 0]@[0, 0, -1/3]": ("bond_006", 1),
            "[0, -2/3, -2/3]@[1/3, 0, 0]": ("bond_007", 1),
            "[0, 2/3, -2/3]@[-1/3, 0, 0]": ("bond_008", 1),
            "[-2/3, 0, -2/3]@[0, 1/3, 0]": ("bond_009", 1),
            "[2/3, 0, -2/3]@[0, -1/3, 0]": ("bond_010", 1),
        },
    },
    "data": {
        "cluster_site": {"S_001": ["site_001"], "S_002": ["site_002", "site_003", "site_004", "site_005"]},
        "cluster_bond": {
            "B_001": ["bond_001", "bond_002", "bond_003", "bond_004"],
            "B_002": ["bond_005", "bond_006", "bond_007", "bond_008", "bond_009", "bond_010"],
        },
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
                (0, 1),
                "[1/3, 1/3, 1/3]",
                "[0, 0, 0];[1/3, 1/3, 1/3]",
            ),
            "bond_002": (
                "[-1/3, -1/3, 1/3]@[-1/6, -1/6, 1/6]",
                [1, 5, 10, 12, 20, 22],
                (0, 2),
                "[-1/3, -1/3, 1/3]",
                "[0, 0, 0];[-1/3, -1/3, 1/3]",
            ),
            "bond_003": (
                "[1/3, -1/3, -1/3]@[1/6, -1/6, -1/6]",
                [2, 6, 11, 14, 18, 23],
                (0, 3),
                "[1/3, -1/3, -1/3]",
                "[0, 0, 0];[1/3, -1/3, -1/3]",
            ),
            "bond_004": (
                "[-1/3, 1/3, -1/3]@[-1/6, 1/6, -1/6]",
                [3, 7, 9, 13, 19, 21],
                (0, 4),
                "[-1/3, 1/3, -1/3]",
                "[0, 0, 0];[-1/3, 1/3, -1/3]",
            ),
            "bond_005": (
                "[-2/3, -2/3, 0]@[0, 0, 1/3]",
                [0, -1, -12, 15],
                (1, 2),
                "[-2/3, -2/3, 0]",
                "[1/3, 1/3, 1/3];[-1/3, -1/3, 1/3]",
            ),
            "bond_006": (
                "[-2/3, 2/3, 0]@[0, 0, -1/3]",
                [2, -3, 18, -21],
                (3, 4),
                "[-2/3, 2/3, 0]",
                "[1/3, -1/3, -1/3];[-1/3, 1/3, -1/3]",
            ),
            "bond_007": (
                "[0, -2/3, -2/3]@[1/3, 0, 0]",
                [4, -6, 16, -23],
                (1, 3),
                "[0, -2/3, -2/3]",
                "[1/3, 1/3, 1/3];[1/3, -1/3, -1/3]",
            ),
            "bond_008": (
                "[0, 2/3, -2/3]@[-1/3, 0, 0]",
                [5, -7, -13, 20],
                (2, 4),
                "[0, 2/3, -2/3]",
                "[-1/3, -1/3, 1/3];[-1/3, 1/3, -1/3]",
            ),
            "bond_009": (
                "[-2/3, 0, -2/3]@[0, 1/3, 0]",
                [8, -9, 17, -19],
                (1, 4),
                "[-2/3, 0, -2/3]",
                "[1/3, 1/3, 1/3];[-1/3, 1/3, -1/3]",
            ),
            "bond_010": (
                "[2/3, 0, -2/3]@[0, -1/3, 0]",
                [10, -11, -14, 22],
                (2, 3),
                "[2/3, 0, -2/3]",
                "[-1/3, -1/3, 1/3];[1/3, -1/3, -1/3]",
            ),
        },
        "cluster_atomic": {
            (0, 0): [(0, 0, "M_001"), (0, 1, "M_002"), (1, 1, "M_003")],
            (1, 1): [(4, 4, "M_001")],
            (2, 2): [(5, 5, "M_001")],
            (3, 3): [(6, 6, "M_001")],
            (4, 4): [(7, 7, "M_001")],
            (0, 1): [(0, 4, "M_001"), (1, 4, "M_004")],
            (0, 2): [(0, 5, "M_001"), (1, 5, "M_004")],
            (0, 3): [(0, 6, "M_001"), (1, 6, "M_004")],
            (0, 4): [(0, 7, "M_001"), (1, 7, "M_004")],
            (1, 2): [(4, 5, "M_001")],
            (3, 4): [(6, 7, "M_001")],
            (1, 3): [(4, 6, "M_001")],
            (2, 4): [(5, 7, "M_001")],
            (1, 4): [(4, 7, "M_001")],
            (2, 3): [(5, 6, "M_001")],
        },
        "atomic_braket": {
            "M_001": (["s"], ["s"]),
            "M_002": (["s"], ["px", "py", "pz"]),
            "M_003": (["px", "py", "pz"], ["px", "py", "pz"]),
            "M_004": (["px", "py", "pz"], ["s"]),
        },
    },
    "detail": {
        "rep_bond_all": {
            "C_H": [{}, {"C:H:1:1": ("[0, 0, 0];[1/3, 1/3, 1/3]", "4a", "D", 1)}],
            "H_H": [{}, {"H:H:1:1": ("[1/3, 1/3, 1/3];[-1/3, -1/3, 1/3]", "6b", "ND", 1)}],
        },
        "max_neighbor": 10,
        "version": "1.1.1",
    },
}
