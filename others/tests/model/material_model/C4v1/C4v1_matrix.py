"""
== SAMB in full matrix form in real space (* only for crystal) ===
- model : model name.
- molecule : molecule or crystal ?
- group : (tag, detailed str)
- dimension : dimension of full matrix
- ket : ket basis list, orbital@site
- version : MultiPie version
- k_point* : representative k points
- k_path* : high-symmetry line in k space
- cell_site : { name_idx(pset): (position, SOs) }
- A* : transform matrix, [a1,a2,a3]
- matrix : { "z_#": {(n1, n2, n3, a, b): matrix element} }
"""
C4v1 = {
    "model": "C4v1",
    "molecule": False,
    "group": ("C4v^1", "space group No. 99 : C4v^1 / P4mm : PG C4v"),
    "dimension": 2,
    "ket": ["(s,U)@A_1", "(s,D)@A_1"],
    "cell_site": {"A_1": ("[0, 0, 0]", "[1,2,3,4,5,6,7,8]")},
    "version": "1.1.14",
    "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
    "k_path": "Γ-X",
    "A": "[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]",
    "matrix": {
        "z_001": {(0, 0, 0, 0, 0): "sqrt(2)/2", (0, 0, 0, 1, 1): "sqrt(2)/2"},
        "z_002": {
            (0, 1, 0, 0, 0): "sqrt(2)/4",
            (0, -1, 0, 0, 0): "sqrt(2)/4",
            (0, 1, 0, 1, 1): "sqrt(2)/4",
            (0, -1, 0, 1, 1): "sqrt(2)/4",
            (1, 0, 0, 0, 0): "sqrt(2)/4",
            (-1, 0, 0, 0, 0): "sqrt(2)/4",
            (1, 0, 0, 1, 1): "sqrt(2)/4",
            (-1, 0, 0, 1, 1): "sqrt(2)/4",
        },
        "z_003": {
            (0, 1, 0, 0, 1): "sqrt(2)*I/4",
            (0, -1, 0, 1, 0): "-sqrt(2)*I/4",
            (0, 1, 0, 1, 0): "sqrt(2)*I/4",
            (0, -1, 0, 0, 1): "-sqrt(2)*I/4",
            (1, 0, 0, 0, 1): "-sqrt(2)/4",
            (-1, 0, 0, 1, 0): "-sqrt(2)/4",
            (1, 0, 0, 1, 0): "sqrt(2)/4",
            (-1, 0, 0, 0, 1): "sqrt(2)/4",
        },
        "z_004": {(0, 0, 1, 0, 0): "1/2", (0, 0, -1, 0, 0): "1/2", (0, 0, 1, 1, 1): "1/2", (0, 0, -1, 1, 1): "1/2"},
    },
}
