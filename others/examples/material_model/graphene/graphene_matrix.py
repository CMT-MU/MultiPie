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
graphene = {
    "model": "graphene",
    "molecule": False,
    "group": ("D6h^1", "space group No. 191 : D6h^1 / P6/mmm : PG D6h"),
    "dimension": 4,
    "ket": ["(pz,U)@A_1", "(pz,D)@A_1", "(pz,U)@A_2", "(pz,D)@A_2"],
    "cell_site": {
        "A_1": ("[1/3, 2/3, 0]", "[1,6,7,8,9,10,14,15,16,17,23,24]"),
        "A_2": ("[2/3, 1/3, 0]", "[2,3,4,5,11,12,13,18,19,20,21,22]"),
    },
    "version": "1.1.14",
    "k_point": {"Γ": "[0, 0, 0]", "M": "[1/2, 0, 0]", "K": "[1/3, 1/3, 0]", "K'": "[-1/3, -1/3, 0]"},
    "k_path": "Γ-M-K-Γ-K'",
    "A": "[[1.0, -0.5, 0.0], [0.0, 0.86602540378444, 0.0], [0.0, 0.0, 1.0]]",
    "matrix": {
        "z_001": {(0, 0, 0, 0, 0): "1/2", (0, 0, 0, 1, 1): "1/2", (0, 0, 0, 2, 2): "1/2", (0, 0, 0, 3, 3): "1/2"},
        "z_002": {
            (0, 1, 0, 2, 0): "sqrt(3)/6",
            (0, -1, 0, 0, 2): "sqrt(3)/6",
            (0, 1, 0, 3, 1): "sqrt(3)/6",
            (0, -1, 0, 1, 3): "sqrt(3)/6",
            (0, 0, 0, 2, 0): "sqrt(3)/6",
            (0, 0, 0, 0, 2): "sqrt(3)/6",
            (0, 0, 0, 3, 1): "sqrt(3)/6",
            (0, 0, 0, 1, 3): "sqrt(3)/6",
            (-1, 0, 0, 2, 0): "sqrt(3)/6",
            (1, 0, 0, 0, 2): "sqrt(3)/6",
            (-1, 0, 0, 3, 1): "sqrt(3)/6",
            (1, 0, 0, 1, 3): "sqrt(3)/6",
        },
    },
}
