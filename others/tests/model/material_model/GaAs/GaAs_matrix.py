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
GaAs = {
    "model": "GaAs",
    "molecule": False,
    "group": ("Td^2", "space group No. 216 : Td^2 / F-43m : PG Td"),
    "dimension": 6,
    "ket": ["px@Ga_1", "py@Ga_1", "pz@Ga_1", "px@As_1", "py@As_1", "pz@As_1"],
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
    "version": "1.1.14",
    "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
    "k_path": "Γ-X",
    "A": "[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]",
    "matrix": {
        "z_001": {(0, 0, 0, 0, 0): "sqrt(3)/3", (0, 0, 0, 1, 1): "sqrt(3)/3", (0, 0, 0, 2, 2): "sqrt(3)/3"},
        "z_002": {(0, 0, 0, 3, 3): "sqrt(3)/3", (0, 0, 0, 4, 4): "sqrt(3)/3", (0, 0, 0, 5, 5): "sqrt(3)/3"},
        "z_003": {
            (0, 0, 0, 3, 0): "sqrt(6)/3",
            (0, 0, 0, 0, 3): "sqrt(6)/3",
            (0, 0, 0, 4, 1): "sqrt(6)/3",
            (0, 0, 0, 1, 4): "sqrt(6)/3",
            (0, 0, 0, 5, 2): "sqrt(6)/3",
            (0, 0, 0, 2, 5): "sqrt(6)/3",
        },
        "z_004": {
            (0, 0, 0, 4, 2): "0",
            (0, 0, 0, 2, 4): "0",
            (0, 0, 0, 5, 1): "0",
            (0, 0, 0, 1, 5): "0",
            (0, 0, 0, 3, 2): "0",
            (0, 0, 0, 2, 3): "0",
            (0, 0, 0, 5, 0): "0",
            (0, 0, 0, 0, 5): "0",
            (0, 0, 0, 3, 1): "0",
            (0, 0, 0, 1, 3): "0",
            (0, 0, 0, 4, 0): "0",
            (0, 0, 0, 0, 4): "0",
        },
    },
}