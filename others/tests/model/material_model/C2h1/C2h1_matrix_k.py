"""
== SAMB in full matrix form (* only for crystal) ===
- model : model name.
- molecule : molecule or crystal ?
- group : (tag, detailed str)
- dimension : dimension of full matrix
- ket : ket basis list, orbital@site
- version : MultiPie version
- k_point* : representative k points
- k_path* : high-symmetry line in k space
- A* : transform matrix, [a1,a2,a3]
- bond* : { "bond_#": "vector" }
- matrix : { "z_#": "matrix" }
"""
C2h1 = {
    "model": "C2h1",
    "molecule": False,
    "group": ("C2h^1", "space group No. 10 : C2h^1 / P2/m (b-axis) : PG C2h"),
    "dimension": 4,
    "ket": ["(s,U)@A_1", "(s,D)@A_1", "(s,U)@B_1", "(s,D)@B_1"],
    "version": "1.1.14",
    "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
    "k_path": "Γ-X",
    "A": "[[1.0, 0.0, 0.0], [0.0, 1.2, 0.0], [0.0, 0.0, 1.0]]",
    "bond": {
        "bond_001": "[1, 0, 0]",
        "bond_002": "[0, 0, 1]",
        "bond_003": "[1, 0, 0]",
        "bond_004": "[0, 0, 1]",
        "bond_005": "[1/2, 1/2, 0]",
        "bond_006": "[-1/2, 1/2, 0]",
        "bond_007": "[-1/2, -1/2, 0]",
        "bond_008": "[1/2, -1/2, 0]",
    },
    "matrix": {
        "z_001": "[[sqrt(2)/2, 0, 0, 0], [0, sqrt(2)/2, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]",
        "z_002": "[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, sqrt(2)/2, 0], [0, 0, 0, sqrt(2)/2]]",
        "z_003": "[[c001, 0, 0, 0], [0, c001, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]",
        "z_004": "[[c002, 0, 0, 0], [0, c002, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]",
        "z_005": "[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, c003, 0], [0, 0, 0, c003]]",
        "z_006": "[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, c004, 0], [0, 0, 0, c004]]",
        "z_007": "[[0, 0, c005/2 + c006/2, 0], [0, 0, 0, c005/2 + c006/2], [c005/2 + c006/2, 0, 0, 0], [0, c005/2 + c006/2, 0, 0]]",
        "z_008": "[[0, 0, 0, -c005/2 - c006/2], [0, 0, c005/2 + c006/2, 0], [0, c005/2 + c006/2, 0, 0], [-c005/2 - c006/2, 0, 0, 0]]",
        "z_009": "[[0, 0, 0, -I*c005/2 + I*c006/2], [0, 0, -I*c005/2 + I*c006/2, 0], [0, I*c005/2 - I*c006/2, 0, 0], [I*c005/2 - I*c006/2, 0, 0, 0]]",
        "z_010": "[[0, 0, I*c005/2 - I*c006/2, 0], [0, 0, 0, -I*c005/2 + I*c006/2], [-I*c005/2 + I*c006/2, 0, 0, 0], [0, I*c005/2 - I*c006/2, 0, 0]]",
    },
}
