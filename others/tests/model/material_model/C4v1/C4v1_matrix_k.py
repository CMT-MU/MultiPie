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
C4v1 = {
    "model": "C4v1",
    "molecule": False,
    "group": ("C4v^1", "space group No. 99 : C4v^1 / P4mm : PG C4v"),
    "dimension": 2,
    "ket": ["(s,U)@A_1", "(s,D)@A_1"],
    "version": "1.1.14",
    "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
    "k_path": "Γ-X",
    "A": "[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]",
    "bond": {"bond_001": "[0, 1, 0]", "bond_002": "[1, 0, 0]", "bond_003": "[0, 0, -1]"},
    "matrix": {
        "z_001": "[[sqrt(2)/2, 0], [0, sqrt(2)/2]]",
        "z_002": "[[sqrt(2)*c001/2 + sqrt(2)*c002/2, 0], [0, sqrt(2)*c001/2 + sqrt(2)*c002/2]]",
        "z_003": "[[0, sqrt(2)*s001/2 + sqrt(2)*I*s002/2], [sqrt(2)*s001/2 - sqrt(2)*I*s002/2, 0]]",
        "z_004": "[[c003, 0], [0, c003]]",
    },
}
