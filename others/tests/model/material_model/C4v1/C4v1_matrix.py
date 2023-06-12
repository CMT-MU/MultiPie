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
    "version": "1.1.7",
    "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
    "k_path": "Γ-X",
    "A": "[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]",
    "bond": {"bond_001": "[0, 0, -1]", "bond_002": "[0, 1, 0]", "bond_003": "[1, 0, 0]"},
    "matrix": {
        "z_001": "[[sqrt(2)/2, 0], [0, sqrt(2)/2]]",
        "z_002": "[[c001, 0], [0, c001]]",
        "z_003": "[[sqrt(2)*c002/2 + sqrt(2)*c003/2, 0], [0, sqrt(2)*c002/2 + sqrt(2)*c003/2]]",
        "z_004": "[[0, sqrt(2)*s002/2 + sqrt(2)*I*s003/2], [sqrt(2)*s002/2 - sqrt(2)*I*s003/2, 0]]",
    },
}
