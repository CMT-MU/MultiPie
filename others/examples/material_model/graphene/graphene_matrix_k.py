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
graphene = {
    "model": "graphene",
    "molecule": False,
    "group": ("D6h^1", "space group No. 191 : D6h^1 / P6/mmm : PG D6h"),
    "dimension": 4,
    "ket": ["(pz,U)@A_1", "(pz,D)@A_1", "(pz,U)@A_2", "(pz,D)@A_2"],
    "version": "1.1.14",
    "k_point": {"Γ": "[0, 0, 0]", "M": "[1/2, 0, 0]", "K": "[1/3, 1/3, 0]", "K'": "[-1/3, -1/3, 0]"},
    "k_path": "Γ-M-K-Γ-K'",
    "A": "[[1.0, -0.5, 0.0], [0.0, 0.86602540378444, 0.0], [0.0, 0.0, 1.0]]",
    "bond": {"bond_001": "[1/3, 2/3, 0]", "bond_002": "[1/3, -1/3, 0]", "bond_003": "[-2/3, -1/3, 0]"},
    "matrix": {
        "z_001": "[[1/2, 0, 0, 0], [0, 1/2, 0, 0], [0, 0, 1/2, 0], [0, 0, 0, 1/2]]",
        "z_002": "[[0, 0, sqrt(3)*c001/6 + sqrt(3)*c002/6 + sqrt(3)*c003/6 + sqrt(3)*I*s001/6 + sqrt(3)*I*s002/6 + sqrt(3)*I*s003/6, 0], [0, 0, 0, sqrt(3)*c001/6 + sqrt(3)*c002/6 + sqrt(3)*c003/6 + sqrt(3)*I*s001/6 + sqrt(3)*I*s002/6 + sqrt(3)*I*s003/6], [sqrt(3)*c001/6 + sqrt(3)*c002/6 + sqrt(3)*c003/6 - sqrt(3)*I*s001/6 - sqrt(3)*I*s002/6 - sqrt(3)*I*s003/6, 0, 0, 0], [0, sqrt(3)*c001/6 + sqrt(3)*c002/6 + sqrt(3)*c003/6 - sqrt(3)*I*s001/6 - sqrt(3)*I*s002/6 - sqrt(3)*I*s003/6, 0, 0]]",
    },
}
