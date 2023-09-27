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
    "dimension": 2,
    "ket": ["pz@A_1", "pz@A_2"],
    "version": "1.1.14",
    "k_point": {"Γ": "[0, 0, 0]", "M": "[1/2, 0, 0]", "K": "[1/3, 1/3, 0]", "K'": "[-1/3, -1/3, 0]"},
    "k_path": "Γ-M-K-Γ-K'",
    "A": "[[2.435, -1.2175, 0.0], [0.0, 2.10877185821511, 0.0], [0.0, 0.0, 10.0]]",
    "bond": {
        "bond_001": "[1/3, 2/3, 0]",
        "bond_002": "[1/3, -1/3, 0]",
        "bond_003": "[-2/3, -1/3, 0]",
        "bond_004": "[0, 1, 0]",
        "bond_005": "[0, 1, 0]",
        "bond_006": "[1, 1, 0]",
        "bond_007": "[1, 0, 0]",
        "bond_008": "[1, 1, 0]",
        "bond_009": "[1, 0, 0]",
    },
    "matrix": {
        "z_001": "[[sqrt(2)/2, 0], [0, sqrt(2)/2]]",
        "z_002": "[[0, sqrt(6)*c001/6 + sqrt(6)*c002/6 + sqrt(6)*c003/6 + sqrt(6)*I*s001/6 + sqrt(6)*I*s002/6 + sqrt(6)*I*s003/6], [sqrt(6)*c001/6 + sqrt(6)*c002/6 + sqrt(6)*c003/6 - sqrt(6)*I*s001/6 - sqrt(6)*I*s002/6 - sqrt(6)*I*s003/6, 0]]",
        "z_003": "[[sqrt(3)*c004/3 + sqrt(3)*c006/3 + sqrt(3)*c007/3, 0], [0, sqrt(3)*c004/3 + sqrt(3)*c006/3 + sqrt(3)*c007/3]]",
    },
}
