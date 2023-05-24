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
grapheneAB = {
    "model": "grapheneAB",
    "molecule": False,
    "group": ("D3h^1", "space group No. 187 : D3h^1 / P-6m2 : PG D3h"),
    "dimension": 3,
    "ket": ["s@A_1", "px@B_1", "py@B_1"],
    "version": "1.1.1",
    "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
    "k_path": "Γ-X",
    "A": "[[2.435, -1.2175, 0.0], [0.0, 2.10877185821511, 0.0], [0.0, 0.0, 10.0]]",
    "bond": {
        "bond_001": "[1/3, 2/3, 0]",
        "bond_002": "[1/3, -1/3, 0]",
        "bond_003": "[-2/3, -1/3, 0]",
        "bond_004": "[0, 1, 0]",
        "bond_005": "[1, 1, 0]",
        "bond_006": "[1, 0, 0]",
        "bond_007": "[1, 0, 0]",
        "bond_008": "[1, 1, 0]",
        "bond_009": "[0, 1, 0]",
    },
    "matrix": {
        "z_001": "[[1, 0, 0], [0, 0, 0], [0, 0, 0]]",
        "z_002": "[[0, 0, 0], [0, sqrt(2)/2, 0], [0, 0, sqrt(2)/2]]",
        "z_003": "[[0, -sqrt(2)*c002/4 + sqrt(2)*c003/4 + sqrt(2)*I*s002/4 - sqrt(2)*I*s003/4, -sqrt(6)*c001/6 + sqrt(6)*c002/12 + sqrt(6)*c003/12 + sqrt(6)*I*s001/6 - sqrt(6)*I*s002/12 - sqrt(6)*I*s003/12], [-sqrt(2)*c002/4 + sqrt(2)*c003/4 - sqrt(2)*I*s002/4 + sqrt(2)*I*s003/4, 0, 0], [-sqrt(6)*c001/6 + sqrt(6)*c002/12 + sqrt(6)*c003/12 - sqrt(6)*I*s001/6 + sqrt(6)*I*s002/12 + sqrt(6)*I*s003/12, 0, 0]]",
        "z_004": "[[sqrt(6)*c004/3 + sqrt(6)*c005/3 + sqrt(6)*c006/3, 0, 0], [0, 0, 0], [0, 0, 0]]",
        "z_005": "[[0, 0, 0], [0, sqrt(3)*c007/3 + sqrt(3)*c008/3 + sqrt(3)*c009/3, 0], [0, 0, sqrt(3)*c007/3 + sqrt(3)*c008/3 + sqrt(3)*c009/3]]",
        "z_006": "[[0, 0, 0], [0, -sqrt(3)*c007/3 + sqrt(3)*c008/6 + sqrt(3)*c009/6, -c008/2 + c009/2], [0, -c008/2 + c009/2, sqrt(3)*c007/3 - sqrt(3)*c008/6 - sqrt(3)*c009/6]]",
        "z_007": "[[0, 0, 0], [0, 0, -sqrt(3)*I*s007/3 + sqrt(3)*I*s008/3 - sqrt(3)*I*s009/3], [0, sqrt(3)*I*s007/3 - sqrt(3)*I*s008/3 + sqrt(3)*I*s009/3, 0]]",
    },
}
