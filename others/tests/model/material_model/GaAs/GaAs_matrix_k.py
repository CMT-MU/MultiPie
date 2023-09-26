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
GaAs = {
    "model": "GaAs",
    "molecule": False,
    "group": ("Td^2", "space group No. 216 : Td^2 / F-43m : PG Td"),
    "dimension": 6,
    "ket": ["px@Ga_1", "py@Ga_1", "pz@Ga_1", "px@As_1", "py@As_1", "pz@As_1"],
    "version": "1.1.14",
    "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
    "k_path": "Γ-X",
    "A": "[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]",
    "bond": {
        "bond_001": "[1/4, 1/4, 1/4]",
        "bond_002": "[-1/4, -1/4, 1/4]",
        "bond_003": "[1/4, -1/4, -1/4]",
        "bond_004": "[-1/4, 1/4, -1/4]",
    },
    "matrix": {
        "z_001": "[[sqrt(3)/3, 0, 0, 0, 0, 0], [0, sqrt(3)/3, 0, 0, 0, 0], [0, 0, sqrt(3)/3, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]",
        "z_002": "[[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, sqrt(3)/3, 0, 0], [0, 0, 0, 0, sqrt(3)/3, 0], [0, 0, 0, 0, 0, sqrt(3)/3]]",
        "z_003": "[[0, 0, 0, sqrt(6)*c001/12 + sqrt(6)*c002/12 + sqrt(6)*c003/12 + sqrt(6)*c004/12 + sqrt(6)*I*s001/12 + sqrt(6)*I*s002/12 + sqrt(6)*I*s003/12 + sqrt(6)*I*s004/12, 0, 0], [0, 0, 0, 0, sqrt(6)*c001/12 + sqrt(6)*c002/12 + sqrt(6)*c003/12 + sqrt(6)*c004/12 + sqrt(6)*I*s001/12 + sqrt(6)*I*s002/12 + sqrt(6)*I*s003/12 + sqrt(6)*I*s004/12, 0], [0, 0, 0, 0, 0, sqrt(6)*c001/12 + sqrt(6)*c002/12 + sqrt(6)*c003/12 + sqrt(6)*c004/12 + sqrt(6)*I*s001/12 + sqrt(6)*I*s002/12 + sqrt(6)*I*s003/12 + sqrt(6)*I*s004/12], [sqrt(6)*c001/12 + sqrt(6)*c002/12 + sqrt(6)*c003/12 + sqrt(6)*c004/12 - sqrt(6)*I*s001/12 - sqrt(6)*I*s002/12 - sqrt(6)*I*s003/12 - sqrt(6)*I*s004/12, 0, 0, 0, 0, 0], [0, sqrt(6)*c001/12 + sqrt(6)*c002/12 + sqrt(6)*c003/12 + sqrt(6)*c004/12 - sqrt(6)*I*s001/12 - sqrt(6)*I*s002/12 - sqrt(6)*I*s003/12 - sqrt(6)*I*s004/12, 0, 0, 0, 0], [0, 0, sqrt(6)*c001/12 + sqrt(6)*c002/12 + sqrt(6)*c003/12 + sqrt(6)*c004/12 - sqrt(6)*I*s001/12 - sqrt(6)*I*s002/12 - sqrt(6)*I*s003/12 - sqrt(6)*I*s004/12, 0, 0, 0]]",
        "z_004": "[[0, 0, 0, 0, sqrt(3)*c001/12 + sqrt(3)*c002/12 - sqrt(3)*c003/12 - sqrt(3)*c004/12 + sqrt(3)*I*s001/12 + sqrt(3)*I*s002/12 - sqrt(3)*I*s003/12 - sqrt(3)*I*s004/12, sqrt(3)*c001/12 - sqrt(3)*c002/12 - sqrt(3)*c003/12 + sqrt(3)*c004/12 + sqrt(3)*I*s001/12 - sqrt(3)*I*s002/12 - sqrt(3)*I*s003/12 + sqrt(3)*I*s004/12], [0, 0, 0, sqrt(3)*c001/12 + sqrt(3)*c002/12 - sqrt(3)*c003/12 - sqrt(3)*c004/12 + sqrt(3)*I*s001/12 + sqrt(3)*I*s002/12 - sqrt(3)*I*s003/12 - sqrt(3)*I*s004/12, 0, sqrt(3)*c001/12 - sqrt(3)*c002/12 + sqrt(3)*c003/12 - sqrt(3)*c004/12 + sqrt(3)*I*s001/12 - sqrt(3)*I*s002/12 + sqrt(3)*I*s003/12 - sqrt(3)*I*s004/12], [0, 0, 0, sqrt(3)*c001/12 - sqrt(3)*c002/12 - sqrt(3)*c003/12 + sqrt(3)*c004/12 + sqrt(3)*I*s001/12 - sqrt(3)*I*s002/12 - sqrt(3)*I*s003/12 + sqrt(3)*I*s004/12, sqrt(3)*c001/12 - sqrt(3)*c002/12 + sqrt(3)*c003/12 - sqrt(3)*c004/12 + sqrt(3)*I*s001/12 - sqrt(3)*I*s002/12 + sqrt(3)*I*s003/12 - sqrt(3)*I*s004/12, 0], [0, sqrt(3)*c001/12 + sqrt(3)*c002/12 - sqrt(3)*c003/12 - sqrt(3)*c004/12 - sqrt(3)*I*s001/12 - sqrt(3)*I*s002/12 + sqrt(3)*I*s003/12 + sqrt(3)*I*s004/12, sqrt(3)*c001/12 - sqrt(3)*c002/12 - sqrt(3)*c003/12 + sqrt(3)*c004/12 - sqrt(3)*I*s001/12 + sqrt(3)*I*s002/12 + sqrt(3)*I*s003/12 - sqrt(3)*I*s004/12, 0, 0, 0], [sqrt(3)*c001/12 + sqrt(3)*c002/12 - sqrt(3)*c003/12 - sqrt(3)*c004/12 - sqrt(3)*I*s001/12 - sqrt(3)*I*s002/12 + sqrt(3)*I*s003/12 + sqrt(3)*I*s004/12, 0, sqrt(3)*c001/12 - sqrt(3)*c002/12 + sqrt(3)*c003/12 - sqrt(3)*c004/12 - sqrt(3)*I*s001/12 + sqrt(3)*I*s002/12 - sqrt(3)*I*s003/12 + sqrt(3)*I*s004/12, 0, 0, 0], [sqrt(3)*c001/12 - sqrt(3)*c002/12 - sqrt(3)*c003/12 + sqrt(3)*c004/12 - sqrt(3)*I*s001/12 + sqrt(3)*I*s002/12 + sqrt(3)*I*s003/12 - sqrt(3)*I*s004/12, sqrt(3)*c001/12 - sqrt(3)*c002/12 + sqrt(3)*c003/12 - sqrt(3)*c004/12 - sqrt(3)*I*s001/12 + sqrt(3)*I*s002/12 - sqrt(3)*I*s003/12 + sqrt(3)*I*s004/12, 0, 0, 0, 0]]",
    },
}
