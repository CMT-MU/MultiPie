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
SrVO3 = {
    "model": "SrVO3",
    "molecule": False,
    "group": ("Oh^1", "space group No. 221 : Oh^1 / Pm-3m : PG Oh"),
    "dimension": 3,
    "ket": ["dyz@V_1", "dzx@V_1", "dxy@V_1"],
    "version": "1.1.7",
    "k_point": {"R": "[1/2,1/2,1/2]", "Γ": "[0,0,0]", "X": "[1/2,0,0]", "M": "[1/2,1/2,0]"},
    "k_path": "R-Γ-X-R-M-Γ",
    "A": "[[3.8409, 0.0, 0.0], [0.0, 3.8409, 0.0], [0.0, 0.0, 3.8409]]",
    "bond": {
        "bond_001": "[0, 0, 1]",
        "bond_002": "[1, 0, 0]",
        "bond_003": "[0, 1, 0]",
        "bond_004": "[0, 1, 1]",
        "bond_005": "[0, 1, -1]",
        "bond_006": "[1, 0, -1]",
        "bond_007": "[1, -1, 0]",
        "bond_008": "[1, 0, 1]",
        "bond_009": "[1, 1, 0]",
    },
    "matrix": {
        "z_001": "[[sqrt(3)/3, 0, 0], [0, sqrt(3)/3, 0], [0, 0, sqrt(3)/3]]",
        "z_002": "[[sqrt(2)*c001/3 + sqrt(2)*c002/3 + sqrt(2)*c003/3, 0, 0], [0, sqrt(2)*c001/3 + sqrt(2)*c002/3 + sqrt(2)*c003/3, 0], [0, 0, sqrt(2)*c001/3 + sqrt(2)*c002/3 + sqrt(2)*c003/3]]",
        "z_003": "[[-c001/3 + 2*c002/3 - c003/3, 0, 0], [0, -c001/3 - c002/3 + 2*c003/3, 0], [0, 0, 2*c001/3 - c002/3 - c003/3]]",
        "z_004": "[[c004/3 + c005/3 + c006/3 + c007/3 + c008/3 + c009/3, 0, 0], [0, c004/3 + c005/3 + c006/3 + c007/3 + c008/3 + c009/3, 0], [0, 0, c004/3 + c005/3 + c006/3 + c007/3 + c008/3 + c009/3]]",
        "z_005": "[[-2*sqrt(5)*c004/15 - 2*sqrt(5)*c005/15 + sqrt(5)*c006/15 + sqrt(5)*c007/15 + sqrt(5)*c008/15 + sqrt(5)*c009/15, -sqrt(10)*c007/10 + sqrt(10)*c009/10, -sqrt(10)*c006/10 + sqrt(10)*c008/10], [-sqrt(10)*c007/10 + sqrt(10)*c009/10, sqrt(5)*c004/15 + sqrt(5)*c005/15 - 2*sqrt(5)*c006/15 + sqrt(5)*c007/15 - 2*sqrt(5)*c008/15 + sqrt(5)*c009/15, sqrt(10)*c004/10 - sqrt(10)*c005/10], [-sqrt(10)*c006/10 + sqrt(10)*c008/10, sqrt(10)*c004/10 - sqrt(10)*c005/10, sqrt(5)*c004/15 + sqrt(5)*c005/15 + sqrt(5)*c006/15 - 2*sqrt(5)*c007/15 + sqrt(5)*c008/15 - 2*sqrt(5)*c009/15]]",
        "z_006": "[[-sqrt(30)*c004/15 - sqrt(30)*c005/15 + sqrt(30)*c006/30 + sqrt(30)*c007/30 + sqrt(30)*c008/30 + sqrt(30)*c009/30, sqrt(15)*c007/15 - sqrt(15)*c009/15, sqrt(15)*c006/15 - sqrt(15)*c008/15], [sqrt(15)*c007/15 - sqrt(15)*c009/15, sqrt(30)*c004/30 + sqrt(30)*c005/30 - sqrt(30)*c006/15 + sqrt(30)*c007/30 - sqrt(30)*c008/15 + sqrt(30)*c009/30, -sqrt(15)*c004/15 + sqrt(15)*c005/15], [sqrt(15)*c006/15 - sqrt(15)*c008/15, -sqrt(15)*c004/15 + sqrt(15)*c005/15, sqrt(30)*c004/30 + sqrt(30)*c005/30 + sqrt(30)*c006/30 - sqrt(30)*c007/15 + sqrt(30)*c008/30 - sqrt(30)*c009/15]]",
    },
}
