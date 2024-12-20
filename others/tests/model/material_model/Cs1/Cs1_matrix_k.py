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
Cs1 = {
    "model": "Cs1",
    "molecule": False,
    "group": ("Cs^1", "space group No. 6 : Cs^1 / Pm (b-axis) : PG Cs"),
    "dimension": 4,
    "ket": ["(px,U)@A_1", "(px,D)@A_1", "(py,U)@A_1", "(py,D)@A_1"],
    "version": "1.1.14",
    "k_point": {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"},
    "k_path": "Γ-X",
    "A": "[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]",
    "bond": {
        "bond_001": "[-1, 0, 0]",
        "bond_002": "[0, 0, -1]",
        "bond_003": "[0, 1, 0]",
        "bond_004": "[-1, 1, 0]",
        "bond_005": "[-1, -1, 0]",
        "bond_006": "[1, 0, 1]",
        "bond_007": "[0, 1, -1]",
        "bond_008": "[0, -1, -1]",
        "bond_009": "[1, 0, -1]",
    },
    "matrix": {
        "z_001": "[[1/2, 0, 0, 0], [0, 1/2, 0, 0], [0, 0, 1/2, 0], [0, 0, 0, 1/2]]",
        "z_002": "[[1/2, 0, 0, 0], [0, 1/2, 0, 0], [0, 0, -1/2, 0], [0, 0, 0, -1/2]]",
        "z_003": "[[0, 0, -I/2, 0], [0, 0, 0, I/2], [I/2, 0, 0, 0], [0, -I/2, 0, 0]]",
        "z_004": "[[0, 0, 0, -I/2], [0, 0, -I/2, 0], [0, I/2, 0, 0], [I/2, 0, 0, 0]]",
        "z_005": "[[sqrt(2)*c001/2, 0, 0, 0], [0, sqrt(2)*c001/2, 0, 0], [0, 0, sqrt(2)*c001/2, 0], [0, 0, 0, sqrt(2)*c001/2]]",
        "z_006": "[[sqrt(2)*c001/2, 0, 0, 0], [0, sqrt(2)*c001/2, 0, 0], [0, 0, -sqrt(2)*c001/2, 0], [0, 0, 0, -sqrt(2)*c001/2]]",
        "z_007": "[[0, 0, -sqrt(2)*I*c001/2, 0], [0, 0, 0, sqrt(2)*I*c001/2], [sqrt(2)*I*c001/2, 0, 0, 0], [0, -sqrt(2)*I*c001/2, 0, 0]]",
        "z_008": "[[0, 0, 0, -sqrt(2)*I*c001/2], [0, 0, -sqrt(2)*I*c001/2, 0], [0, sqrt(2)*I*c001/2, 0, 0], [sqrt(2)*I*c001/2, 0, 0, 0]]",
        "z_009": "[[0, sqrt(38)*I*s001/19, 0, 3*sqrt(38)*s001/38], [-sqrt(38)*I*s001/19, 0, 3*sqrt(38)*s001/38, 0], [0, 3*sqrt(38)*s001/38, 0, -2*sqrt(38)*I*s001/19], [3*sqrt(38)*s001/38, 0, 2*sqrt(38)*I*s001/19, 0]]",
        "z_010": "[[0, -7*sqrt(19)*I*s001/38, 0, -sqrt(19)*s001/38], [7*sqrt(19)*I*s001/38, 0, -sqrt(19)*s001/38, 0], [0, -sqrt(19)*s001/38, 0, -5*sqrt(19)*I*s001/38], [-sqrt(19)*s001/38, 0, 5*sqrt(19)*I*s001/38, 0]]",
        "z_011": "[[0, 0, sqrt(2)*s001/2, 0], [0, 0, 0, -sqrt(2)*s001/2], [sqrt(2)*s001/2, 0, 0, 0], [0, -sqrt(2)*s001/2, 0, 0]]",
        "z_012": "[[0, I*s001/2, 0, -s001/2], [-I*s001/2, 0, -s001/2, 0], [0, -s001/2, 0, -I*s001/2], [-s001/2, 0, I*s001/2, 0]]",
        "z_013": "[[sqrt(2)*c002/2, 0, 0, 0], [0, sqrt(2)*c002/2, 0, 0], [0, 0, sqrt(2)*c002/2, 0], [0, 0, 0, sqrt(2)*c002/2]]",
        "z_014": "[[sqrt(2)*c002/2, 0, 0, 0], [0, sqrt(2)*c002/2, 0, 0], [0, 0, -sqrt(2)*c002/2, 0], [0, 0, 0, -sqrt(2)*c002/2]]",
        "z_015": "[[0, 0, -sqrt(2)*I*c002/2, 0], [0, 0, 0, sqrt(2)*I*c002/2], [sqrt(2)*I*c002/2, 0, 0, 0], [0, -sqrt(2)*I*c002/2, 0, 0]]",
        "z_016": "[[0, 0, 0, -sqrt(2)*I*c002/2], [0, 0, -sqrt(2)*I*c002/2, 0], [0, sqrt(2)*I*c002/2, 0, 0], [sqrt(2)*I*c002/2, 0, 0, 0]]",
        "z_017": "[[0, sqrt(38)*I*s002/19, 0, 3*sqrt(38)*s002/38], [-sqrt(38)*I*s002/19, 0, 3*sqrt(38)*s002/38, 0], [0, 3*sqrt(38)*s002/38, 0, -2*sqrt(38)*I*s002/19], [3*sqrt(38)*s002/38, 0, 2*sqrt(38)*I*s002/19, 0]]",
        "z_018": "[[0, -7*sqrt(19)*I*s002/38, 0, -sqrt(19)*s002/38], [7*sqrt(19)*I*s002/38, 0, -sqrt(19)*s002/38, 0], [0, -sqrt(19)*s002/38, 0, -5*sqrt(19)*I*s002/38], [-sqrt(19)*s002/38, 0, 5*sqrt(19)*I*s002/38, 0]]",
        "z_019": "[[0, 0, sqrt(2)*s002/2, 0], [0, 0, 0, -sqrt(2)*s002/2], [sqrt(2)*s002/2, 0, 0, 0], [0, -sqrt(2)*s002/2, 0, 0]]",
        "z_020": "[[0, I*s002/2, 0, -s002/2], [-I*s002/2, 0, -s002/2, 0], [0, -s002/2, 0, -I*s002/2], [-s002/2, 0, I*s002/2, 0]]",
        "z_021": "[[sqrt(2)*c003/2, 0, 0, 0], [0, sqrt(2)*c003/2, 0, 0], [0, 0, sqrt(2)*c003/2, 0], [0, 0, 0, sqrt(2)*c003/2]]",
        "z_022": "[[sqrt(2)*c003/2, 0, 0, 0], [0, sqrt(2)*c003/2, 0, 0], [0, 0, -sqrt(2)*c003/2, 0], [0, 0, 0, -sqrt(2)*c003/2]]",
        "z_023": "[[0, 0, -sqrt(2)*I*c003/2, 0], [0, 0, 0, sqrt(2)*I*c003/2], [sqrt(2)*I*c003/2, 0, 0, 0], [0, -sqrt(2)*I*c003/2, 0, 0]]",
        "z_024": "[[0, 0, 0, -sqrt(2)*I*c003/2], [0, 0, -sqrt(2)*I*c003/2, 0], [0, sqrt(2)*I*c003/2, 0, 0], [sqrt(2)*I*c003/2, 0, 0, 0]]",
        "z_025": "[[sqrt(2)*s003/2, 0, 0, 0], [0, -sqrt(2)*s003/2, 0, 0], [0, 0, sqrt(2)*s003/2, 0], [0, 0, 0, -sqrt(2)*s003/2]]",
        "z_026": "[[0, 2*sqrt(38)*s003/19, 0, -3*sqrt(38)*I*s003/38], [2*sqrt(38)*s003/19, 0, 3*sqrt(38)*I*s003/38, 0], [0, -3*sqrt(38)*I*s003/38, 0, -sqrt(38)*s003/19], [3*sqrt(38)*I*s003/38, 0, -sqrt(38)*s003/19, 0]]",
        "z_027": "[[0, 5*sqrt(19)*s003/38, 0, sqrt(19)*I*s003/38], [5*sqrt(19)*s003/38, 0, -sqrt(19)*I*s003/38, 0], [0, sqrt(19)*I*s003/38, 0, 7*sqrt(19)*s003/38], [-sqrt(19)*I*s003/38, 0, 7*sqrt(19)*s003/38, 0]]",
        "z_028": "[[-sqrt(2)*s003/2, 0, 0, 0], [0, sqrt(2)*s003/2, 0, 0], [0, 0, sqrt(2)*s003/2, 0], [0, 0, 0, -sqrt(2)*s003/2]]",
        "z_029": "[[0, -s003/2, 0, -I*s003/2], [-s003/2, 0, I*s003/2, 0], [0, -I*s003/2, 0, s003/2], [I*s003/2, 0, s003/2, 0]]",
        "z_030": "[[0, 0, sqrt(2)*I*s003/2, 0], [0, 0, 0, sqrt(2)*I*s003/2], [-sqrt(2)*I*s003/2, 0, 0, 0], [0, -sqrt(2)*I*s003/2, 0, 0]]",
        "z_031": "[[c004/2 + c005/2, 0, 0, 0], [0, c004/2 + c005/2, 0, 0], [0, 0, c004/2 + c005/2, 0], [0, 0, 0, c004/2 + c005/2]]",
        "z_032": "[[c004/2 + c005/2, 0, 0, 0], [0, c004/2 + c005/2, 0, 0], [0, 0, -c004/2 - c005/2, 0], [0, 0, 0, -c004/2 - c005/2]]",
        "z_033": "[[0, 0, c004/2 - c005/2, 0], [0, 0, 0, c004/2 - c005/2], [c004/2 - c005/2, 0, 0, 0], [0, c004/2 - c005/2, 0, 0]]",
        "z_034": "[[0, 0, -I*c004/2 - I*c005/2, 0], [0, 0, 0, I*c004/2 + I*c005/2], [I*c004/2 + I*c005/2, 0, 0, 0], [0, -I*c004/2 - I*c005/2, 0, 0]]",
        "z_035": "[[0, 0, 0, -I*c004/2 - I*c005/2], [0, 0, -I*c004/2 - I*c005/2, 0], [0, I*c004/2 + I*c005/2, 0, 0], [I*c004/2 + I*c005/2, 0, 0, 0]]",
        "z_036": "[[0, 0, 0, -c004/2 + c005/2], [0, 0, c004/2 - c005/2, 0], [0, c004/2 - c005/2, 0, 0], [-c004/2 + c005/2, 0, 0, 0]]",
        "z_037": "[[0, sqrt(19)*I*s004/19 + sqrt(19)*I*s005/19, 0, 3*sqrt(19)*s004/38 + 3*sqrt(19)*s005/38], [-sqrt(19)*I*s004/19 - sqrt(19)*I*s005/19, 0, 3*sqrt(19)*s004/38 + 3*sqrt(19)*s005/38, 0], [0, 3*sqrt(19)*s004/38 + 3*sqrt(19)*s005/38, 0, -2*sqrt(19)*I*s004/19 - 2*sqrt(19)*I*s005/19], [3*sqrt(19)*s004/38 + 3*sqrt(19)*s005/38, 0, 2*sqrt(19)*I*s004/19 + 2*sqrt(19)*I*s005/19, 0]]",
        "z_038": "[[s004/2 - s005/2, 0, 0, 0], [0, -s004/2 + s005/2, 0, 0], [0, 0, s004/2 - s005/2, 0], [0, 0, 0, -s004/2 + s005/2]]",
        "z_039": "[[0, 2*sqrt(19)*s004/19 - 2*sqrt(19)*s005/19, 0, -3*sqrt(19)*I*s004/38 + 3*sqrt(19)*I*s005/38], [2*sqrt(19)*s004/19 - 2*sqrt(19)*s005/19, 0, 3*sqrt(19)*I*s004/38 - 3*sqrt(19)*I*s005/38, 0], [0, -3*sqrt(19)*I*s004/38 + 3*sqrt(19)*I*s005/38, 0, -sqrt(19)*s004/19 + sqrt(19)*s005/19], [3*sqrt(19)*I*s004/38 - 3*sqrt(19)*I*s005/38, 0, -sqrt(19)*s004/19 + sqrt(19)*s005/19, 0]]",
        "z_040": "[[0, -7*sqrt(38)*I*s004/76 - 7*sqrt(38)*I*s005/76, 0, -sqrt(38)*s004/76 - sqrt(38)*s005/76], [7*sqrt(38)*I*s004/76 + 7*sqrt(38)*I*s005/76, 0, -sqrt(38)*s004/76 - sqrt(38)*s005/76, 0], [0, -sqrt(38)*s004/76 - sqrt(38)*s005/76, 0, -5*sqrt(38)*I*s004/76 - 5*sqrt(38)*I*s005/76], [-sqrt(38)*s004/76 - sqrt(38)*s005/76, 0, 5*sqrt(38)*I*s004/76 + 5*sqrt(38)*I*s005/76, 0]]",
        "z_041": "[[0, 5*sqrt(38)*s004/76 - 5*sqrt(38)*s005/76, 0, sqrt(38)*I*s004/76 - sqrt(38)*I*s005/76], [5*sqrt(38)*s004/76 - 5*sqrt(38)*s005/76, 0, -sqrt(38)*I*s004/76 + sqrt(38)*I*s005/76, 0], [0, sqrt(38)*I*s004/76 - sqrt(38)*I*s005/76, 0, 7*sqrt(38)*s004/76 - 7*sqrt(38)*s005/76], [-sqrt(38)*I*s004/76 + sqrt(38)*I*s005/76, 0, 7*sqrt(38)*s004/76 - 7*sqrt(38)*s005/76, 0]]",
        "z_042": "[[0, 0, s004/2 + s005/2, 0], [0, 0, 0, -s004/2 - s005/2], [s004/2 + s005/2, 0, 0, 0], [0, -s004/2 - s005/2, 0, 0]]",
        "z_043": "[[0, sqrt(2)*I*s004/4 + sqrt(2)*I*s005/4, 0, -sqrt(2)*s004/4 - sqrt(2)*s005/4], [-sqrt(2)*I*s004/4 - sqrt(2)*I*s005/4, 0, -sqrt(2)*s004/4 - sqrt(2)*s005/4, 0], [0, -sqrt(2)*s004/4 - sqrt(2)*s005/4, 0, -sqrt(2)*I*s004/4 - sqrt(2)*I*s005/4], [-sqrt(2)*s004/4 - sqrt(2)*s005/4, 0, sqrt(2)*I*s004/4 + sqrt(2)*I*s005/4, 0]]",
        "z_044": "[[-s004/2 + s005/2, 0, 0, 0], [0, s004/2 - s005/2, 0, 0], [0, 0, s004/2 - s005/2, 0], [0, 0, 0, -s004/2 + s005/2]]",
        "z_045": "[[0, -sqrt(2)*s004/4 + sqrt(2)*s005/4, 0, -sqrt(2)*I*s004/4 + sqrt(2)*I*s005/4], [-sqrt(2)*s004/4 + sqrt(2)*s005/4, 0, sqrt(2)*I*s004/4 - sqrt(2)*I*s005/4, 0], [0, -sqrt(2)*I*s004/4 + sqrt(2)*I*s005/4, 0, sqrt(2)*s004/4 - sqrt(2)*s005/4], [sqrt(2)*I*s004/4 - sqrt(2)*I*s005/4, 0, sqrt(2)*s004/4 - sqrt(2)*s005/4, 0]]",
        "z_046": "[[0, 0, I*s004/2 - I*s005/2, 0], [0, 0, 0, I*s004/2 - I*s005/2], [-I*s004/2 + I*s005/2, 0, 0, 0], [0, -I*s004/2 + I*s005/2, 0, 0]]",
        "z_047": "[[sqrt(2)*c006/2, 0, 0, 0], [0, sqrt(2)*c006/2, 0, 0], [0, 0, sqrt(2)*c006/2, 0], [0, 0, 0, sqrt(2)*c006/2]]",
        "z_048": "[[sqrt(2)*c006/2, 0, 0, 0], [0, sqrt(2)*c006/2, 0, 0], [0, 0, -sqrt(2)*c006/2, 0], [0, 0, 0, -sqrt(2)*c006/2]]",
        "z_049": "[[0, 0, -sqrt(2)*I*c006/2, 0], [0, 0, 0, sqrt(2)*I*c006/2], [sqrt(2)*I*c006/2, 0, 0, 0], [0, -sqrt(2)*I*c006/2, 0, 0]]",
        "z_050": "[[0, 0, 0, -sqrt(2)*I*c006/2], [0, 0, -sqrt(2)*I*c006/2, 0], [0, sqrt(2)*I*c006/2, 0, 0], [sqrt(2)*I*c006/2, 0, 0, 0]]",
        "z_051": "[[0, sqrt(38)*I*s006/19, 0, 3*sqrt(38)*s006/38], [-sqrt(38)*I*s006/19, 0, 3*sqrt(38)*s006/38, 0], [0, 3*sqrt(38)*s006/38, 0, -2*sqrt(38)*I*s006/19], [3*sqrt(38)*s006/38, 0, 2*sqrt(38)*I*s006/19, 0]]",
        "z_052": "[[0, -7*sqrt(19)*I*s006/38, 0, -sqrt(19)*s006/38], [7*sqrt(19)*I*s006/38, 0, -sqrt(19)*s006/38, 0], [0, -sqrt(19)*s006/38, 0, -5*sqrt(19)*I*s006/38], [-sqrt(19)*s006/38, 0, 5*sqrt(19)*I*s006/38, 0]]",
        "z_053": "[[0, 0, sqrt(2)*s006/2, 0], [0, 0, 0, -sqrt(2)*s006/2], [sqrt(2)*s006/2, 0, 0, 0], [0, -sqrt(2)*s006/2, 0, 0]]",
        "z_054": "[[0, I*s006/2, 0, -s006/2], [-I*s006/2, 0, -s006/2, 0], [0, -s006/2, 0, -I*s006/2], [-s006/2, 0, I*s006/2, 0]]",
        "z_055": "[[c007/2 + c008/2, 0, 0, 0], [0, c007/2 + c008/2, 0, 0], [0, 0, c007/2 + c008/2, 0], [0, 0, 0, c007/2 + c008/2]]",
        "z_056": "[[c007/2 + c008/2, 0, 0, 0], [0, c007/2 + c008/2, 0, 0], [0, 0, -c007/2 - c008/2, 0], [0, 0, 0, -c007/2 - c008/2]]",
        "z_057": "[[0, 0, c007/2 - c008/2, 0], [0, 0, 0, c007/2 - c008/2], [c007/2 - c008/2, 0, 0, 0], [0, c007/2 - c008/2, 0, 0]]",
        "z_058": "[[0, 0, -I*c007/2 - I*c008/2, 0], [0, 0, 0, I*c007/2 + I*c008/2], [I*c007/2 + I*c008/2, 0, 0, 0], [0, -I*c007/2 - I*c008/2, 0, 0]]",
        "z_059": "[[0, 0, 0, -I*c007/2 - I*c008/2], [0, 0, -I*c007/2 - I*c008/2, 0], [0, I*c007/2 + I*c008/2, 0, 0], [I*c007/2 + I*c008/2, 0, 0, 0]]",
        "z_060": "[[0, 0, 0, -c007/2 + c008/2], [0, 0, c007/2 - c008/2, 0], [0, c007/2 - c008/2, 0, 0], [-c007/2 + c008/2, 0, 0, 0]]",
        "z_061": "[[0, sqrt(19)*I*s007/19 + sqrt(19)*I*s008/19, 0, 3*sqrt(19)*s007/38 + 3*sqrt(19)*s008/38], [-sqrt(19)*I*s007/19 - sqrt(19)*I*s008/19, 0, 3*sqrt(19)*s007/38 + 3*sqrt(19)*s008/38, 0], [0, 3*sqrt(19)*s007/38 + 3*sqrt(19)*s008/38, 0, -2*sqrt(19)*I*s007/19 - 2*sqrt(19)*I*s008/19], [3*sqrt(19)*s007/38 + 3*sqrt(19)*s008/38, 0, 2*sqrt(19)*I*s007/19 + 2*sqrt(19)*I*s008/19, 0]]",
        "z_062": "[[s007/2 - s008/2, 0, 0, 0], [0, -s007/2 + s008/2, 0, 0], [0, 0, s007/2 - s008/2, 0], [0, 0, 0, -s007/2 + s008/2]]",
        "z_063": "[[0, 2*sqrt(19)*s007/19 - 2*sqrt(19)*s008/19, 0, -3*sqrt(19)*I*s007/38 + 3*sqrt(19)*I*s008/38], [2*sqrt(19)*s007/19 - 2*sqrt(19)*s008/19, 0, 3*sqrt(19)*I*s007/38 - 3*sqrt(19)*I*s008/38, 0], [0, -3*sqrt(19)*I*s007/38 + 3*sqrt(19)*I*s008/38, 0, -sqrt(19)*s007/19 + sqrt(19)*s008/19], [3*sqrt(19)*I*s007/38 - 3*sqrt(19)*I*s008/38, 0, -sqrt(19)*s007/19 + sqrt(19)*s008/19, 0]]",
        "z_064": "[[0, -7*sqrt(38)*I*s007/76 - 7*sqrt(38)*I*s008/76, 0, -sqrt(38)*s007/76 - sqrt(38)*s008/76], [7*sqrt(38)*I*s007/76 + 7*sqrt(38)*I*s008/76, 0, -sqrt(38)*s007/76 - sqrt(38)*s008/76, 0], [0, -sqrt(38)*s007/76 - sqrt(38)*s008/76, 0, -5*sqrt(38)*I*s007/76 - 5*sqrt(38)*I*s008/76], [-sqrt(38)*s007/76 - sqrt(38)*s008/76, 0, 5*sqrt(38)*I*s007/76 + 5*sqrt(38)*I*s008/76, 0]]",
        "z_065": "[[0, 5*sqrt(38)*s007/76 - 5*sqrt(38)*s008/76, 0, sqrt(38)*I*s007/76 - sqrt(38)*I*s008/76], [5*sqrt(38)*s007/76 - 5*sqrt(38)*s008/76, 0, -sqrt(38)*I*s007/76 + sqrt(38)*I*s008/76, 0], [0, sqrt(38)*I*s007/76 - sqrt(38)*I*s008/76, 0, 7*sqrt(38)*s007/76 - 7*sqrt(38)*s008/76], [-sqrt(38)*I*s007/76 + sqrt(38)*I*s008/76, 0, 7*sqrt(38)*s007/76 - 7*sqrt(38)*s008/76, 0]]",
        "z_066": "[[0, 0, s007/2 + s008/2, 0], [0, 0, 0, -s007/2 - s008/2], [s007/2 + s008/2, 0, 0, 0], [0, -s007/2 - s008/2, 0, 0]]",
        "z_067": "[[0, sqrt(2)*I*s007/4 + sqrt(2)*I*s008/4, 0, -sqrt(2)*s007/4 - sqrt(2)*s008/4], [-sqrt(2)*I*s007/4 - sqrt(2)*I*s008/4, 0, -sqrt(2)*s007/4 - sqrt(2)*s008/4, 0], [0, -sqrt(2)*s007/4 - sqrt(2)*s008/4, 0, -sqrt(2)*I*s007/4 - sqrt(2)*I*s008/4], [-sqrt(2)*s007/4 - sqrt(2)*s008/4, 0, sqrt(2)*I*s007/4 + sqrt(2)*I*s008/4, 0]]",
        "z_068": "[[-s007/2 + s008/2, 0, 0, 0], [0, s007/2 - s008/2, 0, 0], [0, 0, s007/2 - s008/2, 0], [0, 0, 0, -s007/2 + s008/2]]",
        "z_069": "[[0, -sqrt(2)*s007/4 + sqrt(2)*s008/4, 0, -sqrt(2)*I*s007/4 + sqrt(2)*I*s008/4], [-sqrt(2)*s007/4 + sqrt(2)*s008/4, 0, sqrt(2)*I*s007/4 - sqrt(2)*I*s008/4, 0], [0, -sqrt(2)*I*s007/4 + sqrt(2)*I*s008/4, 0, sqrt(2)*s007/4 - sqrt(2)*s008/4], [sqrt(2)*I*s007/4 - sqrt(2)*I*s008/4, 0, sqrt(2)*s007/4 - sqrt(2)*s008/4, 0]]",
        "z_070": "[[0, 0, I*s007/2 - I*s008/2, 0], [0, 0, 0, I*s007/2 - I*s008/2], [-I*s007/2 + I*s008/2, 0, 0, 0], [0, -I*s007/2 + I*s008/2, 0, 0]]",
        "z_071": "[[sqrt(2)*c009/2, 0, 0, 0], [0, sqrt(2)*c009/2, 0, 0], [0, 0, sqrt(2)*c009/2, 0], [0, 0, 0, sqrt(2)*c009/2]]",
        "z_072": "[[sqrt(2)*c009/2, 0, 0, 0], [0, sqrt(2)*c009/2, 0, 0], [0, 0, -sqrt(2)*c009/2, 0], [0, 0, 0, -sqrt(2)*c009/2]]",
        "z_073": "[[0, 0, -sqrt(2)*I*c009/2, 0], [0, 0, 0, sqrt(2)*I*c009/2], [sqrt(2)*I*c009/2, 0, 0, 0], [0, -sqrt(2)*I*c009/2, 0, 0]]",
        "z_074": "[[0, 0, 0, -sqrt(2)*I*c009/2], [0, 0, -sqrt(2)*I*c009/2, 0], [0, sqrt(2)*I*c009/2, 0, 0], [sqrt(2)*I*c009/2, 0, 0, 0]]",
        "z_075": "[[0, sqrt(38)*I*s009/19, 0, 3*sqrt(38)*s009/38], [-sqrt(38)*I*s009/19, 0, 3*sqrt(38)*s009/38, 0], [0, 3*sqrt(38)*s009/38, 0, -2*sqrt(38)*I*s009/19], [3*sqrt(38)*s009/38, 0, 2*sqrt(38)*I*s009/19, 0]]",
        "z_076": "[[0, -7*sqrt(19)*I*s009/38, 0, -sqrt(19)*s009/38], [7*sqrt(19)*I*s009/38, 0, -sqrt(19)*s009/38, 0], [0, -sqrt(19)*s009/38, 0, -5*sqrt(19)*I*s009/38], [-sqrt(19)*s009/38, 0, 5*sqrt(19)*I*s009/38, 0]]",
        "z_077": "[[0, 0, sqrt(2)*s009/2, 0], [0, 0, 0, -sqrt(2)*s009/2], [sqrt(2)*s009/2, 0, 0, 0], [0, -sqrt(2)*s009/2, 0, 0]]",
        "z_078": "[[0, I*s009/2, 0, -s009/2], [-I*s009/2, 0, -s009/2, 0], [0, -s009/2, 0, -I*s009/2], [-s009/2, 0, I*s009/2, 0]]",
    },
}
