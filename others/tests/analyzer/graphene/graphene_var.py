"""
Correspondence between zj and atomic variable.
- correspondence for each bond cluster, dict[bond_name, dict[zj, expression in terms of atomic variables] ].
- only for SAMB with identity irrep.
"""

graphene_var = {
    "A;A_001_1": {"z2": "sqrt(6)*g_{(pz_A2,pz_A1)}"},
    "A;A_002_1": {"z3": "2*sqrt(3)*g_{(pz_A1,pz_A1)}"},
    "A;A_003_1": {"z4": "sqrt(6)*g_{(pz_A2,pz_A1)}"},
    "A;A_004_1": {"z5": "2*sqrt(3)*g_{(pz_A2,pz_A1)}"},
    "A;A_005_1": {"z6": "2*sqrt(3)*g_{(pz_A1,pz_A1)}"},
    "A;A_006_1": {"z7": "2*sqrt(3)*g_{(pz_A1,pz_A1)}"},
    "A;A_007_1": {"z8": "2*sqrt(3)*g_{(pz_A2,pz_A1)}"},
    "A;A_008_1": {"z9": "sqrt(6)*g_{(pz_A2,pz_A1)}"},
    "A;A_009_1": {"z10": "2*sqrt(3)*g_{(pz_A2,pz_A1)}"},
    "A;A_010_1": {"z11": "2*sqrt(6)*g_{(pz_A1,pz_A1)}"},
}
