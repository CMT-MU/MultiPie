"""
Correspondence between zj and atomic variable.
- correspondence for each bond cluster, dict[bond_name, dict[zj, expression in terms of atomic variables] ].
- only for SAMB with identity irrep.
"""

graphene_var = {"C;C_001_1": {"z2": "sqrt(6)*g_{(pz_C1,pz_C2)}"}, "C;C_002_1": {"z3": "sqrt(3)*g_{(pz_C1,pz_C1)}"}}
