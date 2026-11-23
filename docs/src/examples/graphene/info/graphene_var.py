"""
Correspondence between zj and atomic variable.
- correspondence for each bond cluster, dict[bond_name, dict["bond" or zj, (head site, tail site) or expression in terms of atomic variables] ].
- "bond": head and tail sites of bond 1, (atom(sublattice), atom(sublattice)).
- only for SAMB with identity irrep.
- the matrix element of bond 1 between the head orbital a and the tail orbital b is H_{a,b} = g_{a,b} + i h_{a,b}.
- orbital: as it is (spinless), e.g., pxu = (px,u) (lgs), p(3/2,1/2) = (3/2,1/2) of p (jml).
"""

graphene_var = {"C;C_001_1": {"bond": ("C(2)", "C(1)"), "z2": "sqrt(6)*g_{pz,pz}"}, "C;C_002_1": {"bond": ("C(1)", "C(1)"), "z3": "2*sqrt(3)*g_{pz,pz}"}}
