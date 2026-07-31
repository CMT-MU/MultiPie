"""
Correspondence between zj and atomic variable.
- correspondence for each bond cluster, dict[bond_name, dict[zj, expression in terms of atomic variables] ].
- only for SAMB with identity irrep.
- atomic variable of (m,n) component at bond 1 is given by g(m,n) +i h(m,n).
"""

Te_var = {
    "A;A_001_1": {
        "z5": "sqrt(2)*(g_{px_A3,px_A2}+g_{py_A3,py_A2}+g_{pz_A3,pz_A2})",
        "z6": "-g_{px_A3,px_A2}-g_{py_A3,py_A2}+2*g_{pz_A3,pz_A2}",
        "z7": "sqrt(3)*(g_{px_A3,px_A2}-g_{py_A3,py_A2})",
        "z8": "g_{px_A3,py_A2}-sqrt(2)*g_{px_A3,pz_A2}-g_{py_A3,px_A2}+sqrt(2)*g_{pz_A3,px_A2}",
        "z9": "sqrt(3)*(-g_{py_A3,pz_A2}-g_{pz_A3,py_A2})",
        "z10": "sqrt(2)*g_{px_A3,py_A2}+g_{px_A3,pz_A2}-sqrt(2)*g_{py_A3,px_A2}-g_{pz_A3,px_A2}",
    }
}
