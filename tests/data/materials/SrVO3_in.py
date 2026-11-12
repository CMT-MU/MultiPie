# SrVO3_in.py: cubic perovskite (Pm-3m), V t2g and O p orbitals (Sr is omitted).
SrVO3_in = {
    "model": "SrVO3",
    "group": 221,  # Pm-3m.
    "cell": {"a": 3.84},
    #
    "site": {
        "V": ("[1/2,1/2,1/2]", ["dyz", "dxz", "dxy"]),  # 1b, t2g.
        "O": ("[1/2,1/2,0]", "p"),  # 3c.
    },
    "bond": [("V", "O", [1]), ("O", "O", [1]), ("V", "V", [1])],
    #
    "spinful": False,
    "SAMB_select": {"Gamma": "IR"},  # identity irrep, for the Hamiltonian.
}
