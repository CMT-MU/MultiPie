# MoS2_in.py: monolayer 2H-MoS2 (P-6m2, D3h), Mo d and S p orbitals.
MoS2_in = {
    "model": "MoS2",
    "group": 187,  # P-6m2.
    "cell": {"a": 3.16, "c": 20.0},  # vacuum along c.
    #
    "site": {
        "Mo": ("[0,0,0]", "d"),  # 1a.
        "S": ("[1/3,2/3,0.078]", "p"),  # 2h, z = 1.56/20.
    },
    "bond": [("Mo", "S", [1]), ("Mo", "Mo", [1]), ("S", "S", [1, 2])],
    #
    "spinful": False,
    "SAMB_select": {"Gamma": "IR"},  # identity irrep, for the Hamiltonian.
}
