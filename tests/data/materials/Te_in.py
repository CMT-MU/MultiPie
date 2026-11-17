# Te_in.py: trigonal Te (P3_121, chiral D3), spinful p orbitals.
Te_in = {
    "model": "Te",
    "group": 152,  # P3_121.
    "cell": {"a": 4.45, "c": 5.93},
    #
    "site": {"Te": ("[0.2636,0,1/3]", "p")},  # 3a.
    "bond": [("Te", "Te", [1, 2])],  # intra- and inter-chain.
    #
    "spinful": True,
    "SAMB_select": {"Gamma": "IR"},  # identity irrep, for the Hamiltonian (Q and G).
}
