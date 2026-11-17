# Te_ctrl.py: illustrative weights (not fitted), with the spin-orbit coupling.
Te_ctrl = {
    "samb": {
        "model": "Te",
        "parameter": {
            "z6": 1.0,  # on-site spin-orbit coupling.
            "z10": -2.0,  # intra-chain hopping.
            "z11": 1.0,  # intra-chain hopping.
            "z31": -0.5,  # inter-chain hopping.
        },
    },
    "output": {
        "dispersion": {
            "k_path": "Γ-M-K-Γ-A-L-H-A",
            "k_point": {
                "Γ": "[0,0,0]",
                "M": "[1/2,0,0]",
                "K": "[1/3,1/3,0]",
                "A": "[0,0,1/2]",
                "L": "[1/2,0,1/2]",
                "H": "[1/3,1/3,1/2]",
            },
            "local": ["Sz"],
        },
    },
}
