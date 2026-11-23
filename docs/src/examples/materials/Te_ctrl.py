# Te_ctrl.py: illustrative weights (not fitted), with the symmetric and antisymmetric spin-orbit couplings.
Te_ctrl = {
    "samb": {
        "model": "Te",
        "parameter": {
            "z9": 1.0,  # on-site (symmetric) spin-orbit coupling, l.s.
            "z10": -2.0,  # intra-chain hopping.
            "z13": 1.0,  # intra-chain hopping.
            "z16": 0.3,  # intra-chain antisymmetric spin-orbit coupling (sigma_z).
            "z17": 0.3,  # intra-chain antisymmetric spin-orbit coupling (sigma_x, sigma_y).
            "z31": -0.5,  # inter-chain hopping.
        },
    },
    "output": {
        "dispersion": {
            "k_path": "Γ-M-K-Γ|A'-Γ-A-L-H-A",  # A'-Γ-A: k_z and -k_z.
            "k_point": {
                "Γ": "[0,0,0]",
                "A'": "[0,0,-1/2]",
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
