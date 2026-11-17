# MoS2_ctrl.py: illustrative weights (not fitted).
MoS2_ctrl = {
    "samb": {
        "model": "MoS2",
        "parameter": {
            "z2": 1.0,  # Mo d crystal field.
            "z4": -6.0,  # S p on-site.
            "z15": 3.0,  # Mo-S hopping.
            "z16": 1.0,  # Mo-S hopping.
            "z6": -0.5,  # Mo-Mo hopping.
        },
    },
    "output": {"dispersion": {"k_path": "Γ-M-K-Γ", "k_point": {"Γ": "[0,0,0]", "M": "[1/2,0,0]", "K": "[1/3,1/3,0]"}}},
}
