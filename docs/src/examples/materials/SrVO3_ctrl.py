# SrVO3_ctrl.py: illustrative weights (not fitted), on-site energies and nearest-neighbor V-O hopping.
SrVO3_ctrl = {
    "samb": {
        "model": "SrVO3",
        "parameter": {
            "z1": -12.0,  # O p on-site (uniform).
            "z3": 0.0,  # V t2g on-site.
            "z8": 3.0,  # V-O hopping (t2g-p pi bond).
            "z4": 0.5,  # O-O nearest-neighbor hopping.
        },
    },
    "output": {
        "dispersion": {
            "k_path": "Γ-X-M-Γ-R-X",
            "k_point": {"Γ": "[0,0,0]", "X": "[1/2,0,0]", "M": "[1/2,1/2,0]", "R": "[1/2,1/2,1/2]"},
        }
    },
}
