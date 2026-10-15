# graphene_post_ctrl.py
graphene_post_ctrl = {
    "samb": {
        "model": "graphene",
        "parameter": {
            "z2": 1.0,
        },
    },
    "grid": (120, 120, 1),  # k grid for DOS, Fermi energy, and callback_k_full (2D: N3 = 1).
    "output": {
        "dos": "auto",  # partial DOS for each site cluster and orbital.
        "n_electron": 1,  # half filling (spinless: one state per band, two bands).
        "save_k_path": True,  # save Ek, Uk along k path in .npz.
        "callback_k_path": "graphene_callback.py",
        "callback_k_full": "graphene_callback.py",
        "dispersion": {
            "k_path": "Γ-M-K-Γ",
            "k_point": {"Γ": "[0,0,0]", "K": "[1/3,1/3,0]", "M": "[1/2,0,0]"},
            "weight": [["pz@C(1)"]],  # fat band of pz at C(1).
        },
    },
}
