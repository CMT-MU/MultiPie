Te = {
    "model": "Te",
    "group": 152,
    "cell": {"a": 4.4580000, "b": 4.4580000, "c": 5.925},
    #
    "site": {"A": ("[0.274,0,1/3]", ["px", "py", "pz"])},
    "bond": [("A", "A", [1])],
    #
    "spinful": False,
    #
    "k_point": {
        "Γ": "[0, 0, 0]",
        "A": "[0, 0, 1/2]",
        "M": "[1/2, 0, 0]",
        "K": "[1/3, 1/3, 0]",
        "H": "[1/3, 1/3, 1/2]",
        "L": "[1/2, 0, 1/2]",
    },
    "k_path": "A-Γ-H-A-L-H-K-Γ-M-K",
    #
    "generate": {"model_type": "phonon", "fourier_transform": False},
    "detail": {"max_neighbor": 15, "cell_range": (-5, 5, -5, 5, -5, 5)},
}
