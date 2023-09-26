graphene = {
    "model": "graphene",
    "group": 191,
    #
    "site": {"A": ("[1/3,2/3,0]", "pz")},
    "bond": [("A", "A", [1])],
    #
    "spinful": True,
    #
    "k_point": {"Γ": "[0, 0, 0]", "M": "[1/2, 0, 0]", "K": "[1/3, 1/3, 0]", "K'": "[-1/3, -1/3, 0]"},
    "k_path": "Γ-M-K-Γ-K'",
    #
    "option": {"view": [0, 0, 1]},
    "generate": {"fourier_transform": False},
    "detail": {"max_neighbor": 15, "cell_range": (-5, 5, -5, 5, -5, 5)},
}
