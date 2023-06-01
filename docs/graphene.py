"""
sample input file for graphene.
"""
graphene = {
    "model": "graphene",  #  name of model.
    "group": 191,  # No. of space group.
    "cell": {"c": 4},  # set large enough interlayer distance.
    #
    "site": {"C": ("[1/3,2/3,0]", "pz")},  # positions of C site and its orbital.
    "bond": [("C", "C", [1, 2, 3, 4, 5, 6])],  # C-C bonds up to 6th neighbors.
    #
    "spinful": False,  # spinless.
    #
    "k_point": {"Γ": "[0, 0, 0]", "M": "[1/2, 0, 0]", "K": "[1/3, 1/3, 0]"},  # def. of k points.
    "k_path": "Γ-K-M-Γ",  # high-symmetry line.
}
