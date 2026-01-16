# graphene.py
"""
input file for graphene.
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
}
